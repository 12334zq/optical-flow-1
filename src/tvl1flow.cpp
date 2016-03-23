
// This program is free software: you can use, modify and/or redistribute it
// under the terms of the simplified BSD License. You should have received a
// copy of this license along this program. If not, see
// <http://www.opensource.org/licenses/bsd-license.html>.
//
// Copyright (C) 2011, Javier Sánchez Pérez <jsanchez@dis.ulpgc.es>
// All rights reserved.


#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

#include "operators.h"
#include "bicubic_interpolation.h"
#include "zoom.h"
#include "utils.h"
#include "tvl1flow.h"

#define MAX_ITERATIONS 300
#define PRESMOOTHING_SIGMA 0.8
#define GRAD_IS_ZERO 1E-10

/**
 * Implementation of the Zach, Pock and Bischof dual TV-L1 optic flow method
 *
 * see reference:
 *  [1] C. Zach, T. Pock and H. Bischof, "A Duality Based Approach for Realtime
 *      TV-L1 Optical Flow", In Proceedings of Pattern Recognition (DAGM),
 *      Heidelberg, Germany, pp. 214-223, 2007
 *
 *
 * Details on the total variation minimization scheme can be found in:
 *  [2] A. Chambolle, "An Algorithm for Total Variation Minimization and
 *      Applications", Journal of Mathematical Imaging and Vision, 20: 89-97, 2004
 **/


/**
 *
 * Function to compute the optical flow in one scale
 *
 **/
void
Dual_TVL1_optic_flow(ofpix_t *I0, // source image
                     ofpix_t *I1, // target image
                     ofpix_t *u1, // x component of the optical flow
                     ofpix_t *u2, // y component of the optical flow
                     const int nx, // image width
                     const int ny, // image height
                     const double tau, // time step
                     const double lambda, // weight parameter for the data term
                     const double theta, // weight parameter for (u - v)²
                     const int warps, // number of warpings per scale
                     const double epsilon, // tolerance for numerical convergence
                     const bool verbose // enable/disable the verbose mode
                     )
{
    const int size = nx * ny;
    const double l_t = lambda * theta;
    ofpix_t *I1x    = new ofpix_t[size];
    ofpix_t *I1y    = new ofpix_t[size];
    ofpix_t *I1w    = new ofpix_t[size];
    ofpix_t *I1wx   = new ofpix_t[size];
    ofpix_t *I1wy   = new ofpix_t[size];
    ofpix_t *rho_c  = new ofpix_t[size];
    ofpix_t *v1     = new ofpix_t[size];
    ofpix_t *v2     = new ofpix_t[size];
    ofpix_t *p11    = new ofpix_t[size];
    ofpix_t *p12    = new ofpix_t[size];
    ofpix_t *p21    = new ofpix_t[size];
    ofpix_t *p22    = new ofpix_t[size];
    ofpix_t *div    = new ofpix_t[size];
    ofpix_t *grad   = new ofpix_t[size];
    ofpix_t *div_p1 = new ofpix_t[size];
    ofpix_t *div_p2 = new ofpix_t[size];
    ofpix_t *u1x    = new ofpix_t[size];
    ofpix_t *u1y    = new ofpix_t[size];
    ofpix_t *u2x    = new ofpix_t[size];
    ofpix_t *u2y    = new ofpix_t[size];

    centered_gradient(I1, I1x, I1y, nx, ny);

    // initialization of p
    for (int i = 0; i < size; i++) {
        p11[i] = p12[i] = 0.0;
        p21[i] = p22[i] = 0.0;
    }

    for (int warpings = 0; warpings < warps; warpings++) {
        // compute the warping of the target image and its derivatives
        bicubic_interpolation_warp(I1,  u1, u2, I1w,  nx, ny, true);
        bicubic_interpolation_warp(I1x, u1, u2, I1wx, nx, ny, true);
        bicubic_interpolation_warp(I1y, u1, u2, I1wy, nx, ny, true);

#pragma omp parallel for
        for (int i = 0; i < size; i++) {
            const double Ix2 = I1wx[i] * I1wx[i];
            const double Iy2 = I1wy[i] * I1wy[i];

            // store the |Grad(I1)|^2
            grad[i] = (Ix2 + Iy2);

            // compute the constant part of the rho function
            rho_c[i] = (I1w[i] - I1wx[i] * u1[i]
                        - I1wy[i] * u2[i] - I0[i]);
        }

        int n = 0;
        double error = INFINITY;
        while (error > epsilon * epsilon && n < MAX_ITERATIONS) {
            n++;
            // estimate the values of the variable (v1, v2)
            // (thresholding opterator TH)
#pragma omp parallel for
            for (int i = 0; i < size; i++) {
                const double rho = rho_c[i]
                                  + (I1wx[i] * u1[i] + I1wy[i] * u2[i]);
                double d1, d2;

                if (rho < -l_t * grad[i]) {
                    d1 = l_t * I1wx[i];
                    d2 = l_t * I1wy[i];
                } else   {
                    if (rho > l_t * grad[i]) {
                        d1 = -l_t * I1wx[i];
                        d2 = -l_t * I1wy[i];
                    } else   {
                        if (grad[i] < GRAD_IS_ZERO) {
                            d1 = d2 = 0;
                        } else {
                            double fi = -rho / grad[i];
                            d1 = fi * I1wx[i];
                            d2 = fi * I1wy[i];
                        }
                    }
                }

                v1[i] = u1[i] + d1;
                v2[i] = u2[i] + d2;
            }

            // compute the divergence of the dual variable (p1, p2)
            divergence(p11, p12, div_p1, nx, ny);
            divergence(p21, p22, div_p2, nx, ny);

            // estimate the values of the optical flow (u1, u2)
            error = 0.0;
#pragma omp parallel for reduction(+:error)
            for (int i = 0; i < size; i++) {
                const double u1k = u1[i];
                const double u2k = u2[i];

                u1[i] = v1[i] + theta * div_p1[i];
                u2[i] = v2[i] + theta * div_p2[i];

                error += (u1[i] - u1k) * (u1[i] - u1k) +
                         (u2[i] - u2k) * (u2[i] - u2k);
            }
            error /= size;

            // compute the gradient of the optical flow (Du1, Du2)
            forward_gradient(u1, u1x, u1y, nx, ny);
            forward_gradient(u2, u2x, u2y, nx, ny);

            // estimate the values of the dual variable (p1, p2)
#pragma omp parallel for
            for (int i = 0; i < size; i++) {
                const double taut = tau / theta;
                const double g1   = hypot(u1x[i], u1y[i]);
                const double g2   = hypot(u2x[i], u2y[i]);
                const double ng1  = 1.0 + taut * g1;
                const double ng2  = 1.0 + taut * g2;

                p11[i] = (p11[i] + taut * u1x[i]) / ng1;
                p12[i] = (p12[i] + taut * u1y[i]) / ng1;
                p21[i] = (p21[i] + taut * u2x[i]) / ng2;
                p22[i] = (p22[i] + taut * u2y[i]) / ng2;
            }
        }

        if (verbose) {
            fprintf(stderr, "Warping: %d, "
                    "Iterations: %d, "
                    "Error: %f\n", warpings, n, error);
        }
    }

    // delete allocated memory
    delete [] I1x;
    delete [] I1y;
    delete [] I1w;
    delete [] I1wx;
    delete [] I1wy;
    delete [] rho_c;
    delete [] v1;
    delete [] v2;
    delete [] p11;
    delete [] p12;
    delete [] p21;
    delete [] p22;
    delete [] div;
    delete [] grad;
    delete [] div_p1;
    delete [] div_p2;
    delete [] u1x;
    delete [] u1y;
    delete [] u2x;
    delete [] u2y;
} // Dual_TVL1_optic_flow

/**
 *
 * Function to compute the optical flow using multiple scales
 *
 **/
void
Dual_TVL1_optic_flow_multiscale(ofpix_t *I0, // source image
                                ofpix_t *I1, // target image
                                ofpix_t *u1, // x component of the optical flow
                                ofpix_t *u2, // y component of the optical flow
                                const int nxx, // image width
                                const int nyy, // image height
                                const double tau, // time step
                                const double lambda, // weight parameter for the data term
                                const double theta, // weight parameter for (u - v)²
                                const int nscales, // number of scales
                                const double zfactor, // factor for building the image piramid
                                const int warps, // number of warpings per scale
                                const double epsilon, // tolerance for numerical convergence
                                const bool verbose // enable/disable the verbose mode
                                )
{
    int size = nxx * nyy;

    // allocate memory for the pyramid structure
    ofpix_t **I0s = new ofpix_t* [nscales];
    ofpix_t **I1s = new ofpix_t* [nscales];
    ofpix_t **u1s = new ofpix_t* [nscales];
    ofpix_t **u2s = new ofpix_t* [nscales];
    int    *nx  = new int [nscales];
    int    *ny  = new int [nscales];

    I0s[0] = new ofpix_t[size];
    I1s[0] = new ofpix_t[size];

    u1s[0] = u1;
    u2s[0] = u2;
    nx [0] = nxx;
    ny [0] = nyy;

    // normalize the images between 0 and 255
    image_normalization_2(I0, I1, I0s[0], I1s[0], size);

    // pre-smooth the original images
    gaussian(I0s[0], nx[0], ny[0], PRESMOOTHING_SIGMA);
    gaussian(I1s[0], nx[0], ny[0], PRESMOOTHING_SIGMA);

    // create the scales
    for (int s = 1; s < nscales; s++) {
        zoom_size(nx[s - 1], ny[s - 1], &nx[s], &ny[s], zfactor);
        const int sizes = nx[s] * ny[s];

        // allocate memory
        I0s[s] = new ofpix_t[sizes];
        I1s[s] = new ofpix_t[sizes];
        u1s[s] = new ofpix_t[sizes];
        u2s[s] = new ofpix_t[sizes];

        // zoom in the images to create the pyramidal structure
        zoom_out(I0s[s - 1], I0s[s], nx[s - 1], ny[s - 1], zfactor);
        zoom_out(I1s[s - 1], I1s[s], nx[s - 1], ny[s - 1], zfactor);
    }

    // initialize the flow at the coarsest scale
    for (int i = 0; i < nx[nscales - 1] * ny[nscales - 1]; i++) {
        u1s[nscales - 1][i] = u2s[nscales - 1][i] = 0.0;
    }

    // pyramidal structure for computing the optical flow
    for (int s = nscales - 1; s >= 0; s--) {
        if (verbose) {
            fprintf(stderr, "Scale %d: %dx%d\n", s, nx[s], ny[s]);
        }

        // compute the optical flow at the current scale
        Dual_TVL1_optic_flow(
            I0s[s], I1s[s], u1s[s], u2s[s], nx[s], ny[s],
            tau, lambda, theta, warps, epsilon, verbose
            );

        // if this was the last scale, finish now
        if (!s) {
            break;
        }

        // otherwise, upsample the optical flow

        // zoom the optical flow for the next finer scale
        zoom_in(u1s[s], u1s[s - 1], nx[s], ny[s], nx[s - 1], ny[s - 1]);
        zoom_in(u2s[s], u2s[s - 1], nx[s], ny[s], nx[s - 1], ny[s - 1]);

        // scale the optical flow with the appropriate zoom factor
        for (int i = 0; i < nx[s - 1] * ny[s - 1]; i++) {
            u1s[s - 1][i] *= 1.0 / zfactor;
            u2s[s - 1][i] *= 1.0 / zfactor;
        }
    }

    // delete allocated memory
    for (int i = 1; i < nscales; i++) {
        delete [] I0s[i];
        delete [] I1s[i];
        delete [] u1s[i];
        delete [] u2s[i];
    }
    delete [] I0s[0];
    delete [] I1s[0];

    delete [] I0s;
    delete [] I1s;
    delete [] u1s;
    delete [] u2s;
    delete [] nx;
    delete [] ny;
} // Dual_TVL1_optic_flow_multiscale

