// This program is free software: you can use, modify and/or redistribute it
// under the terms of the simplified BSD License. You should have received a
// copy of this license along this program. If not, see
// <http://www.opensource.org/licenses/bsd-license.html>.
//
// Copyright (C) 2012, Javier Sánchez Pérez <jsanchez@dis.ulpgc.es>
// All rights reserved.

#include "brox_optic_flow.h"

//#include <omp.h>
#include <vector>
#include <iostream>
#include <cmath>
#include <algorithm>

#include "operators.h"
#include "zoom.h"
#include "bicubic_interpolation.h"
#include "brox_spatial_mask.h"
#include "utils.h"

#define EPSILON 0.001
#define MAXITER 300
#define SOR_PARAMETER 1.9
#define GAUSSIAN_SIGMA 0.8

/**
 *
 * Compute the coefficients of the robust functional (data term)
 *
 **/
static
void
psi_data(const ofpix_t *I1, //first image
         const ofpix_t *I2, //second image
         const ofpix_t *I2x, //gradient of the second image
         const ofpix_t *I2y, //gradient of the second image
         const ofpix_t *du, //motion increment
         const ofpix_t *dv, //motion increment
         ofpix_t *psip, //output coefficients
         const int nx, //image width
         const int ny //image height
         )
{
    const int size = nx * ny;

    //compute 1/(sqrt((I2-I1+I2x*du+I2y*dv)²+e²) in each pixel
    //(equation (5) in the article)
#pragma omp parallel for
    for (int i = 0; i < size; i++) {
        const double dI  = I2[i] - I1[i] + I2x[i] * du[i] + I2y[i] * dv[i];
        const double dI2 = dI * dI;

        psip[i] = 1. / sqrt(dI2 + EPSILON * EPSILON);
    }
}

/**
 *
 * Compute the coefficients of the robust functional (gradient term)
 *
 **/
static
void
psi_gradient(const ofpix_t *I1x, //gradient of the first image
             const ofpix_t *I1y, //gradient of the first image
             const ofpix_t *I2x, //gradient of the second image
             const ofpix_t *I2y, //gradient of the second image
             const ofpix_t *I2xx, //second derivatives of the second image
             const ofpix_t *I2xy, //second derivatives of the second image
             const ofpix_t *I2yy, //second derivatives of the second image
             const ofpix_t *du, //motion increment
             const ofpix_t *dv, //motion increment
             ofpix_t *psip, //output coefficients
             const int nx, //image width
             const int ny //image height
             )
{
    const int size = nx * ny;

    //compute 1/(sqrt(|DI2-DI1+HI2*(du,dv)|²+e²) in each pixel
    //(equation (5) in the article)
#pragma omp parallel for
    for (int i = 0; i < size; i++) {
        const double dIx = I2x[i] - I1x[i] + I2xx[i] * du[i] + I2xy[i] * dv[i];
        const double dIy = I2y[i] - I1y[i] + I2xy[i] * du[i] + I2yy[i] * dv[i];
        const double dI2 = dIx * dIx + dIy * dIy;

        psip[i] = 1. / sqrt(dI2 + EPSILON * EPSILON);
    }
}

/**
 *
 * Compute the coefficients of the robust functional (smoothness term)
 *
 **/
static
void
psi_smooth(const ofpix_t *ux, //gradient of x component of the optical flow
           const ofpix_t *uy, //gradient of x component of the optical flow
           const ofpix_t *vx, //gradient of y component of the optical flow
           const ofpix_t *vy, //gradient of y component of the optical flow
           ofpix_t *psi, //output coefficients
           const int nx, //image width
           const int ny //image height
           )
{
    const int size = nx * ny;

    //compute 1/(sqrt(ux²+uy²+vx²+vy²+e²) in each pixel
    //(equation (5) in the article)
#pragma omp parallel for
    for (int i = 0; i < size; i++) {
        const double du  = ux[i] * ux[i] + uy[i] * uy[i];
        const double dv  = vx[i] * vx[i] + vy[i] * vy[i];
        const double d2  = du + dv;

        psi[i] = 1. / sqrt(d2 + EPSILON * EPSILON);
    }
}

/**
 *
 *  SOR iteration in one position
 *
 */
static
inline double
sor_iteration(const ofpix_t *Au, //constant part of the numerator of u
              const ofpix_t *Av, //constant part of the numerator of v
              const ofpix_t *Du, //denominator of u
              const ofpix_t *Dv, //denominator of v
              const ofpix_t *D, //constant part of the numerator
              ofpix_t       *du,//x component of the motion increment
              ofpix_t       *dv,//y component of the motion increment
              const double alpha, //alpha smoothness parameter
              const ofpix_t *psi1, //coefficients of the divergence
              const ofpix_t *psi2,
              const ofpix_t *psi3,
              const ofpix_t *psi4,
              const int i, //current row
              const int i0, //previous row
              const int i1, //following row
              const int j, //current column
              const int nx, //number of columns
              const int j0, //previous column
              const int j1 //following column
              )
{
    //set the SOR extrapolation parameter
    const double w = SOR_PARAMETER;

    //calculate the position in the array
    const int k = i * nx + j;

    //compute the divergence part of the numerator (equation (10))
    const double div_du = psi1[k] * du[k + i1] + psi2[k] * du[k - i0] +
                          psi3[k] * du[k + j1] + psi4[k] * du[k - j0];
    const double div_dv = psi1[k] * dv[k + i1] + psi2[k] * dv[k - i0] +
                          psi3[k] * dv[k + j1] + psi4[k] * dv[k - j0];
    const double duk = du[k];
    const double dvk = dv[k];

    //update the motion increment (equation (12))
    du[k] = (1. - w) * du[k] + w * (Au[k] - D[k] * dv[k] + alpha * div_du) / Du[k];
    dv[k] = (1. - w) * dv[k] + w * (Av[k] - D[k] * du[k] + alpha * div_dv) / Dv[k];

    //return the covergence error in this position (equation (13))
    return (du[k] - duk) * (du[k] - duk) + (dv[k] - dvk) * (dv[k] - dvk);
}

/**
 *
 * Compute the optic flow with the Brox spatial method
 *
 **/
static
void
brox_optic_flow(   const ofpix_t *I1,//first image
                   const ofpix_t *I2, //second image
                   ofpix_t *u, //x component of the optical flow
                   ofpix_t *v, //y component of the optical flow
                   const int nx, //image width
                   const int ny, //image height
                   const double alpha, //smoothness parameter
                   const double gamma, //gradient term parameter
                   const double TOL, //stopping criterion threshold
                   const int inner_iter, //number of inner iterations
                   const int outer_iter, //number of outer iterations
                   const int number_of_threads, // number of threads for the parallel code
                   const bool verbose //switch on messages
                   )
{
    const int size = nx * ny;

    //allocate memory
    ofpix_t *du    = new ofpix_t[size];
    ofpix_t *dv    = new ofpix_t[size];
    ofpix_t *ux    = new ofpix_t[size];
    ofpix_t *uy    = new ofpix_t[size];
    ofpix_t *vx    = new ofpix_t[size];
    ofpix_t *vy    = new ofpix_t[size];
    ofpix_t *I1x   = new ofpix_t[size];
    ofpix_t *I1y   = new ofpix_t[size];
    ofpix_t *I2x   = new ofpix_t[size];
    ofpix_t *I2y   = new ofpix_t[size];
    ofpix_t *I2w   = new ofpix_t[size];
    ofpix_t *I2wx  = new ofpix_t[size];
    ofpix_t *I2wy  = new ofpix_t[size];
    ofpix_t *I2xx  = new ofpix_t[size];
    ofpix_t *I2yy  = new ofpix_t[size];
    ofpix_t *I2xy  = new ofpix_t[size];
    ofpix_t *I2wxx = new ofpix_t[size];
    ofpix_t *I2wyy = new ofpix_t[size];
    ofpix_t *I2wxy = new ofpix_t[size];
    ofpix_t *div_u = new ofpix_t[size];
    ofpix_t *div_v = new ofpix_t[size];
    ofpix_t *div_d = new ofpix_t[size];
    ofpix_t *Au    = new ofpix_t[size];
    ofpix_t *Av    = new ofpix_t[size];
    ofpix_t *Du    = new ofpix_t[size];
    ofpix_t *Dv    = new ofpix_t[size];
    ofpix_t *D     = new ofpix_t[size];
    ofpix_t *psid  = new ofpix_t[size];
    ofpix_t *psig  = new ofpix_t[size];
    ofpix_t *psis  = new ofpix_t[size];
    ofpix_t *psi1  = new ofpix_t[size];
    ofpix_t *psi2  = new ofpix_t[size];
    ofpix_t *psi3  = new ofpix_t[size];
    ofpix_t *psi4  = new ofpix_t[size];

    //compute the gradient of the images
    centered_gradient(I1, I1x, I1y, nx, ny, 1);
    centered_gradient(I2, I2x, I2y, nx, ny, 1);

    //compute second order derivatives
    Dxx(I2, I2xx, nx, ny, 1);
    Dyy(I2, I2yy, nx, ny, 1);
    Dxy(I2, I2xy, nx, ny, 1);

    //outer iterations loop
    for (int no = 0; no < outer_iter; no++) {
        //warp the second image and its derivatives
        bicubic_interpolation_warp(I2,   u, v, I2w,   nx, ny, true);
        bicubic_interpolation_warp(I2x,  u, v, I2wx,  nx, ny, true);
        bicubic_interpolation_warp(I2y,  u, v, I2wy,  nx, ny, true);
        bicubic_interpolation_warp(I2xx, u, v, I2wxx, nx, ny, true);
        bicubic_interpolation_warp(I2xy, u, v, I2wxy, nx, ny, true);
        bicubic_interpolation_warp(I2yy, u, v, I2wyy, nx, ny, true);

        //compute the flow gradient
        centered_gradient(u, ux, uy, nx, ny, 1);
        centered_gradient(v, vx, vy, nx, ny, 1);

        //compute robust function Psi for the smoothness term
        psi_smooth(ux, uy, vx, vy, psis, nx, ny);

        //compute coefficients of Psi functions in divergence
        brox_spatial_psi_divergence(psis, psi1, psi2, psi3, psi4, nx, ny);

        //compute the divergence for the gradient of w (equation (8))
        brox_spatial_divergence_u(u, v, psi1, psi2, psi3, psi4, div_u, div_v, nx, ny);

#pragma omp parallel for
        for (int i = 0; i < size; i++) {
            //compute the coefficents of dw[i] in the smoothness term
            //(equation (10))
            div_d[i] = alpha * (psi1[i] + psi2[i] + psi3[i] + psi4[i]);

            //initialize the motion increment
            du[i] = dv[i] = 0;
        }

        //inner iterations loop
        for (int ni = 0; ni < inner_iter; ni++) {
            //compute robust function Psi for the data and gradient terms
            psi_data(I1, I2w, I2wx, I2wy, du, dv,  psid, nx, ny);
            psi_gradient(I1x, I1y, I2wx, I2wy, I2wxx, I2wxy, I2wyy, du, dv, psig, nx, ny);

            //store constant parts of the numerical scheme (equation (11))
            for (int i = 0; i < size; i++) {
                const ofpix_t p = psid[i];
                const double g = gamma * psig[i];

                //brightness constancy term
                const double dif = I2w[i] - I1[i];
                const double BNu = -p * dif * I2wx[i];
                const double BNv = -p * dif * I2wy[i];
                const double BDu = p * I2wx[i] * I2wx[i];
                const double BDv = p * I2wy[i] * I2wy[i];

                //gradient constancy term
                const double dx  = (I2wx[i] - I1x[i]);
                const double dy  = (I2wy[i] - I1y[i]);
                const double GNu = -g * (dx * I2wxx[i] + dy * I2wxy[i]);
                const double GNv = -g * (dx * I2wxy[i] + dy * I2wyy[i]);
                const double GDu =  g * (I2wxx[i] * I2wxx[i] + I2wxy[i] * I2wxy[i]);
                const double GDv =  g * (I2wyy[i] * I2wyy[i] + I2wxy[i] * I2wxy[i]);
                const double DI  = (I2wxx[i] + I2wyy[i]) * I2wxy[i];
                const double Duv =  p * I2wy[i] * I2wx[i] + g * DI;

                Au[i] = BNu + GNu + alpha * div_u[i];
                Av[i] = BNv + GNv + alpha * div_v[i];
                Du[i] = BDu + GDu + div_d[i];
                Dv[i] = BDv + GDv + div_d[i];
                D [i] = Duv;
            }

            //sor iterations loop
            double error = 1000;
            int nsor = 0;

            while (error > TOL && nsor < MAXITER) {
                error = 0;
                nsor++;

                //update the motion increment in the center of the images
#pragma omp parallel for reduction(+:error) num_threads((number_of_threads < (ny-3)) ? number_of_threads : (ny-3))
                for (int i = 1; i < ny - 1; i++) {
                    for (int j = 1; j < nx - 1; j++) {
                        error += sor_iteration(
                            Au, Av, Du, Dv, D, du, dv, alpha,
                            psi1, psi2, psi3, psi4,
                            i, nx, nx, j, nx, 1, 1
                            );
                    }
                }

                //update the motion increment in the first and last rows
                for (int j = 1; j < nx - 1; j++) {
                    error += sor_iteration(
                        Au, Av, Du, Dv, D, du, dv, alpha,
                        psi1, psi2, psi3, psi4,
                        0, 0, nx, j, nx, 1, 1
                        );

                    error += sor_iteration(
                        Au, Av, Du, Dv, D, du, dv, alpha,
                        psi1, psi2, psi3, psi4,
                        ny - 1, nx, 0, j, nx, 1, 1
                        );
                }

                //update the motion increment in the first and last columns
                for (int i = 1; i < ny - 1; i++) {
                    error += sor_iteration(
                        Au, Av, Du, Dv, D, du, dv, alpha,
                        psi1, psi2, psi3, psi4,
                        i, nx, nx, 0, nx, 0, 1
                        );

                    error += sor_iteration(
                        Au, Av, Du, Dv, D, du, dv, alpha,
                        psi1, psi2, psi3, psi4,
                        i, nx, nx, nx - 1, nx, 1, 0
                        );
                }

                //process the top-left corner (0,0)
                error += sor_iteration(
                    Au, Av, Du, Dv, D, du, dv, alpha,
                    psi1, psi2, psi3, psi4,
                    0, 0, nx, 0, nx, 0, 1
                    );

                //process the top-right corner (0,nx-1)
                error += sor_iteration(
                    Au, Av, Du, Dv, D, du, dv, alpha,
                    psi1, psi2, psi3, psi4,
                    0, 0, nx, nx - 1, nx, 1, 0
                    );

                //process the bottom-left corner (ny-1,0)
                error += sor_iteration(
                    Au, Av, Du, Dv, D, du, dv, alpha,
                    psi1, psi2, psi3, psi4,
                    ny - 1, nx, 0, 0, nx, 0, 1
                    );

                //process the bottom-right corner (ny-1,nx-1)
                error += sor_iteration(
                    Au, Av, Du, Dv, D, du, dv, alpha,
                    psi1, psi2, psi3, psi4,
                    ny - 1, nx, 0, nx - 1, nx, 1, 0
                    );

                error = sqrt(error / size);
            }

            if (verbose) {
                std::cout << "Iterations: " << nsor << std::endl;
            }
        }

        //update the flow with the estimated motion increment
        for (int i = 0; i < size; i++) {
            u[i] += du[i];
            v[i] += dv[i];
        }
    }

    //delete allocated memory
    delete []du;
    delete []dv;

    delete []ux;
    delete []uy;
    delete []vx;
    delete []vy;

    delete []I1x;
    delete []I1y;
    delete []I2x;
    delete []I2y;
    delete []I2w;
    delete []I2wx;
    delete []I2wy;
    delete []I2xx;
    delete []I2yy;
    delete []I2xy;
    delete []I2wxx;
    delete []I2wyy;
    delete []I2wxy;

    delete []div_u;
    delete []div_v;
    delete []div_d;

    delete []Au;
    delete []Av;
    delete []Du;
    delete []Dv;
    delete []D;

    delete []psid;
    delete []psig;
    delete []psis;
    delete []psi1;
    delete []psi2;
    delete []psi3;
    delete []psi4;
} // brox_optic_flow

/**
 *
 *  Multiscale approach for computing the optical flow
 *
 **/
void
brox_optic_flow_spatial(const ofpix_t *I1, //first image
                        const ofpix_t *I2, //second image
                        ofpix_t *u, //x component of the optical flow
                        ofpix_t *v, //y component of the optical flow
                        const int nxx, //image width
                        const int nyy, //image height
                        const double alpha, //smoothness parameter
                        const double gamma, //gradient term parameter
                        const int nscales, //number of scales
                        const double nu, //downsampling factor
                        const double TOL, //stopping criterion threshold
                        const int inner_iter, //number of inner iterations
                        const int outer_iter, //number of outer iterations
                        const bool verbose //switch on messages
                        )
{
    int size = nxx * nyy;
    std::vector<ofpix_t *> I1s(nscales);
    std::vector<ofpix_t *> I2s(nscales);
    std::vector<ofpix_t *> us (nscales);
    std::vector<ofpix_t *> vs (nscales);
    std::vector<int> nx(nscales);
    std::vector<int> ny(nscales);

    I1s[0] = new ofpix_t[size];
    I2s[0] = new ofpix_t[size];

    //normalize the input images between 0 and 255
    image_normalization_2(I1, I2, I1s[0], I2s[0], size);

    //presmoothing the finest scale images
    gaussian(I1s[0], nxx, nyy, GAUSSIAN_SIGMA);
    gaussian(I2s[0], nxx, nyy, GAUSSIAN_SIGMA);

    us [0] = u;
    vs [0] = v;
    nx [0] = nxx;
    ny [0] = nyy;

    //create the scales
    for (int s = 1; s < nscales; s++) {
        zoom_size(nx[s - 1], ny[s - 1], &nx[s], &ny[s], nu);
        const int sizes = nx[s] * ny[s];

        I1s[s] = new ofpix_t[sizes];
        I2s[s] = new ofpix_t[sizes];
        us[s]  = new ofpix_t[sizes];
        vs[s]  = new ofpix_t[sizes];

        //compute the zoom from the previous scale
        zoom_out(I1s[s - 1], I1s[s], nx[s - 1], ny[s - 1], nu);
        zoom_out(I2s[s - 1], I2s[s], nx[s - 1], ny[s - 1], nu);
    }

    //initialization of the optical flow at the coarsest scale
    for (int i = 0; i < nx[nscales - 1] * ny[nscales - 1]; i++) {
        us[nscales - 1][i] = vs[nscales - 1][i] = 0.0;
    }

    int number_of_threads = 0;
#pragma omp parallel reduction(+:number_of_threads)
    number_of_threads += 1;

    //pyramidal approach for computing the optical flow
    for (int s = nscales - 1; s >= 0; s--) {
        if (verbose) {
            std::cout << "Scale: " << s << std::endl;
        }

        //compute the optical flow for the current scale
        brox_optic_flow(
            I1s[s], I2s[s], us[s], vs[s], nx[s], ny[s],
            alpha, gamma, TOL, inner_iter, outer_iter, number_of_threads, verbose
            );

        //if it is not the finer scale, then upsample the optical flow and adapt it conveniently
        if (s) {
            zoom_in(us[s], us[s - 1], nx[s], ny[s], nx[s - 1], ny[s - 1]);
            zoom_in(vs[s], vs[s - 1], nx[s], ny[s], nx[s - 1], ny[s - 1]);

            for (int i = 0; i < nx[s - 1] * ny[s - 1]; i++) {
                us[s - 1][i] *= 1.0 / nu;
                vs[s - 1][i] *= 1.0 / nu;
            }
        }
    }

    //delete allocated memory
    delete [] I1s[0];
    delete [] I2s[0];

    for (int i = 1; i < nscales; i++) {
        delete [] I1s[i];
        delete [] I2s[i];
        delete [] us [i];
        delete [] vs [i];
    }
} // brox_optic_flow_spatial

