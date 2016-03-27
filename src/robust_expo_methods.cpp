// This program is free software: you can use, modify and/or redistribute it
// under the terms of the simplified BSD License. You should have received a
// copy of this license along this program. If not, see
// <http://www.opensource.org/licenses/bsd-license.html>.
//
// Copyright (C) 2012, Javier Sánchez Pérez <jsanchez@dis.ulpgc.es>
// Copyright (C) 2014, Agustín Salgado de la Nuez <asalgado@dis.ulpgc.es>
// Copyright (C) 2014, 2015 Nelson Monzón López <nmonzon@ctim.es>
// All rights reserved.

#include "robust_expo_methods.h"

#include <cmath>
#include <vector>
#include <iostream>

#include "operators.h"
#include "zoom.h"
#include "utils.h"
#include "bicubic_interpolation.h"
#include "robust_expo_generic_tensor.h"
#include "robust_expo_smoothness.h"

#define MAXITER 300
#define SOR_PARAMETER 1.9
#define GAUSSIAN_SIGMA 0.8

using namespace std;


/**
 *
 * Compute the coefficients of the robust functional (data term)
 *
 **/
void
psi_data (const ofpix_t *I1, //first image
          const ofpix_t *I2, //second image
          const ofpix_t *I2x, //gradient of the second image
          const ofpix_t *I2y, //gradient of the second image
          const ofpix_t *du, //motion increment
          const ofpix_t *dv, //motion increment
          ofpix_t *psip, //output coefficients
          const int nx, //image width
          const int ny, //image height
          const int nz  //image channels
          )
{
    const int size = nx * ny;

    //compute 1/(sqrt(Sum(I2-I1+I2x*du+I2y*dv)²+e²) in each pixel
    #pragma omp parallel for
    for (int i = 0; i < size; i++) {
        double dI2 = 0;
        for (int k = 0; k < nz; k++) {
            int real_index = i * nz + k;
            const double dI  = I2[real_index] + I2x[real_index] * du[i] + I2y[real_index] * dv[i] - I1[real_index];

            dI2 += dI * dI;
        }

        psip[i] = ( 1. / sqrt (dI2 + ROBUST_EXPO_EPSILON * ROBUST_EXPO_EPSILON) );
    }
}

/**
 *
 * Compute the coefficients of the robust functional (gradient term)
 *
 **/

void
psi_gradient (const ofpix_t *I1x, //gradient of the first image
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
              const int ny, //image height
              const int nz // nº image channels
              )
{
    const int size = nx * ny;

    //compute 1/(sqrt(|DI2-DI1=HI2*(du,dv)|²+e²) in each pixel
    #pragma omp parallel for
    for (int i = 0; i < size; i++) {
        double dI2 = 0;

        for (int k = 0; k < nz; k++) {
            int real_index = i * nz + k;
            const double dIx = I2x[real_index] + I2xx[real_index] * du[i] + I2xy[real_index] * dv[i] - I1x[real_index];
            const double dIy = I2y[real_index] + I2xy[real_index] * du[i] + I2yy[real_index] * dv[i] - I1y[real_index];

            dI2 += dIx * dIx + dIy * dIy;
        }

        psip[i] = ( 1. / sqrt (dI2 + ROBUST_EXPO_EPSILON * ROBUST_EXPO_EPSILON) );
    }
}

/**
 *
 *  SOR iteration in one position
 *
 */
inline double
sor_iteration(const ofpix_t *Au, //constant part of the numerator of u
              const ofpix_t *Av, //constant part of the numerator of vPriscila Estévez
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
              const int j0, //previous column
              const int j1, //following column
              const int nx //number of columns
              )
{
    //set the SOR extrapolation parameter
    const double w = SOR_PARAMETER;

    //calculate the position in the array
    const int k = i * nx + j;

    //compute the divergence part of the numerator
    const double div_du = psi1[k] * du[k + j1] + psi2[k] * du[k - j0] +
                          psi3[k] * du[k + i1] + psi4[k] * du[k - i0];
    const double div_dv = psi1[k] * dv[k + j1] + psi2[k] * dv[k - j0] +
                          psi3[k] * dv[k + i1] + psi4[k] * dv[k - i0];
    const double duk = du[k];
    const double dvk = dv[k];

    //update the motion increment
    du[k] = (1. - w) * du[k] + w * (Au[k] - D[k] * dv[k] + alpha * div_du) / Du[k];
    dv[k] = (1. - w) * dv[k] + w * (Av[k] - D[k] * du[k] + alpha * div_dv) / Dv[k];

    //return the covergence error in this position
    return (du[k] - duk) * (du[k] - duk) + (dv[k] - dvk) * (dv[k] - dvk);
}

/**
 *
 * Compute the optic flow with the matrix
 *
 **/
void
robust_expo_methods(const ofpix_t *I1, //first image
                    const ofpix_t *I2, //second image
                    ofpix_t *u, //x component of the optical flow
                    ofpix_t *v, //y component of the optical flow
                    const int nx, //image width
                    const int ny, //image height
                    const int nz, // number of color channels in the image
                    const int method_type, // choose the diffusion strategy
                    const double alpha, // smoothness parameter
                    const double gamma, // gradient term parameter
                    const double lambda, // coefficient parameter for the decreasing function (if needed)
                    const double TOL, // stopping criterion threshold
                    const int inner_iter, // number of inner iterations
                    const int outer_iter, // number of outer iterations
                    const int number_of_threads, // number of threads for the parallel code
                    const bool verbose // switch on messages
                    )
{
    const int size_flow = nx * ny;
    const int size      = size_flow * nz;

    //allocate memory
    ofpix_t *du    = new ofpix_t[size_flow];
    ofpix_t *dv    = new ofpix_t[size_flow];
    ofpix_t *ux    = new ofpix_t[size_flow];
    ofpix_t *uy    = new ofpix_t[size_flow];
    ofpix_t *vx    = new ofpix_t[size_flow];
    ofpix_t *vy    = new ofpix_t[size_flow];
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
    ofpix_t *div_u = new ofpix_t[size_flow];
    ofpix_t *div_v = new ofpix_t[size_flow];
    ofpix_t *div_d = new ofpix_t[size_flow];
    ofpix_t *Au    = new ofpix_t[size_flow];
    ofpix_t *Av    = new ofpix_t[size_flow];
    ofpix_t *Du    = new ofpix_t[size_flow];
    ofpix_t *Dv    = new ofpix_t[size_flow];
    ofpix_t *D     = new ofpix_t[size_flow];
    ofpix_t *expo  = new ofpix_t[size_flow];
    ofpix_t *psid  = new ofpix_t[size_flow];
    ofpix_t *psig  = new ofpix_t[size_flow];
    ofpix_t *psis  = new ofpix_t[size_flow];
    ofpix_t *psi1  = new ofpix_t[size_flow];
    ofpix_t *psi2  = new ofpix_t[size_flow];
    ofpix_t *psi3  = new ofpix_t[size_flow];
    ofpix_t *psi4  = new ofpix_t[size_flow];

    //compute the gradient of the images
    centered_gradient(I1, I1x, I1y, nx, ny, nz);
    centered_gradient(I2, I2x, I2y, nx, ny, nz);

    //compute second order derivatives
    Dxx(I2, I2xx, nx, ny, nz);
    Dyy(I2, I2yy, nx, ny, nz);
    Dxy(I2, I2xy, nx, ny, nz);


    //compute the smoothness_coefficients including the robust function Phi for the smoothness term
    robust_expo_exponential_calculation (I1x, I1y, size_flow, size, nz, alpha, lambda, method_type, expo);

    //outer iterations loop
    for (int no = 0; no < outer_iter; no++) {
        // Warp the second image and its derivatives
        bicubic_interpolation_warp_color (I2, u, v, I2w, nx, ny, nz, true);
        bicubic_interpolation_warp_color (I2x, u, v, I2wx, nx, ny, nz, true);
        bicubic_interpolation_warp_color (I2y, u, v, I2wy, nx, ny, nz, true);
        bicubic_interpolation_warp_color (I2xx, u, v, I2wxx, nx, ny, nz, true);
        bicubic_interpolation_warp_color (I2xy, u, v, I2wxy, nx, ny, nz, true);
        bicubic_interpolation_warp_color (I2yy, u, v, I2wyy, nx, ny, nz, true);

        // Compute the flow gradient
        centered_gradient (u, ux, uy, nx, ny, 1);
        centered_gradient (v, vx, vy, nx, ny, 1);

        //compute robust function Phi for the smoothness term
        robust_expo_psi_smooth(ux, uy, vx, vy, expo, size_flow, psis);

        //compute coefficients of Phi functions in divergence
        robust_expo_psi_divergence (psi1, psi2, psi3, psi4, psis, nx, ny);

        //Calculate the divergence
        robust_expo_divergence (u, psi1, psi2, psi3, psi4, nx, ny, div_u);
        robust_expo_divergence (v, psi1, psi2, psi3, psi4, nx, ny, div_v);

        #pragma omp parallel for
        for (int i = 0; i < size_flow; i++) {
            //compute the coefficents of dw[i] in the smoothness term using gradient exponent
            div_d[i] = alpha * (psi1[i] + psi2[i] + psi3[i] + psi4[i]);

            //initialize the motion increment
            du[i] = dv[i] = 0;
        }

        //inner iterations loop
        for (int ni = 0; ni < inner_iter; ni++) {
            //compute robust function Psi for the data and gradient terms
            psi_data (I1, I2w, I2wx, I2wy, du, dv, psid, nx, ny, nz);
            psi_gradient (I1x, I1y, I2wx, I2wy, I2wxx, I2wxy, I2wyy, du, dv, psig, nx, ny, nz);

            //store constant parts of the numerical scheme (equation (11))
            double BNu, dif, BNv, BDu, BDv, dx, dy, GNu, GNv, GDu, GDv, DI_Gradient, DI_Data;
            int index_flow = 0;

            for (int index_image = 0; index_image < size; index_image += nz) {
                BNu = BNv = BDu = BDv = dx = dy = GNu = GNv = GDu = GDv = DI_Gradient = DI_Data = 0;

                for (int index_multichannel = 0; index_multichannel < nz; index_multichannel++) {
                    const int real_index = index_image + index_multichannel;

                    //brightness constancy term
                    dif  = I2w[real_index] - I1[real_index];
                    BNu += dif * I2wx[real_index];
                    BNv += dif * I2wy[real_index];
                    BDu += I2wx[real_index] * I2wx[real_index];
                    BDv += I2wy[real_index] * I2wy[real_index];
                    DI_Data += (I2wy [real_index] * I2wx [real_index]);

                    //gradient constancy term
                    dx   = (I2wx[real_index] - I1x[real_index]);
                    dy   = (I2wy[real_index] - I1y[real_index]);
                    GNu += (dx * I2wxx[real_index] + dy * I2wxy[real_index]);
                    GNv += (dx * I2wxy[real_index] + dy * I2wyy[real_index]);
                    GDu += (I2wxx[real_index] * I2wxx[real_index] + I2wxy[real_index] * I2wxy[real_index]);
                    GDv += (I2wyy[real_index] * I2wyy[real_index] + I2wxy[real_index] * I2wxy[real_index]);
                    DI_Gradient  += (I2wxx[real_index] + I2wyy[real_index]) * I2wxy[real_index];
                }

                const double g = gamma * psig[index_flow];

                BNu  = -psid[index_flow] * BNu;
                BNv  = -psid[index_flow] * BNv;
                BDu  =  psid[index_flow] * BDu;
                BDv  =  psid[index_flow] * BDv;
                GNu  = -g * GNu;
                GNv  = -g * GNv;
                GDu  =  g * GDu;
                GDv  =  g * GDv;

                Au[index_flow] = BNu + GNu + alpha * div_u[index_flow];
                Av[index_flow] = BNv + GNv + alpha * div_v[index_flow];
                Du[index_flow] = BDu + GDu + div_d[index_flow];
                Dv[index_flow] = BDv + GDv + div_d[index_flow];
                D [index_flow] = psid[index_flow] * DI_Data + g * DI_Gradient;

                index_flow++;
            } // end image loop


            //sor iterations loop
            double error = 1000;
            int nsor = 0;

            while (error > TOL && nsor < MAXITER) {
                error = 0;
                nsor++;

                //#pragma omp parallel for reduction(+:error)
                #pragma omp parallel for reduction(+:error) num_threads((number_of_threads < (ny-3)) ? number_of_threads : (ny-3))
                //update the motion increment in the center of the images
                for (int i = 1; i < ny - 1; i++) {
                    for (int j = 1; j < nx - 1; j++) {
                        error += sor_iteration(
                            Au, Av, Du, Dv, D, du, dv, alpha,
                            psi1, psi2, psi3, psi4,
                            i, nx, nx, j, 1, 1, nx
                            );
                    }
                }

                //update the motion increment in the first and last rows
                for (int j = 1; j < nx - 1; j++) {
                    error += sor_iteration(
                        Au, Av, Du, Dv, D, du, dv, alpha,
                        psi1, psi2, psi3, psi4,
                        0, 0, nx, j, 1, 1, nx
                        );

                    error += sor_iteration(
                        Au, Av, Du, Dv, D, du, dv, alpha,
                        psi1, psi2, psi3, psi4,
                        ny - 1, nx, 0, j, 1, 1, nx
                        );
                }

                //update the motion increment in the first and last columns
                for (int i = 1; i < ny - 1; i++) {
                    error += sor_iteration(
                        Au, Av, Du, Dv, D, du, dv, alpha,
                        psi1, psi2, psi3, psi4,
                        i, nx, nx, 0, 0, 1, nx
                        );

                    error += sor_iteration(
                        Au, Av, Du, Dv, D, du, dv, alpha,
                        psi1, psi2, psi3, psi4,
                        i, nx, nx, nx - 1, 1, 0, nx
                        );
                }

                //process the top-left corner (0,0)
                error += sor_iteration(
                    Au, Av, Du, Dv, D, du, dv, alpha,
                    psi1, psi2, psi3, psi4,
                    0, 0, nx, 0, 0, 1, nx
                    );

                //process the top-right corner (0,nx-1)
                error += sor_iteration(
                    Au, Av, Du, Dv, D, du, dv, alpha,
                    psi1, psi2, psi3, psi4,
                    0, 0, nx, nx - 1, 1, 0, nx
                    );

                //process the bottom-left corner (ny-1,0)
                error += sor_iteration(
                    Au, Av, Du, Dv, D, du, dv, alpha,
                    psi1, psi2, psi3, psi4,
                    ny - 1, nx, 0, 0, 0, 1, nx
                    );

                //process the bottom-right corner (ny-1,nx-1)
                error += sor_iteration(
                    Au, Av, Du, Dv, D, du, dv, alpha,
                    psi1, psi2, psi3, psi4,
                    ny - 1, nx, 0, nx - 1, 1, 0, nx
                    );

                error = sqrt(error / size);
            }
            if (verbose) {
                std::cout << "Iterations: " << nsor << " Error: " << error << std::endl;
            }
        }

        //update the flow with the estimated motion increment
        for (int i = 0; i < size_flow; i++) {
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

    delete []expo;

    delete []psid;
    delete []psig;
    delete []psis;
    delete []psi1;
    delete []psi2;
    delete []psi3;
    delete []psi4;
} // robust_expo_methods

/**
 *
 *  Multiscale approach for computing the optical flow
 *
 **/
void
robust_expo_methods(const ofpix_t *I1, // first image
                    const ofpix_t *I2, // second image
                    ofpix_t *u, // x component of the optical flow
                    ofpix_t *v, // y component of the optical flow
                    const int nxx, // image width
                    const int nyy, // image height
                    const int nzz, // number of color channels in image
                    const int method_type, // choose the diffusion strategy
                    const double alpha, // smoothness parameter
                    const double gamma, // gradient term parameter
                    const double lambda, // coefficient parameter for the decreasing function (if needed)
                    const int nscales, // number of scales
                    const double nu, // downsampling factor
                    const double TOL, // stopping criterion threshold
                    const int inner_iter, // number of inner iterations
                    const int outer_iter, // number of outer iterations
                    const bool verbose // switch on messages
                    )
{
    int size = nxx * nyy * nzz;
    std::vector<ofpix_t *> I1s(nscales);
    std::vector<ofpix_t *> I2s(nscales);
    std::vector<ofpix_t *> us (nscales);
    std::vector<ofpix_t *> vs (nscales);
    std::vector<int> nx(nscales);
    std::vector<int> ny(nscales);

    I1s[0] = new ofpix_t[size];
    I2s[0] = new ofpix_t[size];

    //normalize the input images between 0 and 255
    image_normalization_2_color(I1, I2, I1s[0], I2s[0], size, nzz);

    //presmoothing the finest scale images
    gaussian(I1s[0], nxx, nyy, nzz, GAUSSIAN_SIGMA);
    gaussian(I2s[0], nxx, nyy, nzz, GAUSSIAN_SIGMA);

    us [0] = u;
    vs [0] = v;
    nx [0] = nxx;
    ny [0] = nyy;

    //create the scales
    for (int s = 1; s < nscales; s++) {
        zoom_size(nx[s - 1], ny[s - 1], &nx[s], &ny[s], nu);
        const int sizes_flow = nx[s] * ny[s];
        const int sizes = sizes_flow * nzz;

        I1s[s] = new ofpix_t[sizes];
        I2s[s] = new ofpix_t[sizes];
        us[s]  = new ofpix_t[sizes_flow];
        vs[s]  = new ofpix_t[sizes_flow];

        //compute the zoom from the previous scale
        zoom_out_color(I1s[s - 1], I1s[s], nx[s - 1], ny[s - 1], nzz, nu);
        zoom_out_color(I2s[s - 1], I2s[s], nx[s - 1], ny[s - 1], nzz, nu);
    }

    //initialization of the optical flow at the coarsest scale
    for (int i = 0; i < nx[nscales - 1] * ny[nscales - 1]; i++) {
        us[nscales - 1][i] = vs[nscales - 1][i] = 0.0;
    }

    // readapt alpha for multichannel information
    int alpha_adapted_for_nchannels = alpha * nzz;
    int number_of_threads = 0;
    #pragma omp parallel reduction(+:number_of_threads)
    number_of_threads += 1;

    //pyramidal approach for computing the optical flow
    for (int s = nscales - 1; s >= 0; s--) {
        if (verbose) {
            std::cout << "Scale: " << s << std::endl;
        }

        robust_expo_methods(
            I1s[s], I2s[s], us[s], vs[s], nx[s], ny[s], nzz,
            method_type, alpha_adapted_for_nchannels, gamma, lambda,
            TOL, inner_iter, outer_iter, number_of_threads, verbose
            );

        //if it is not the finer scale, then upsample the optical flow and adapt it conveniently
        if (s) {
            zoom_in (us[s], us[s - 1], nx[s], ny[s], nx[s - 1], ny[s - 1]);
            zoom_in (vs[s], vs[s - 1], nx[s], ny[s], nx[s - 1], ny[s - 1]);

            for (int i = 0; i < nx[s - 1] * ny[s - 1]; i++) {
                us[s - 1][i] *= 1.0 / nu;
                vs[s - 1][i] *= 1.0 / nu;
            }
        }
    }

    //delete allocated memory
    delete []I1s[0];
    delete []I2s[0];

    for (int i = 1; i < nscales; i++) {
        delete []I1s[i];
        delete []I2s[i];
        delete []us [i];
        delete []vs [i];
    }
} // robust_expo_methods

