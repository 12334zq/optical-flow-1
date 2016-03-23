// This program is free software: you can use, modify and/or redistribute it
// under the terms of the simplified BSD License. You should have received a
// copy of this license along this program. If not, see
// <http://www.opensource.org/licenses/bsd-license.html>.
//
// Copyright (C) 2012, Coloma Ballester <coloma.ballester@upf.edu>
// Copyright (C) 2013-2014 J. F. Garamendi <jf.garamendi@upf.edu>

#include "tvl1occflow_solvers.h"

#include <cmath>
#include "of.h"
#include "operators.h"
#include "tvl1occflow_tv_rof_box.h"
#include "tvl1occflow_constants.h"

using namespace std;
/*
 * Function name: project
 * Author:
 * Parameters:
 *              eta1 : (in/out) vectorized matrix.
 *              eta2 : (in/out) vectorized matrix.
 *              nx   : (in) Scalar integer. image width
 *		ny   : (in) Scalar integer. image height
 * Output:
 *              VOID
 * Description: for each i component of eta1 and eta2 performs eta1/sqrt(eta1[i]^2 + eta2[i]^2)
 *              and   eta2/sqrt(eta1[i]^2 + eta2[i]^2)
 *
 */
inline void project(ofpix_t *eta1, ofpix_t *eta2, const int nx, const int ny);
inline void
project(ofpix_t *eta1,
        ofpix_t *eta2,
        const int nx,
        const int ny)
{
    const int size = nx * ny;
    double norm2, norm;

    for (int j = 0; j < size; j++) {
        norm2 = eta1[j] * eta1[j] + eta2[j] * eta2[j];
        if (norm2 < IS_ZERO) {
            eta1[j] = 0.0;
            eta2[j] = 0.0;
        } else {
            norm = sqrt(norm2);
            eta1[j] = eta1[j] / norm;
            eta2[j] = eta2[j] / norm;
        }
    }
}

void
Solver_wrt_v(ofpix_t *u1,
             ofpix_t *u2,
             ofpix_t *v1,
             ofpix_t *v2,
             ofpix_t *chi,
             const ofpix_t *I1wx,
             const ofpix_t *I1wy,
             const ofpix_t *I_1wx,
             const ofpix_t *I_1wy,
             const ofpix_t *rho1_c,
             const ofpix_t *rho3_c,
             ofpix_t *Vfwd_1,
             ofpix_t *Vfwd_2,
             ofpix_t *Vbck_1,
             ofpix_t *Vbck_2,
             const ofpix_t *grad1,
             const ofpix_t *grad3,
             const double alpha,
             const double theta,
             const double lambda,
             const int nx,
             const int ny)
{
    const int size = nx * ny;
    const double l_t = lambda * theta;
    const double _1pat = 1. + alpha * theta;
    const double at_d_1pat = alpha * theta / _1pat;
    const double lt_d_1pat = 2. * lambda * theta / _1pat;

#pragma omp parallel for
    for (int i = 0; i < size; i++) {
        double d1 = 0, d2 = 0;

        // We stored the data in V1 for using them when chi==0
        double rho1 = rho1_c[i] + (I1wx[i] * u1[i] + I1wy[i] * u2[i]);
        if (rho1 < -l_t * grad1[i]) {
            d1 = l_t * I1wx[i];
            d2 = l_t * I1wy[i];
        } else {
            if (rho1 > l_t * grad1[i]) {
                d1 = -l_t * I1wx[i];
                d2 = -l_t * I1wy[i];
            } else {
                if (grad1[i] < IS_ZERO) {
                    d1 = d2 = 0;
                } else {
                    d1 = -rho1 * I1wx[i] / grad1[i];
                    d2 = -rho1 * I1wy[i] / grad1[i];
                }
            }
        }
        Vfwd_1[i] = u1[i] + d1;
        Vfwd_2[i] = u2[i] + d2;

        //  we stored the data in V3 for using them when chi==1
        double rho3 = rho3_c[i] - (I_1wx[i] * u1[i] + I_1wy[i] * u2[i]);
        double A = rho3 + at_d_1pat * (I_1wx[i] * u1[i] + I_1wy[i] * u2[i]);

        if (A < -lt_d_1pat * grad3[i]) {
            d1 = -lt_d_1pat * I_1wx[i];
            d2 = -lt_d_1pat * I_1wy[i];

            Vbck_1[i] = (u1[i] / _1pat) + d1;
            Vbck_2[i] = (u2[i] / _1pat) + d2;
        } else {
            if (A > lt_d_1pat * grad3[i]) {
                d1 = lt_d_1pat * I_1wx[i];
                d2 = lt_d_1pat * I_1wy[i];

                Vbck_1[i] = (u1[i] / _1pat) + d1;
                Vbck_2[i] = (u2[i] / _1pat) + d2;
            } else {
                if (grad3[i] < IS_ZERO) {
                    d1 = d2 = 0;
                } else {
                    d1 = rho3 * I_1wx[i] / grad3[i];
                    d2 = rho3 * I_1wy[i] / grad3[i];
                }
                Vbck_1[i] = u1[i] + d1;
                Vbck_2[i] = u2[i] + d2;
            }
        }

        if (chi[i] < THR_CHI) {
            v1[i] = Vfwd_1[i];
            v2[i] = Vfwd_2[i];
        } else {
            v1[i] = Vbck_1[i];
            v2[i] = Vbck_2[i];
        }
    } //end for i.
} // Solver_wrt_v

void
Solver_wrt_u(ofpix_t *u1,
             ofpix_t *u2,
             const ofpix_t *v1,
             const ofpix_t *v2,
             const ofpix_t *chi,
             const ofpix_t *g,
             const double theta,
             const double beta,
             const int nx,
             const int ny)
{
    // f1,f2, chix, chiy and p variables are static because Solver_wrt_u will be called inside a loop
    // and we don't want to store memory at each iteration (malloc is a very computtional expensive function)
#warning "non-reentrant code"
    static ofpix_t *p11 = 0, *p12 = 0, *p21 = 0, *p22 = 0;
    static int old_nx = 0;  //, old_ny=0;

    if (old_nx != nx) {
        old_nx = nx;
        delete [] p11;
        delete [] p12;
        delete [] p21;
        delete [] p22;

        p11 = new ofpix_t[nx * ny];
        p12 = new ofpix_t[nx * ny];
        p21 = new ofpix_t[nx * ny];
        p22 = new ofpix_t[nx * ny];


        for (int i = 0; i < nx * ny; i++) {
            p11[i] = 0.0;
            p12[i] = 0.0;
            p21[i] = 0.0;
            p22[i] = 0.0;
        }
    }

    ofpix_t *chix = new ofpix_t[nx * ny];
    ofpix_t *chiy = new ofpix_t[nx * ny];
    forward_gradient(chi, chix, chiy, nx, ny);

    ofpix_t *f1 = new ofpix_t[nx * ny];
    ofpix_t *f2 = new ofpix_t[nx * ny];

    for (int i = 0; i < nx * ny; i++) {
        f1[i] = v1[i] / theta + beta * chix[i];
        f2[i] = v2[i] / theta + beta * chiy[i];

        //For initial approximation we consider the dual variable as zero,
        //so div(p) term is omitted
        u1[i] = v1[i] + theta * beta * chix[i];
        u2[i] = v2[i] + theta * beta * chiy[i];
    }
    delete [] chix;
    delete [] chiy;

    //For each component of the optical flow, we perform the minimize the  Rudin-Osher-Fatemi (modified) problem
    Scalar_ROF_BoxCellCentered(u1, f1, p11, p12, g, theta, OMEGA, nx, ny,
                               MAX_ITERATIONS_U);
    Scalar_ROF_BoxCellCentered(u2, f2, p21, p22, g, theta, OMEGA, nx, ny,
                               MAX_ITERATIONS_U);

    delete [] f1;
    delete [] f2;
} // Solver_wrt_u

void
Solver_wrt_chi(const ofpix_t *u1,
               const ofpix_t *u2,
               ofpix_t *chi,
               const ofpix_t *I1wx,
               const ofpix_t *I1wy,
               const ofpix_t *I_1wx,
               const ofpix_t *I_1wy,
               const ofpix_t *rho1_c,
               const ofpix_t *rho3_c,
               const ofpix_t *Vfwd_1,
               const ofpix_t *Vfwd_2,
               const ofpix_t *Vbck_1,
               const ofpix_t *Vbck_2,
               const ofpix_t *g,
               const double lambda,
               const double theta,
               const double alpha,
               const double beta,
               const double tau_chi,
               const double tau_eta,
               const int nx,
               const int ny)
{
#warning "non-reentrant code"
    static int old_nx = 0;      //, old_ny=0;
    static ofpix_t *eta1 = 0, *eta2;

    if (old_nx != nx) {
        old_nx = nx;

        delete [] eta1;
        delete [] eta2;

        eta1 = new ofpix_t[nx * ny];
        eta2 = new ofpix_t[nx * ny];
    }

    //Calculus of eta
    for (int n_chi = 0; n_chi < MAX_ITERATIONS_CHI; n_chi++) {
        //calculate the gradient of chi
        ofpix_t *chix = new ofpix_t[nx * ny];
        ofpix_t *chiy = new ofpix_t[nx * ny];
        forward_gradient(chi, chix, chiy, nx, ny);
        for (int i = 0; i < nx * ny; i++) {
#warning "eta1 and eta2 are used uninitialized"
            eta1[i] = eta1[i] + tau_eta * g[i] * chix[i];
            eta2[i] = eta2[i] + tau_eta * g[i] * chiy[i];
        }
        delete [] chix;
        delete [] chiy;

        project(eta1, eta2, nx, ny);
        //end of caluculus of eta

        //Caluculus of chi
        //apply the g function to eta
        ofpix_t *geta1 = new ofpix_t[nx * ny];
        ofpix_t *geta2 = new ofpix_t[nx * ny];
        for (int i = 0; i < nx * ny; i++) {
            geta1[i] = g[i] * eta1[i];
            geta2[i] = g[i] * eta2[i];
        }

        //Divergence of geta=(eta1, eta2)
        ofpix_t *div_eta = new ofpix_t[nx * ny];
        divergence(geta1, geta2, div_eta, nx, ny);
        delete [] geta1;
        delete [] geta2;

        //Divergence of u
        ofpix_t *div_u = new ofpix_t[nx * ny];
        divergence(u1, u2, div_u, nx, ny);

        for (int i = 0; i < nx * ny; i++) {
            /*
               //using one variable for v
               double rho1 = rho1_c[i] + (I1wx[i]*v1[i]+I1wy[i]*v2[i]);
               double abs_rho1 = (rho1< 0) ? -rho1 : rho1;

               double rho3 = rho3_c[i] - (I_1wx[i]*v1[i]+I_1wy[i]*v2[i]);
               double abs_rho3 = (rho3< 0) ? -rho3 : rho3;


               double F = lambda*(abs_rho3-abs_rho1);
               double G = alpha*0.5*(v1[i]*v1[i]+v2[i]*v2[i]);
             */

            //Using two variables for v
            double rho1 = rho1_c[i] + (I1wx[i] * Vfwd_1[i] + I1wy[i] * Vfwd_2[i]);
            double abs_rho1 = (rho1 < 0.) ? -rho1 : rho1;
            double rho3 = rho3_c[i] - (I_1wx[i] * Vbck_1[i] + I_1wy[i] * Vbck_2[i]);
            double abs_rho3 = (rho3 < 0.) ? -rho3 : rho3;
            double F = 0;
            double G = 0;
            if (chi[i] < 0.5) {
                F = -lambda * abs_rho1;
                G = -(0.5 / theta)
                    * ( (Vfwd_1[i] - u1[i]) * (Vfwd_1[i] - u1[i])
                        + (Vfwd_2[i] - u2[i]) * (Vfwd_2[i] - u2[i]) );
            } else {
                F = lambda * abs_rho3;
                G = (0.5 / theta)
                    * ( (Vbck_1[i] - u1[i]) * (Vbck_1[i] - u1[i])
                        + (Vbck_2[i] - u2[i]) * (Vbck_2[i] - u2[i]) )
                    + alpha * theta * (Vbck_1[i] * Vbck_1[i] + Vbck_2[i] * Vbck_2[i]);
            }

            chi[i] = chi[i] + tau_chi * (div_eta[i] - F - G - beta * div_u[i]);

            //Proyecci into [0,1]
            if (chi[i] > 1.) {
                chi[i] = 1.;
            } else if (chi[i] < 0.) {
                chi[i] = 0.;
            }
        } //end for i
        delete [] div_eta;
        delete [] div_u;
    } //end for n_chi
} // Solver_wrt_chi

