// This program is free software: you can use, modify and/or redistribute it
// under the terms of the simplified BSD License. You should have received a
// copy of this license along this program. If not, see
// <http://www.opensource.org/licenses/bsd-license.html>.
//
// Copyright (C) 2012, Javier Sánchez Pérez <jsanchez@dis.ulpgc.es>
// All rights reserved.

#include "brox_temporal_mask.h"

#include "of.h"

/**
 *
 * Compute the coefficients of the divergence term
 *
 */
void
brox_temporal_psi_divergence(const ofpix_t *psi, //robust functional
                             ofpix_t *psi1, //coefficients of divergence
                             ofpix_t *psi2, //coefficients of divergence
                             ofpix_t *psi3, //coefficients of divergence
                             ofpix_t *psi4, //coefficients of divergence
                             ofpix_t *psi5, //coefficients of divergence
                             ofpix_t *psi6, //coefficients of divergence
                             const int nx, //image width
                             const int ny, //image height
                             const int nz //image depth
                             )
{
    const int df = nx * ny;

    //compute the spatial psi functions for all frames
    for (int f = 0; f < nz; f++) {
        #pragma omp parallel for
        for (int i = 1; i < ny - 1; i++) {
            for (int j = 1; j < nx - 1; j++) {
                const int k = f * df + i * nx + j;

                psi1[k] = 0.5 * (psi[k + nx] + psi[k]);
                psi2[k] = 0.5 * (psi[k - nx] + psi[k]);
                psi3[k] = 0.5 * (psi[k + 1]  + psi[k]);
                psi4[k] = 0.5 * (psi[k - 1]  + psi[k]);
            }
        }

        //calculate coefficients in the first and last rows
        #pragma omp parallel for
        for (int j = 1; j < nx - 1; j++) {
            int k = f * df + j;

            psi1[k] = 0.5 * (psi[k + nx] + psi[k]);
            psi2[k] = 0;
            psi3[k] = 0.5 * (psi[k + 1]  + psi[k]);
            psi4[k] = 0.5 * (psi[k - 1]  + psi[k]);

            k = f * df + (ny - 1) * nx + j;

            psi1[k] = 0;
            psi2[k] = 0.5 * (psi[k - nx] + psi[k]);
            psi3[k] = 0.5 * (psi[k + 1]  + psi[k]);
            psi4[k] = 0.5 * (psi[k - 1]  + psi[k]);
        }

        //calculate coefficients in the first and last columns
        #pragma omp parallel for
        for (int i = 1; i < ny - 1; i++) {
            int k = f * df + i * nx;

            psi1[k] = 0.5 * (psi[k + nx] + psi[k]);
            psi2[k] = 0.5 * (psi[k - nx] + psi[k]);
            psi3[k] = 0.5 * (psi[k + 1] + psi[k]);
            psi4[k] = 0;

            k = f * df + (i + 1) * nx - 1;

            psi1[k] = 0.5 * (psi[k + nx] + psi[k]);
            psi2[k] = 0.5 * (psi[k - nx] + psi[k]);
            psi3[k] = 0;
            psi4[k] = 0.5 * (psi[k - 1]  + psi[k]);
        }

        //up-left corner (0,0)
        int k = f * df;
        psi1[k] = 0.5 * (psi[k + nx] + psi[k]);
        psi3[k] = 0.5 * (psi[k + 1]  + psi[k]);
        psi2[k] = psi4[k] = 0;

        //up-right corner (nx,0)
        k = f * df + nx - 1;
        psi1[k] = 0.5 * (psi[k + nx] + psi[k]);
        psi4[k] = 0.5 * (psi[k - 1]  + psi[k]);
        psi2[k] = psi3[k] = 0;

        //bottom-left corner (0,ny)
        k = f * df + (ny - 1) * nx;
        psi2[k] = 0.5 * (psi[k - nx] + psi[k]);
        psi3[k] = 0.5 * (psi[k + 1]  + psi[k]);
        psi1[k] = psi4[k] = 0;

        //bottom-right corner (nx,ny)
        k = f * df + ny * nx - 1;
        psi2[k] = 0.5 * (psi[k - nx] + psi[k]);
        psi4[k] = 0.5 * (psi[k - 1]  + psi[k]);
        psi1[k] = psi3[k] = 0;
    }

    if (nz > 1) {
        //compute the coefficients in the interior frames
        for (int f = 1; f < nz - 1; f++) {
            for (int i = 0; i < df; i++) {
                const int k = f * df + i;
                psi5[k] = 0.5 * (psi[k - df] + psi[k]);
                psi6[k] = 0.5 * (psi[k + df] + psi[k]);
            }
        }

        //compute the coefficients in the first and last frames
        for (int i = 0; i < df; i++) {
            int k = i;
            psi5[k] = 0;
            psi6[k] = 0.5 * (psi[k + df] + psi[k]);

            k = (nz - 1) * df + i;
            psi5[k] = 0.5 * (psi[k - df] + psi[k]);
            psi6[k] = 0;
        }
    } else {
        for (int i = 0; i < df; i++) {
            psi5[i] = psi6[i] = 0;
        }
    }
} // brox_temporal_psi_divergence

/**
 *
 * Compute the divergence of the optical flow
 *
 */
void
brox_temporal_divergence_u(const ofpix_t *u, //x component of optical flow
                           const ofpix_t *v, //y component of optical flow
                           const ofpix_t *psi1, //coefficients of divergence
                           const ofpix_t *psi2, //coefficients of divergence
                           const ofpix_t *psi3, //coefficients of divergence
                           const ofpix_t *psi4, //coefficients of divergence
                           const ofpix_t *psi5, //coefficients of divergence
                           const ofpix_t *psi6, //coefficients of divergence
                           ofpix_t *div_u, //computed divergence for u
                           ofpix_t *div_v, //computed divergence for v
                           const int nx, //image width
                           const int ny, //image height
                           const int nz //image depth
                           )
{
    const int df = nx * ny;

    //calculate the divergence in the center body of the image
    for (int f = 0; f < nz; f++) {
        #pragma omp parallel for
        for (int i = 1; i < ny - 1; i++) {
            for (int j = 1; j < nx - 1; j++) {
                const int k = f * df + i * nx + j;

                div_u[k] = psi1[k] * (u[k + nx] - u[k]) + psi2[k] * (u[k - nx] - u[k]) +
                           psi3[k] * (u[k + 1]  - u[k]) + psi4[k] * (u[k - 1]  - u[k]);
                div_v[k] = psi1[k] * (v[k + nx] - v[k]) + psi2[k] * (v[k - nx] - v[k]) +
                           psi3[k] * (v[k + 1]  - v[k]) + psi4[k] * (v[k - 1]  - v[k]);
            }
        }

        //calculate the divergence in the first and last rows
        #pragma omp parallel for
        for (int j = 1; j < nx - 1; j++) {
            int k = f * df + j;
            div_u[k] = psi1[k] * (u[k + nx] - u[k]) + psi3[k] * (u[k + 1] - u[k]) + psi4[k] * (u[k - 1] - u[k]);
            div_v[k] = psi1[k] * (v[k + nx] - v[k]) + psi3[k] * (v[k + 1] - v[k]) + psi4[k] * (v[k - 1] - v[k]);

            k = f * df + (ny - 1) * nx + j;
            div_u[k] = psi2[k] * (u[k - nx] - u[k]) + psi3[k] * (u[k + 1] - u[k]) + psi4[k] * (u[k - 1] - u[k]);
            div_v[k] = psi2[k] * (v[k - nx] - v[k]) + psi3[k] * (v[k + 1] - v[k]) + psi4[k] * (v[k - 1] - v[k]);
        }

        //calculate the divergence in the first and last columns
        #pragma omp parallel for
        for (int i = 1; i < ny - 1; i++) {
            int k = f * df + i * nx;
            div_u[k] = psi1[k] * (u[k + nx] - u[k]) + psi2[k] * (u[k - nx] - u[k]) + psi3[k] * (u[k + 1] - u[k]);
            div_v[k] = psi1[k] * (v[k + nx] - v[k]) + psi2[k] * (v[k - nx] - v[k]) + psi3[k] * (v[k + 1] - v[k]);

            k = f * df + (i + 1) * nx - 1;
            div_u[k] = psi1[k] * (u[k + nx] - u[k]) + psi2[k] * (u[k - nx] - u[k]) + psi4[k] * (u[k - 1] - u[k]);
            div_v[k] = psi1[k] * (v[k + nx] - v[k]) + psi2[k] * (v[k - nx] - v[k]) + psi4[k] * (v[k - 1] - v[k]);
        }

        //up-left corner (0,0)
        int k = f * df;
        div_u[k] = psi1[k] * (u[k + nx] - u[k]) + psi3[k] * (u[k + 1] - u[k]);
        div_v[k] = psi1[k] * (v[k + nx] - v[k]) + psi3[k] * (v[k + 1] - v[k]);

        //up-right corner (nx,0)
        k = f * df + nx - 1;
        div_u[k] = psi1[k] * (u[k + nx] - u[k]) + psi4[k] * (u[k - 1] - u[k]);
        div_v[k] = psi1[k] * (v[k + nx] - v[k]) + psi4[k] * (v[k - 1] - v[k]);

        //bottom-left corner (0,ny)
        k = f * df + (ny - 1) * nx;
        div_u[k] = psi2[k] * (u[k - nx] - u[k]) + psi3[k] * (u[k + 1] - u[k]);
        div_v[k] = psi2[k] * (v[k - nx] - v[k]) + psi3[k] * (v[k + 1] - v[k]);

        //bottom-right corner (nx,ny)
        k = f * df + ny * nx - 1;
        div_u[k] = psi2[k] * (u[k - nx] - u[k]) + psi4[k] * (u[k - 1] - u[k]);
        div_v[k] = psi2[k] * (v[k - nx] - v[k]) + psi4[k] * (v[k - 1] - v[k]);
    }

    if (nz > 1) {
        //calculate the divergence in the interior frames
        for (int f = 1; f < nz - 1; f++) {
            for (int i = 0; i < df; i++) {
                const int k = f * df + i;
                div_u[k] += psi5[k] * (u[k - df] - u[k]) + psi6[k] * (u[k + df] - u[k]);
                div_v[k] += psi5[k] * (v[k - df] - v[k]) + psi6[k] * (v[k + df] - v[k]);
            }
        }

        //calculate the divergence in the first and last frames
        for (int i = 0; i < df; i++) {
            int k = i;
            div_u[k] += psi6[k] * (u[k + df] - u[k]);
            div_v[k] += psi6[k] * (v[k + df] - v[k]);

            k = (nz - 1) * df + i;
            div_u[k] += psi5[k] * (u[k - df] - u[k]);
            div_v[k] += psi5[k] * (v[k - df] - v[k]);
        }
    }
} // brox_temporal_divergence_u

