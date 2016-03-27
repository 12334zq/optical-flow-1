// This program is free software: you can use, modify and/or redistribute it
// under the terms of the simplified BSD License. You should have received a
// copy of this license along this program. If not, see
// <http://www.opensource.org/licenses/bsd-license.html>.
//
// Copyright (C) 2012, Javier Sánchez Pérez <jsanchez@dis.ulpgc.es>
// All rights reserved.

#include "brox_spatial_mask.h"

/**
 *
 * Compute the coefficients of the divergence term
 *
 */
void
brox_spatial_psi_divergence(const ofpix_t *psi, //robust functional
                            ofpix_t *psi1, //coefficients of divergence
                            ofpix_t *psi2, //coefficients of divergence
                            ofpix_t *psi3, //coefficients of divergence
                            ofpix_t *psi4, //coefficients of divergence
                            const int nx, //image width
                            const int ny //image height
                            )
{
    //calculate coefficients in the center body of the image
    #pragma omp parallel for
    for (int i = 1; i < ny - 1; i++) {
        for (int j = 1; j < nx - 1; j++) {
            const int k = i * nx + j;

            psi1[k] = 0.5 * (psi[k + nx] + psi[k]);
            psi2[k] = 0.5 * (psi[k - nx] + psi[k]);
            psi3[k] = 0.5 * (psi[k + 1] + psi[k]);
            psi4[k] = 0.5 * (psi[k - 1] + psi[k]);
        }
    }

    //calculate coefficients in the first and last rows
    #pragma omp parallel for
    for (int j = 1; j < nx - 1; j++) {
        psi1[j] = 0.5 * (psi[j + nx] + psi[j]);
        psi2[j] = 0;
        psi3[j] = 0.5 * (psi[j + 1] + psi[j]);
        psi4[j] = 0.5 * (psi[j - 1] + psi[j]);

        const int k  = (ny - 1) * nx + j;

        psi1[k] = 0;
        psi2[k] = 0.5 * (psi[k - nx] + psi[k]);
        psi3[k] = 0.5 * (psi[k + 1] + psi[k]);
        psi4[k] = 0.5 * (psi[k - 1] + psi[k]);
    }

    //calculate coefficients in the first and last columns
    #pragma omp parallel for
    for (int i = 1; i < ny - 1; i++) {
        const int k  = i * nx;

        psi1[k] = 0.5 * (psi[k + nx] + psi[k]);
        psi2[k] = 0.5 * (psi[k - nx] + psi[k]);
        psi3[k] = 0.5 * (psi[k + 1] + psi[k]);
        psi4[k] = 0;

        const int j  = (i + 1) * nx - 1;

        psi1[j] = 0.5 * (psi[j + nx] + psi[j]);
        psi2[j] = 0.5 * (psi[j - nx] + psi[j]);
        psi3[j] = 0;
        psi4[j] = 0.5 * (psi[j - 1] + psi[j]);
    }

    //up-left corner (0,0)
    psi1[0] = 0.5 * (psi[nx] + psi[0]);
    psi3[0] = 0.5 * (psi[1] + psi[0]);
    psi2[0] = psi4[0] = 0;

    //up-right corner (nx,0)
    psi1[nx - 1] = 0.5 * (psi[nx - 1 + nx] + psi[nx - 1]);
    psi4[nx - 1] = 0.5 * (psi[nx - 2] + psi[nx - 1]);
    psi2[nx - 1] = psi3[nx - 1] = 0;

    //bottom-left corner (0,ny)
    psi2[(ny - 1) * nx] = 0.5 * (psi[(ny - 2) * nx] + psi[(ny - 1) * nx]);
    psi3[(ny - 1) * nx] = 0.5 * (psi[(ny - 1) * nx + 1] + psi[(ny - 1) * nx]);
    psi1[(ny - 1) * nx] = psi4[(ny - 1) * nx] = 0;


    //bottom-right corner (nx,ny)
    psi2[ny * nx - 1] = 0.5 * (psi[ny * nx - 1 - nx] + psi[ny * nx - 1]);
    psi4[ny * nx - 1] = 0.5 * (psi[ny * nx - 2] + psi[ny * nx - 1]);
    psi1[ny * nx - 1] = psi3[ny * nx - 1] = 0;
} // brox_spatial_psi_divergence

/**
 *
 * Compute the divergence of the optical flow
 *
 */
void
brox_spatial_divergence_u(const ofpix_t *u, //x component of optical flow
                          const ofpix_t *v, //y component of optical flow
                          const ofpix_t *psi1, //coefficients of divergence
                          const ofpix_t *psi2, //coefficients of divergence
                          const ofpix_t *psi3, //coefficients of divergence
                          const ofpix_t *psi4, //coefficients of divergence
                          ofpix_t *div_u, //computed divergence for u
                          ofpix_t *div_v, //computed divergence for v
                          const int nx, //image width
                          const int ny //image height
                          )
{
    //calculate the divergence in the center body of the image
    #pragma omp parallel for
    for (int i = 1; i < ny - 1; i++) {
        for (int j = 1; j < nx - 1; j++) {
            const int k = i * nx + j;

            div_u[k] = psi1[k] * (u[k + nx] - u[k]) + psi2[k] * (u[k - nx] - u[k]) +
                       psi3[k] * (u[k + 1]  - u[k]) + psi4[k] * (u[k - 1]  - u[k]);
            div_v[k] = psi1[k] * (v[k + nx] - v[k]) + psi2[k] * (v[k - nx] - v[k]) +
                       psi3[k] * (v[k + 1]  - v[k]) + psi4[k] * (v[k - 1]  - v[k]);
        }
    }

    //calculate the divergence in the first and last rows
    #pragma omp parallel for
    for (int j = 1; j < nx - 1; j++) {
        div_u[j] = psi1[j] * (u[j + nx] - u[j]) + psi3[j] * (u[j + 1]  - u[j]) + psi4[j] * (u[j - 1] - u[j]);
        div_v[j] = psi1[j] * (v[j + nx] - v[j]) + psi3[j] * (v[j + 1]  - v[j]) + psi4[j] * (v[j - 1] - v[j]);

        const int k  = (ny - 1) * nx + j;

        div_u[k] = psi2[k] * (u[k - nx] - u[k]) + psi3[k] * (u[k + 1] - u[k]) + psi4[k] * (u[k - 1] - u[k]);
        div_v[k] = psi2[k] * (v[k - nx] - v[k]) + psi3[k] * (v[k + 1] - v[k]) + psi4[k] * (v[k - 1] - v[k]);
    }

    //calculate the divergence in the first and last columns
    #pragma omp parallel for
    for (int i = 1; i < ny - 1; i++) {
        const int k  = i * nx;

        div_u[k] = psi1[k] * (u[k + nx] - u[k]) + psi2[k] * (u[k - nx] - u[k]) + psi3[k] * (u[k + 1] - u[k]);
        div_v[k] = psi1[k] * (v[k + nx] - v[k]) + psi2[k] * (v[k - nx] - v[k]) + psi3[k] * (v[k + 1] - v[k]);

        const int j  = (i + 1) * nx - 1;

        div_u[j] = psi1[j] * (u[j + nx] - u[j]) + psi2[j] * (u[j - nx] - u[j]) + psi4[j] * (u[j - 1] - u[j]);
        div_v[j] = psi1[j] * (v[j + nx] - v[j]) + psi2[j] * (v[j - nx] - v[j]) + psi4[j] * (v[j - 1] - v[j]);
    }

    //up-left corner (0,0)
    div_u[0] = psi1[0] * (u[nx] - u[0]) + psi3[0] * (u[1] - u[0]);
    div_v[0] = psi1[0] * (v[nx] - v[0]) + psi3[0] * (v[1] - v[0]);

    //up-right corner (nx,0)
    div_u[nx - 1] = psi1[nx - 1] * (u[nx - 1 + nx] - u[nx - 1]) + psi4[nx - 1] * (u[nx - 2]  - u[nx - 1]);
    div_v[nx - 1] = psi1[nx - 1] * (v[nx - 1 + nx] - v[nx - 1]) + psi4[nx - 1] * (v[nx - 2]  - v[nx - 1]);

    //bottom-left corner (0,ny)
    div_u[(ny - 1) * nx] = psi2[(ny - 1) * nx] * (u[(ny - 2) * nx]     - u[(ny - 1) * nx]) +
                           psi3[(ny - 1) * nx] * (u[(ny - 1) * nx + 1] - u[(ny - 1) * nx]);
    div_v[(ny - 1) * nx] = psi2[(ny - 1) * nx] * (v[(ny - 2) * nx]     - v[(ny - 1) * nx]) +
                           psi3[(ny - 1) * nx] * (v[(ny - 1) * nx + 1] - v[(ny - 1) * nx]);

    //bottom-right corner (nx,ny)
    div_u[ny * nx - 1] = psi2[ny * nx - 1] * (u[ny * nx - 1 - nx] - u[ny * nx - 1]) +
                         psi4[ny * nx - 1] * (u[ny * nx - 2]      - u[ny * nx - 1]);
    div_v[ny * nx - 1] = psi2[ny * nx - 1] * (v[ny * nx - 1 - nx] - v[ny * nx - 1]) +
                         psi4[ny * nx - 1] * (v[ny * nx - 2]      - v[ny * nx - 1]);
} // brox_spatial_divergence_u

