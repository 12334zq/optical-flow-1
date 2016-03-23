// This program is free software: you can use, modify and/or redistribute it
// under the terms of the simplified BSD License. You should have received a
// copy of this license along this program. If not, see
// <http://www.opensource.org/licenses/bsd-license.html>.
//
// Copyright (C) 2012, Javier Sánchez Pérez <jsanchez@dis.ulpgc.es>
// All rights reserved.

#ifndef BROX_SPATIAL_MASK_H
#define BROX_SPATIAL_MASK_H

#include "of.h"

/**
 *
 * Compute the second order X derivative
 *
 */
void brox_spatial_Dxx(const ofpix_t *I, //input image
                      ofpix_t *Ixx, //oputput derivative
                      const int nx, //image width
                      const int ny //image height
                      );


/**
 *
 * Compute the second order Y derivative
 *
 */
void brox_spatial_Dyy(const ofpix_t *I, //input image
                      ofpix_t *Iyy, //oputput derivative
                      const int nx, //image width
                      const int ny //image height
                      );


/**
 *
 * Compute the second order XY derivative
 *
 */
void brox_spatial_Dxy(const ofpix_t *I, //input image
                      ofpix_t *Ixy, //oputput derivative
                      const int nx, //image width
                      const int ny //image height
                      );


/**
 *
 * Compute the coefficients of the divergence term
 *
 */
void brox_spatial_psi_divergence(const ofpix_t *psi, //robust functional
                                 ofpix_t *psi1, //coefficients of divergence
                                 ofpix_t *psi2, //coefficients of divergence
                                 ofpix_t *psi3, //coefficients of divergence
                                 ofpix_t *psi4, //coefficients of divergence
                                 const int nx, //image width
                                 const int ny //image height
                                 );


/**
 *
 * Compute the divergence of the optical flow
 *
 */
void brox_spatial_divergence_u(const ofpix_t *u,    //x component of optical flow
                               const ofpix_t *v, //y component of optical flow
                               const ofpix_t *psi1, //coefficients of divergence
                               const ofpix_t *psi2, //coefficients of divergence
                               const ofpix_t *psi3, //coefficients of divergence
                               const ofpix_t *psi4, //coefficients of divergence
                               ofpix_t *div_u, //computed divergence for u
                               ofpix_t *div_v, //computed divergence for v
                               const int nx, //image width
                               const int ny //image height
                               );

#endif // ifndef BROX_SPATIAL_MASK_H
