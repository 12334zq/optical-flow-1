// This program is free software: you can use, modify and/or redistribute it
// under the terms of the simplified BSD License. You should have received a
// copy of this license along this program. If not, see
// <http://www.opensource.org/licenses/bsd-license.html>.
//
// Copyright (C) 2012, Javier Sánchez Pérez <jsanchez@dis.ulpgc.es>
// All rights reserved.

#ifndef BROX_TEMPORAL_MASK_H
#define BROX_TEMPORAL_MASK_H

#include "of.h"

/**
 *
 * Compute the coefficients of the divergence term
 *
 */
void brox_temporal_psi_divergence(
    const ofpix_t *psi, //robust functional
    ofpix_t *psi1,      //coefficients of divergence
    ofpix_t *psi2,      //coefficients of divergence
    ofpix_t *psi3,      //coefficients of divergence
    ofpix_t *psi4,      //coefficients of divergence    
    ofpix_t *psi5,      //coefficients of divergence 
    ofpix_t *psi6,      //coefficients of divergence 
    const int nx,     //image width
    const int ny,     //image height
    const int nz      //image depth
                                  );


/**
 *
 * Compute the divergence of the optical flow
 *
 */
void brox_temporal_divergence_u(
    const ofpix_t *u,    //x component of optical flow
    const ofpix_t *v,    //y component of optical flow
    const ofpix_t *psi1, //coefficients of divergence
    const ofpix_t *psi2, //coefficients of divergence
    const ofpix_t *psi3, //coefficients of divergence
    const ofpix_t *psi4, //coefficients of divergence
    const ofpix_t *psi5, //coefficients of divergence
    const ofpix_t *psi6, //coefficients of divergence
    ofpix_t *div_u,      //computed divergence for u
    ofpix_t *div_v,      //computed divergence for v
    const int nx,      //image width 
    const int ny,      //image height
    const int nz       //image depth
                                );


#endif
