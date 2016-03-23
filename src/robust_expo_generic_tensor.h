// This program is free software: you can use, modify and/or redistribute it
// under the terms of the simplified BSD License. You should have received a
// copy of this license along this program. If not, see
// <http://www.opensource.org/licenses/bsd-license.html>.
//
// Copyright (C) 2012, Javier Sánchez Pérez <jsanchez@dis.ulpgc.es>
// Copyright (C) 2014, Nelson Monzón López <nmonzon@ctim.es>
// Copyright (C) 2014, Agustín Salgado de la Nuez <asalgado@dis.ulpgc.es>
// All rights reserved.

#ifndef OF_ROBUST_EXPO_GENERIC_TENSOR_H
#define OF_ROBUST_EXPO_GENERIC_TENSOR_H

#include "of.h"

/**
 *
 * Compute the coefficients of the divergence term
 *
 */
void robust_expo_psi_divergence(ofpix_t *psi1,      //coefficients of divergence
                                ofpix_t *psi2,
                                ofpix_t *psi3,
                                ofpix_t *psi4,
                                const ofpix_t *psi, //robust functional
                                const int nx, //image width
                                const int ny //image height
                                );

/**
 *
 * Compute the divergence of the optical flow
 *
 */
void robust_expo_divergence(const ofpix_t *u,    //component of optical flow
                            const ofpix_t *psi1, //coefficients of divergence
                            const ofpix_t *psi2,
                            const ofpix_t *psi3,
                            const ofpix_t *psi4,
                            const int nx, //image width
                            const int ny, //image height
                            ofpix_t *div //computed divergence for u
                            );


#endif
