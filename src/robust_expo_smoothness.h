// This program is free software: you can use, modify and/or redistribute it
// under the terms of the simplified BSD License. You should have received a
// copy of this license along this program. If not, see
// <http://www.opensource.org/licenses/bsd-license.html>.
//
// Copyright (C) 2012, Javier Sánchez Pérez <jsanchez@dis.ulpgc.es>
// Copyright (C) 2014, Agustín Salgado de la Nuez <asalgado@dis.ulpgc.es>
// Copyright (C) 2014, 2015 Nelson Monzón López <nmonzon@ctim.es>
// All rights reserved.

#ifndef OF_ROBUST_EXPO_SMOOTHNESS_H
#define OF_ROBUST_EXPO_SMOOTHNESS_H

#include "of.h"

#define ROBUST_EXPO_EPSILON 0.001

/**
 *
 * Compute the coefficients of the robust functional (smoothness term)
 *
 **/
void robust_expo_psi_smooth(const ofpix_t *ux, //gradient of x component of the optical flow
                            const ofpix_t *uy, //gradient of x component of the optical flow
                            const ofpix_t *vx, //gradient of y component of the optical flow
                            const ofpix_t *vy, //gradient of y component of the optical flow
                            const ofpix_t *expo, //exponential smoothing factor
                            const int size_flow,
                            ofpix_t       *psi       //output coefficients
                            );

/**
**  Calculate the exponential values.
**
**/
void robust_expo_exponential_calculation(const ofpix_t *Ix,      // Computed Image 1 derivative in x
                                         const ofpix_t *Iy, // Computed Image 1 derivative in y
                                         const int size_flow, // size of the flow field
                                         const int size, // size of the multi-channel image
                                         const int nz, // nº of image channels
                                         const double alpha, // smoothness weight
                                         const double lambda, // Coeffient for decreasing function
                                         const int method_type, // (1 = DF, 2 = DF_BETA, 3 = DF_AUTO)
                                         ofpix_t       *expo// e^(lambda * DI)
                                         );

#endif
