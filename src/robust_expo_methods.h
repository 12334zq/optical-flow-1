// This program is free software: you can use, modify and/or redistribute it
// under the terms of the simplified BSD License. You should have received a
// copy of this license along this program. If not, see
// <http://www.opensource.org/licenses/bsd-license.html>.
//
// Copyright (C) 2012, Javier Sánchez Pérez <jsanchez@dis.ulpgc.es>
// Copyright (C) 2014, Agustín Salgado de la Nuez <asalgado@dis.ulpgc.es>
// Copyright (C) 2014, 2015 Nelson Monzón López <nmonzon@ctim.es>
// All rights reserved.

#ifndef OF_ROBUST_EXPO_METHODS_H
#define OF_ROBUST_EXPO_METHODS_H

#include "of.h"

/**
 *
 *  Multiscale approach for computing the optical flow
 *
 **/
void robust_expo_methods(const ofpix_t *I1,          // first image
                         const ofpix_t *I2, // second image
                         ofpix_t *u, // x component of the optical flow
                         ofpix_t *v, // y component of the optical flow
                         const int nxx, // image width
                         const int nyy, // image height
                         const int nzz, // number of color channels in image
                         const int method_type, // choose the diffusion strategy
                         const float alpha, // smoothness parameter
                         const float gamma, // gradient term parameter
                         const float lambda, // coefficient parameter for the decreasing function (if needed)
                         const int nscales, // number of scales
                         const float nu, // downsampling factor
                         const float TOL, // stopping criterion threshold
                         const int inner_iter, // number of inner iterations
                         const int outer_iter, // number of outer iterations
                         const bool verbose // switch on messages
                         );

#endif
