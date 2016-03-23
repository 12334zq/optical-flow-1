
// This program is free software: you can use, modify and/or redistribute it
// under the terms of the simplified BSD License. You should have received a
// copy of this license along this program. If not, see
// <http://www.opensource.org/licenses/bsd-license.html>.
//
// Copyright (C) 2011, Javier Sánchez Pérez <jsanchez@dis.ulpgc.es>
// All rights reserved.


#ifndef DUAL_TVL1_OPTIC_FLOW_H
#define DUAL_TVL1_OPTIC_FLOW_H

#include "of.h"

/**
 * Implementation of the Zach, Pock and Bischof dual TV-L1 optic flow method
 *
 * see reference:
 *  [1] C. Zach, T. Pock and H. Bischof, "A Duality Based Approach for Realtime
 *      TV-L1 Optical Flow", In Proceedings of Pattern Recognition (DAGM),
 *      Heidelberg, Germany, pp. 214-223, 2007
 *
 *
 * Details on the total variation minimization scheme can be found in:
 *  [2] A. Chambolle, "An Algorithm for Total Variation Minimization and
 *      Applications", Journal of Mathematical Imaging and Vision, 20: 89-97, 2004
 **/


/**
 *
 * Function to compute the optical flow in one scale
 *
 **/
void Dual_TVL1_optic_flow(ofpix_t *I0,           // source image
                          ofpix_t *I1, // target image
                          ofpix_t *u1, // x component of the optical flow
                          ofpix_t *u2, // y component of the optical flow
                          const int nx, // image width
                          const int ny, // image height
                          const double tau, // time step
                          const double lambda, // weight parameter for the data term
                          const double theta, // weight parameter for (u - v)²
                          const int warps, // number of warpings per scale
                          const double epsilon, // tolerance for numerical convergence
                          const bool verbose // enable/disable the verbose mode
                          );


/**
 *
 * Function to compute the optical flow using multiple scales
 *
 **/
void Dual_TVL1_optic_flow_multiscale(ofpix_t *I0,           // source image
                                     ofpix_t *I1, // target image
                                     ofpix_t *u1, // x component of the optical flow
                                     ofpix_t *u2, // y component of the optical flow
                                     const int nxx, // image width
                                     const int nyy, // image height
                                     const double tau, // time step
                                     const double lambda, // weight parameter for the data term
                                     const double theta, // weight parameter for (u - v)²
                                     const int nscales, // number of scales
                                     const double zfactor, // factor for building the image piramid
                                     const int warps, // number of warpings per scale
                                     const double epsilon, // tolerance for numerical convergence
                                     const bool verbose // enable/disable the verbose mode
                                     );

#endif//DUAL_TVL1_OPTIC_FLOW_H
