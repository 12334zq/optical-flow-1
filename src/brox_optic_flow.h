// This program is free software: you can use, modify and/or redistribute it
// under the terms of the simplified BSD License. You should have received a
// copy of this license along this program. If not, see
// <http://www.opensource.org/licenses/bsd-license.html>.
//
// Copyright (C) 2012, Javier Sánchez Pérez <jsanchez@dis.ulpgc.es>
// All rights reserved.

#ifndef BROX_OPTIC_FLOW_H
#define BROX_OPTIC_FLOW_H

#include "of.h"

/**
 *
 *  Multiscale approach for computing the optical flow
 *
 **/
void brox_optic_flow_spatial(const ofpix_t *I1,         //first image
                             const ofpix_t *I2, //second image
                             ofpix_t *u, //x component of the optical flow
                             ofpix_t *v, //y component of the optical flow
                             const int nxx, //image width
                             const int nyy, //image height
                             const double alpha, //smoothness parameter
                             const double gamma, //gradient term parameter
                             const int nscales, //number of scales
                             const double nu, //downsampling factor
                             const double TOL, //stopping criterion threshold
                             const int inner_iter, //number of inner iterations
                             const int outer_iter, //number of outer iterations
                             const bool verbose //switch on messages
                             );


/**
 *
 *  Multiscale approach for computing the optical flow
 *
 **/
void brox_optic_flow_temporal(const ofpix_t *I,          //sequence of images
                              ofpix_t *u, //x component of the optical flow
                              ofpix_t *v, //y component of the optical flow
                              const int nxx, //image width
                              const int nyy, //image height
                              const int frames, //number of frames
                              const double alpha, //smoothness parameter
                              const double gamma, //gradient term parameter
                              const int nscales, //number of scales
                              const double nu, //downsampling factor
                              const double TOL, //stopping criterion threshold
                              const int inner_iter, //number of inner iterations
                              const int outer_iter, //number of outer iterations
                              const bool verbose //switch on messages
                              );

#endif
