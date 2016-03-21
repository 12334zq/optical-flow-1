// This program is free software: you can use, modify and/or redistribute it
// under the terms of the simplified BSD License. You should have received a
// copy of this license along this program. If not, see
// <http://www.opensource.org/licenses/bsd-license.html>.
//
// Copyright (C) 2012, Javier Sánchez Pérez <jsanchez@dis.ulpgc.es>
// All rights reserved.

#ifndef BROX_OPTIC_FLOW_H
#define BROX_OPTIC_FLOW_H

#include <stdbool.h>

/**
  *
  *  Multiscale approach for computing the optical flow
  *
**/
void brox_optic_flow_spatial(
    const float *I1,         //first image
    const float *I2,         //second image
    float *u, 		      //x component of the optical flow
    float *v, 		      //y component of the optical flow
    const int    nxx,        //image width
    const int    nyy,        //image height
    const float  alpha,      //smoothness parameter
    const float  gamma,      //gradient term parameter
    const int    nscales,    //number of scales
    const float  nu,         //downsampling factor
    const float  TOL,        //stopping criterion threshold
    const int    inner_iter, //number of inner iterations
    const int    outer_iter, //number of outer iterations
    const bool   verbose     //switch on messages
                     );


/**
  *
  *  Multiscale approach for computing the optical flow
  *
**/
void brox_optic_flow_temporal(
    const float *I,          //sequence of images
    float *u,                //x component of the optical flow
    float *v,                //y component of the optical flow
    const int    nxx,        //image width
    const int    nyy,        //image height
    const int    frames,     //number of frames
    const float  alpha,      //smoothness parameter
    const float  gamma,      //gradient term parameter
    const int    nscales,    //number of scales
    const float  nu,         //downsampling factor
    const float  TOL,        //stopping criterion threshold
    const int    inner_iter, //number of inner iterations
    const int    outer_iter, //number of outer iterations
    const bool   verbose     //switch on messages
                              );

#endif
