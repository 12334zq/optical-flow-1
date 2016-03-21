// This program is free software: you can use, modify and/or redistribute it
// under the terms of the simplified BSD License. You should have received a
// copy of this license along this program. If not, see
// <http://www.opensource.org/licenses/bsd-license.html>.
//
// Copyright (C) 2012, Javier Sánchez Pérez <jsanchez@dis.ulpgc.es>
// All rights reserved.


#ifndef GAUSSIAN_H
#define GAUSSIAN_H


#define BOUNDARY_CONDITION_DIRICHLET 0
#define BOUNDARY_CONDITION_REFLECTING 1
#define BOUNDARY_CONDITION_PERIODIC 2

#define DEFAULT_GAUSSIAN_WINDOW_SIZE 5
#define DEFAULT_BOUNDARY_CONDITION BOUNDARY_CONDITION_REFLECTING

/**
 *
 * Convolution with a Gaussian
 *
 */
 void gaussian(
    float *I,             //input/output image
    const int xdim,       //image width
    const int ydim,       //image height
    const double sigma,   //Gaussian sigma
    const int bc=1,       //boundary condition
    const int precision=5 //defines the size of the window
               );

#endif

