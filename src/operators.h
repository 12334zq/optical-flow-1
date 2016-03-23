// This program is free software: you can use, modify and/or redistribute it
// under the terms of the simplified BSD License. You should have received a
// copy of this license along this program. If not, see
// <http://www.opensource.org/licenses/bsd-license.html>.
//
// Copyright (C) 2012, Javier Sánchez Pérez <jsanchez@dis.ulpgc.es>
// All rights reserved.

#ifndef OF_OPERATORS_H
#define OF_OPERATORS_H

#include "of.h"

/**
 *
 * Details on how to compute the divergence and the grad(u) can be found in:
 * [2] A. Chambolle, "An Algorithm for Total Variation Minimization and
 * Applications", Journal of Mathematical Imaging and Vision, 20: 89-97, 2004
 *
 **/


/**
 *
 * Function to compute the divergence with backward differences
 * (see [2] for details)
 *
 **/
void divergence(const ofpix_t *v1, // x component of the vector field
                const ofpix_t *v2, // y component of the vector field
                ofpix_t *div, // output divergence
                const int nx, // image width
                const int ny // image height
                );

/**
 *
 * Function to compute the gradient with forward differences
 * (see [2] for details)
 *
 **/
void forward_gradient(const ofpix_t *f, //input image
                      ofpix_t *fx, //computed x derivative
                      ofpix_t *fy, //computed y derivative
                      const int nx, //image width
                      const int ny //image height
                      );


/**
 *
 * Compute the second order X derivative
 *
 */
void Dxx(const ofpix_t *I, //input image
         ofpix_t *Ixx, //oputput derivative
         const int nx, //image width
         const int ny,      //image height
         const int nz               //number of color channels in the image
         );


/**
 *
 * Compute the second order Y derivative
 *
 */
void Dyy(const ofpix_t *I, //input image
         ofpix_t *Iyy, //oputput derivative
         const int nx, //image width
         const int ny,      //image height
         const int nz               //number of color channels in the image
         );


/**
 *
 * Compute the second order XY derivative
 *
 */
void Dxy(const ofpix_t *I, //input image
         ofpix_t *Ixy, //oputput derivative
         const int nx, //image width
         const int ny,      //image height
         const int nz               //number of color channels in the image
         );

/**
 *
 * Function to compute the gradient with centered differences
 *
 */
void centered_gradient(const ofpix_t *input, //input image
                       ofpix_t *dx, //computed x derivative
                       ofpix_t *dy, //computed y derivative
                       const int nx, //image width
                       const int ny,        //image height
                       const int nz               //number of color channels in the image
                       );


/**
 *
 * Function to compute the 3D gradient with centered differences
 *
 */
void centered_gradient3(const ofpix_t *input, //input image
                        ofpix_t *dx, //x derivative
                        ofpix_t *dy, //y derivative
                        ofpix_t *dz, //z derivative
                        const int nx, //image width
                        const int ny, //image height
                        const int nz //image depth
                        );

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
void gaussian(ofpix_t *I,   //input/output image
              const int xdim, //image width
              const int ydim, //image height
              const double sigma, //Gaussian sigma
              const int bc = DEFAULT_BOUNDARY_CONDITION, //boundary condition
              const int precision = DEFAULT_GAUSSIAN_WINDOW_SIZE //defines the size of the window
              );


#endif // ifndef OF_OPERATORS_H
