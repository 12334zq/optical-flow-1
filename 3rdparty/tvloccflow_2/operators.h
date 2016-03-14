/*
 * operators.h
 *
 *  Created on: 03/12/2013
 *      Author: juanfran
 */

#ifndef OPERATORS_H_
#define OPERATORS_H_

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define BOUNDARY_CONDITION_DIRICHLET 0
#define BOUNDARY_CONDITION_REFLECTING 1
#define BOUNDARY_CONDITION_PERIODIC 2

#define DEFAULT_GAUSSIAN_WINDOW_SIZE 5
#define DEFAULT_BOUNDARY_CONDITION BOUNDARY_CONDITION_REFLECTING

void forward_gradient(const double *f, //input image
		double *fx,      //computed x derivative
		double *fy,      //computed y derivative
		const int nx,   //image width
		const int ny    //image height
);

void centered_gradient(double *input, double *dx, double *dy, const int nx,
		const int ny);
void gaussian(double *I,             // input/output image
		const int xdim,       // image width
		const int ydim,       // image height
		const double sigma    // Gaussian sigma
);
void divergence(const double *v1, const double *v2, double *div, const int nx,
		const int ny);
#endif /* OPERATORS_H_ */
