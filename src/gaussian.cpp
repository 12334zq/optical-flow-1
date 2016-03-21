// This program is free software: you can use, modify and/or redistribute it
// under the terms of the simplified BSD License. You should have received a
// copy of this license along this program. If not, see
// <http://www.opensource.org/licenses/bsd-license.html>.
//
// Copyright (C) 2012, Javier Sánchez Pérez <jsanchez@dis.ulpgc.es>
// All rights reserved.

#include "gaussian.h"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include "xmalloc.h"


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
    const int boundary_condition,       //boundary condition
    const int window_size //defines the size of the window
)
{
	const double den  = 2*sigma*sigma;
	const int   size = (int) (window_size * sigma) + 1 ;
	const int   bdx  = xdim + size;
	const int   bdy  = ydim + size;

	if (boundary_condition && size > xdim) {
		fprintf(stderr, "GaussianSmooth: sigma too large\n");
		abort();
	}

	// compute the coefficients of the 1D convolution kernel
	double B[size];
	for(int i = 0; i < size; i++)
		B[i] = 1 / (sigma * sqrt(2.0 * 3.1415926)) * exp(-i * i / den);

	// normalize the 1D convolution kernel
	double norm = 0;
	for(int i = 0; i < size; i++)
		norm += B[i];
	norm *= 2;
	norm -= B[0];
	for(int i = 0; i < size; i++)
		B[i] /= norm;

	// convolution of each line of the input image
	double *R = (double*)xmalloc((size + xdim + size)*sizeof*R);

	for (int k = 0; k < ydim; k++)
	{
		int i, j;
		for (i = size; i < bdx; i++)
			R[i] = I[k * xdim + i - size];

		switch (boundary_condition)
		{
		case BOUNDARY_CONDITION_DIRICHLET:
			for(i = 0, j = bdx; i < size; i++, j++)
				R[i] = R[j] = 0;
			break;

		case BOUNDARY_CONDITION_REFLECTING:
			for(i = 0, j = bdx; i < size; i++, j++) {
				R[i] = I[k * xdim + size-i];
				R[j] = I[k * xdim + xdim-i-1];
			}
			break;

		case BOUNDARY_CONDITION_PERIODIC:
			for(i = 0, j = bdx; i < size; i++, j++) {
				R[i] = I[k * xdim + xdim-size+i];
				R[j] = I[k * xdim + i];
			}
			break;
		}

		for (i = size; i < bdx; i++)
		{
			double sum = B[0] * R[i];
			for (j = 1; j < size; j++ )
				sum += B[j] * ( R[i-j] + R[i+j] );
			I[k * xdim + i - size] = sum;
		}
	}

	// convolution of each column of the input image
	double *T = (double*)xmalloc((size + ydim + size)*sizeof*T);

	for (int k = 0; k < xdim; k++)
	{
		int i, j;
		for (i = size; i < bdy; i++)
			T[i] = I[(i - size) * xdim + k];

		switch (boundary_condition)
		{
		case BOUNDARY_CONDITION_DIRICHLET:
			for (i = 0, j = bdy; i < size; i++, j++)
				T[i] = T[j] = 0;
			break;

		case BOUNDARY_CONDITION_REFLECTING:
			for (i = 0, j = bdy; i < size; i++, j++) {
				T[i] = I[(size-i) * xdim + k];
				T[j] = I[(ydim-i-1) * xdim + k];
			}
			break;

		case BOUNDARY_CONDITION_PERIODIC:
			for( i = 0, j = bdx; i < size; i++, j++) {
				T[i] = I[(ydim-size+i) * xdim + k];
				T[j] = I[i * xdim + k];
			}
			break;
		}

		for (i = size; i < bdy; i++)
		{
			double sum = B[0] * T[i];
			for (j = 1; j < size; j++ )
				sum += B[j] * (T[i-j] + T[i+j]);
			I[(i - size) * xdim + k] = sum;
		}
	}

	free(R);
	free(T);
}



#if 0
/**
 *
 * Convolution with a Gaussian
 *
 */
inline void gaussian(
    const float *in,      //input image
    float *out,           //output image
    const int xdim,       //image width
    const int ydim,       //image height
    const double sigma,   //Gaussian sigma
    const int bc=1,       //boundary condition
    const int precision=5 //defines the size of the window
)
{
    for(int i = 0; i < xdim*ydim; i++) out[i] = in[i];

    gaussian(out, xdim, ydim, sigma, bc, precision);
}
#endif

