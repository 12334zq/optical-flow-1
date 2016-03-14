/*
 * bicubic_iterpolation.h
 */

#ifndef BICUBIC_ITERPOLATION_H_
#define BICUBIC_ITERPOLATION_H_
#include <stdbool.h>

int neumann_bc(int x, int nx, bool *out);

/**
 *
 * Periodic boundary condition test
 *
 **/
int periodic_bc(int x, int nx, bool *out);

/**
 *
 * Symmetric boundary condition test
 *
 **/
int symmetric_bc(int x, int nx, bool *out);

/**
 *
 * Cubic interpolation in one dimension
 *
 **/
double cubic_interpolation_cell(double v[4],  //interpolation points
		double x      //point to be interpolated
		);

/**
 *
 * Bicubic interpolation in two dimensions
 *
 **/
double bicubic_interpolation_cell(double p[4][4], //array containing the interpolation points
		double x,       //x position to be interpolated
		double y        //y position to be interpolated
		);

/**
 *
 * Compute the bicubic interpolation of a point in an image.
 * Detect if the point goes outside the image domain.
 *
 **/
double bicubic_interpolation_at(const double *input, //image to be interpolated
		const double uu,    //x component of the vector field
		const double vv,    //y component of the vector field
		const int nx,    //image width
		const int ny,    //image height
		bool border_out //if true, return zero outside the region
		);

/**
 *
 * Compute the bicubic interpolation of an image.
 *
 **/
void bicubic_interpolation_warp(const double *input,     // image to be warped
		const double *u,         // x component of the vector field
		const double *v,         // y component of the vector field
		double *output,    // image warped with bicubic interpolation
		const int nx,        // image width
		const int ny,        // image height
		bool border_out // if true, put zeros outside the region
		);

double me_interpolate_bilinear(const double *in, int ncol, double x, double y);

double me_interpolate_bicubic(const double *in, int nx, int ny, double x,
		double y);

#endif /* BICUBIC_ITERPOLATION_H_ */
