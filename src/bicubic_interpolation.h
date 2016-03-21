/*
 * bicubic_iterpolation.h
 */

#ifndef BICUBIC_ITERPOLATION_H_
#define BICUBIC_ITERPOLATION_H_
#include <stdbool.h>

/**
 *
 * Compute the bicubic interpolation of a point in an image.
 * Detect if the point goes outside the image domain.
 *
 **/
double
bicubic_interpolation_at (const float *input,	//image to be interpolated
                          const double uu,		//x component of the vector field
                          const double vv,		//y component of the vector field
                          const int nx,		//width of the image
                          const int ny,		//height of the image
                          const bool border_out=false	//if true, put zeros outside the region
                          );

/**
 *
 * Compute the bicubic interpolation of a point in an image.
 * Detect if the point goes outside the image domain.
 *
 **/
double
bicubic_interpolation_at_color (const float *input,	//image to be interpolated
                          const double uu,		//x component of the vector field
                          const double vv,		//y component of the vector field
                          const int nx,		//width of the image
                          const int ny,		//height of the image
                          const int nz,            //number of channels of the image
                          const int k,  		//actual channel
                          const bool border_out=false	//if true, put zeros outside the region
                          );

/**
 *
 * Compute the bicubic interpolation of an image.
 *
 **/
void
bicubic_interpolation_warp (const float *input,	//image to be warped
                            const float *u,		//x component of the vector field
                            const float *v,		//y component of the vector field
                            float *output,		//warped output image with bicubic interpolation
                            const int nx,		//width of the image
                            const int ny,		//height of the image
                            bool border_out = false	//if true, put zeros outside the region
                            );

/**
 *
 * Compute the bicubic interpolation of an image.
 *
 **/
void
bicubic_interpolation_warp_color (const float *input,	//image to be warped
                            const float *u,		//x component of the vector field
                            const float *v,		//y component of the vector field
                            float *output,		//warped output image with bicubic interpolation
                            const int nx,		//width of the image
                            const int ny,		//height of the image
                            const int nz,		//number of channels of the image	
                            bool border_out = false	//if true, put zeros outside the region
                            );

double me_interpolate_bilinear(const float *in, int ncol, double x, double y);

double me_interpolate_bicubic(const float *in, int nx, int ny, double x,
		double y);

void me_image_restriction(const float *in, float *out, int ncol, int nrow,
                          int new_ncol, int new_nrow);

#endif /* BICUBIC_ITERPOLATION_H_ */
