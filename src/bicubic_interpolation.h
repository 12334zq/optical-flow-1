/*
 * bicubic_iterpolation.h
 */

#ifndef OF_BICUBIC_ITERPOLATION_H
#define OF_BICUBIC_ITERPOLATION_H

#include "of.h"

/**
 *
 * Compute the bicubic interpolation of a point in an image.
 * Detect if the point goes outside the image domain.
 *
 **/
double bicubic_interpolation_at (const ofpix_t *input, //image to be interpolated
                                 const double uu, //x component of the vector field
                                 const double vv, //y component of the vector field
                                 const int nx, //width of the image
                                 const int ny, //height of the image
                                 const bool border_out = false //if true, put zeros outside the region
                                 );

/**
 *
 * Compute the bicubic interpolation of a point in an image.
 * Detect if the point goes outside the image domain.
 *
 **/
double bicubic_interpolation_at_color (const ofpix_t *input,   //image to be interpolated
                                       const double uu, //x component of the vector field
                                       const double vv, //y component of the vector field
                                       const int nx, //width of the image
                                       const int ny, //height of the image
                                       const int nz, //number of channels of the image
                                       const int k, //actual channel
                                       const bool border_out = false //if true, put zeros outside the region
                                       );

/**
 *
 * Compute the bicubic interpolation of an image.
 *
 **/
void bicubic_interpolation_warp (const ofpix_t *input,   //image to be warped
                                 const ofpix_t *u,  //x component of the vector field
                                 const ofpix_t *v,  //y component of the vector field
                                 ofpix_t *output,   //warped output image with bicubic interpolation
                                 const int nx,  //width of the image
                                 const int ny,  //height of the image
                                 bool border_out = false //if true, put zeros outside the region
                                 );

/**
 *
 * Compute the bicubic interpolation of an image.
 *
 **/
void bicubic_interpolation_warp_color (const ofpix_t *input, //image to be warped
                                       const ofpix_t *u, //x component of the vector field
                                       const ofpix_t *v, //y component of the vector field
                                       ofpix_t *output, //warped output image with bicubic interpolation
                                       const int nx, //width of the image
                                       const int ny, //height of the image
                                       const int nz, //number of channels of the image
                                       bool border_out = false //if true, put zeros outside the region
                                       );

double me_interpolate_bilinear(const ofpix_t *in, int ncol, double x, double y);

double me_interpolate_bicubic(const ofpix_t *in, int nx, int ny, double x,
                              double y);

void me_image_restriction(const ofpix_t *in, ofpix_t *out, int ncol, int nrow,
                          int new_ncol, int new_nrow);

#endif /* BICUBIC_ITERPOLATION_H_ */
