#ifndef HORN_SCHUNK_H
#define HORN_SCHUNK_H

#include "of.h"

// run n iterations of the classical Horn-Schunck method
void hs(ofpix_t *u, ofpix_t *v, ofpix_t *a, ofpix_t *b, int w, int h,
        int n, double alpha);

/**
 *
 *  Horn & Schunck method for optical flow estimation at a single scale
 *
 */
void horn_schunck_optical_flow(const ofpix_t *I1,             // source image
                               const ofpix_t *I2, // target image
                               ofpix_t       *u,// x component of optical flow
                               ofpix_t       *v,// y component of optical flow
                               const int nx, // image width
                               const int ny, // image height
                               const double alpha, // smoothing parameter
                               const int warps, // number of warpings per scale
                               const double TOL, // stopping criterion threshold
                               const int maxiter, // maximum number of iterations
                               const bool verbose // switch on messages
                               );

/**
 *
 *  Procedure to handle the pyramidal approach.
 *  This procedure relies on the previous functions to calculate
 *  large optical flow fields using a pyramidal scheme.
 *
 */
void horn_schunck_pyramidal(const ofpix_t *I1, // source image
                            const ofpix_t *I2, // target image
                            ofpix_t       *u,// x component of optical flow
                            ofpix_t       *v,// y component of optical flow
                            const int nx, // image width
                            const int ny, // image height
                            const double alpha, // smoothing weight
                            const int nscales, // number of scales
                            const double zfactor, // zoom factor
                            const int warps, // number of warpings per scale
                            const double TOL, // stopping criterion threshold
                            const int maxiter, // maximum number of iterations
                            const bool verbose // switch on messages
                            );

#endif
