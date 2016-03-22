#ifndef HORN_SCHUNK_H
#define HORN_SCHUNK_H

// run n iterations of the classical Horn-Schunck method
void hs(float *u, float *v, float *a, float *b, int w, int h,
        int n, float alpha);

/**
 *
 *  Horn & Schunck method for optical flow estimation at a single scale
 *
 */
void horn_schunck_optical_flow(
	const float *I1,             // source image
	const float *I2,             // target image
	float       *u,              // x component of optical flow
	float       *v,              // y component of optical flow
	const int    nx,             // image width
	const int    ny,             // image height
	const float  alpha,          // smoothing parameter
	const int    warps,          // number of warpings per scale
	const float  TOL,            // stopping criterion threshold
	const int    maxiter,        // maximum number of iterations
	const bool   verbose         // switch on messages
                               );

/**
 *
 *  Procedure to handle the pyramidal approach.
 *  This procedure relies on the previous functions to calculate
 *  large optical flow fields using a pyramidal scheme.
 *
 */
 void horn_schunck_pyramidal(
	const float *I1,              // source image
	const float *I2,              // target image
	float       *u,               // x component of optical flow
	float       *v,               // y component of optical flow
	const int    nx,              // image width
	const int    ny,              // image height
	const float  alpha,           // smoothing weight
	const int    nscales,         // number of scales
	const float  zfactor,         // zoom factor
	const int    warps,           // number of warpings per scale
	const float  TOL,             // stopping criterion threshold
	const int    maxiter,         // maximum number of iterations
	const bool   verbose          // switch on messages
                             );

#endif
