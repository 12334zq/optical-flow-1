#ifndef OF_UTILS_H
#define OF_UTILS_H

#include "of.h"

ofpix_t *me_sgauss(double std, int n);
void me_sepconvol(const ofpix_t *in, ofpix_t *out, int nx, int ny, const ofpix_t *filter_x,
		const ofpix_t *filter_y, int size_x, int size_y);

void me_median_filtering(ofpix_t *in, int nx, int ny, int wsize);

/**
 *
 * Function to normalize the images between 0 and 255
 *
 **/
void image_normalization_1(
                           const ofpix_t *I,    //input images
                           float       *In,   //normalized output images
                           int          size  //size of the images
);

/**
 *
 * Function to normalize the images between 0 and 255
 *
 **/
void image_normalization_2(
                           const ofpix_t *I1,   //input image 1
                           const ofpix_t *I2,   //input image 2
                           ofpix_t       *I1n,  //normalized output image 1
                           ofpix_t       *I2n,  //normalized output image 2
                           int          size  //size of the image
);

/*
 * Function name: image_normalization_3
 * Author:
 * Parameters:
 * 		I0 : (in/out) vectorized matrix.
 * 		I1 : (in/out) vectorized matrix.
 * 		I3 : (in/out) vectorized matrix.
 * 		filtI0 : (in) vectorized matrix.
 * 		size   : (in) Integer. size of the vectorized matrix
 * Output:
 * 		VOID
 * Description: Computes the normalization to [0,255] of I0, I1, I2 and returns the normalized vector in the same
 * 				variables. For normalization, divides each component of  I0, I1, I2 by the maximum
 * 			of the three matrix. The number 3 comes from the number of images.
 *
 */
void image_normalization_3(ofpix_t *I0, ofpix_t *I1, ofpix_t *I2, int size);

/*
 * Function name: image_normalization_4
 * Author:
 * Parameters:
 * 		I_1 : (in) vectorized matrix.
 * 		I0 : (in) vectorized matrix.
 * 		I1 : (in) vectorized matrix.
 * 		filtI0 : (in) vectorized matrix.
 * 		I_1n : (out) Normalized vectorized matrix.
 * 		I0n : (out) Normalized vectorized matrix.
 * 		I1n : (out) Normalized vectorized matrix.
 * 		filtI0n : (out) Normalized vectorized matrix.
 * 		size   : (in) Integer. size of the vectorized matrix
 * Output:
 * 		VOID
 * Description: Computes the normalization to [0,255] of I_1, I0, I1, filtI0 and returns the normalized vector in
 * 			I_1n, I0n, I1n, filtI0n. For normalization, divides each component of  I_1, I0, I1, filtI0 by the maximum
 * 			of the four matrix. The number 4 comes from the number of images.
 *
 */
void image_normalization_4(const ofpix_t *I_1, const ofpix_t *I0,
		const ofpix_t *I1, const ofpix_t *filtI0, ofpix_t *I_1n, ofpix_t *I0n,
		ofpix_t *I1n, ofpix_t *filtI0n, int size);

/*
 * Function name: L2error
 * Author:
 * Parameters:
 * 		u1prev  : (in/out) Vectorized matrix corresponding to the first component of uprev.
 * 		u2prev  : (in/out) Vectorized matrix corresponding to the second component of uprev.
 * 		u1:     : (in) Vectorized matrix corresponding to the first component of u.
 * 		u2:     : (in) Vectorized matrix corresponding to the first component of u.
 * 		size    : Number of elements of the matrix. It shoud be the same for  u1, u2, u1prev and u2prev
 * Output:
 * 		error   : mean (|u-uprev|^2);
 * Description: Let u=(u1, u2) and uprev=(u1prev, u2prev) be two different flow fields.
 *              L2error performs the mean(|u-uprev|^2) where |á| is the euclidean norm.
 *              After computing the L2error, uprev is replaced by u.
 *
 */

/*
 * Function name: getminmax
 * Author:
 * Parameters:
 * 		min : (out) Scalar.
 * 		max : (out) Scalar.
 * 		x   : (in) vector
 *		n   : (in) Scalar integer corresponding to the length of x.
 * Output:
 * 		VOID
 * Description: Returns in min, max the minimum and maximum of x.
 *
 */
void getminmax(ofpix_t *min, ofpix_t *max, const ofpix_t *x, int n);

#endif
