#ifndef UTILS_H
#define UTILS_H

double *me_sgauss(double std, int n);
void me_sepconvol(double *in, double *out, int nx, int ny, double *filter_x,
		double *filter_y, int size_x, int size_y);

void me_save_image(double *in, int nx, int ny);
int me_median_compare(const void *i, const void *j);
void me_median_filtering(double *in, int nx, int ny, int wsize);

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
void image_normalization_3(double *I0, double *I1, double *I2, int size);

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
void image_normalization_4(const double *I_1, const double *I0,
		const double *I1, const double *filtI0, double *I_1n, double *I0n,
		double *I1n, double *filtI0n, int size);

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
void getminmax(double *min, double *max, const double *x, int n);

#endif
