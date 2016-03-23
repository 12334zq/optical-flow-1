// This program is free software: you can use, modify and/or redistribute it
// under the terms of the simplified BSD License. You should have received a
// copy of this license along this program. If not, see
// <http://www.opensource.org/licenses/bsd-license.html>.
//
// Copyright (C) 2012, Coloma Ballester <coloma.ballester@upf.edu>
// Copyright (C) 2013-2014 J. F. Garamendi <jf.garamendi@upf.edu>
// All rights reserved.

#include "utils.h"

#include <cmath>
#include <cstdlib>

ofpix_t *me_sgauss(double std, int n) {
	int i, shift;
	double sum, v;
	ofpix_t *out;

	out = new ofpix_t[n];

	shift = -0.5 * (double) (n - 1);

	if (n == 1) {
		out[0] = 1.0;
	} else {
		/* store Gaussian signal */
		for (i = (n + 1) / 2; i--;) {
			v = ((double) i + (double) shift) / (double) std;
			out[i] = out[n - 1 - i] = exp(-0.5 * v * v);
		}
		/* normalize to get unit mass */
		for (sum = 0.0, i = n; i--;)
			sum += (double) out[i];
		for (i = n; i--;)
			out[i] /= sum;
	}

	return (out);
}

void me_sepconvol(const ofpix_t *in, ofpix_t *out, int nx, int ny, const ofpix_t *filter_x,
		const ofpix_t *filter_y, int size_x, int size_y) {
	ofpix_t *tmp;
	int nx1, ny1, x, y, org, i, s;

	double sum;

	nx1 = nx - 1;
	ny1 = ny - 1;

	/* Initialize temporal image */

	tmp = new ofpix_t[nx * ny];

	/* convolution along x axis */

	org = (size_x - 1) >> 1;

	for (y = 0; y <= ny1; y++) {
		for (x = 0; x <= nx1; x++) {
			sum = 0.0;

			for (i = size_x; i--;) {
				s = x - (i - org);

				while ((s < 0) || (s > nx1)) {
					if (s < 0)
						s = 0 - s - 1;

					if (s > nx1)
						s = nx1 - (s - nx1) + 1;
				}

				sum += filter_x[i] * in[y * nx + s];
			}

			tmp[y * nx + x] = sum;
		}
	}

	/* convolution along y axis */

	org = (size_y - 1) >> 1;

	for (y = 0; y <= ny1; y++) {
		for (x = 0; x <= nx1; x++) {
			sum = 0.0;

			for (i = size_y; i--;) {
				s = y - (i - org);

				while ((s < 0) || (s > ny1)) {
					if (s < 0)
						s = 0 - s - 1;

					if (s > ny1)
						s = ny1 - (s - ny1) + 1;
				}

				sum += filter_y[i] * tmp[s * nx + x];
			}

			out[y * nx + x] = sum;
		}
	}

	/* Free memory */

	delete [] tmp;
}

static
int me_median_compare(const void *i, const void *j) {
	ofpix_t *val1, *val2;

	val1 = (ofpix_t *) i;
	val2 = (ofpix_t *) j;

	if (*val1 < *val2)
		return -1;

	if (*val1 > *val2)
		return 1;

	return 0;
}

void me_median_filtering(ofpix_t *in, int nx, int ny, int wsize)
{
    ofpix_t *median_vector = new ofpix_t[wsize * wsize];
    ofpix_t *out;
	int i, pixels, nrow, ncol, border;
	int x, y, xx, yy, xx0, yy0;

	border = wsize >> 1;
	out = new ofpix_t[nx * ny];

	nrow = ny;
	ncol = nx;

	for (x = 0; x < ncol; x++) {
		for (y = 0; y < nrow; y++) {
			i = 0;

			for (yy = y - border; yy <= y + border; yy++)
				for (xx = x - border; xx <= x + border; xx++) {
					xx0 = xx;
					yy0 = yy;

					/* Symmetry */

					if (xx0 < 0)
						xx0 = -xx0 - 1;

					if (xx0 >= ncol)
						xx0 = 2 * ncol - xx0 - 1;

					if (yy0 < 0)
						yy0 = -yy0 - 1;

					if (yy0 >= nrow)
						yy0 = 2 * nrow - yy0 - 1;

					/* Access data */

					median_vector[i++] = in[yy0 * ncol + xx0];
				}

            std::qsort(median_vector, i, sizeof(ofpix_t), me_median_compare);
			out[y * ncol + x] = median_vector[i / 2];
		}
	}

	pixels = nx * ny;

	for (i = 0; i < pixels; i++)
		in[i] = out[i];

    delete [] median_vector;
	delete [] out;
}

/* BORRAR SIN COMPARARLO CON NADA
 void warping(
 const double *input, const double *u, const double *v,
 double *output, const int nx, const int ny
 )
 {
 #pragma omp parallel
 {
 #pragma omp for schedule(dynamic) nowait
 for(int j = 0; j < ny; j++)
 for(int i = 0; i < nx; i++)
 {

 const double uu = (double) (i + u[i + nx * j]);
 const double vv = (double) (j + v[i + nx * j]);

 if ((uu < 0) || (uu > (nx - 1)) ||
 (vv < 0) || (vv > (ny - 1)))
 {
 output[i + nx * j] = 0;
 }
 else
 {
 output[i + nx * j] = me_interpolate_bicubic(input, nx, ny, uu, vv);
 }
 }
 }
 }
 */


/**
 *
 * Function to normalize the images between 0 and 255
 *
 **/
void image_normalization_1(
                         const ofpix_t *I,    //input images
                         float       *In,   //normalized output images
                         int          size  //size of the images
)
{
    //find the max and min values of the images
    ofpix_t max, min;

    getminmax(&min, &max, I, size);
    const ofpix_t den = max - min;

    if(den > 0)
        //normalize the images between 0 and 255
        for(int i = 0; i < size; i++)
            In[i] = 255.0 * (I[i] - min) / den;

    else
        //copy the original data
        for(int i = 0; i < size; i++)
            In[i] = I[i];
}

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
)
{
    //compute the max and min values of the images
    ofpix_t max1, min1, max2, min2, max, min;

    getminmax(&min1, &max1, I1, size);
    getminmax(&min2, &max2, I2, size);

    //compute the global max and min
    max = max1;
    min = min1;

    if (max2 > max)
        max = max2;

    if (min2 < min)
        min = min2;

    const ofpix_t den = max - min;

    if(den > 0)
        //normalize the images between 0 and 255
#pragma omp parallel for
        for(int i = 0; i < size; i++)
        {
            I1n[i] = 255.0 * (I1[i] - min) / den;
            I2n[i] = 255.0 * (I2[i] - min) / den;
        }
    
    else
        //copy the original data
#pragma omp parallel for
        for(int i = 0; i < size; i++)
        {
            I1n[i] = I1[i];
            I2n[i] = I2[i];
        }
}

/**
 *
 * Function to normalize the images between 0 and 255
 *
 **/
void image_normalization_3(ofpix_t *I0, ofpix_t *I1, ofpix_t *I2, int size) {
	ofpix_t max0, min0, max1, min1, max2, min2, max, min;

	getminmax(&min0, &max0, I0, size);
	getminmax(&min1, &max1, I1, size);
	getminmax(&min2, &max2, I2, size);

	max = max0;
	min = min0;

	if (max1 > max)
		max = max1;

	if (min1 < min)
		min = min1;

	if (max2 > max)
		max = max2;

	if (min2 < min)
		min = min2;

	const ofpix_t den = max - min;

#pragma omp parallel for
	for (int i = 0; i < size; i++) {
		I0[i] = 255.0 * (I0[i] - min) / den;
		I1[i] = 255.0 * (I1[i] - min) / den;
		I2[i] = 255.0 * (I2[i] - min) / den;
	}
}

void image_normalization_4(const ofpix_t *I_1, // input image-1
		const ofpix_t *I0,  // input image0
		const ofpix_t *I1,  // input image1
		const ofpix_t *filtI0, //Smooth vefsion of I0
		ofpix_t *I_1n,		// normalized output image -1
		ofpix_t *I0n,       // normalized output image0
		ofpix_t *I1n,       // normalized output image1
		ofpix_t *filtI0n,   // normalized output image filtI0
		int size          // size of the image
		) {
	ofpix_t max_1, max0, max1, maxf0, min_1, min0, minf0, min1;

	// obtain the max and min of each image
	getminmax(&min_1, &max_1, I_1, size);
	getminmax(&min0, &max0, I0, size);
	getminmax(&min1, &max1, I1, size);
	getminmax(&minf0, &maxf0, filtI0, size);

	// obtain the max and min of images
	double max = (max_1 > max0) ? max_1 : max0;
	max = (max > max1) ? max : max1;
	max = (max > maxf0) ? max : maxf0;

	double min = (min_1 < min0) ? min_1 : min0;
	min = (min < min1) ? min : min1;
	min = (min < minf0) ? min : minf0;

	const double den = max - min;

	if (den > 0)
		// normalize both images
		for (int i = 0; i < size; i++) {
			I_1n[i] = 255.0 * (I_1[i] - min) / den;
			I0n[i] = 255.0 * (I0[i] - min) / den;
			I1n[i] = 255.0 * (I1[i] - min) / den;
			filtI0n[i] = 255.0 * (filtI0[i] - min) / den;
		}

	else
		// copy the original images
		for (int i = 0; i < size; i++) {
			I_1n[i] = I_1[i];
			I0n[i] = I0[i];
			I1n[i] = I1[i];
			filtI0n[i] = filtI0[i];
		}
}

/**
 *
 * Compute the max and min of an array
 *
 **/
void getminmax(ofpix_t *min,     // output min
		ofpix_t *max,     // output max
		const ofpix_t *x, // input array
		int n           // array size
		) {
	*min = *max = x[0];
	for (int i = 1; i < n; i++) {
		if (x[i] < *min)
			*min = x[i];
		if (x[i] > *max)
			*max = x[i];
	}
}
