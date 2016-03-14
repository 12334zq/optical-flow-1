// This program is free software: you can use, modify and/or redistribute it
// under the terms of the simplified BSD License. You should have received a
// copy of this license along this program. If not, see
// <http://www.opensource.org/licenses/bsd-license.html>.
//
// Copyright (C) 2011, Javier Sˆnchez PŽrez <jsanchez@dis.ulpgc.es>
// Copyright (C) 2012, Coloma Ballester <coloma.ballester@upf.edu>
// Copyright (C) 2013-2014 J. F. Garamendi <jf.garamendi@upf.edu>

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

#include "operators.h"
#include "solvers.h"
#include "bicubic_interpolation.h"
#include "pyramid.h"
#include "warp.h"
#include "utils.h"
#include "tvl1OccFlow_lib.h"
#include "xmalloc.h"

#define malloc(m) xmalloc(m)

/**
 * Implementation of the Ballester, Garrido, Lazcano and Caselles  TV-L1 optic flow method with Occlussions
 *
 * see reference:
 *  [1] C. Ballester, L. Garrido, V. Lazcano, V. Caselles.
 *      "A TV-L1 Optical Flow Method with Occlusion Detection".
 *      Lecture Notes in Computer Science. 7476:31-40, 2012
 **/

/////////////////////////////////////////////////////////////////////////////////////////////////////
// Declaration of private (to this file) functions. (Public functions are declarated in .h header file)
/*
 * Function name: Threshold
 * Author:
 * Parameters:
 * 		vector         : (in/out) Pointer to the vector to threshold
 * 		thresh_vector  : (in) Threshold
 * 		size:          : (in) Size of the vector
 * Output:
 * 		NONE
 * Description: Given a vector, performs a thresholding of every component of the vector.
 */
inline void Threshold(double *vector, const double threshold, const int size);

/*
 * Function name: choosed_g
 * Author:
 * Parameters:
 * 		g_fun   : (out) Vectorized matrix corresponding to apply the function number "choice" to I
 * 		choice  : (in) Scalar integer corresponding to function to apply:
 * 		               1: Identity
 * 		               2: (1+g_factor|grad(I)|)^(-1)
 * 		*I:     : (in) Vectorized matrix corresponding to the image
 * 		g_factor: (in) Scalar. g_factor value
 * 		nx      : (in) Scalar Integer. image width
 *		ny      : (in) Scalar Integer. image height
 * Output:
 * 		VOID
 * Description: Depending on choice, Apply a function to the image I.
 *
 */
inline void choosed_g(double *g_fun, const int choice, double *I,
		const double g_factor, const int nx, const int ny);

inline double L2error(double *u1prev, double *u2prev, const double *u1,
		const double *u2, const int size);

// End of declaration
/////////////////////////////////////////////////////////////////////////////////////////////////////

inline double L2error(double *u1prev, double *u2prev, const double *u1,
		const double *u2, const int size) {
	double error = 0.0;

	for (int i = 0; i < size; i++) {
		error += (u1[i] - u1prev[i]) * (u1[i] - u1prev[i])
				+ (u2[i] - u2prev[i]) * (u2[i] - u2prev[i]);

		u1prev[i] = u1[i];
		u2prev[i] = u2[i];
	}
	error /= size;

	return (error);
}

inline void Threshold(double *vector, const double threshold, const int size) {
#pragma omp parallel for
	for (int j = 0; j < size; j++) {
		vector[j] = (vector[j] > threshold);
	}

}

inline void choosed_g(double *g_fun, const int choice, double *I,
		const double g_factor, const int nx, const int ny) {
	const int size = nx * ny;
	double *Ix = (double *) malloc(sizeof(double) * size);
	double *Iy = (double *) malloc(sizeof(double) * size);

	switch (choice) {
	case 1: //Identity

		for (int i = 0; i < size; i++)
			g_fun[i] = 1.;
		break;

	case 2:

		centered_gradient(I, Ix, Iy, nx, ny);
#pragma omp parallel for
		for (int i = 0; i < size; i++) {
			double gggrad = sqrt(Ix[i] * Ix[i] + Iy[i] * Iy[i]);
			double aux = 1. + g_factor * gggrad;
			g_fun[i] = 1. / aux;
		}

		break;
	}
}

/**
 *
 * Function to compute the optical flow in one scale
 *
 **/
void Dual_TVL1_optic_flow(double *I_1,         // Previous frame to source Image
		double *I0,           // source image
		double *I1,           // target image
		double *filtI0,		  //Image used for computing the g function
		double *u1,           // x component of the optical flow
		double *u2,           // y component of the optical flow
		double *chi,          // Occlusion map
		const int nx,      // image width
		const int ny,      // image height
		const double lambda,  // weight parameter for the data term
		const double alpha,   // weight paramenter for chi|v|^2
		const double beta,    // weight parameter for chi*div(u)
		const double theta,   // weight parameter for (u - v)Â²
		const int warps,   // number of warpings per scale
		const double epsilon, // tolerance for numerical convergence
		const bool verbose  // enable/disable the verbose mode
		) {
	const int size = nx * ny;

	size_t sf = sizeof(double);
	double *I1x = (double *) malloc(size * sf);
	double *I1y = (double *) malloc(size * sf);
	double *I1w = (double *) malloc(size * sf);
	double *I1wx = (double *) malloc(size * sf);
	double *I1wy = (double *) malloc(size * sf);
	double *I_1x = (double *) malloc(size * sf);
	double *I_1y = (double *) malloc(size * sf);
	double *I_1w = (double *) malloc(size * sf);
	double *I_1wx = (double *) malloc(size * sf);
	double *I_1wy = (double *) malloc(size * sf);
	double *rho_c = (double *) malloc(size * sf);
	double *rho1_c = (double *) malloc(size * sf);
	double *rho3_c = (double *) malloc(size * sf);
	double *v1 = (double *) malloc(size * sf);
	double *v2 = (double *) malloc(size * sf);
	double *v11 = (double *) malloc(size * sf);
	double *v12 = (double *) malloc(size * sf);
	double *v31 = (double *) malloc(size * sf);
	double *v32 = (double *) malloc(size * sf);
	double *gp1 = (double *) malloc(size * sf);
	double *gp2 = (double *) malloc(size * sf);
	double *div = (double *) malloc(size * sf);
	double *grad = (double *) malloc(size * sf);
	double *grad1 = (double *) malloc(size * sf);
	double *grad3 = (double *) malloc(size * sf);
	double *g_function = (double *) malloc(size * sf);
	double *u1prev = (double *) malloc(size * sf);
	double *u2prev = (double *) malloc(size * sf);

	if (verbose)
		printf("verbose\n");

	//store g_function, used in the regularizing term (for the optical flow and for the occlusion map chi)
	choosed_g(g_function, G_CHOICE, filtI0, G_FACTOR, nx, ny);

	centered_gradient(I1, I1x, I1y, nx, ny);
	centered_gradient(I_1, I_1x, I_1y, nx, ny);

	// initialization of variables
	for (int i = 0; i < size; i++) {
		u1prev[i] = u1[i];
		u2prev[i] = u2[i];

		v1[i] = v2[i] = 0.0;

		v11[i] = v12[i] = 0.0;
		v31[i] = v32[i] = 0.0;

		gp1[i] = gp2[i] = 0.0;

		grad1[i] = grad3[i] = 0.0;
	}

	for (int warpings = 0; warpings < warps; warpings++) {

		//***********WARPING COMPUTATION
		// compute the warping of the target image and its derivatives
		//less final error  than bicubic_interpolation warp

		// I1(x+u)
		warping(I1, u1, u2, I1w, nx, ny);
		warping(I1x, u1, u2, I1wx, nx, ny);
		warping(I1y, u1, u2, I1wy, nx, ny);

		//Im1(x-u)
#pragma omp parallel for
		for (int i = 0; i < size; i++) {
			gp1[i] = -u1[i];
			gp2[i] = -u2[i];
		}

		warping(I_1, gp1, gp2, I_1w, nx, ny);
		warping(I_1x, gp1, gp2, I_1wx, nx, ny);
		warping(I_1y, gp1, gp2, I_1wy, nx, ny);

		//*********** END OF WARPING COMPUTATION

#pragma omp parallel for
		for (int i = 0; i < size; i++) {
			double Ix2 = 0.0;
			double Iy2 = 0.0;

			// store the |Grad(I1)|^2
			Ix2 = I1wx[i] * I1wx[i];
			Iy2 = I1wy[i] * I1wy[i];
			grad1[i] = (Ix2 + Iy2);

			// store the |Grad(I_{-1})|^2
			Ix2 = I_1wx[i] * I_1wx[i];
			Iy2 = I_1wy[i] * I_1wy[i];

			grad3[i] = (Ix2 + Iy2);

			// compute the constant part of the rho function
			rho1_c[i] = (I1w[i] - I1wx[i] * u1[i] - I1wy[i] * u2[i] - I0[i]); //rho1
			rho3_c[i] = (I_1w[i] + I_1wx[i] * u1[i] + I_1wy[i] * u2[i] - I0[i]); //rho_{-1}
		} //for i

		// Iterative Process
		int n = 0;
		double error = INFINITY;
		n = epsilon - epsilon;
		n = 0;
		//while (error > epsilon  && n < EXT_MAX_ITERATIONS) //comentado para pruebas, porque para probar hay que tener un numero fijo de iteraciones
		while (n < EXT_MAX_ITERATIONS) {
			n++;

			//Relaxation steps

			//1.- Minimization with respect to (wrt) v, the auxiliary variable
			Solver_wrt_v(u1, u2, v1, v2, chi, I1wx, I1wy, I_1wx, I_1wy, rho1_c,
					rho3_c, v11, v12, v31, v32, grad1, grad3, alpha, theta,
					lambda, nx, ny);

			//2.- Minimization with respect to (wrt) u, the optical flow map
			Solver_wrt_u(u1, u2, v1, v2, chi, g_function, theta, beta, nx, ny);
			// This median filtering improves the results (See "An Improved algorithm for TV-L1 Optical Flow"), but I don't like it.
			me_median_filtering(u1, nx, ny, 3);
			me_median_filtering(u2, nx, ny, 3);

			//3.- Minimization with respect to (wrt) chi, the occlusion map
			Solver_wrt_chi(u1, u2, chi, I1wx, I1wy, I_1wx, I_1wy, rho1_c,
					rho3_c, v11, v12, v31, v32, g_function, lambda, theta,
					alpha, beta, TAU_CHI, TAU_ETA, nx, ny);

			// difference between two consecutive iterations
			error = L2error(u1prev, u2prev, u1, u2, size);

		} //While error and ext_max_iterations

		if (verbose) {

			fprintf(stderr, "Warping: %d, "
					"Iterations: %d, "
					"Error: %e\n", warpings, n, error);

		} //if verbose
	} //for warpings

	// delete allocated memory
	free(I1x);
	free(I1y);
	free(I1w);
	free(I1wx);
	free(I1wy);
	free(I_1x);
	free(I_1y);
	free(I_1w);
	free(I_1wx);
	free(I_1wy);
	free(rho_c);
	free(rho1_c);
	free(rho3_c);
	free(v1);
	free(v2);
	free(v11);
	free(v12);
	free(v31);
	free(v32);
	free(gp1);
	free(gp2);
	free(div);
	free(grad);
	free(grad1);
	free(grad3);
	free(g_function);
	free(u1prev);
	free(u2prev);

}

/**
 *
 * Function to compute the optical flow using multiple scales
 *
 **/
void Dual_TVL1_optic_flow_multiscale(double *I_1, // Previous frame to source image
		double *I0,           // source image
		double *I1,           // target image
		double *filtI0,   // Smoothed version of I0 for computing the g function
		double *u1,           // x component of the optical flow
		double *u2,           // y component of the optical flow
		double *chi, 		 // Occlusion map
		const int nxx,     // image width
		const int nyy,     // image height
		const double lambda,  // weight parameter for the data term
		const double alpha,   // weight paramenter for chi|v|^2
		const double beta,    // weight parameter for chi*div(u)
		const double theta,   // weight parameter for (u - v)Â²
		const int nscales, // number of scales
		const double zfactor, // factor for building the image piramid
		const int warps,   // number of warpings per scale
		const double epsilon, // tolerance for numerical convergence
		const bool verbose  // enable/disable the verbose mode
		) {
	int size = nxx * nyy;

	// allocate memory for the pyramid structure
	double **I_1s = (double **) malloc(nscales * sizeof(double*));
	double **I0s = (double **) malloc(nscales * sizeof(double*));
	double **I1s = (double **) malloc(nscales * sizeof(double*));
	double **filtI0s = (double **) malloc(nscales * sizeof(double*));
	double **u1s = (double **) malloc(nscales * sizeof(double*));
	double **u2s = (double **) malloc(nscales * sizeof(double*));
	double **chis = (double **) malloc(nscales * sizeof(double*));
	int *nx = (int *) malloc(nscales * sizeof(int));
	int *ny = (int *) malloc(nscales * sizeof(int));

	I_1s[0] = (double *) malloc(size * sizeof(double));
	I0s[0] = (double *) malloc(size * sizeof(double));
	I1s[0] = (double *) malloc(size * sizeof(double));
	filtI0s[0] = (double *) malloc(size * sizeof(double));

	u1s[0] = u1;
	u2s[0] = u2;

	chis[0] = chi;

	nx[0] = nxx;
	ny[0] = nyy;

	// normalize the images between 0 and 255
	image_normalization_4(I_1, I0, I1, filtI0, I_1s[0], I0s[0], I1s[0],
			filtI0s[0], size);

	//Initialize the images, the flow and the occlussion map at the finest scale
	for (int i = 0; i < size; i++) {
		I_1s[0][i] = I_1[i];
		I0s[0][i] = I0[i];
		I1s[0][i] = I1[i];
		filtI0s[0][i] = filtI0[i];

		u1s[0][i] = 0.0;
		u2s[0][i] = 0.0;

		chis[0][i] = 0.0;
	}

	// pre-smooth the original images
	gaussian(I_1s[0], nx[0], ny[0], PRESMOOTHING_SIGMA);
	gaussian(I0s[0], nx[0], ny[0], PRESMOOTHING_SIGMA);
	gaussian(I1s[0], nx[0], ny[0], PRESMOOTHING_SIGMA);
	gaussian(filtI0s[0], nx[0], ny[0], PRESMOOTHING_SIGMA);

	// create the scales. 0 --> The finest
	for (int s = 1; s < nscales; s++) {
		zoom_size(nx[s - 1], ny[s - 1], &nx[s], &ny[s], zfactor);
		const int sizes = nx[s] * ny[s];

		// allocate memory
		I_1s[s] = (double *) malloc(sizes * sizeof(double));
		I0s[s] = (double *) malloc(sizes * sizeof(double));
		I1s[s] = (double *) malloc(sizes * sizeof(double));
		filtI0s[s] = (double *) malloc(sizes * sizeof(double));
		u1s[s] = (double *) malloc(sizes * sizeof(double));
		u2s[s] = (double *) malloc(sizes * sizeof(double));
		chis[s] = (double *) malloc(sizes * sizeof(double));

		if (s == nscales - 1) {
			//initialize the flow and the occlussion map at the coarsest scale
			for (int i = 0; i < sizes; i++) {
				u1s[s][i] = u2s[s][i] = chis[s][i] = 0.0;

			}
		}

		//Down-Sampling the images to create the pyramidal structure
		DownSampling(I_1s[s - 1], I_1s[s], nx[s - 1], ny[s - 1], zfactor);
		DownSampling(I0s[s - 1], I0s[s], nx[s - 1], ny[s - 1], zfactor);
		DownSampling(filtI0s[s - 1], filtI0s[s], nx[s - 1], ny[s - 1], zfactor);
		DownSampling(I1s[s - 1], I1s[s], nx[s - 1], ny[s - 1], zfactor);

	}

	// pyramidal structure for computing the optical flow
	for (int s = nscales - 1; s >= 0; s--) {

		// compute the optical flow at the current scale
		Dual_TVL1_optic_flow(I_1s[s], I0s[s], I1s[s],      // frames
				filtI0s[s],                   // Filtered Image (for g function)
				u1s[s], u2s[s],               // Optical flow map
				chis[s],					  // Occlusion map
				nx[s], ny[s],			      // Size of the images
				lambda, alpha, beta, theta,  // Energy functional parameters
				warps, epsilon, verbose   // Scheme parameters
				);

		// if it is not  the last scale,  upsample the optical flow
		if (s) {
			// zoom the optical flow for the next finer scale
			UpSampling(u1s[s], u1s[s - 1], nx[s], ny[s], nx[s - 1], ny[s - 1]);
			UpSampling(u2s[s], u2s[s - 1], nx[s], ny[s], nx[s - 1], ny[s - 1]);
			UpSampling(chis[s], chis[s - 1], nx[s], ny[s], nx[s - 1],
					ny[s - 1]);

			// scale the optical flow with the appropriate zoom factor
			for (int i = 0; i < nx[s - 1] * ny[s - 1]; i++) {
				u1s[s - 1][i] *= (double) 1.0 / zfactor;
				u2s[s - 1][i] *= (double) 1.0 / zfactor;
			}
		} else //Otherwise, threshold the occlusion map and finish
		{
			Threshold(chis[s], THR_CHI, nx[s] * ny[s]);
		}
	}
	// delete allocated memory
	for (int i = 1; i < nscales; i++) {
		free(I_1s[i]);
		free(I0s[i]);
		free(I1s[i]);
		free(filtI0s[i]);
		free(u1s[i]);
		free(u2s[i]);
		free(chis[i]);
	}
	free(nx);
	free(ny);
}
