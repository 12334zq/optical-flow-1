// This file is free software: you can use, modify and/or
// redistribute it under the terms of the GNU General Public
// License as published by the Free Software Foundation, either
// version 3 of the License, or (at your option) any later
// version. You should have received a copy of this license along
// this program. If not, see <http://www.gnu.org/licenses/>.

// Copyright (C) 2013 J. F. Garamendi <jf.garamendi@upf.edu>
// All rights reserved.

#include "TV_ROF_Box.h"

void Scalar_ROF_BoxCellCentered(double *u, const double *f, double *initialP1,
		double *initialP2, const double *g_function, const double lambda,
		const double omega, //Omega-relaxation
		const int nx, const int ny, const int nIter) {

	//Variable declaration
	const int size = nx * ny;
	const int nx_stg = nx * 2 + 1;
	const int ny_stg = ny * 2 + 1;
	const int size_stg = nx_stg * ny_stg;

	int i_v = 0; //i index for the vectorized staggered grid
	// Vectorized indices for the Staggered Grid.
	// m means "minus", p means "plus", so ip1_jm2_v means position (i+1,j-2)
	int im2_v = 0;
	int ip2_v = 0;
	int im1_jm2_v = 0;
	int ip1_jm2_v = 0;
	int im2_jp1_v = 0;
	int ip1_jp2_v = 0;
	int ip2_jp1_v = 0;
	int ip2_jm1_v = 0;
	int im1_jp2_v = 0;
	int im2_jm1_v = 0;
	int im3_v = 0;
	int ip3_v = 0;
	int jm3_v = 0;
	int jp3_v = 0;
	int im1_v = 0;
	int ip1_v = 0;
	int jm1_v = 0;
	int jp1_v = 0;
	double *stgGrid_F = NULL;
	double *stgGrid_P = NULL; //where it will be placed p11 and p12
	double *stgGrid_alfa = NULL; //where it will be placed |grad(u1)|
	double *ux = NULL;
	double *uy = NULL;
	//Beta is from the Garamendi's paper, don't from coloma's paper
	double *beta = NULL; //Vector for storing the beta values of the 4x4 general matrix
	double W = 0.0; //Variables W,N,S,E for the free term of the algebraic system of equations
	double N = 0.0;
	double S = 0.0;
	double E = 0.0;
	double den = 1.0;
	double a = 0.0;
	double b = 0.0;
	double alf = 0.0;
	double gam = 0.0;
	double x = 0.0;
	double y = 0.0;
	double c = 0.0;
	int i_ns = 0; //i index for the non-staggered image
	//END OF VARIABLE DECLARATION

	//if (stgGrid_F==NULL)
	{
		stgGrid_F = (double *) malloc(size_stg * sizeof(double));
	}
	//if (stgGrid_P==NULL)
	{
		stgGrid_P = (double *) malloc(size_stg * sizeof(double));
	}
	//if (stgGrid_alfa==NULL)
	{
		stgGrid_alfa = (double *) malloc(size_stg * sizeof(double));
	}

	//if (ux==NULL)
	{
		ux = (double *) malloc(size * sizeof(double));
	}
	//if (uy==NULL)
	{
		uy = (double *) malloc(size * sizeof(double));
	}

	//if (beta==NULL)
	{
		beta = (double *) malloc(4 * sizeof(double));

	}

	i_v = 0; //i index for the vectorized staggered grid

	//Initialitation for stgGrid_F, stgGrid_P and stgGrid_alfa
	for (int i = 0; i < size_stg; i++) {
		stgGrid_F[i] = 0.0;
		stgGrid_P[i] = 0.0;
		stgGrid_alfa[i] = 0.0;
	}

	//initialization of beta.
	for (int i = 0; i < 4; i++) {
		beta[i] = 0.0;
	}

	//Initialization of the center of the cells
	i_ns = 0; //First row
	i_v = 0;
	for (int i = 1; i <= ny_stg - 2; i = i + 2) {
		for (int j = 1; j <= nx_stg - 2; j = j + 2) {
			i_v = i * nx_stg + j;
			jp1_v = i * nx_stg + j + 1;
			ip1_v = (i + 1) * nx_stg + j;

			stgGrid_F[i_v] = f[i_ns];
			stgGrid_P[ip1_v] = initialP1[i_ns];
			stgGrid_P[jp1_v] = initialP2[i_ns];

			i_ns++;
		}
	}

	//Derivative of  of the inner x nodes on upper and bottom  sides of the cells.
	//for the derivative of f
	//It is the gradient of the center nodes and it is placed at the sides of the cell
	i_v = 0;
	im1_v = 0;
	ip1_v = 0;

	//(i,j) is the coordinate of node (NOT: the center of the cell),
	//placed at the upper and bottom side of the cell
	for (int i = 2; i <= ny_stg - 3; i = i + 2) {
		for (int j = 1; j <= nx_stg - 2; j = j + 2) {
			i_v = i * nx_stg + j;
			im1_v = (i - 1) * nx_stg + j;
			ip1_v = (i + 1) * nx_stg + j;

			stgGrid_F[i_v] = (stgGrid_F[ip1_v] - stgGrid_F[im1_v]);
		}
	}

	//Initialization of the inner o nodes on the sides of the cells. It is the gradient of the center nodes
	// placed at the sides of the cell
	i_v = 0;
	jm1_v = 0;
	jp1_v = 0;

	//(i,j) is the coordinate of node . (NOT: the center of the cell)
	for (int i = 1; i <= ny_stg - 2; i = i + 2) {
		for (int j = 2; j <= nx_stg - 3; j = j + 2) {
			i_v = i * nx_stg + j;
			jm1_v = i * nx_stg + (j - 1);
			jp1_v = i * nx_stg + (j + 1);

			stgGrid_F[i_v] = (stgGrid_F[jp1_v] - stgGrid_F[jm1_v]);
		}
	}
	//END of initialization

	for (int iter = 1; iter <= nIter; iter++) {

		// Compute the alfa values. This value corresponds to |grad u| / g

		i_v = 0;
		ip1_v = 0; //i vector index for the previous value in j direction
		jp1_v = 0; //i vector index for the next value in j direction

		i_ns = 0;

		//In these nested loop (i,j) is representing the center of the cell
		forward_gradient(u, ux, uy, nx, ny);

		for (int i = 1; i <= ny_stg - 2; i = i + 2) {
			for (int j = 1; j <= nx_stg - 2; j = j + 2) {
				i_v = i * nx_stg + j;
				ip1_v = (i + 1) * nx_stg + j;
				jp1_v = i * nx_stg + (j + 1);

				stgGrid_alfa[ip1_v] = hypot(ux[i_ns], uy[i_ns])
						/ (lambda * g_function[i_ns]);

				//Direct interpolation. We should try others kind of interpolation
				stgGrid_alfa[jp1_v] = stgGrid_alfa[ip1_v];

				i_ns++;
			}
		}

		//again, (i,j) corresponds to the center of cells
		////////////Boundary conditions//////////////
		// North-West corner
		int ii = 1;
		int jj = 1;

		ip1_v = (ii + 1) * nx_stg + jj;
		ip1_jp2_v = (ii + 1) * nx_stg + (jj + 2);
		ip2_jp1_v = (ii + 2) * nx_stg + (jj + 1);

		ip2_jm1_v = (ii + 2) * nx_stg + (jj - 1);
		im1_jp2_v = (ii - 1) * nx_stg + (jj + 2);
		jp3_v = ii * nx_stg + (jj + 3);
		ip3_v = (ii + 3) * nx_stg + jj;
		jp1_v = ii * nx_stg + (jj + 1);

		beta[0] = 0.0;
		beta[1] = 0.0;

		beta[2] = -2 - stgGrid_alfa[ip1_v];
		beta[3] = -2 - stgGrid_alfa[jp1_v];

		S = -stgGrid_P[ip3_v] - stgGrid_P[ip2_jp1_v] + stgGrid_P[ip2_jm1_v]
				- stgGrid_F[ip1_v];

		E = -stgGrid_P[jp3_v] - stgGrid_P[ip1_jp2_v] + stgGrid_P[im1_jp2_v]
				- stgGrid_F[jp1_v];

		den = beta[2] * beta[3] - 1;

		stgGrid_P[ip1_v] = (1 - omega) * stgGrid_P[ip1_v]
				+ omega * (S * beta[3] + E) / den;
		stgGrid_P[jp1_v] = (1 - omega) * stgGrid_P[jp1_v]
				+ omega * (E * beta[2] + S) / den;

		//North side. Remember, (i,j) is the center of the cell.
		beta[1] = 0;
		ii = 1;

		for (int j = 3; j <= nx_stg - 4; j = j + 2) {
			im1_jp2_v = (ii - 1) * nx_stg + (j + 2);
			ip1_v = (ii + 1) * nx_stg + j;
			ip1_jp2_v = (ii + 1) * nx_stg + (j + 2);
			ip2_jp1_v = (ii + 2) * nx_stg + (j + 1);
			ip2_jm1_v = (ii + 2) * nx_stg + (j - 1);
			ip1_jm2_v = (ii + 1) * nx_stg + (j - 2);
			im1_jm2_v = (ii - 1) * nx_stg + (j - 2);
			jp3_v = ii * nx_stg + (j + 3);
			ip3_v = (ii + 3) * nx_stg + j;
			jp1_v = ii * nx_stg + (j + 1);
			jm1_v = ii * nx_stg + (j - 1);
			jm3_v = ii * nx_stg + (j - 3);

			beta[0] = -2 - stgGrid_alfa[jm1_v];
			beta[2] = -2 - stgGrid_alfa[ip1_v];
			beta[3] = -2 - stgGrid_alfa[jp1_v];

			W = -stgGrid_F[jm1_v] - stgGrid_P[jm3_v] + stgGrid_P[ip1_jm2_v]
					- stgGrid_P[im1_jm2_v];
			S = -stgGrid_F[ip1_v] - stgGrid_P[ip3_v] - stgGrid_P[ip2_jp1_v]
					+ stgGrid_P[ip2_jm1_v];
			E = -stgGrid_F[jp1_v] - stgGrid_P[jp3_v] - stgGrid_P[ip1_jp2_v]
					+ stgGrid_P[im1_jp2_v];

			den = beta[0] * beta[2] * beta[3] - beta[0] - beta[2] - beta[3] - 2;

			stgGrid_P[jm1_v] = (1 - omega) * stgGrid_P[jm1_v]
					+ omega
							* (W * beta[2] * beta[3] - E * beta[2] - S * beta[3]
									- W - E - S) / den;

			stgGrid_P[ip1_v] = (1 - omega) * stgGrid_P[ip1_v]
					+ omega
							* (S * beta[0] * beta[3] - W * beta[3] + E * beta[0]
									- W + E - S) / den;

			stgGrid_P[jp1_v] = (1 - omega) * stgGrid_P[jp1_v]
					+ omega
							* (E * beta[0] * beta[2] - W * beta[2] + S * beta[0]
									- W - E + S) / den;

		}

		//East-North Corner
		ii = 1;
		jj = nx_stg - 2;

		ip1_v = (ii + 1) * nx_stg + jj;
		jm1_v = ii * nx_stg + (jj - 1);
		ip2_jp1_v = (ii + 2) * nx_stg + (jj + 1);
		ip1_jm2_v = (ii + 1) * nx_stg + (jj - 2);
		im1_jm2_v = (ii - 1) * nx_stg + (jj - 2);
		ip2_jm1_v = (ii + 2) * nx_stg + (jj - 1);
		ip3_v = (ii + 3) * nx_stg + jj;
		jm3_v = ii * nx_stg + (jj - 3);

		beta[1] = 0;
		beta[3] = 0;

		beta[0] = -2 - stgGrid_alfa[jm1_v];
		beta[2] = -2 - stgGrid_alfa[ip1_v];

		den = beta[0] * beta[2] - 1;

		W = -stgGrid_P[jm3_v] + stgGrid_P[ip1_jm2_v] - stgGrid_P[im1_jm2_v]
				- stgGrid_F[jm1_v];

		S = -stgGrid_P[ip3_v] - stgGrid_P[ip2_jp1_v] + stgGrid_P[ip2_jm1_v]
				- stgGrid_F[ip1_v];

		stgGrid_P[jm1_v] = (1 - omega) * stgGrid_P[jm1_v]
				+ omega * (W * beta[2] - S) / den;
		stgGrid_P[ip1_v] = (1 - omega) * stgGrid_P[ip1_v]
				+ omega * (S * beta[0] - W) / den;

		//inner cells
		for (int i = 3; i <= ny_stg - 4; i = i + 2) {
			//west side
			jj = 1;

			im1_v = (i - 1) * nx_stg + jj;
			im3_v = (i - 3) * nx_stg + jj;

			im1_jp2_v = (i - 1) * nx_stg + (jj + 2);
			im2_jp1_v = (i - 2) * nx_stg + (jj + 1);
			im2_jm1_v = (i - 2) * nx_stg + (jj - 1);

			ip1_v = (i + 1) * nx_stg + jj;
			ip1_jp2_v = (i + 1) * nx_stg + (jj + 2);
			ip2_jp1_v = (i + 2) * nx_stg + (jj + 1);

			ip2_jm1_v = (i + 2) * nx_stg + (jj - 1);
			jp3_v = i * nx_stg + (jj + 3);
			ip3_v = (i + 3) * nx_stg + jj;
			jp1_v = i * nx_stg + (jj + 1);

			beta[0] = 0;

			beta[1] = -2 - stgGrid_alfa[im1_v];
			beta[2] = -2 - stgGrid_alfa[ip1_v];
			beta[3] = -2 - stgGrid_alfa[jp1_v];

			den = beta[1] * beta[2] * beta[3] - (beta[1] + beta[2] + beta[3])
					- 2;

			N = -stgGrid_P[im3_v] + stgGrid_P[im2_jp1_v] - stgGrid_P[im2_jm1_v]
					- stgGrid_F[im1_v];
			S = -stgGrid_P[ip3_v] - stgGrid_P[ip2_jp1_v] + stgGrid_P[ip2_jm1_v]
					- stgGrid_F[ip1_v];

			E = -stgGrid_P[jp3_v] - stgGrid_P[ip1_jp2_v] + stgGrid_P[im1_jp2_v]
					- stgGrid_F[jp1_v];

			stgGrid_P[im1_v] = (1 - omega) * stgGrid_P[im1_v]
					+ omega
							* (beta[2] * beta[3] * N - E * beta[2] - S * beta[3]
									- N - S - E) / den;
			stgGrid_P[ip1_v] = (1 - omega) * stgGrid_P[ip1_v]
					+ omega
							* (beta[1] * beta[3] * S + E * beta[1] - N * beta[3]
									- N - S + E) / den;
			stgGrid_P[jp1_v] = (1 - omega) * stgGrid_P[jp1_v]
					+ omega
							* (beta[1] * beta[2] * E - N * beta[2] + S * beta[1]
									- N + S - E) / den;

			// %%%%%%%%%%%%%%%% INNER POINT %%%%%%%%%%%%%%%%%%%%

			for (int j = 3; j <= nx_stg - 4; j = j + 2) {
				i_v = i * nx_stg + j;
				im1_v = i_v - nx_stg; 					// (i-1) * nx_stg + j;
				im2_v = im1_v - nx_stg;					// (i-2) * nx_stg + j
				im3_v = im2_v - nx_stg; 				// (i-3) * nx_stg + j;
				im1_jp2_v = im1_v + 2;				// (i-1) * nx_stg + (j+2);
				im2_jp1_v = im2_v + 1;				// (i-2) * nx_stg + (j+1);
				im2_jm1_v = im2_v - 1; 				// (i-2) * nx_stg + (j-1);
				ip1_v = i_v + nx_stg;					// (i+1) * nx_stg + j;
				ip2_v = ip1_v + nx_stg;
				ip3_v = ip2_v + nx_stg;					// (i+3) * nx_stg + j;
				ip1_jp2_v = ip1_v + 2;				// (i+1) * nx_stg + (j+2);
				ip2_jp1_v = ip2_v + 1;				// (i+2) * nx_stg + (j+1);
				ip2_jm1_v = ip2_v - 1;				// (i+2) * nx_stg + (j-1);
				jp3_v = i_v + 3; 						// i * nx_stg + (j+3);
				jp1_v = i_v + 1;						// i * nx_stg + (j+1);
				jm1_v = i_v - 1;						// i * nx_stg + (j-1);
				jm3_v = i_v - 3; 						// i * nx_stg + (j-3);
				ip1_jm2_v = ip1_v - 2;				// (i+1) * nx_stg + (j-2);
				im1_jm2_v = im1_v - 2;				// (i-1) * nx_stg + (j-2);

				//P1
				beta[0] = -2 - stgGrid_alfa[jm1_v];
				beta[1] = -2 - stgGrid_alfa[im1_v];
				beta[2] = -2 - stgGrid_alfa[ip1_v];
				beta[3] = -2 - stgGrid_alfa[jp1_v];

				//Cramer's Rule
				den = beta[0]
						* (beta[1] * beta[2] * beta[3] - beta[1] - beta[2]
								- beta[3] - 2)
						+ (-beta[2] * beta[3] - beta[2] - beta[3] - 1)
						+ (-beta[1] * beta[3] - beta[1] - beta[3] - 1)
						- (beta[1] * beta[2] + beta[1] + beta[2] + 1);

				W = -stgGrid_P[jm3_v] + stgGrid_P[ip1_jm2_v]
						- stgGrid_P[im1_jm2_v] - stgGrid_F[jm1_v];
				N = -stgGrid_P[im3_v] + stgGrid_P[im2_jp1_v]
						- stgGrid_P[im2_jm1_v] - stgGrid_F[im1_v];
				S = -stgGrid_P[ip3_v] - stgGrid_P[ip2_jp1_v]
						+ stgGrid_P[ip2_jm1_v] - stgGrid_F[ip1_v];
				E = -stgGrid_P[jp3_v] - stgGrid_P[ip1_jp2_v]
						+ stgGrid_P[im1_jp2_v] - stgGrid_F[jp1_v];

				/*
				 // Cramer's Rule
				 stgGrid_P[jm1_v] = (1-omega)*stgGrid_P[jm1_v]+omega*(W*(beta[1]*beta[2]*beta[3]-beta[1] -beta[2] -beta[3] -2)
				 +(N*beta[2]*beta[3] -E -S -E*beta[2] -N -S*beta[3])
				 +(N*beta[3] - E*beta[1] +S -E +N -S*beta[1]*beta[3])
				 -( -N + E*beta[1]*beta[2] +S -E -N*beta[2] + S*beta[1]))
				 / den;
				 stgGrid_P[im1_v] = (1-omega)*stgGrid_P[im1_v]+omega*(beta[0]*(N*beta[2]*beta[3] -E -S -E*beta[2] -N -S*beta[3])
				 - W*(-beta[2]*beta[3] -beta[2] -beta[3] -1)
				 + (-S*beta[3] +E -N -S -N*beta[3] -E)
				 - (S + N*beta[2] +E -S +N +E*beta[2]))
				 / den;
				 stgGrid_P1[ip1_v] = (1-omega)*stgGrid_P1[ip1_v]+omega*(beta[0]*(S*beta[1]*beta[3] -N +E -S +E*beta[1] -N*beta[3])
				 +(-S*beta[3]-N-S -N*beta[3])
				 +W*(-beta[3]-beta[1]-1-beta[1]*beta[3])
				 -(-E +S*beta[1] +S -E*beta[1]))
				 / den;
				 stgGrid_P1[jp1_v] = (1-omega)*stgGrid_P1[jp1_v]+omega*(beta[0]*(E*beta[1]*beta[2] +S -N -N*beta[2] +S*beta[1] -E)
				 +(-E*beta[2]+S-N-N*beta[2] -S -E)
				 +(-E +S*beta[1] +S -E*beta[1])
				 -W*(beta[1]*beta[2] +1 +beta[2] +beta[1]))
				 /den;
				 */

				// Gauss Elimination
				a = 1 / beta[0];
				b = -(beta[0] + 1) / (beta[0] * beta[1] - 1);
				alf = 1 + a;
				gam = -a + b * alf;
				x = N + a * W;
				y = -a * W + b * x;
				c = (1 - gam) / (beta[2] + gam);

				stgGrid_P[jp1_v] = (1 - omega) * stgGrid_P[jp1_v]
						+ omega * (E + y + c * (S + y))
								/ (beta[3] + gam + c * (gam - 1));
				stgGrid_P[ip1_v] = (1 - omega) * stgGrid_P[ip1_v]
						+ omega * (S + y + stgGrid_P[jp1_v] * (1 - gam))
								/ (beta[2] + gam);
				stgGrid_P[im1_v] = (1 - omega) * stgGrid_P[im1_v]
						+ omega
								* (x
										- alf
												* (stgGrid_P[jp1_v]
														+ stgGrid_P[ip1_v]))
								/ (beta[1] - a);
				stgGrid_P[jm1_v] = (1 - omega) * stgGrid_P[jm1_v]
						+ omega
								* (W + stgGrid_P[im1_v] - stgGrid_P[ip1_v]
										- stgGrid_P[jp1_v]) / (beta[0]);

			}

			//%%%%%%%%%%%%%%%%% END OF INNER CELLS%%%%%%%%%%%%%%

			//EAST SIDE

			jj = nx_stg - 2;
			im1_v = (i - 1) * nx_stg + jj;
			im3_v = (i - 3) * nx_stg + jj;
			im2_jp1_v = (i - 2) * nx_stg + (jj + 1);
			im2_jm1_v = (i - 2) * nx_stg + (jj - 1);
			ip1_v = (i + 1) * nx_stg + jj;
			ip2_jp1_v = (i + 2) * nx_stg + (jj + 1);
			ip2_jm1_v = (i + 2) * nx_stg + (jj - 1);
			ip3_v = (i + 3) * nx_stg + jj;
			jm1_v = i * nx_stg + (jj - 1);
			jm3_v = i * nx_stg + (jj - 3);
			ip1_jm2_v = (i + 1) * nx_stg + (jj - 2);
			im1_jm2_v = (i - 1) * nx_stg + (jj - 2);

			beta[0] = -2 - stgGrid_alfa[jm1_v];
			beta[1] = -2 - stgGrid_alfa[im1_v];
			beta[2] = -2 - stgGrid_alfa[ip1_v];
			beta[3] = 0;

			den = (beta[0] * beta[1] * beta[2])
					+ (-beta[0] - beta[1] - beta[2] - 2);

			W = -stgGrid_P[jm3_v] + stgGrid_P[ip1_jm2_v] - stgGrid_P[im1_jm2_v]
					- stgGrid_F[jm1_v];
			N = -stgGrid_P[im3_v] + stgGrid_P[im2_jp1_v] - stgGrid_P[im2_jm1_v]
					- stgGrid_F[im1_v];
			S = -stgGrid_P[ip3_v] - stgGrid_P[ip2_jp1_v] + stgGrid_P[ip2_jm1_v]
					- stgGrid_F[ip1_v];

			stgGrid_P[jm1_v] = (1 - omega) * stgGrid_P[jm1_v]
					+ omega
							* (W * beta[1] * beta[2] - S + N - S * beta[1] - W
									+ N * beta[2]) / den;
			stgGrid_P[im1_v] = (1 - omega) * stgGrid_P[im1_v]
					+ omega
							* (N * beta[0] * beta[2] + W - S - N - S * beta[0]
									+ W * beta[2]) / den;
			stgGrid_P[ip1_v] = (1 - omega) * stgGrid_P[ip1_v]
					+ omega
							* (S * beta[0] * beta[1] - N - W - W * beta[1]
									- N * beta[0] - S) / den;

		}

		//South-West Corner
		ii = ny_stg - 2;
		jj = 1;

		im1_v = (ii - 1) * nx_stg + jj;
		im3_v = (ii - 3) * nx_stg + jj;
		im1_jp2_v = (ii - 1) * nx_stg + (jj + 2);
		im2_jp1_v = (ii - 2) * nx_stg + (jj + 1);
		im2_jm1_v = (ii - 2) * nx_stg + (jj - 1);
		ip1_jp2_v = (ii + 1) * nx_stg + (jj + 2);
		jp3_v = ii * nx_stg + (jj + 3);
		jp1_v = ii * nx_stg + (jj + 1);

		beta[0] = 0;
		beta[1] = -2 - stgGrid_alfa[im1_v];
		beta[2] = 0;
		beta[3] = -2 - stgGrid_alfa[jp1_v];

		N = -stgGrid_P[im3_v] + stgGrid_P[im2_jp1_v] - stgGrid_P[im2_jm1_v]
				- stgGrid_F[im1_v];

		E = -stgGrid_P[jp3_v] - stgGrid_P[ip1_jp2_v] + stgGrid_P[im1_jp2_v]
				- stgGrid_F[jp1_v];

		den = beta[3] * beta[1] - 1;

		stgGrid_P[im1_v] = (1 - omega) * stgGrid_P[im1_v]
				+ omega * (beta[3] * N - E) / den;
		stgGrid_P[jp1_v] = (1 - omega) * stgGrid_P[jp1_v]
				+ omega * (beta[1] * E - N) / den;

		//South Side
		ii = ny_stg - 2;
		beta[2] = 0;
		for (int j = 3; j <= nx_stg - 4; j = j + 2) {
			im1_v = (ii - 1) * nx_stg + j;
			im3_v = (ii - 3) * nx_stg + j;
			im1_jp2_v = (ii - 1) * nx_stg + (j + 2);
			im2_jp1_v = (ii - 2) * nx_stg + (j + 1);
			im2_jm1_v = (ii - 2) * nx_stg + (j - 1);
			ip1_jp2_v = (ii + 1) * nx_stg + (j + 2);
			jp3_v = ii * nx_stg + (j + 3);
			jp1_v = ii * nx_stg + (j + 1);
			jm1_v = ii * nx_stg + (j - 1);
			jm3_v = ii * nx_stg + (j - 3);
			ip1_jm2_v = (ii + 1) * nx_stg + (j - 2);
			im1_jm2_v = (ii - 1) * nx_stg + (j - 2);

			beta[0] = -2 - stgGrid_alfa[jm1_v];
			beta[1] = -2 - stgGrid_alfa[im1_v];
			beta[3] = -2 - stgGrid_alfa[jp1_v];

			den = beta[0] * beta[1] * beta[3] - beta[0] - beta[1] - beta[3] - 2;

			W = -stgGrid_P[jm3_v] + stgGrid_P[ip1_jm2_v] - stgGrid_P[im1_jm2_v]
					- stgGrid_F[jm1_v];
			N = -stgGrid_P[im3_v] + stgGrid_P[im2_jp1_v] - stgGrid_P[im2_jm1_v]
					- stgGrid_F[im1_v];
			E = -stgGrid_P[jp3_v] - stgGrid_P[ip1_jp2_v] + stgGrid_P[im1_jp2_v]
					- stgGrid_F[jp1_v];

			stgGrid_P[jm1_v] = (1 - omega) * stgGrid_P[jm1_v]
					+ omega
							* (W * beta[1] * beta[3] - E + N - E * beta[1] - W
									+ N * beta[3]) / den;
			stgGrid_P[im1_v] = (1 - omega) * stgGrid_P[im1_v]
					+ omega
							* (N * beta[0] * beta[3] + W - E - N - E * beta[0]
									+ W * beta[3]) / den;
			stgGrid_P[jp1_v] = (1 - omega) * stgGrid_P[jp1_v]
					+ omega
							* (E * beta[0] * beta[1] - N - W - W * beta[1]
									- N * beta[0] - E) / den;

		}

		//SOUTH-EAST CORNER
		ii = ny_stg - 2;
		jj = nx_stg - 2;

		im1_v = (ii - 1) * nx_stg + jj;
		im3_v = (ii - 3) * nx_stg + jj;
		im2_jp1_v = (ii - 2) * nx_stg + (jj + 1);
		im2_jm1_v = (ii - 2) * nx_stg + (jj - 1);
		jm1_v = ii * nx_stg + (jj - 1);
		jm3_v = ii * nx_stg + (jj - 3);
		ip1_jm2_v = (ii + 1) * nx_stg + (jj - 2);
		im1_jm2_v = (ii - 1) * nx_stg + (jj - 2);

		beta[2] = 0;
		beta[3] = 0;

		beta[0] = -2 - stgGrid_alfa[jm1_v];
		beta[1] = -2 - stgGrid_alfa[im1_v];

		den = beta[0] * beta[1] - 1;

		W = -stgGrid_P[jm3_v] + stgGrid_P[ip1_jm2_v] - stgGrid_P[im1_jm2_v]
				- stgGrid_F[jm1_v];
		N = -stgGrid_P[im3_v] + stgGrid_P[im2_jp1_v] - stgGrid_P[im2_jm1_v]
				- stgGrid_F[im1_v];

		stgGrid_P[jm1_v] = (1 - omega) * stgGrid_P[jm1_v]
				+ omega * (W * beta[1] + N) / den;
		stgGrid_P[im1_v] = (1 - omega) * stgGrid_P[im1_v]
				+ omega * (N * beta[0] + W) / den;

		// Computing u
		i_v = 0;
		//error = 0.0;
		const double uk = 0.0;
		for (int i = 1; i <= ny_stg - 2; i = i + 2)
			for (int j = 1; j <= nx_stg - 2; j = j + 2) {
				im1_v = (i - 1) * nx_stg + j;
				ip1_v = (i + 1) * nx_stg + j;

				jm1_v = i * nx_stg + (j - 1);
				jp1_v = i * nx_stg + (j + 1);

				//uk = u[i_v];

				u[i_v] = lambda * f[i_v]
						+ lambda
								* (stgGrid_P[ip1_v] - stgGrid_P[im1_v]
										+ stgGrid_P[jp1_v] - stgGrid_P[jm1_v]);

				//error += (u[i_v]-uk)*(u[i_v]-uk);

				initialP1[i_v] = stgGrid_P[ip1_v];
				initialP2[i_v] = stgGrid_P[jp1_v];

				i_v++;

			}
		//error /= size;

	}					//for iter
	free(stgGrid_P);
	free(stgGrid_F);
	free(stgGrid_alfa);
	free(uy);
	free(ux);
	free(beta);
}

