// This program is free software: you can use, modify and/or redistribute it
// under the terms of the simplified BSD License. You should have received a
// copy of this license along this program. If not, see
// <http://www.opensource.org/licenses/bsd-license.html>.
//
// Copyright (C) 2012, Javier Sánchez Pérez <jsanchez@dis.ulpgc.es>
// All rights reserved.

#include "bicubic_interpolation.h"

#include <stdlib.h>
#include <math.h>

#define BOUNDARY_CONDITION 0
//0 Neumann
//1 Periodic
//2 Symmetric

/**
  *
  * Neumann boundary condition test
  *
**/
static
int
neumann_bc (int x, int nx, bool *out)
{
  if (x < 0)
    {
      x = 0;
      *out = true;
    }
  else if (x >= nx)
    {
      x = nx - 1;
      *out = true;
    }

  return x;
}

/**
  *
  * Periodic boundary condition test
  *
**/
static
int
periodic_bc (int x, int nx, bool *out)
{
  if (x < 0)
    {
      const int n = 1 - (int) (x / (nx + 1));
      const int ixx = x + n * nx;

      x = ixx % nx;
      *out = true;
    }
  else if (x >= nx)
    {
      x = x % nx;
      *out = true;
    }

  return x;
}


/**
  *
  * Symmetric boundary condition test
  *
**/
static
int
symmetric_bc (int x, int nx, bool *out)
{
  if (x < 0)
    {
      const int border = nx - 1;
      const int xx = -x;
      const int n = (int) (xx / border) % 2;

      if (n)
	x = border - (xx % border);
      else
	x = xx % border;
      *out = true;
    }
  else if (x >= nx)
    {
      const int border = nx - 1;
      const int n = (int) (x / border) % 2;

      if (n)
	x = border - (x % border);
      else
	x = x % border;
      *out = true;
    }

  return x;
}


/**
 *
 * Cubic interpolation in one dimension
 *
 **/
static
double cubic_interpolation_cell(double v[4],  //interpolation points
		double x      //point to be interpolated
) {
	return v[1]
	         + 0.5 * x
	         * (v[2] - v[0]
	                     + x
	                     * (2.0 * v[0] - 5.0 * v[1] + 4.0 * v[2]
	                                                          - v[3]
	                                                              + x
	                                                              * (3.0 * (v[1] - v[2])
	                                                            		  + v[3] - v[0])));
}

/**
 *
 * Bicubic interpolation in two dimensions
 *
 **/
static
double bicubic_interpolation_cell(double p[4][4], //array containing the interpolation points
		double x,       //x position to be interpolated
		double y        //y position to be interpolated
) {
	double v[4];
	v[0] = cubic_interpolation_cell(p[0], y);
	v[1] = cubic_interpolation_cell(p[1], y);
	v[2] = cubic_interpolation_cell(p[2], y);
	v[3] = cubic_interpolation_cell(p[3], y);
	return cubic_interpolation_cell(v, x);
}


/**
 *
 * Compute the bicubic interpolation of a point in an image.
 * Detect if the point goes outside the image domain.
 *
 **/
double bicubic_interpolation_at(const ofpix_t *input, //image to be interpolated
		const double uu,    //x component of the vector field
		const double vv,    //y component of the vector field
		const int nx,    //image width
		const int ny,    //image height
		bool border_out //if true, return zero outside the region
) {
	const int sx = (uu < 0) ? -1 : 1;
	const int sy = (vv < 0) ? -1 : 1;

	int x, y, mx, my, dx, dy, ddx, ddy;
	bool out = false;

	//apply the corresponding boundary conditions
	switch (BOUNDARY_CONDITION) {

	case 0:
		x = neumann_bc((int) uu, nx, &out);
		y = neumann_bc((int) vv, ny, &out);
		mx = neumann_bc((int) uu - sx, nx, &out);
		my = neumann_bc((int) vv - sx, ny, &out);
		dx = neumann_bc((int) uu + sx, nx, &out);
		dy = neumann_bc((int) vv + sy, ny, &out);
		ddx = neumann_bc((int) uu + 2 * sx, nx, &out);
		ddy = neumann_bc((int) vv + 2 * sy, ny, &out);
		break;

	case 1:
		x = periodic_bc((int) uu, nx, &out);
		y = periodic_bc((int) vv, ny, &out);
		mx = periodic_bc((int) uu - sx, nx, &out);
		my = periodic_bc((int) vv - sx, ny, &out);
		dx = periodic_bc((int) uu + sx, nx, &out);
		dy = periodic_bc((int) vv + sy, ny, &out);
		ddx = periodic_bc((int) uu + 2 * sx, nx, &out);
		ddy = periodic_bc((int) vv + 2 * sy, ny, &out);
		break;

	case 2:
		x = symmetric_bc((int) uu, nx, &out);
		y = symmetric_bc((int) vv, ny, &out);
		mx = symmetric_bc((int) uu - sx, nx, &out);
		my = symmetric_bc((int) vv - sx, ny, &out);
		dx = symmetric_bc((int) uu + sx, nx, &out);
		dy = symmetric_bc((int) vv + sy, ny, &out);
		ddx = symmetric_bc((int) uu + 2 * sx, nx, &out);
		ddy = symmetric_bc((int) vv + 2 * sy, ny, &out);
		break;

	default:
		x = neumann_bc((int) uu, nx, &out);
		y = neumann_bc((int) vv, ny, &out);
		mx = neumann_bc((int) uu - sx, nx, &out);
		my = neumann_bc((int) vv - sx, ny, &out);
		dx = neumann_bc((int) uu + sx, nx, &out);
		dy = neumann_bc((int) vv + sy, ny, &out);
		ddx = neumann_bc((int) uu + 2 * sx, nx, &out);
		ddy = neumann_bc((int) vv + 2 * sy, ny, &out);
		break;
	}

	if (out && border_out)
		return 0.0;

	else {
		//obtain the interpolation points of the image
		const double p11 = input[mx + nx * my];
		const double p12 = input[x + nx * my];
		const double p13 = input[dx + nx * my];
		const double p14 = input[ddx + nx * my];

		const double p21 = input[mx + nx * y];
		const double p22 = input[x + nx * y];
		const double p23 = input[dx + nx * y];
		const double p24 = input[ddx + nx * y];

		const double p31 = input[mx + nx * dy];
		const double p32 = input[x + nx * dy];
		const double p33 = input[dx + nx * dy];
		const double p34 = input[ddx + nx * dy];

		const double p41 = input[mx + nx * ddy];
		const double p42 = input[x + nx * ddy];
		const double p43 = input[dx + nx * ddy];
		const double p44 = input[ddx + nx * ddy];

		//create array
		double pol[4][4] = { { p11, p21, p31, p41 }, { p12, p22, p32, p42 }, {
				p13, p23, p33, p43 }, { p14, p24, p34, p44 } };

		//return interpolation
		return bicubic_interpolation_cell(pol, uu - x, vv - y);
	}
}

/**
  *
  * Compute the bicubic interpolation of a point in an image. 
  * Detects if the point goes outside the image domain
  *
**/
double
bicubic_interpolation_at_color (const ofpix_t *input,	//image to be interpolated
		       const double uu,		//x component of the vector field
		       const double vv,		//y component of the vector field
		       const int nx,		//width of the image
		       const int ny,		//height of the image
		       const int nz,            //number of channels of the image
		       const int k,  		//actual channel
		       const bool border_out	//if true, put zeros outside the region
  )
{
  const int sx = (uu < 0) ? -1 : 1;
  const int sy = (vv < 0) ? -1 : 1;

  int x, y, mx, my, dx, dy, ddx, ddy;
  bool out = false;

  switch (BOUNDARY_CONDITION)
    {

    case 0:
      x = neumann_bc ((int) uu, nx, &out);
      y = neumann_bc ((int) vv, ny, &out);
      mx = neumann_bc ((int) uu - sx, nx, &out);
      my = neumann_bc ((int) vv - sx, ny, &out);
      dx = neumann_bc ((int) uu + sx, nx, &out);
      dy = neumann_bc ((int) vv + sy, ny, &out);
      ddx = neumann_bc ((int) uu + 2 * sx, nx, &out);
      ddy = neumann_bc ((int) vv + 2 * sy, ny, &out);
      break;

    case 1:
      x = periodic_bc ((int) uu, nx, &out);
      y = periodic_bc ((int) vv, ny, &out);
      mx = periodic_bc ((int) uu - sx, nx, &out);
      my = periodic_bc ((int) vv - sx, ny, &out);
      dx = periodic_bc ((int) uu + sx, nx, &out);
      dy = periodic_bc ((int) vv + sy, ny, &out);
      ddx = periodic_bc ((int) uu + 2 * sx, nx, &out);
      ddy = periodic_bc ((int) vv + 2 * sy, ny, &out);
      break;

    case 2:
      x = symmetric_bc ((int) uu, nx, &out);
      y = symmetric_bc ((int) vv, ny, &out);
      mx = symmetric_bc ((int) uu - sx, nx, &out);
      my = symmetric_bc ((int) vv - sx, ny, &out);
      dx = symmetric_bc ((int) uu + sx, nx, &out);
      dy = symmetric_bc ((int) vv + sy, ny, &out);
      ddx = symmetric_bc ((int) uu + 2 * sx, nx, &out);
      ddy = symmetric_bc ((int) vv + 2 * sy, ny, &out);
      break;

    default:
      x = neumann_bc ((int) uu, nx, &out);
      y = neumann_bc ((int) vv, ny, &out);
      mx = neumann_bc ((int) uu - sx, nx, &out);
      my = neumann_bc ((int) vv - sx, ny, &out);
      dx = neumann_bc ((int) uu + sx, nx, &out);
      dy = neumann_bc ((int) vv + sy, ny, &out);
      ddx = neumann_bc ((int) uu + 2 * sx, nx, &out);
      ddy = neumann_bc ((int) vv + 2 * sy, ny, &out);
      break;
    }

  if (out && border_out)

    return 0.0;

  else
    {
      //obtain the interpolation points of the image
      const float p11 = input[(mx  + nx * my) * nz + k];
      const float p12 = input[(x   + nx * my) * nz + k];
      const float p13 = input[(dx  + nx * my) * nz + k];
      const float p14 = input[(ddx + nx * my) * nz + k];

      const float p21 = input[(mx  + nx * y) * nz + k];
      const float p22 = input[(x   + nx * y) * nz + k];
      const float p23 = input[(dx  + nx * y) * nz + k];
      const float p24 = input[(ddx + nx * y) * nz + k];

      const float p31 = input[(mx  + nx * dy) * nz + k];
      const float p32 = input[(x   + nx * dy) * nz + k];
      const float p33 = input[(dx  + nx * dy) * nz + k];
      const float p34 = input[(ddx + nx * dy) * nz + k];

      const float p41 = input[(mx  + nx * ddy) * nz + k];
      const float p42 = input[(x   + nx * ddy) * nz + k];
      const float p43 = input[(dx  + nx * ddy) * nz + k];
      const float p44 = input[(ddx + nx * ddy) * nz + k];

      //create array
      double pol[4][4] = { 
	{p11, p21, p31, p41}, {p12, p22, p32, p42},
	{p13, p23, p33, p43}, {p14, p24, p34, p44}
      };

      //return interpolation
      return bicubic_interpolation_cell (pol, (float) uu - x, (float) vv - y);
    }
}


/**
 *
 * Compute the bicubic interpolation of an image.
 *
 **/
void bicubic_interpolation_warp(const ofpix_t *input,     // image to be warped
		const ofpix_t *u,         // x component of the vector field
		const ofpix_t *v,         // y component of the vector field
		ofpix_t *output,    // image warped with bicubic interpolation
		const int nx,        // image width
		const int ny,        // image height
		bool border_out // if true, put zeros outside the region
) {
#pragma omp parallel for
	for (int i = 0; i < ny; i++)
		for (int j = 0; j < nx; j++) {
			const int p = i * nx + j;
			const double uu = (double) (j + u[p]);
			const double vv = (double) (i + v[p]);

			// obtain the bicubic interpolation at position (uu, vv)
			output[p] = bicubic_interpolation_at(input, uu, vv, nx, ny,
					border_out);
		}
}

/**
  *
  * Compute the bicubic interpolation of an image.
  *
**/
void
bicubic_interpolation_warp_color (const ofpix_t *input,	//image to be warped
		       const ofpix_t *u,		//x component of the vector field
		       const ofpix_t *v,		//y component of the vector field
		       ofpix_t *output,		//warped output image with bicubic interpolation
		       const int nx,		//width of the image
		       const int ny,		//height of the image
		       const int nz,		//number of channels of the image	
		       bool border_out  	//if true, put zeros outside the region
  )
{
  for(int k = 0; k < nz; k++){
    #pragma omp parallel for
    for (int i = 0; i < ny; i++)
       for (int j = 0; j < nx; j++){

	      const int p = i * nx + j;
	      const double uu = (double) (j + u[p]);
	      const double vv = (double) (i + v[p]);

	      //obtain the bicubic interpolation at position (uu, vv)
	      output[p * nz + k] = bicubic_interpolation_at_color (input, uu, vv, nx, ny, nz, k, border_out);
      
      }
  } // end multi-channel loop
}

double me_interpolate_bilinear(const ofpix_t *in, int ncol, double x, double y) {
	int l, k, offset;
	const ofpix_t *x0, *x1, *x2, *x3;
	double a, b, b_1, a_1, phi;

	l = floor(x);
	k = floor(y);

	a = x - l;
	b = y - k;

	offset = k * ncol + l;

	a_1 = 1.0 - a;
	b_1 = 1.0 - b;

	x0 = x1 = in + offset;
	x1++;
	x2 = x3 = x0 + ncol;
	x3++;

	if ((!a) || (!b)) {
		if ((!a) && (!b))
			phi = *x0;
		else if (!a)
			phi = b_1 * (*x0) + b * (*x2);
		else
			phi = a_1 * (*x0) + a * (*x1);
	} else
		phi = b_1 * (a_1 * (*x0) + a * (*x1)) + b * (a_1 * (*x2) + a * (*x3));

	return (phi);
}

double me_interpolate_bicubic(const ofpix_t *in, int nx, int ny, double x,
		double y) {
	int N, M, ncol, nrow, j, k, offset;
	const ofpix_t *p_double, *p_double1;
	double phi, sh, sv;
	double mem_c[16], *mem_pc[4], **c;
	double mem_u[8], *uh, *uv;
	double mem_conv[4], *tmp_conv;
	double sh2, sh3, sv2, sv3;

	/* x,y   : the coordinates of the point to be interpolated
	 uh, uv: the horizontal and vertical components of the convolution kernel.
	 its values depend on sh, sv respectively.
	 mem_u : the space of memory where uv, uh are stored
	 sh    : the distance in the horizontal axis between the point to be interpolated
	 and the interpolation node.
	 sv    : the distance in the vertical axis between the point to be interpolated
	 and the interpolation node.
	 s     : auxiliar variable which contains the computed distances between the point
	 and the interpolation node.
	 phi   : the interpolated value of point (x,y).
	 c     : the image values of the interpolation nodes.
	 mem_c : the memory where the interpolation node values are stored.
	 mem_pc: it is interpreted as a vector of pointers to double. This memory is
	 used to create the matrix where the interpolation nodes are stored.
	 */

	/* uh, uv has range uv[-1..2] and uh[-1..2] */

	uh = mem_u + 1;
	uv = mem_u + 5;

	/* c has range c[-1..2][-1..2] */

	c = mem_pc + 1;

	c[-1] = mem_c + 1;
	c[0] = mem_c + 5;
	c[1] = mem_c + 9;
	c[2] = mem_c + 13;

	/* Memory for tmp_conv[-1..2] */

	tmp_conv = mem_conv + 1;

	/* Initialize variables */

	j = x;
	k = y;

	/* Get nrow and ncol of the image */

	nrow = ny;
	ncol = nx;

	/* Initialize M, N */

	M = nrow - 1;
	N = ncol - 1;

	if ((j <= 0) || (j >= (N - 1))) {
		phi = me_interpolate_bilinear(in, nx, x, y);
		return (phi);
	}

	if ((k <= 0) || (k >= (M - 1))) {
		phi = me_interpolate_bilinear(in, nx, x, y);
		return (phi);
	}

	/* Compute distance of (x,y) to (j,k) */

	sh = x - j;
	sv = y - k;

	sv2 = sv * sv;
	sv3 = sv2 * sv;

	sh2 = sh * sh;
	sh3 = sh2 * sh;

	/* Initialze weights associated to the interpolation nodes */

	uh[-1] = (-sh3 + 2.0 * sh2 - sh) * 0.5;
	uv[-1] = (-sv3 + 2.0 * sv2 - sv) * 0.5;

	uh[0] = (3.0 * sh3 - 5.0 * sh2 + 2.0) * 0.5;
	uv[0] = (3.0 * sv3 - 5.0 * sv2 + 2.0) * 0.5;

	uh[1] = (-3.0 * sh3 + 4.0 * sh2 + sh) * 0.5;
	uv[1] = (-3.0 * sv3 + 4.0 * sv2 + sv) * 0.5;

	uh[2] = (sh3 - sh2) * 0.5;
	uv[2] = (sv3 - sv2) * 0.5;

	/* We are sure that these values cannot fall outside the image */

	offset = k * ncol + j;
	p_double = in + offset;
	p_double1 = p_double + ncol;

	c[0][0] = *p_double;
	c[1][0] = *(p_double + 1);
	c[0][1] = *(p_double1);
	c[1][1] = *(p_double1 + 1);

	/* For the remaining values we have to check if they fall outside the image */

	if (j == 0) {
		c[-1][0] = 3 * in[k * ncol] - 3 * in[k * ncol + 1] + in[k * ncol + 2];
		c[-1][1] = 3 * in[(k + 1) * ncol] - 3 * in[(k + 1) * ncol + 1]
		                                           + in[(k + 1) * ncol + 2];
	} else {
		c[-1][0] = in[k * ncol + (j - 1)];
		c[-1][1] = in[(k + 1) * ncol + (j - 1)];
	}

	if (j == (N - 1)) {
		c[2][0] = 3 * in[k * ncol + j] - 3 * in[k * ncol + (j - 1)]
		                                        + in[k * ncol + (j - 2)];
		c[2][1] = 3 * in[(k + 1) * ncol + j] - 3 * in[(k + 1) * ncol + (j - 1)]
		                                              + in[(k + 1) * ncol + (j - 2)];
	} else {
		c[2][0] = in[k * ncol + (j + 2)];
		c[2][1] = in[(k + 1) * ncol + (j + 2)];
	}

	if (k == 0) {
		c[0][-1] = 3 * in[j] - 3 * in[ncol + j] + in[2 * ncol + j];
		c[1][-1] = 3 * in[j + 1] - 3 * in[ncol + j + 1] + in[2 * ncol + j + 1];
	} else {
		c[0][-1] = in[(k - 1) * ncol + j];
		c[1][-1] = in[(k - 1) * ncol + (j + 1)];
	}

	if (k == (M - 1)) {
		c[0][2] = 3 * in[k * ncol + j] - 3 * in[(k - 1) * ncol + j]
		                                        + in[(k - 2) * ncol + j];
		c[1][2] = 3 * in[k * ncol + (j + 1)] - 3 * in[(k - 1) * ncol + (j + 1)]
		                                              + in[(k - 2) * ncol + (j + 1)];
	} else {
		c[0][2] = in[(k + 2) * ncol + j];
		c[1][2] = in[(k + 2) * ncol + (j + 1)];
	}

	if ((j == 0) && (k == 0))
		c[-1][-1] = 3 * c[0][-1] - 3 * c[1][-1] + c[2][-1];
	else
		c[-1][-1] = in[(k - 1) * ncol + (j - 1)];

	if ((j == (N - 1)) && (k == 0))
		c[2][-1] = 3 * c[1][-1] - 3 * c[0][-1] + c[-1][-1];
	else
		c[2][-1] = in[(k - 1) * ncol + (j + 2)];

	if ((j == 0) && (k == (M - 1)))
		c[-1][2] = 3 * c[0][2] - 3 * c[1][2] + c[2][2];
	else
		c[-1][2] = in[(k - 1) * ncol + (j + 2)];

	if ((j == (N - 1)) && (k == (M - 1)))
		c[2][2] = 3 * c[1][2] - 3 * c[0][2] + c[0][2];
	else
		c[2][2] = in[(k + 2) * ncol + (j + 2)];

	/* Ok, now compute convolution */

	for (j = -1; j <= 2; j++) {
		phi = 0;

		for (k = -1; k <= 2; k++)
			phi += c[j][k] * uv[k];

		tmp_conv[j] = phi;
	}

	phi = 0;

	for (j = -1; j <= 2; j++)
		phi += tmp_conv[j] * uh[j];

	/*
	 phi = 0;

	 for(j = -1; j <= 2; j++)
	 for(k = -1; k <= 2; k++)
	 phi += c[j][k] * uh[j] * uv[k];
	 */

	return (phi);
}

void me_image_restriction(const ofpix_t *in, ofpix_t *out, int ncol, int nrow,
		int new_ncol, int new_nrow) {
	int i, j;
	double coord_x, coord_y;
	double min_x, min_y;
	double value, gamma_x, gamma_y;

	/* compute gamma factors */

	gamma_x = (double) ncol / (double) new_ncol;
	gamma_y = (double) nrow / (double) new_nrow;

	/* min_x and min_y */

	min_x = gamma_x / 2.0 - 0.5;
	min_y = gamma_y / 2.0 - 0.5;

	/* Set samples in output image. We use here bilinear interpolation. */

	for (i = 0; i < new_nrow; i++) {
		coord_y = min_y + i * gamma_y;

		for (j = 0; j < new_ncol; j++) {
			coord_x = min_x + j * gamma_x;

			value = me_interpolate_bilinear(in, ncol, coord_x, coord_y);
			out[i * (int) new_ncol + j] = value;
		}
	}
}

