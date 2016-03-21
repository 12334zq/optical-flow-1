// This program is free software: you can use, modify and/or redistribute it
// under the terms of the simplified BSD License. You should have received a
// copy of this license along this program. If not, see
// <http://www.opensource.org/licenses/bsd-license.html>.
//
// Copyright (C) 2012, Javier Sánchez Pérez <jsanchez@dis.ulpgc.es>
// All rights reserved.

#include "mask.h"


/**
 *
 * Details on how to compute the divergence and the grad(u) can be found in:
 * [2] A. Chambolle, "An Algorithm for Total Variation Minimization and
 * Applications", Journal of Mathematical Imaging and Vision, 20: 89-97, 2004
 *
 **/


/**
 *
 * Function to compute the divergence with backward differences
 * (see [2] for details)
 *
 **/
void divergence(
		const float *v1, // x component of the vector field
		const float *v2, // y component of the vector field
		float *div,      // output divergence
		const int nx,    // image width
		const int ny     // image height
	       )
{
	// compute the divergence on the central body of the image
#pragma omp parallel for schedule(dynamic)
	for (int i = 1; i < ny-1; i++)
	{
		for(int j = 1; j < nx-1; j++)
		{
			const int p  = i * nx + j;
			const int p1 = p - 1;
			const int p2 = p - nx;

			const float v1x = v1[p] - v1[p1];
			const float v2y = v2[p] - v2[p2];

			div[p] = v1x + v2y;
		}
	}

	// compute the divergence on the first and last rows
	for (int j = 1; j < nx-1; j++)
	{
		const int p = (ny-1) * nx + j;

		div[j] = v1[j] - v1[j-1] + v2[j];
		div[p] = v1[p] - v1[p-1] - v2[p-nx];
	}

	// compute the divergence on the first and last columns
	for (int i = 1; i < ny-1; i++)
	{
		const int p1 = i * nx;
		const int p2 = (i+1) * nx - 1;

		div[p1] =  v1[p1]   + v2[p1] - v2[p1 - nx];
		div[p2] = -v1[p2-1] + v2[p2] - v2[p2 - nx];

	}

	div[0]         =  v1[0] + v2[0];
	div[nx-1]      = -v1[nx - 2] + v2[nx - 1];
	div[(ny-1)*nx] =  v1[(ny-1)*nx] - v2[(ny-2)*nx];
	div[ny*nx-1]   = -v1[ny*nx - 2] - v2[(ny-1)*nx - 1];
}


/**
 *
 * Function to compute the gradient with forward differences
 * (see [2] for details)
 *
 **/
void forward_gradient(
		const float *f, //input image
		float *fx,      //computed x derivative
		float *fy,      //computed y derivative
		const int nx,   //image width
		const int ny    //image height
		)
{
	// compute the gradient on the central body of the image
#pragma omp parallel for schedule(dynamic)
	for (int i = 0; i < ny-1; i++)
	{
		for(int j = 0; j < nx-1; j++)
		{
			const int p  = i * nx + j;
			const int p1 = p + 1;
			const int p2 = p + nx;

			fx[p] = f[p1] - f[p];
			fy[p] = f[p2] - f[p];
		}
	}

	// compute the gradient on the last row
	for (int j = 0; j < nx-1; j++)
	{
		const int p = (ny-1) * nx + j;

		fx[p] = f[p+1] - f[p];
		fy[p] = 0;
	}

	// compute the gradient on the last column
	for (int i = 1; i < ny; i++)
	{
		const int p = i * nx-1;

		fx[p] = 0;
		fy[p] = f[p+nx] - f[p];
	}

	fx[ny * nx - 1] = 0;
	fy[ny * nx - 1] = 0;
}

/**
 *
 * Function to apply a 3x3 mask to an image
 *
 */
static
void mask3x3(
    const float *input, //input image
    float *output,      //output image
    const int nx,       //image width
    const int ny,       //image height
    const float *mask   //mask to be applied
)
{
    //apply the mask to the center body of the image
#pragma omp parallel for
    for(int i = 1; i < ny-1; i++)
    {
	for(int j = 1; j < nx-1; j++)
	{
	    double sum = 0;
	    for(int l = 0; l < 3; l++)
	    {
		for(int m = 0; m < 3; m++)
		{
		    int p = (i + l -1) * nx + j + m -1;
		    sum += input[p] * mask[l * 3 + m];
		}
	    }
	    int k = i * nx + j;
	    output[k] = sum;
	}
    }

    //apply the mask to the first and last rows
#pragma omp parallel for
    for(int j = 1; j < nx-1; j++)
    {
	double sum = 0;
	sum += input[j-1] * (mask[0] + mask[3]);
	sum += input[ j ] * (mask[1] + mask[4]);
	sum += input[j+1] * (mask[2] + mask[5]);

	sum += input[nx + j-1] * mask[6];
	sum += input[nx +  j ] * mask[7];
	sum += input[nx + j+1] * mask[8];

	output[j] = sum;

	sum = 0;
	sum += input[(ny-2)*nx+j-1] * mask[0];
	sum += input[(ny-2)*nx+j  ] * mask[1];
	sum += input[(ny-2)*nx+j+1] * mask[2];

	sum += input[(ny-1)*nx+j-1] * (mask[6] + mask[3]);
	sum += input[(ny-1)*nx+j  ] * (mask[7] + mask[4]);
	sum += input[(ny-1)*nx+j+1] * (mask[8] + mask[5]);

	output[(ny-1)*nx + j] = sum;
    }

    //apply the mask to the first and last columns
#pragma omp parallel for
    for(int i = 1; i < ny-1; i++)
    {
	double sum = 0;
	sum += input[(i - 1)*nx]   * (mask[0] + mask[1]);
	sum += input[(i - 1)*nx+1] * mask[2];

	sum += input[i * nx]   * (mask[3] + mask[4]);
	sum += input[i * nx+1] * mask[5];

	sum += input[(i + 1)*nx]   * (mask[6] + mask[7]);
	sum += input[(i + 1)*nx+1] * mask[8];

	output[i*nx] = sum;

	sum = 0;
	sum += input[i * nx-2] * mask[0];
	sum += input[i * nx-1] * (mask[1] + mask[2]);

	sum += input[(i + 1)*nx-2] * mask[3];
	sum += input[(i + 1)*nx-1] * (mask[4] + mask[5]);

	sum += input[(i + 2)*nx-2] * mask[6];
	sum += input[(i + 2)*nx-1] * (mask[7] + mask[8]);

	output[i*nx + nx -1] = sum;
    }

    //apply the mask to the four corners
    output[0] = input[0]    * (mask[0] + mask[1] + mask[3] + mask[4]) +
	input[1]    * (mask[2] + mask[5]) +
	input[nx]   * (mask[6] + mask[7]) +
	input[nx+1] * mask[8];

    output[nx-1] =
	input[nx-2]   * (mask[0] + mask[3]) +
	input[nx-1]   * (mask[1] + mask[2] + mask[4] + mask[5]) +
	input[2*nx-2] * mask[6] +
	input[2*nx-1] * (mask[7] + mask[8]);

    output[(ny-1)*nx] =
	input[(ny-2)*nx]   * (mask[0] + mask[1]) +
	input[(ny-2)*nx+1] *  mask[2] +
	input[(ny-1)*nx]   * (mask[3] + mask[4] + mask[6] + mask[7]) +
	input[(ny-1)*nx+1] * (mask[5] + mask[8]);

    output[ny*nx-1] =
	input[(ny-1)*nx-2] * mask[0] +
	input[(ny-1)*nx-1] * (mask[1] + mask[2]) +
	input[ny*nx-2] * (mask[3] + mask[6]) +
	input[ny*nx-1] * (mask[4] + mask[5] + mask[7] + mask[8]);
}

/**
 *
 * Compute the second order X derivative
 *
 */
void Dxx(
    const float *I, //input image
    float *Ixx,     //oputput derivative
    const int nx,   //image width
    const int ny    //image height
)
{
    //mask of second derivative
    float M[]  = {0., 0., 0.,
		  1.,-2., 1.,
		  0., 0., 0.};

    //computing the second derivative
    mask3x3(I, Ixx, nx, ny, M);
}


/**
 *
 * Compute the second order Y derivative
 *
 */
void Dyy(
    const float *I, //input image
    float *Iyy,     //oputput derivative
    const int nx,   //image width
    const int ny    //image height
)
{
    //mask of second derivative
    float M[]  = {0., 1., 0.,
		  0.,-2., 0.,
		  0., 1., 0.};

    //computing the second derivative
    mask3x3(I, Iyy, nx, ny, M);
}


/**
 *
 * Compute the second order XY derivative
 *
 */
void Dxy(
    const float *I, //input image
    float *Ixy,     //oputput derivative
    const int nx,   //image width
    const int ny    //image height
)
{
    //mask of second derivative
    float M[]  = {1./4., 0.,-1./4.,
		  0.,    0., 0.,
		  -1./4., 0., 1./4.};

    //computing the second derivative
    mask3x3(I, Ixy, nx, ny, M);
}

/**
 *
 * Function to compute the gradient with centered differences
 *
 */
void centered_gradient(
    const float *input, //input image
    float *dx,          //computed x derivative
    float *dy,          //computed y derivative
    const int nx,       //image width
    const int ny        //image height
)
{
    //gradient in the center body of the image
#pragma omp parallel for schedule(dynamic)
    for(int i = 1; i < ny-1; i++)
    {
	for(int j = 1; j < nx-1; j++)
	{
	    const int k = i * nx + j;
	    dx[k] = 0.5*(input[k+1] - input[k-1]);
	    dy[k] = 0.5*(input[k+nx] - input[k-nx]);
	}
    }

    //gradient in the first and last rows
#pragma omp parallel for
    for(int j = 1; j < nx-1; j++)
    {
	dx[j] = 0.5*(input[j+1] - input[j-1]);
	dy[j] = 0.5*(input[j+nx] - input[j]);

	const int k = (ny - 1) * nx + j;

	dx[k] = 0.5*(input[k+1] - input[k-1]);
	dy[k] = 0.5*(input[k] - input[k-nx]);
    }

    //gradient in the first and last columns
#pragma omp parallel for
    for(int i = 1; i < ny-1; i++)
    {
	const int p = i * nx;
	dx[p] = 0.5*(input[p+1] - input[p]);
	dy[p] = 0.5*(input[p+nx] - input[p-nx]);

	const int k = (i+1) * nx - 1;

	dx[k] = 0.5*(input[k] - input[k-1]);
	dy[k] = 0.5*(input[k+nx] - input[k-nx]);
    }

    //calculate the gradient in the corners
    dx[0] = 0.5*(input[1] - input[0]);
    dy[0] = 0.5*(input[nx] - input[0]);

    dx[nx-1] = 0.5*(input[nx-1] - input[nx-2]);
    dy[nx-1] = 0.5*(input[2*nx-1] - input[nx-1]);

    dx[(ny-1)*nx] = 0.5*(input[(ny-1)*nx + 1] - input[(ny-1)*nx]);
    dy[(ny-1)*nx] = 0.5*(input[(ny-1)*nx] - input[(ny-2)*nx]);

    dx[ny*nx-1] = 0.5*(input[ny*nx-1] - input[ny*nx-1-1]);
    dy[ny*nx-1] = 0.5*(input[ny*nx-1] - input[(ny-1)*nx-1]);
}



/**
 *
 * Compute the 3D gradient with central differences
 *
 */
void centered_gradient3(
    const float *input, //input image
    float *dx,          //x derivative
    float *dy,          //y derivative
    float *dz,          //z derivative
    const int nx,       //image width
    const int ny,       //image height
    const int nz        //image depth
)
{
    const int df = nx * ny;

    //compute the x and y derivatives for all frames
    for(int f = 0; f < nz; f++)
    {
	for(int i = 1; i < ny-1; i++)
	{
	    for(int j = 1; j < nx-1; j++)
	    {
		const int k = f * df + i * nx + j;
		dx[k] = 0.5 * (input[k+1]  - input[k-1]);
		dy[k] = 0.5 * (input[k+nx] - input[k-nx]);
	    }
	}

	//gradient in the first and last rows
	for(int j = 1; j < nx-1; j++)
	{
	    int k = f * df + j;
	    dx[k] = 0.5 * (input[k+1]  - input[k-1]);
	    dy[k] = 0.5 * (input[k+nx] - input[k]);

	    k = f * df + (ny - 1) * nx + j;

	    dx[k] = 0.5 * (input[k+1] - input[k-1]);
	    dy[k] = 0.5 * (input[k]   - input[k-nx]);
	}

	//gradient in the first and last columns
	for(int i = 1; i < ny-1; i++)
	{
	    int k = f * df + i * nx;
	    dx[k] = 0.5 * (input[k+1]  - input[k]);
	    dy[k] = 0.5 * (input[k+nx] - input[k-nx]);

	    k = f * df + (i+1) * nx - 1;

	    dx[k] = 0.5 * (input[k] - input[k-1]);
	    dy[k] = 0.5 * (input[k+nx] - input[k-nx]);
	}

	//calculate the gradient in the corners
	int k = f * df;
	dx[k] = 0.5 * (input[k+1]  - input[k]);
	dy[k] = 0.5 * (input[k+nx] - input[k]);

	k = f * df + nx - 1;
	dx[k] = 0.5 * (input[k] - input[k-1]);
	dy[k] = 0.5 * (input[k+nx] - input[k]);

	k = f * df + (ny - 1) * nx;
	dx[k] = 0.5 * (input[k + 1] - input[k]);
	dy[k] = 0.5 * (input[k] - input[k-nx]);

	k = f * df + ny * nx - 1;
	dx[k] = 0.5 * (input[k] - input[k-1]);
	dy[k] = 0.5 * (input[k] - input[k-nx]);
    }

    if(nz > 1)
    {
	//compute the z derivative for the interior frames
	for(int f = 1; f < nz-1; f++)
	{
	    for(int i = 0; i < df; i++)
	    {
		const int k = f * df + i;
		dz[k] = 0.5 * (input[k+df] - input[k-df]);
	    }
	}

	//compute the z derivative for the first and last frames
	for(int i = 0; i < df; i++)
	{
	    int k = i;
	    dz[k] = 0.5 * (input[k+df] - input[k]);

	    k = (nz-1) * df + i;
	    dz[k] = 0.5 * (input[k] - input[k-df]);
	}
    }
    else
      for(int i = 0; i < df; i++) 
	  dz[i] = 0;

}
