// This program is free software: you can use, modify and/or redistribute it
// under the terms of the simplified BSD License. You should have received a
// copy of this license along this program. If not, see
// <http://www.opensource.org/licenses/bsd-license.html>.
//
// Copyright (C) 2012, Javier Sánchez Pérez <jsanchez@dis.ulpgc.es>
// All rights reserved.

#include "operators.h"

#include <cmath>
#include <stdexcept>

#include "of.h"

//#include "stdio.h"

using namespace std;

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
void
divergence(const ofpix_t *v1, // x component of the vector field
           const ofpix_t *v2, // y component of the vector field
           ofpix_t *div,   // output divergence
           const int nx, // image width
           const int ny  // image height
           )
{
    // compute the divergence on the central body of the image
#pragma omp parallel for schedule(dynamic)
    for (int i = 1; i < ny - 1; i++) {
        for (int j = 1; j < nx - 1; j++) {
            const int p  = i * nx + j;
            const int p1 = p - 1;
            const int p2 = p - nx;
            const double v1x = v1[p] - v1[p1];
            const double v2y = v2[p] - v2[p2];

            div[p] = v1x + v2y;
        }
    }

    // compute the divergence on the first and last rows
    for (int j = 1; j < nx - 1; j++) {
        const int p = (ny - 1) * nx + j;

        div[j] = v1[j] - v1[j - 1] + v2[j];
        div[p] = v1[p] - v1[p - 1] - v2[p - nx];
    }

    // compute the divergence on the first and last columns
    for (int i = 1; i < ny - 1; i++) {
        const int p1 = i * nx;
        const int p2 = (i + 1) * nx - 1;

        div[p1] =  v1[p1]   + v2[p1] - v2[p1 - nx];
        div[p2] = -v1[p2 - 1] + v2[p2] - v2[p2 - nx];
    }

    div[0]         =  v1[0] + v2[0];
    div[nx - 1]      = -v1[nx - 2] + v2[nx - 1];
    div[(ny - 1) * nx] =  v1[(ny - 1) * nx] - v2[(ny - 2) * nx];
    div[ny * nx - 1]   = -v1[ny * nx - 2] - v2[(ny - 1) * nx - 1];
}

/**
 *
 * Function to compute the gradient with forward differences
 * (see [2] for details)
 *
 **/
void
forward_gradient(const ofpix_t *f, //input image
                 ofpix_t *fx, //computed x derivative
                 ofpix_t *fy, //computed y derivative
                 const int nx, //image width
                 const int ny //image height
                 )
{
    // compute the gradient on the central body of the image
#pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < ny - 1; i++) {
        for (int j = 0; j < nx - 1; j++) {
            const int p  = i * nx + j;
            const int p1 = p + 1;
            const int p2 = p + nx;

            fx[p] = f[p1] - f[p];
            fy[p] = f[p2] - f[p];
        }
    }

    // compute the gradient on the last row
    for (int j = 0; j < nx - 1; j++) {
        const int p = (ny - 1) * nx + j;

        fx[p] = f[p + 1] - f[p];
        fy[p] = 0;
    }

    // compute the gradient on the last column
    for (int i = 1; i < ny; i++) {
        const int p = i * nx - 1;

        fx[p] = 0;
        fy[p] = f[p + nx] - f[p];
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
void
mask3x3(const ofpix_t *input, //input image
        ofpix_t *output,  //output image
        const int nx,   //image width
        const int ny,   //image height
        const int nz,          // number of color channels in the image
        const ofpix_t *mask //mask to be applied
        )
{
    int nx_multichannel = nx * nz;

    for (int index_multichannel = 0; index_multichannel < nz; index_multichannel++) {
        //apply the mask to the center body of the image

#pragma omp parallel for
        for (int i = 1; i < ny - 1; i++) {
            for (int j = 1; j < nx - 1; j++) {
                int k = (i * nx + j) * nz + index_multichannel;
                double sum = 0;
                for (int l = 0; l < 3; l++) {
                    for (int m = 0; m < 3; m++) {
                        int p = ( (i + l - 1) * nx + j + m - 1 ) * nz + index_multichannel;
                        sum += input[p] * mask[l * 3 + m];
                    }
                }
                output[k] = sum;
            }
        }

        //apply the mask to the first and last rows
#pragma omp parallel for
        for (int j = 1; j < nx - 1; j++) {
            int index = j * nz + index_multichannel;
            double sum = 0;

            sum += input[index - nz] * (mask[0] + mask[3]);
            sum += input[index] * (mask[1] + mask[4]);
            sum += input[index + nz] * (mask[2] + mask[5]);

            sum += input[nx_multichannel + j - nz] * mask[6];
            sum += input[nx_multichannel + j] * mask[7];
            sum += input[nx_multichannel + j + nz] * mask[8];

            output[j] = sum;

            index = ( (ny - 2) * nx + j ) * nz + index_multichannel;

            sum = 0;
            sum += input[index - nz] * mask[0];
            sum += input[index] * mask[1];
            sum += input[index + nz] * mask[2];

            index = ( (ny - 1) * nx + j ) * nz + index_multichannel;

            sum += input[index - nz] * (mask[6] + mask[3]);
            sum += input[index] * (mask[7] + mask[4]);
            sum += input[index + 1] * (mask[8] + mask[5]);

            output[index] = sum;
        }

        //apply the mask to the first and last columns
#pragma omp parallel for
        for (int i = 1; i < ny - 1; i++) {
            int index = i * nx_multichannel + index_multichannel;
            double sum = 0;
            int index_row = (i - 1) * nx_multichannel + index_multichannel;

            sum += input[index_row] * (mask[0] + mask[1]);
            sum += input[index_row + nz] * mask[2];

            sum += input[index] * (mask[3] + mask[4]);
            sum += input[index + nz] * mask[5];

            index_row = (i + 1) * nx_multichannel + index_multichannel;

            sum += input[index_row] * (mask[6] + mask[7]);
            sum += input[index_row + nz] * mask[8];

            output[index] = sum;

            sum = 0;
            sum += input[index - 2 * nz] * mask[0];
            sum += input[index - nz] * (mask[1] + mask[2]);

            index_row = (i + 1) * nx_multichannel + index_multichannel;

            sum += input[index_row - 2 * nz] * mask[3];
            sum += input[index_row - nz] * (mask[4] + mask[5]);

            index_row = (i + 2) * nx_multichannel + index_multichannel;

            sum += input[index_row - 2 * nz] * mask[6];
            sum += input[index_row - nz] * (mask[7] + mask[8]);

            output[(i * nx + nx - 1) * nz + index_multichannel] = sum;
        }

        //apply the mask to the four corners
        output[index_multichannel] =
            input[index_multichannel] * (mask[0] + mask[1] + mask[3] + mask[4]) +
            input[index_multichannel + nz] * (mask[2] + mask[5]) +
            input[nx_multichannel + index_multichannel] * (mask[6] + mask[7]) +
            input[nx_multichannel + index_multichannel + nz] * mask[8];

        output[nx_multichannel - nz + index_multichannel] =
            input[(nx - 2) * nz + index_multichannel] * (mask[0] + mask[3]) +
            input[(nx - 1) * nz + index_multichannel] * (mask[1] + mask[2] + mask[4] + mask[5]) +
            input[(2 * nx - 2) * nz + index_multichannel] * mask[6] +
            input[(2 * nx - 1) * nz + index_multichannel] * (mask[7] + mask[8]);

        output[(ny - 1) * nx_multichannel + index_multichannel] =
            input[(ny - 2) * nx_multichannel + index_multichannel] * (mask[0] + mask[1]) +
            input[( (ny - 2) * nx + 1 ) * nz + index_multichannel] * mask[2] +
            input[(ny - 1) * nx_multichannel + index_multichannel] * (mask[3] + mask[4] + mask[6] + mask[7]) +
            input[( (ny - 1) * nx + 1 ) * nz + index_multichannel] * (mask[5] + mask[8]);

        output[(ny * nx - 1) * nz + index_multichannel] =
            input[( (ny - 1) * nx - 2 ) * nz + index_multichannel] * mask[0] +
            input[( (ny - 1) * nx - 1 ) * nz + index_multichannel] * (mask[1] + mask[2]) +
            input[(ny * nx - 2) * nz + index_multichannel] * (mask[3] + mask[6]) +
            input[(ny * nx - 1) * nz + index_multichannel] * (mask[4] + mask[5] + mask[7] + mask[8]);
    } // end loop for channels information
} // mask3x3

/**
 *
 * Compute the second order X derivative
 *
 */
void
Dxx(const ofpix_t *I, //input image
    ofpix_t *Ixx,     //oputput derivative
    const int nx,   //image width
    const int ny,    //image height
    const int nz               //number of color channels in the image
    )
{
    //mask of second derivative
    ofpix_t M[]  = {
        0., 0., 0.,
        1., -2., 1.,
        0., 0., 0.
    };

    //computing the second derivative
    mask3x3 (I, Ixx, nx, ny, nz, M);
}

/**
 *
 * Compute the second order Y derivative
 *
 */
void
Dyy(const ofpix_t *I, //input image
    ofpix_t *Iyy,     //oputput derivative
    const int nx,   //image width
    const int ny,    //image height
    const int nz               //number of color channels in the image
    )
{
    //mask of second derivative
    ofpix_t M[]  = {
        0., 1., 0.,
        0., -2., 0.,
        0., 1., 0.
    };

    //computing the second derivative
    mask3x3 (I, Iyy, nx, ny, nz, M);
}

/**
 *
 * Compute the second order XY derivative
 *
 */
void
Dxy(const ofpix_t *I, //input image
    ofpix_t *Ixy,     //oputput derivative
    const int nx,   //image width
    const int ny,    //image height
    const int nz               //number of color channels in the image
    )
{
    //mask of second derivative
    ofpix_t M[]  = {
        1. / 4., 0., -1. / 4.,
        0.,    0., 0.,
        -1. / 4., 0., 1. / 4.
    };

    //computing the second derivative
    mask3x3(I, Ixy, nx, ny, nz, M);
}

/**
 *
 * Function to compute the gradient with centered differences
 *
 */
void
centered_gradient(const ofpix_t *input, //input image
                  ofpix_t *dx, //computed x derivative
                  ofpix_t *dy, //computed y derivative
                  const int nx, //image width
                  const int ny,     //image height
                  const int nz          //number of color channels in the image
                  )
{
    const int nx_multichannel = nx * nz;

    for (int index_multichannel = 0; index_multichannel < nz; index_multichannel++) {
        //gradient in the center body of the image
#pragma omp parallel for
        for (int i = 1; i < ny - 1; i++) {
            for (int j = 1; j < nx - 1; j++) {
                const int k = (i * nx + j) * nz + index_multichannel;

                dx[k] = 0.5 * (input[k + nz] - input[k - nz]);
                dy[k] = 0.5 * (input[k + nx_multichannel] - input[k - nx_multichannel]);
            }
        }

        //gradient in the first and last rows
#pragma omp parallel for
        for (int j = 1; j < nx - 1; j++) {
            const int index = j * nz + index_multichannel;

            dx[j] = 0.5 * (input[index + nz] - input[index - nz]);
            dy[j] = 0.5 * (input[index + nx_multichannel] - input[index]);

            const int k = ( (ny - 1) * nx + j ) * nz + index_multichannel;

            dx[k] = 0.5 * (input[k + nz] - input[k - nz]);
            dy[k] = 0.5 * (input[k] - input[k - nx_multichannel]);
        }

        //gradient in the first and last columns
#pragma omp parallel for
        for (int i = 1; i < ny - 1; i++) {
            const int p = (i * nx_multichannel) + index_multichannel;

            dx[p] = 0.5 * (input[p + nz] - input[p]);
            dy[p] = 0.5 * (input[p + nx_multichannel] - input[p - nx_multichannel]);

            const int k = ( (i + 1) * nx - 1 ) * nz + index_multichannel;

            dx[k] = 0.5 * (input[k] - input[k - nz]);
            dy[k] = 0.5 * (input[k + nx_multichannel] - input[k - nx_multichannel]);
        }


        //calculate the gradient in the corners
        dx[index_multichannel] = 0.5 * (input[index_multichannel + nz] - input[index_multichannel]);
        dy[index_multichannel] = 0.5 * (input[nx_multichannel + index_multichannel] - input[index_multichannel]);

        const int corner_up_right = (nx - 1) * nz + index_multichannel;

        dx[corner_up_right] = 0.5 * (input[corner_up_right] - input[corner_up_right - nz]);
        dy[corner_up_right] = 0.5 * (input[(2 * nx_multichannel) + index_multichannel - nz] - input[corner_up_right]);

        const int corner_down_left = ( (ny - 1) * nx ) * nz + index_multichannel;

        dx[corner_down_left] = 0.5 * (input[corner_down_left + nz] - input[corner_down_left]);
        dy[corner_down_left] = 0.5 * (input[corner_down_left] - input[(ny - 2) * nx_multichannel + index_multichannel]);

        const int corner_down_right = ny * nx_multichannel - nz + index_multichannel;

        dx[corner_down_right] = 0.5 * (input[corner_down_right] - input[corner_down_right - nz]);
        dy[corner_down_right] = 0.5 * (input[corner_down_right] - input[(ny - 1) * nx_multichannel - nz + index_multichannel]);
    } // end loop for multi-channel
} // centered_gradient

/**
 *
 * Compute the 3D gradient with central differences
 *
 */
void
centered_gradient3(const ofpix_t *input, //input image
                   ofpix_t *dx, //x derivative
                   ofpix_t *dy, //y derivative
                   ofpix_t *dz, //z derivative
                   const int nx, //image width
                   const int ny, //image height
                   const int nz //image depth
                   )
{
    const int df = nx * ny;

    //compute the x and y derivatives for all frames
    for (int f = 0; f < nz; f++) {
        for (int i = 1; i < ny - 1; i++) {
            for (int j = 1; j < nx - 1; j++) {
                const int k = f * df + i * nx + j;
                dx[k] = 0.5 * (input[k + 1]  - input[k - 1]);
                dy[k] = 0.5 * (input[k + nx] - input[k - nx]);
            }
        }

        //gradient in the first and last rows
        for (int j = 1; j < nx - 1; j++) {
            int k = f * df + j;
            dx[k] = 0.5 * (input[k + 1]  - input[k - 1]);
            dy[k] = 0.5 * (input[k + nx] - input[k]);

            k = f * df + (ny - 1) * nx + j;

            dx[k] = 0.5 * (input[k + 1] - input[k - 1]);
            dy[k] = 0.5 * (input[k]   - input[k - nx]);
        }

        //gradient in the first and last columns
        for (int i = 1; i < ny - 1; i++) {
            int k = f * df + i * nx;
            dx[k] = 0.5 * (input[k + 1]  - input[k]);
            dy[k] = 0.5 * (input[k + nx] - input[k - nx]);

            k = f * df + (i + 1) * nx - 1;

            dx[k] = 0.5 * (input[k] - input[k - 1]);
            dy[k] = 0.5 * (input[k + nx] - input[k - nx]);
        }

        //calculate the gradient in the corners
        int k = f * df;
        dx[k] = 0.5 * (input[k + 1]  - input[k]);
        dy[k] = 0.5 * (input[k + nx] - input[k]);

        k = f * df + nx - 1;
        dx[k] = 0.5 * (input[k] - input[k - 1]);
        dy[k] = 0.5 * (input[k + nx] - input[k]);

        k = f * df + (ny - 1) * nx;
        dx[k] = 0.5 * (input[k + 1] - input[k]);
        dy[k] = 0.5 * (input[k] - input[k - nx]);

        k = f * df + ny * nx - 1;
        dx[k] = 0.5 * (input[k] - input[k - 1]);
        dy[k] = 0.5 * (input[k] - input[k - nx]);
    }

    if (nz > 1) {
        //compute the z derivative for the interior frames
        for (int f = 1; f < nz - 1; f++) {
            for (int i = 0; i < df; i++) {
                const int k = f * df + i;
                dz[k] = 0.5 * (input[k + df] - input[k - df]);
            }
        }

        //compute the z derivative for the first and last frames
        for (int i = 0; i < df; i++) {
            int k = i;
            dz[k] = 0.5 * (input[k + df] - input[k]);

            k = (nz - 1) * df + i;
            dz[k] = 0.5 * (input[k] - input[k - df]);
        }
    } else {
        for (int i = 0; i < df; i++) {
            dz[i] = 0;
        }
    }
} // centered_gradient3

/**
 *
 * In-place Gaussian smoothing of an image
 *
 */
void
gaussian(ofpix_t *I,        //input/output image
         const int xdim,  //image width
         const int ydim,  //image height
         const double sigma, //Gaussian sigma
         const int boundary_condition,  //boundary condition
         const int window_size //defines the size of the window
         )
{
    const double den  = 2 * sigma * sigma;
    const int size = (int) (window_size * sigma) + 1;
    const int bdx  = xdim + size;
    const int bdy  = ydim + size;

    if ( boundary_condition && (size > xdim) ) {
        throw std::runtime_error("GaussianSmooth: sigma too large");
    }

    // compute the coefficients of the 1D convolution kernel
    double *B = new double[size];
    for (int i = 0; i < size; i++) {
        B[i] = 1 / ( sigma * sqrt(2.0 * 3.1415926) ) * exp(-i * i / den);
    }

    // normalize the 1D convolution kernel
    double norm = 0;
    for (int i = 0; i < size; i++) {
        norm += B[i];
    }
    norm *= 2;
    norm -= B[0];
    for (int i = 0; i < size; i++) {
        B[i] /= norm;
    }

    // convolution of each line of the input image
    ofpix_t *R = new ofpix_t[size + xdim + size];

    for (int k = 0; k < ydim; k++) {
        int i, j;
        for (i = size; i < bdx; i++) {
            R[i] = I[k * xdim + i - size];
        }

        switch (boundary_condition) {
        case BOUNDARY_CONDITION_DIRICHLET:
            for (i = 0, j = bdx; i < size; i++, j++) {
                R[i] = R[j] = 0;
            }
            break;

        case BOUNDARY_CONDITION_REFLECTING:
            for (i = 0, j = bdx; i < size; i++, j++) {
                R[i] = I[k * xdim + size - i];
                R[j] = I[k * xdim + xdim - i - 1];
            }
            break;

        case BOUNDARY_CONDITION_PERIODIC:
            for (i = 0, j = bdx; i < size; i++, j++) {
                R[i] = I[k * xdim + xdim - size + i];
                R[j] = I[k * xdim + i];
            }
            break;
        }

        for (i = size; i < bdx; i++) {
            double sum = B[0] * R[i];
            for (j = 1; j < size; j++) {
                sum += B[j] * ( R[i - j] + R[i + j] );
            }
            I[k * xdim + i - size] = sum;
        }
    }

    // convolution of each column of the input image
    ofpix_t *T = new ofpix_t[size + ydim + size];

    for (int k = 0; k < xdim; k++) {
        int i, j;
        for (i = size; i < bdy; i++) {
            T[i] = I[(i - size) * xdim + k];
        }

        switch (boundary_condition) {
        case BOUNDARY_CONDITION_DIRICHLET:
            for (i = 0, j = bdy; i < size; i++, j++) {
                T[i] = T[j] = 0;
            }
            break;

        case BOUNDARY_CONDITION_REFLECTING:
            for (i = 0, j = bdy; i < size; i++, j++) {
                T[i] = I[(size - i) * xdim + k];
                T[j] = I[(ydim - i - 1) * xdim + k];
            }
            break;

        case BOUNDARY_CONDITION_PERIODIC:
            for (i = 0, j = bdx; i < size; i++, j++) {
                T[i] = I[(ydim - size + i) * xdim + k];
                T[j] = I[i * xdim + k];
            }
            break;
        }

        for (i = size; i < bdy; i++) {
            double sum = B[0] * T[i];
            for (j = 1; j < size; j++) {
                sum += B[j] * (T[i - j] + T[i + j]);
            }
            I[(i - size) * xdim + k] = sum;
        }
    }

    delete [] B;
    delete [] R;
    delete [] T;
} // gaussian

