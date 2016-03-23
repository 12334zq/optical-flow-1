
// This program is free software: you can use, modify and/or redistribute it
// under the terms of the simplified BSD License. You should have received a
// copy of this license along this program. If not, see
// <http://www.opensource.org/licenses/bsd-license.html>.
//
// Copyright (C) 2011, Javier Sánchez Pérez <jsanchez@dis.ulpgc.es>
// Copyright (C) 2012, Coloma Ballester <coloma.ballester@upf.edu>
// Copyright (C) 2013-2014 J. F. Garamendi <jf.garamendi@upf.edu>
// All rights reserved.

#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#ifdef _OPENMP
#include <omp.h>
#endif


#include "iio.h"
#include "tvl1occflow.h"

#include "tvl1occflow_constants.h"


#define MIN(x, y) ( (x) < (y) ? (x) : (y) )


/**
 *
 *  Function to read images using the iio library
 *  It always returns an allocated the image.
 *
 */
static
ofpix_t *
read_image(const char *filename,
           int *w,
           int *h)
{
#ifdef OFPIX_DOUBLE
    ofpix_t *f = iio_read_image_double(filename, w, h);
#else
    ofpix_t *f = iio_read_image_float(filename, w, h);
#endif

    if (!f) {
        fprintf(stderr, "ERROR: could not read image from file "
                "\"%s\"\n", filename);
    }

    return f;
}


/**
 *
 *  Main program:
 *   This program reads the following parameters from the console and
 *   then computes the optical flow:
 *   I_1                  Previous image to I0
 *   I0          first image
 *   I1          second image
 *   I0_Smoothed Image for using with function g
 *   out         name of the output flow field
 *   outOcc      name of the output occlusion map
 *   nprocs      number of threads to use (OpenMP library)
 *   tauEta      Time step in the primal-dual scheme for eta variable
 *   tauChi      Time step in the primal-dual scheme for chi variable
 *   lambda      Data term weight parameter
 *   alpha       Length term weight parameter (in the occlusion region)
 *   beta		  Negative divergence data Term
 *   theta       tightness parameter
 *   nscales     number of scales in the pyramidal structure
 *   zfactor     downsampling factor for creating the scales
 *   nwarps      number of warps per scales
 *   epsilon     stopping criterion threshold for the iterative process
 *   verbose     switch on/off messages
 *
 */
int
main(int argc,
     char *argv[])
{
    if (argc < 3) {
        fprintf(stderr, "Usage: %s I_1 I0 I1 [I0_Smoothed out "
                //              0   1   2  3      4        5
                "outOcc nproc lambda alpha beta theta nscales zfactor nwarps epsilon "
                // 6      7      8   9    10    11        12     13     14     15
                "verbose  ]\n", *argv);

        // 16
        return EXIT_FAILURE;
    }
    // Variable Declaration
    ofpix_t *u1 = 0, *u2 = 0; //field
    ofpix_t *chi = 0;          //Occlussion map
    ofpix_t *I_1 = 0, *I0 = 0, *I1 = 0;  // Previous (I_1), current (I0) and next image (I1)
    ofpix_t *filtI0 = 0;                   //Filtered image used in function g
    int nx_1, ny_1, nx, ny, nx1, ny1, nxf, nyf; //Image sizes

    //read the parameters
    int i = 1;
    const char* image_1_name = argv[i]; i++; //1
    const char* image1_name  = argv[i]; i++; //2
    const char* image2_name  = argv[i]; i++; //3
    const char* image1_Smooth_name = (argc > i) ? argv[i] : argv[2]; i++; //4 If there is no I0_Smoothed, then it will be I0
    const char* outfile = (argc > i) ? argv[i] : PAR_DEFAULT_OUTFLOW;       i++; //5
    const char* outOccFile = (argc > i) ? argv[i] : PAR_DEFAULT_OUT_OCC;       i++; //6
    int nproc   = (argc > i) ? atoi(argv[i]) : PAR_DEFAULT_NPROC;   i++; //7
    double lambda  = (argc > i) ? atof(argv[i]) : PAR_DEFAULT_LAMBDA;  i++; //8
    double alpha   = (argc > i) ? atof(argv[i]) : PAR_DEFAULT_ALPHA;   i++; //9
    double betaW   = (argc > i) ? atof(argv[i]) : PAR_DEFAULT_BETA;    i++; //10
    double theta   = (argc > i) ? atof(argv[i]) : PAR_DEFAULT_THETA;   i++; //11
    int nscales = (argc > i) ? atoi(argv[i]) : PAR_DEFAULT_NSCALES; i++; //12
    double zfactor = (argc > i) ? atof(argv[i]) : PAR_DEFAULT_ZFACTOR; i++; //13
    int nwarps  = (argc > i) ? atoi(argv[i]) : PAR_DEFAULT_NWARPS;  i++; //14
    double epsilon = (argc > i) ? atof(argv[i]) : PAR_DEFAULT_EPSILON; i++; //15
    int verbose = (argc > i) ? atoi(argv[i]) : PAR_DEFAULT_VERBOSE; i++; //16

    //check parameters
    if (nproc < 0) {
        nproc = PAR_DEFAULT_NPROC;
        fprintf(stderr, "warning: "
                "nproc changed to %d\n", nproc);
    }
    if (lambda <= 0) {
        lambda = PAR_DEFAULT_LAMBDA;
        fprintf(stderr, "warning: "
                "lambda changed to %g\n", lambda);
    }
    if (alpha <= 0) {
        alpha = PAR_DEFAULT_ALPHA;
        fprintf(stderr, "warning: "
                "alpha changed to %g\n", alpha);
    }
    if (betaW <= 0) {
        betaW = PAR_DEFAULT_BETA;
        fprintf(stderr, "warning: "
                "beta changed to %g\n", betaW);
    }
    if (theta <= 0) {
        theta = PAR_DEFAULT_THETA;
        if (verbose) {
            fprintf(stderr, "warning: "
                    "theta changed to %g\n", theta);
        }
    }
    if (nscales <= 0) {
        nscales = PAR_DEFAULT_NSCALES;
        fprintf(stderr, "warning: "
                "nscales changed to %d\n", nscales);
    }
    if ( (zfactor <= 0) || (zfactor >= 1) ) {
        zfactor = PAR_DEFAULT_ZFACTOR;
        fprintf(stderr, "warning: "
                "zfactor changed to %g\n", zfactor);
    }
    if (nwarps <= 0) {
        nwarps = PAR_DEFAULT_NWARPS;
        fprintf(stderr, "warning: "
                "nwarps changed to %d\n", nwarps);
    }
    if (epsilon <= 0) {
        epsilon = PAR_DEFAULT_EPSILON;
        fprintf(stderr, "warning: "
                "epsilon changed to %f\n", epsilon);
    }


#ifdef _OPENMP
    if (nproc > 0) {
        omp_set_num_threads(nproc);
    }
#endif

    // read the input images
    I_1    = read_image(image_1_name, &nx_1, &ny_1);
    I0     = read_image(image1_name, &nx, &ny);
    I1     = read_image(image2_name, &nx1, &ny1);
    filtI0 = read_image(image1_Smooth_name, &nxf, &nyf);

    if ( (nx == nx_1) && (nx == nx1) && (nx == nxf) &&
         ( ny == ny_1) && ( ny == ny1) && ( ny == nyf) ) {
        //Set the number of scales according to the size of the
        //images.  The value N is computed to assure that the smaller
        //images of the pyramid don't have a size smaller than 16x16

        const int N = floor( log( (float)MIN(nx, ny) / 16.0 ) / log(1. / zfactor) ) + 1;

        if (N < nscales) {
            nscales = N;
        }

        if (verbose) {
            fprintf(stderr,
                    " nproc=%d   \n lambda=%f \n alpha=%f \n"
                    " beta=%f \n theta=%f \n nscales=%d \n zfactor=%f\n nwarps=%d \n epsilon=%g\n",
                    nproc, lambda, alpha, betaW, theta, nscales,
                    zfactor, nwarps, epsilon);
        }

        //allocate memory for the flow
        u1  = new ofpix_t [nx * ny];
        u2  = new ofpix_t [nx * ny];

        //and the occlusion map
        chi = new ofpix_t [nx * ny];


        for (int i = 0; i < nx * ny; i++) {
            chi[i] = 0.0;
            u1[i]  = 0.0;
            u2[i]  = 0.0;
        }


        //compute the optical flow
        Dual_TVL1_optic_flow_multiscale(
            I_1, I0, I1, filtI0, u1, u2, chi, nx, ny, lambda, alpha, betaW,  theta,
            nscales, zfactor, nwarps, epsilon, verbose);

        //save the optical flow

        float *f = new float[nx * ny * 2];
        for (int i = 0; i < nx * ny; i++) {
            f[2 * i] = u1[i];
            f[2 * i + 1] = u2[i];
        }
        iio_save_image_float_vec( (char *)outfile, f, nx, ny, 2 );

        free(f);

        //save the occlusions
/*		int iv=0;
        FILE * fid=fopen(outOccFile, "w");
        for (int i=0; i<ny; i++)
        {
            for (int j=0; j<nx; j++)
            {
                fprintf(fid, " %f", chi[iv]);
                iv++;
            }
            fprintf(fid, " \n");
        }
        fclose(fid);*/
        //iio_save_image_double((char *)outOccFile, chi, nx, ny);

        float *fOcc = (float *)malloc(sizeof(float) * nx * ny );
        for (int i = 0; i < nx * ny; i++) {
            fOcc[i] = (float)chi[i] * 255;   //Avoid the cast!
        }
        iio_save_image_float( (char *)outOccFile, fOcc, nx, ny );

        free(fOcc);
    }
    //delete allocated memory
    free(I_1);
    free(I0);
    free(I1);
    delete [] u1;
    delete [] u2;
    free(filtI0);
    delete [] chi;

    return EXIT_SUCCESS;
} // main

