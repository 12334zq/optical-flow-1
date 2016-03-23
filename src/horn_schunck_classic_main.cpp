// This program is free software: you can use, modify and/or redistribute it
// under the terms of the simplified BSD License. You should have received a
// copy of this license along this program. If not, see
// <http://www.opensource.org/licenses/bsd-license.html>.
//
// Copyright 2012, Enric Meinhardt Llopis <enric.meinhardt@cmla.ens-cachan.fr>
// All rights reserved.


#include <stdarg.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "iio.h"
#include "horn_schunck.h"


// main function for testing the Horn-Schunck method from the command line
int
main(int argc,
     char *argv[])
{
    if ( (argc != 6) && (argc != 7) ) {
        return fprintf(stderr, "usage:\n\t%s niter alpha a b f\n", *argv);
    }
    //                                        0 1     2     3 4 5
    int niter = atoi(argv[1]);
    double alpha = atof(argv[2]);
    char *filename_a = argv[3];
    char *filename_b = argv[4];
    char *filename_f = argv[5];
    int w, h, ww, hh;
#ifdef OFPIX_DOUBLE
    ofpix_t *a = iio_read_image_double(filename_a, &w, &h);
    ofpix_t *b = iio_read_image_double(filename_b, &ww, &hh);
#else
    ofpix_t *a = iio_read_image_float(filename_a, &w, &h);
    ofpix_t *b = iio_read_image_float(filename_b, &ww, &hh);
#endif
    if ( (w != ww) || (h != hh) ) {
        return fprintf(stderr, "input images size mismatch\n");
    }
    ofpix_t *u = new ofpix_t [2 * w * h];
    ofpix_t *v = u + w * h;
    hs(u, v, a, b, w, h, niter, alpha);
#ifdef OFPIX_DOUBLE
    iio_save_image_double_split(filename_f, u, w, h, 2);
#else
    iio_save_image_float_split(filename_f, u, w, h, 2);
#endif
    delete [] u;

    return EXIT_SUCCESS;
}

