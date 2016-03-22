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
#include "xmalloc.h"
#include "horn_schunck.h"


// main function for testing the Horn-Schunck method from the command line
int main(int argc, char *argv[])
{
	if (argc != 6 && argc != 7)
		return fprintf(stderr,"usage:\n\t%s niter alpha a b f\n",*argv);
	//                                        0 1     2     3 4 5
	int niter = atoi(argv[1]);
	float alpha = atof(argv[2]);
	char *filename_a = argv[3];
	char *filename_b = argv[4];
	char *filename_f = argv[5];
	int w, h, ww, hh;
	float *a = iio_read_image_float(filename_a, &w, &h);
	float *b = iio_read_image_float(filename_b, &ww, &hh);
	if (w != ww || h != hh)
		return fprintf(stderr, "input images size mismatch\n");
	float *u = (float*)xmalloc(2 * w * h * sizeof(float));
	float *v = u + w*h;
	hs(u, v, a, b, w, h, niter, alpha);
	iio_save_image_float_split(filename_f, u, w, h, 2);
	return EXIT_SUCCESS;
}
