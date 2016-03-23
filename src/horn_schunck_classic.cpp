// This program is free software: you can use, modify and/or redistribute it
// under the terms of the simplified BSD License. You should have received a
// copy of this license along this program. If not, see
// <http://www.opensource.org/licenses/bsd-license.html>.
//
// Copyright 2012, Enric Meinhardt Llopis <enric.meinhardt@cmla.ens-cachan.fr>
// All rights reserved.

#include "horn_schunck.h"

#include <stdarg.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "of.h"

// the type of the "getpixel" function
typedef ofpix_t (*extension_operator_float)(ofpix_t *, int, int, int, int);

// getpixel, with neumann boundary conditions
static
ofpix_t
extend_float_image_constant(ofpix_t *x,
                            int w,
                            int h,
                            int i,
                            int j)
{
    if (i < 0) {
        i = 0;
    }
    if (j < 0) {
        j = 0;
    }
    if (i >= w) {
        i = w - 1;
    }
    if (j >= h) {
        j = h - 1;
    }

    return x[j * w + i];
}

// compute the gradient and temporal derivative of the input image pair
static
void
compute_input_derivatives(ofpix_t *Ex,
                          ofpix_t *Ey,
                          ofpix_t *Et,
                          ofpix_t *a,
                          ofpix_t *b,
                          int w,
                          int h)
{
    extension_operator_float p = extend_float_image_constant;

    for (int j = 0; j < h; j++) {
        for (int i = 0; i < w; i++) {
            Ey[j * w + i] = (1.0 / 4) * ( p(a, w, h, i, j + 1) - p(a, w, h, i, j)
                                          + p(a, w, h, i + 1, j + 1) - p(a, w, h, i + 1, j)
                                          + p(b, w, h, i, j + 1) - p(b, w, h, i, j)
                                          + p(b, w, h, i + 1, j + 1) - p(b, w, h, i + 1, j) );
            Ex[j * w + i] = (1.0 / 4) * ( p(a, w, h, i + 1, j) - p(a, w, h, i, j)
                                          + p(a, w, h, i + 1, j + 1) - p(a, w, h, i, j + 1)
                                          + p(b, w, h, i + 1, j) - p(b, w, h, i, j)
                                          + p(b, w, h, i + 1, j + 1) - p(b, w, h, i, j + 1) );
            Et[j * w + i] = (1.0 / 4) * ( p(b, w, h, i, j) - p(a, w, h, i, j)
                                          + p(b, w, h, i + 1, j) - p(a, w, h, i + 1, j)
                                          + p(b, w, h, i, j + 1) - p(a, w, h, i, j + 1)
                                          + p(b, w, h, i + 1, j + 1) - p(a, w, h, i + 1, j + 1) );
        }
    }
}

// compute a local average of a function "u"
static
void
compute_bar(ofpix_t *ubar,
            ofpix_t *u,
            int w,
            int h)
{
    extension_operator_float p = extend_float_image_constant;

    for (int j = 0; j < h; j++) {
        for (int i = 0; i < w; i++) {
            ubar[j * w + i] = (1.0 / 6) * ( p(u, w, h, i - 1, j) + p(u, w, h, i + 1, j)
                                            + p(u, w, h, i, j - 1) + p(u, w, h, i, j + 1) )
                              + (1.0 / 12) * ( p(u, w, h, i - 1, j - 1) + p(u, w, h, i + 1, j - 1)
                                               + p(u, w, h, i - 1, j + 1) + p(u, w, h, i + 1, j + 1) );
        }
    }
}

// compute a sigle iteration of the classical Horn-Schunck method
static
void
hs_iteration(ofpix_t *u,
             ofpix_t *v,
             ofpix_t *Ex,
             ofpix_t *Ey,
             ofpix_t *Et,
             int w,
             int h,
             double alpha)
{
    ofpix_t *ubar = new ofpix_t[w * h];
    ofpix_t *vbar = new ofpix_t[w * h];

    compute_bar(ubar, u, w, h);
    compute_bar(vbar, v, w, h);
    for (int i = 0; i < w * h; i++) {
        ofpix_t t = Ex[i] * ubar[i] + Ey[i] * vbar[i] + Et[i];
        t /= alpha * alpha + Ex[i] * Ex[i] + Ey[i] * Ey[i];
        u[i] = ubar[i] - Ex[i] * t;
        v[i] = vbar[i] - Ey[i] * t;
    }
    delete [] ubar;
    delete [] vbar;
}

// run n iterations of the classical Horn-Schunck method
void
hs(ofpix_t *u,
   ofpix_t *v,
   ofpix_t *a,
   ofpix_t *b,
   int w,
   int h,
   int n,
   double alpha)
{
    ofpix_t *gx = new ofpix_t[w * h];
    ofpix_t *gy = new ofpix_t[w * h];
    ofpix_t *gt = new ofpix_t[w * h];

    compute_input_derivatives(gx, gy, gt, a, b, w, h);
    for (int i = 0; i < w * h; i++) {
        u[i] = v[i] = 0;
    }
    for (int i = 0; i < n; i++) {
        hs_iteration(u, v, gx, gy, gt, w, h, alpha);
    }
    delete [] gx;
    delete [] gy;
    delete [] gt;
}

