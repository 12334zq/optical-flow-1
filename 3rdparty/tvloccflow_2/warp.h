#ifndef WARP_H
#define WARP_H

#include <omp.h>
#include "bicubic_interpolation.h"


void warping(
    const double *input, const double *u, const double *v,
    double *output, const int nx, const int ny
    );
#endif
