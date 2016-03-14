// This program is free software: you can use, modify and/or redistribute it
// under the terms of the simplified BSD License. You should have received a
// copy of this license along this program. If not, see
// <http://www.opensource.org/licenses/bsd-license.html>.
//
// Copyright (C) 2012, Coloma Ballester <coloma.ballester@upf.edu>
// Copyright (C) 2013-2014 J. F. Garamendi <jf.garamendi@upf.edu>
// All rights reserved.

#include "warp.h"

void warping(const double *input, const double *u, const double *v,
		double *output, const int nx, const int ny) {

	for (int j = 0; j < ny; j++)
		for (int i = 0; i < nx; i++) {
			const double uu = (double) (i + u[i + nx * j]);
			const double vv = (double) (j + v[i + nx * j]);

			//if it is outside the limits the flow, it wil be zero
			if ((uu < 0) || (uu > (nx - 1)) || (vv < 0) || (vv > (ny - 1))) {
				output[i + nx * j] = 0;
			} else {
				//We interpolate
				output[i + nx * j] = me_interpolate_bicubic(input, nx, ny, uu,
						vv);
			}
		}
}
