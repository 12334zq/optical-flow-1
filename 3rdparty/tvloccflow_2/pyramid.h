/*
 * zoom.h
 *
 *  Created on: 24/01/2014
 */

#ifndef ZOOM_H_
#define ZOOM_H_

void zoom_size(int nx,      // width of the orignal image
		int ny,      // height of the orignal image
		int *nxx,    // width of the zoomed image
		int *nyy,    // height of the zoomed image
		double factor // zoom factor between 0 and 1
);

void DownSampling(const double *I,    // input image
		double *Iout,       // output image
		const int nx,      // image width
		const int ny,      // image height
		const double factor // zoom factor between 0 and 1
);

void UpSampling(const double *I, // input image
		double *Iout,    // output image
		int nx,         // width of the original image
		int ny,         // height of the original image
		int nxx,        // width of the zoomed image
		int nyy         // height of the zoomed image
);

#endif /* ZOOM_H_ */
