/*
 * constants.h
 *
 *  Created on: 28/01/2014
 *      Author: juanfran
 */

#ifndef CONSTANTS_H_
#define CONSTANTS_H_

//Default program arguments
#define PAR_DEFAULT_OUTFLOW "flow.flo"
#define PAR_DEFAULT_OUT_OCC "occlusions.png"
#define PAR_DEFAULT_NPROC   1
#define PAR_DEFAULT_LAMBDA  0.15
#define PAR_DEFAULT_ALPHA  0.01
#define PAR_DEFAULT_BETA  0.15
#define PAR_DEFAULT_THETA   0.3
#define PAR_DEFAULT_NSCALES 100
#define PAR_DEFAULT_ZFACTOR 0.5
#define PAR_DEFAULT_NWARPS  2
#define PAR_DEFAULT_EPSILON 0.01
#define PAR_DEFAULT_VERBOSE 0

#define EXT_MAX_ITERATIONS 20
#define OMEGA 1.25

#define IS_ZERO 1E-10
#define THR_CHI 0.75
#define MAX_ITERATIONS_CHI 100
#define PRESMOOTHING_SIGMA 0.8

#define MAX_ITERATIONS 20
#define G_CHOICE 2
#define G_FACTOR 0.05

//Numerical parameters
#define MAX_ITERATIONS_U 10
#define TAU_ETA 0.15
#define TAU_CHI 0.15

#endif /* CONSTANTS_H_ */
