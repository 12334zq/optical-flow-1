/*
 * solvers.h
 *
 *  Created on: 03/12/2013
 *      Author: juanfran
 */

#ifndef SOLVERS_H_
#define SOLVERS_H_

#include "operators.h"
#include "TV_ROF_Box.h"
#include "constants.h"



/*
 * Function name: Solver_wrt_v
 * Author : Coloma Ballester 2011
 * Modified by: J.F. Garamendi 2013
 * Parameters
	u1     : (in)  Vectorized matrix. x component of the optical flow. The initial (input) values are used as initialitation
	u2     : (in)  Vectorized Matrix. y component of the optical flow. The initial (input) values are used as initialitation
	v1     : (out) Vectorized Matrix. x component of the  variable v.
	v2     : (out) Vectorized Matrix. y component of the  variable v.
	chi    : (in)  Vectorized Matrix. Occlusion map
	I1wx   : (in)  Vectorized Matrix.
	I1wy   : (in)  Vectorized Matrix.
	I_1wx  :
	I_1wy  :
	rho1_c : (in)  Vectorized matrix. Constant part of the rho function corresponding to the unoccluded pixels
	rho3_c : (in)  Vectorized Matrix. Constant part of the rho function corresponding to the occluded pixels
	Vfwd_1    : (out) Vectorized Matrix. x component of the auxilary variable corresponding with the values to use at the unoccluded pixels.
	Vfwd_2    : (out) Vectorized Matrix. y component of the auxilary variable corresponding with the values to use at the unoccluded pixels.
	Vbck_1    : (out) Vectorized Matrix. x component of the auxilary variable corresponding with the values to use at the occluded pixels.
	Vbck_2    : (out) Vectorized Matrix. y component of the auxilary variable corresponding with the values to use at the occluded pixels.
	grad1  : (in)  Vectorized Matrix.
	grad3  : (in)  Vectorized Matrix.
	alpha  : (in)  Scalar corresponding to alpha parameter.
	theta  : (in)  Scalar corresponding to theta parameter.
	lambda : (in)  Scalar corresponding to lambda parameter.
	nx     : (in)  Scalar integer. image width
	ny     : (in)  Scalar integer. image height
 * Output:
 * 		VOID
 * Description: This function performs the minimization with respect to the  variable v. The minimization is done by thresholding.
 * There exist to variabes v1 (with components Vfwd_1,Vfwd_2)  and v3 (with components v13, v14) used to carry the information of the
 * variable v. v1 will be used  in those pixels withoout occlusion.
 * v3 will be used in the occluded pixels.
 *
 * 		Reference:
 *  	[1] C. Ballester, L. Garrido, V. Lazcano, V. Caselles.
 *      "A TV-L1 Optical Flow Method with Occlusion Detection".
 *      Lecture Notes in Computer Science. 7476:31-40, 2012
 *
 */
void Solver_wrt_v(double *u1, double *u2,
	    double *v1, double *v2,
	    double *chi,
	    const double *I1wx, const double *I1wy,
	    const double *I_1wx, const double *I_1wy,
	    const double *rho1_c, const double *rho3_c,
	    double *Vfwd_1, double *Vfwd_2, double *Vbck_1, double *Vbck_2,
	    const double *grad1, const double *grad3,
		const double alpha, const double theta, const double lambda,
		const int nx, const int ny);


/*
 * Function name: Solver_wrt_u
 * Author: J.F. Garamendi 2013
 * Parameters
	u1     : (in/out) Vectorized matrix. x component of the optical flow. The initial (input) values are used as initialitation
	u2     : (in/out) Vectorized Matrix. y component of the optical flow. The initial (input) values are used as initialitation
	v1     : (in) Vectorized Matrix. x component of the auxilary variable.
	v2     : (in) Vectorized Matrix. y component of the auxilary variable.
	chi    : (in) Vectorized Matrix. Occlusion map
	g      : (in) Vectorized Matrix. Values of function g applied to a filtered version (given by command line) of I0
	theta  : (in) Scalar corresponding to parameter theta
	beta   : (in) Scalar corresponding to parameter beta
	nx     : (in) Scalar integer. image width
	ny     : (in) Scalar integer. image height
	nIter  :
 * Output:
 * 		VOID
 * Description: This function performs the minimization with respect to the  variable flow field u. It prepares the data for calling
 * 	the function that implements the minimization of the (modified) Rudin-Osher-Fatemi problem.
 *
 * 		Reference:
 *  	[1] C. Ballester, L. Garrido, V. Lazcano, V. Caselles.
 *      "A TV-L1 Optical Flow Method with Occlusion Detection".
 *      Lecture Notes in Computer Science. 7476:31-40, 2012
 */
void Solver_wrt_u(double *u1, double *u2,
	    const double *v1, const double *v2,
	    const double *chi,
	    const double *g,
		const double theta, const double beta,
		const int nx, const int ny);


/*
 * Function name: Solver_wrt_v
 * Author : Coloma Ballester 2011
 * Modified by: J.F. Garamendi 2013
 * Parameters
	u1     : (in)  Vectorized matrix. x component of the optical flow. The initial (input) values are used as initialitation
	u2     : (in)  Vectorized Matrix. y component of the optical flow. The initial (input) values are used as initialitation
	chi    : (in)  Vectorized Matrix. Occlusion map
	I1wx   : (in)  Vectorized Matrix.
	I1wy   : (in)  Vectorized Matrix.
	I_1wy  : (in)  Vectorized Matrix.
	rho1_c : (in)  Vectorized matrix. Constant part of the rho function corresponding to the unoccluded pixels
	rho3_c : (in)  Vectorized Matrix. Constant part of the rho function corresponding to the occluded pixels
	Vfwd_1    : (out) Vectorized Matrix. x component of the auxilary variable corresponding with the values to use at the unoccluded pixels.
	Vfwd_2    : (out) Vectorized Matrix. y component of the auxilary variable corresponding with the values to use at the unoccluded pixels.
	Vbck_1    : (out) Vectorized Matrix. x component of the auxilary variable corresponding with the values to use at the occluded pixels.
	Vbck_2    : (out) Vectorized Matrix. y component of the auxilary variable corresponding with the values to use at the occluded pixels.
	g      : (in)  Vectorized Matrix. Values of function g applied to a filtered version (given by command line) of I0
	lambda : (in)  Scalar>0 corresponding to lambda parameter.
	theta  : (in)  Scalar>0 corresponding to theta parameter.
	alpha  : (in)  Scalar>0 corresponding to alpha parameter.
	beta   : (in)  Scalar>0 corresponding to beta parameter
	grad1  : (in)  Vectorized Matrix.
	grad3  : (in)  Vectorized Matrix.
	tau_chi: (in)  Scalar >0. Step size for solving w.r.t chi in the primal-dual scheme
	tau_eta: (in)  Scalar >0. Step size for solving w.r.t chi in the primal-dual scheme
	nx     : (in)  Scalar integer. image width
	ny     : (in)  Scalar integer. image height
 * Output:
 * 		VOID
 * Description: This function performs the minimization with respect to the  variable chi. The minimization is done by a primal-dual scheme,
 * using a v1 (with components Vfwd_1,Vfwd_2)  and v3 (with components v13, v14) respectively for the non-occluded and occluded pixels,
 *
 * 		Reference:
 *  	[1] C. Ballester, L. Garrido, V. Lazcano, V. Caselles.
 *      "A TV-L1 Optical Flow Method with Occlusion Detection".
 *      Lecture Notes in Computer Science. 7476:31-40, 2012
 *
 */
void Solver_wrt_chi(const double *u1, const double *u2,
		double *chi,
		const double *I1wx, const double *I1wy,
		const double *I_1wx, const double *I_1wy,
		const double *rho1_c, const double *rho3_c,
		const double *Vfwd_1, const double *Vfwd_2, const double *Vbck_1, const double *Vbck_2,
		const double *g,
		const double lambda,
		const double theta,
		const double alpha,
		const double beta,
		const double tau_chi,
		const double tau_eta,
		const int nx, const int ny);

#endif /* SOLVERS_H_ */
