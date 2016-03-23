/*
 * solvers.h
 *
 *  Created on: 03/12/2013
 *      Author: juanfran
 */

#ifndef TVL1OCCFLOW_SOLVERS_H
#define TVL1OCCFLOW_SOLVERS_H

#include "of.h"
#include "operators.h"
#include "tvl1occflow_tv_rof_box.h"
#include "tvl1occflow_constants.h"



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
void Solver_wrt_v(ofpix_t *u1, ofpix_t *u2,
	    ofpix_t *v1, ofpix_t *v2,
	    ofpix_t *chi,
	    const ofpix_t *I1wx, const ofpix_t *I1wy,
	    const ofpix_t *I_1wx, const ofpix_t *I_1wy,
	    const ofpix_t *rho1_c, const ofpix_t *rho3_c,
	    ofpix_t *Vfwd_1, ofpix_t *Vfwd_2, ofpix_t *Vbck_1, ofpix_t *Vbck_2,
	    const ofpix_t *grad1, const ofpix_t *grad3,
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
void Solver_wrt_u(ofpix_t *u1, ofpix_t *u2,
	    const ofpix_t *v1, const ofpix_t *v2,
	    const ofpix_t *chi,
	    const ofpix_t *g,
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
void Solver_wrt_chi(const ofpix_t *u1, const ofpix_t *u2,
		ofpix_t *chi,
		const ofpix_t *I1wx, const ofpix_t *I1wy,
		const ofpix_t *I_1wx, const ofpix_t *I_1wy,
		const ofpix_t *rho1_c, const ofpix_t *rho3_c,
		const ofpix_t *Vfwd_1, const ofpix_t *Vfwd_2, const ofpix_t *Vbck_1, const ofpix_t *Vbck_2,
		const ofpix_t *g,
		const double lambda,
		const double theta,
		const double alpha,
		const double beta,
		const double tau_chi,
		const double tau_eta,
		const int nx, const int ny);

#endif /* SOLVERS_H_ */
