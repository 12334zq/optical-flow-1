/*
 * TV_ROF_Box.h
 *
 *  Created on: 03/12/2013
 *      Author: juanfran
 */

#ifndef OF_TVL1OFFFLOW_TV_ROF_BOX_H
#define OF_TVL1OFFFLOW_TV_ROF_BOX_H

#include "of.h"

/*
 * Function name: Scalar_ROF_BoxCellCentered
 * Author: J.F. Garamendi
 * Parameters:
 * 	u         : (in/out) As input, the initial seed, and as ouput the Denoised image
 * 	f         : (in)  Initial image
 * 	initialP1 : (in/out) initial seed for estimating P1 (the first component of the dual variable P)
 initialP2 : (in/out) initial seed for estimating P2 (the second component of the dual variable P)
 g_function: (in) Vectorized matrix representing the function g(x) to apply to the Total Variation Term
 omega:    : (in)Value of omega for omega-relaxation.
 nx        : (in) image width
 ny        : (in) image height
 nIter     : (in) number of iterations
 * Output:
 * 		VOID
 * Description: Implements a slight modified version of the Staggered Box Cell Centered numerical scheme
 *              of [2] for the minimization of the following energy functional:
 *
 *              int_dom{g(x)|grad(u(x))|}dx + int_dom{|u(x)-f(x)|^2}dx
 *
 *              where
 *              int_dom{.}dx is the integral over the square domain of the image
 *              |.| is the euclidean norm
 *              grad(.) is the gradient
 *
 *              The modification of the numerical scheme consist in the implementation of the over-relaxation and
 *              the inclusion of the function g.
 *              Notice that if g(i)==1 and omega==1, the numerical scheme exactly coincides
 *              with BCC in [2]
 *
 *
 *  	[2] J.F. Garamendi, F.J. Gaspar, N. Malpica, E. Schiavi.
 *      "Box Relaxation Schemes in Staggered Discretizations for the Dual Formulation
 *      of Total Variation Minimization"
 *      IEEE Transactions on Image Processing. 22(5):2030-2043, 2013.
 *
 *
 */
void Scalar_ROF_BoxCellCentered(ofpix_t *u, const ofpix_t *f, ofpix_t *initialP1,
		ofpix_t *initialP2, const ofpix_t *g_function, const double lambda,
		const double omega, //Omega-relaxation
		const int nx, const int ny, const int nIter);

#endif /* OF_TVL1OFFFLOW_TV_ROF_BOX_H */
