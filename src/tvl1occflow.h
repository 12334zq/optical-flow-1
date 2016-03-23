/*
 * tvl1OccFlow.h
 *
 *  Created on: 28/11/2013
 *      Author: juanfran
 */

#ifndef OF_TVL1OCCFLOW_H
#define OF_TVL1OCCFLOW_H

#include "of.h"

/*
 * Function name: Dual_TVL1_optic_flow
 * Author: ??
 * Modified by: J.F. Garamendi 2013
 * Parameters
   I_1    : (in) Previous frame to source Image
   I0     : (in) source image
   I1     :  (in)target image
   filtI0 : (in) Image used for computing the g function
   u1     : (in/out) x component of the optical flow. The initial values are used as initialitation
   u2     : (in/out) y component of the optical flow. The initial values are used as initialitation
   chi    : (in/out) Occlusion map
   nx     : (in) image width
   ny     : (in) image height
   tau_eta: (in) time step for the primal-dual minimization on Chi,
   tau_chi: (in) time step for the primal-dual minimization over Chi,
   lambda : (in) weight parameter for the data term
   alpha  : (in) weight paramenter for chi|v|^2
   beta   : (in) weight parameter for chi*div(u)
   theta  : (in) weight parameter for (u - v)^2
   warps  : (in) number of warpings per scale
   epsilon: (in) tolerance for numerical convergence
   verbose: (in) enable/disable the verbose mode
   max_iterations_u:  //Quitar xq no es del modelo de optic flow
   ext_max_iterations:   //Quitar xq no es del modelo de optic flow
 * Output:
 *              VOID
 * Description: This function computes the optical flow [1] between Images I0 (source image) and I1 (target image)
 *              using the previous image I_1 (previous to I0) as image for estimate the values in occluded areas. Also
 *              it computes the occlussion map.
 *              The minimization of functional is done using an alternating scheme for u,v, and the occlussion map. When
 *              is minimized respect to
 *                      + v (artifical variable): Threshold
 *                      + u (flow field)        : Minimization of ROF functional using [2]
 *                      + Chi (Occlusion map)   : Primal-Dual
 *
 *              Reference:
 *      [1] C. Ballester, L. Garrido, V. Lazcano, V. Caselles.
 *      "A TV-L1 Optical Flow Method with Occlusion Detection".
 *      Lecture Notes in Computer Science. 7476:31-40, 2012
 *
 *
 *              Details on the total variation minimization scheme can be found in:
 *      [2] J.F. Garamendi, F.J. Gaspar, N. Malpica, E. Schiavi.
 *      "Box Relaxation Schemes in Staggered Discretizations for the Dual Formulation
 *      of Total Variation Minimization"
 *      IEEE Transactions on Image Processing. 22(5):2030-2043, 2013.
 *
 *
 */
void Dual_TVL1_optic_flow(ofpix_t *I_1,         // Previous frame to source Image
                          ofpix_t *I0, // source image
                          ofpix_t *I1, // target image
                          ofpix_t *filtI0, //Image used for computing the g function
                          ofpix_t *u1, // x component of the optical flow
                          ofpix_t *u2, // y component of the optical flow
                          ofpix_t *chi, // Occlusion map
                          const int nx, // image width
                          const int ny, // image height
                          const double lambda, // weight parameter for the data term
                          const double alpha, // weight paramenter for chi|v|^2
                          const double beta, // weight parameter for chi*div(u)
                          const double theta, // weight parameter for (u - v)²
                          const int warps, // number of warpings per scale
                          const double epsilon, // tolerance for numerical convergence
                          const bool verbose // enable/disable the verbose mode
                          );

/*
 * Function name: Dual_TVL1_optic_flow_multiscale
 * Author: ??
 * Modified by: J.F. Garamendi 2013
 * Parameters
   I_1    : (in) Previous frame to source Image
   I0     : (in) source image
   I1     :  (in)target image
   filtI0 : (in) Image used for computing the g function
   u1     : (in/out) x component of the optical flow. The initial values are used as initialitation
   u2     : (in/out) y component of the optical flow. The initial values are used as initialitation
   chi    : (in/out) Occlusion map
   nxx     : (in) image width
   nyy     : (in) image height
   tau_eta: (in) time step for the primal-dual minimization on Chi (primal equation)
   tau_chi: (in) time step for the primal-dual minimization on Chi (dual equation)
   lambda : (in) weight parameter for the data term
   alpha  : (in) weight paramenter for chi|v|^2
   beta   : (in) weight parameter for chi*div(u)
   theta  : (in) weight parameter for (u - v)^2
   nscales: (in) number of scales in the pyramid
   zfactor: (in) factor for building the image pyramid
   warps  : (in) number of warpings per scale
   epsilon: (in) tolerance for numerical convergence
   verbose: (in) enable/disable the verbose mode
 * Output:
 *              VOID
 * Description: This function builds a multiscale pyramid of images and it calls for each pyramid level to Dual_TVL1_optic_flow.
 *                              The number of used scales is given by the nscales parameter, and the reduction factor is given by zfactor.
 */
void Dual_TVL1_optic_flow_multiscale(ofpix_t *I_1, // Previous frame to source image
                                     ofpix_t *I0, // source image
                                     ofpix_t *I1, // target image
                                     ofpix_t *filtI0, // Smoothed version of I0 for computing the g function
                                     ofpix_t *u1, // x component of the optical flow
                                     ofpix_t *u2, // y component of the optical flow
                                     ofpix_t *chi, // Occlusion map
                                     const int nxx, // image width
                                     const int nyy, // image height
                                     const double lambda, // weight parameter for the data term
                                     const double alpha, // weight paramenter for chi|v|^2
                                     const double beta, // weight parameter for chi*div(u)
                                     const double theta, // weight parameter for (u - v)²
                                     const int nscales, // number of scales
                                     const double zfactor, // factor for building the image piramid
                                     const int warps, // number of warpings per scale
                                     const double epsilon, // tolerance for numerical convergence
                                     const bool verbose // enable/disable the verbose mode
                                     );

#endif /* TVL1OCCFLOW_H_ */
