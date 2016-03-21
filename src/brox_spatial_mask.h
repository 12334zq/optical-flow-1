// This program is free software: you can use, modify and/or redistribute it
// under the terms of the simplified BSD License. You should have received a
// copy of this license along this program. If not, see
// <http://www.opensource.org/licenses/bsd-license.html>.
//
// Copyright (C) 2012, Javier Sánchez Pérez <jsanchez@dis.ulpgc.es>
// All rights reserved.

#ifndef BROX_SPATIAL_MASK_H
#define BROX_SPATIAL_MASK_H

/**
 *
 * Compute the second order X derivative
 *
 */
void brox_spatial_Dxx(
    const float *I, //input image
    float *Ixx,     //oputput derivative
    const int nx,   //image width
    const int ny    //image height
                 );


/**
 *
 * Compute the second order Y derivative
 *
 */
void brox_spatial_Dyy(
    const float *I, //input image
    float *Iyy,     //oputput derivative
    const int nx,   //image width
    const int ny    //image height
                 );


/**
 *
 * Compute the second order XY derivative
 *
 */
void brox_spatial_Dxy(
    const float *I, //input image
    float *Ixy,     //oputput derivative
    const int nx,   //image width
    const int ny    //image height
                 );



/**
 *
 * Compute the coefficients of the divergence term
 *
 */
void brox_spatial_psi_divergence(
    const float *psi, //robust functional
    float *psi1,      //coefficients of divergence
    float *psi2,      //coefficients of divergence
    float *psi3,      //coefficients of divergence
    float *psi4,      //coefficients of divergence
    const int nx,     //image width
    const int ny      //image height
                                 );



/**
 *
 * Compute the divergence of the optical flow
 *
 */
void brox_spatial_divergence_u(
    const float *u,    //x component of optical flow
    const float *v,    //y component of optical flow
    const float *psi1, //coefficients of divergence
    const float *psi2, //coefficients of divergence
    const float *psi3, //coefficients of divergence
    const float *psi4, //coefficients of divergence
    float *div_u,      //computed divergence for u
    float *div_v,      //computed divergence for v
    const int nx,      //image width 
    const int ny       //image height
                               );

#endif
