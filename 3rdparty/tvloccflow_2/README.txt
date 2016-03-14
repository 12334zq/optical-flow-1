TVL1FLOW

A program for optical flow and occlusion estimation based on total variation, the L1 norm and
divergence of the optical flow.

This program is part of an IPOL publication:


Juan Francisco Garamendi <jf.garamendi@upf.edu>. Departament de Tecnologies de la Informació i les Comunicacions Universitat Pompeu Fabra
Coloma Ballester <coloma.ballester@upf.edu>. Departament de Tecnologies de la Informació i les Comunicacions Universitat Pompeu Fabra
Lluis Garrido
Vanel Lazcano Departament de Tecnologies de la Informació i les Comunicacions Universitat Pompeu Fabra
Vicent Caselles Departament de Tecnologies de la Informació i les Comunicacions Universitat Pompeu Fabra


Version 1, released on

This software is distributed under the terms of the GNU license (see file
license.txt)

Required libraries: libpng, libtiff

Compilation instructions: run "make" to produce an executable named "tvl1Occflow"

Usage instructions:

./tvl1Occflow I_1.png I0.png I1.png [I0_Smoothed.png out.flo outOcc.png nproc lambda alpha beta theta nscales zfactor nwarps epsilon verbose  ]

where the parameters between brackets are optional and
I_1.png		Previous image to I0.png
I0.png		first input image
I1.png		second input image
I0_Smoothed	Image for using with function g (default: I0)
out.flo		Filename for the output optical flow (default: flow.flo)
outOcc.png      Filename for the output occlusion map (default: occlusions.png)
nprocs		number of threads to use (OpenMP library) (NPROCS=0, all processors available)
lambda		Data term weight parameter (default: 0.15)
alpha		Length term weight parameter (in the occlusion region) (default 0.01)
beta		weight for the Negative divergence data Term (default: 0.15)
theta		tightness parameter (default: 0.3)
nscales		number of scales in the pyramidal structure (default, maximum number)
zfactor		downsampling factor for creating the scales (default, 0.5)
nwarps 		number of warps per scales (default, 5)
epsilon		stopping criterion threshold for the iterative process (default, 0.01)
verbose		switch on/off (1/0) messages (default, 0, i.e. no verbose)

(note for developers, the default values are located in the constants.h file)


Simple example:

./tvl1Occflow frame09_160x120.png frame10_160x120.png frame11_160x120.png

Example with reasonable default parameters:

./tvl1Occflow ./frame09.png ./frame10.png ./frame11.png ./frame10_filtered.png ./flow.flo ./occlusions.png 1  0.8  0.01 1 0.3 5 0.5 5 1e-09 1

