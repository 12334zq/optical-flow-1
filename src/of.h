#ifndef OF_OF_H
#define OF_OF_H

#define OFPIX_DOUBLE // uncomment to use arrays of double instead of float everywhere

#ifdef OFPIX_DOUBLE
typedef double ofpix_t;
#else
typedef float ofpix_t;
#endif

#endif
