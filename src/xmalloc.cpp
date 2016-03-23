// This program is free software: you can use, modify and/or redistribute it
// under the terms of the simplified BSD License. You should have received a
// copy of this license along this program. If not, see
// <http://www.opensource.org/licenses/bsd-license.html>.
//
// Copyright 2012, Enric Meinhardt Llopis <enric.meinhardt@cmla.ens-cachan.fr>
// All rights reserved.

#include <stdio.h>
#include <stdlib.h>
#include "xmalloc.h"

//J.F. Garamendi (7-5-2014) We should avoid the use of this function because if for some
//reason, there is no memory available, the program aborts without any chance
// of recovering the status.

/* this function is like "malloc", but it returns always a valid pointer */
void *
xmalloc(size_t size)
{
    //Locate memory
    void *p = malloc(size);

    if (!p) {
        //If there is no memory, abort
        exit( fprintf(stderr, "out of memory\n") );
    }

    return p;
}

