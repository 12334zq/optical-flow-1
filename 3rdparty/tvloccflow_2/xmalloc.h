#ifndef XMALLOC_H_
#define XMALLOC_H_
#include <stdio.h>
#include <stdlib.h>

//J.F. Garamendi (7-5-2014) We should avoid the use of this function because if for some
//reason, there is no memory available, the program aborts without any chance
// of recovering the status.

/* this function is like "malloc", but it returns always a valid pointer */
void *xmalloc(size_t size);

#endif /* XMALLOC_H_ */
