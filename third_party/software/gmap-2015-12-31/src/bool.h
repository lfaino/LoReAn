/* $Id: bool.h 155282 2014-12-12 19:42:54Z twu $ */
#ifndef BOOL_INCLUDED
#define BOOL_INCLUDED

/* typedef enum{false,true} bool; */

typedef unsigned char bool;
#ifdef USE_MPI
#define MPI_BOOL_T MPI_UNSIGNED_CHAR
#endif

#define false 0
#define true 1



#endif
