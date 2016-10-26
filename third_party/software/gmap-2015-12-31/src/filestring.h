/* $Id: filestring.h 159426 2015-02-25 00:35:16Z twu $ */
#ifndef FILESTRING_INCLUDED
#define FILESTRING_INCLUDED

#ifdef USE_MPI
#include <mpi.h>
#include "mpidebug.h"
#endif

#include <stdio.h>
#include "samflags.h"

#define FPRINTF Filestring_put
#define PUTC Filestring_putc


#define T Filestring_T
typedef struct T *T;

extern int
Filestring_id (T this);
extern void
Filestring_set_split_output (T this, int split_output);
extern SAM_split_output_type
Filestring_split_output (T this);
extern T
Filestring_new (int id);
extern void
Filestring_free (T *old);
extern void
Filestring_stringify (T this);
extern void
Filestring_print (
#ifdef USE_MPI
		  MPI_File fp,
#else
		  FILE *fp,
#endif
		  T this);
extern char *
Filestring_get (int *strlength, T this);
extern void
Filestring_put (T this, const char *format, ...);
extern void
Filestring_putc (char c, T this);
extern void
Filestring_puts (T this, char *string, int strlength);

#ifdef USE_MPI
extern char *
Filestring_extract (int *strlength, T this);
extern void
Filestring_send (T this, int dest, int tag, MPI_Comm comm);
extern char *
Filestring_recv (int *strlength, int source, int tag, MPI_Comm comm);
#endif


#undef T
#endif


