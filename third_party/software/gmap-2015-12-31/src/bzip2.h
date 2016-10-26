/* $Id: bzip2.h 157221 2015-01-22 18:38:57Z twu $ */
#ifndef BZIP2_INCLUDED
#define BZIP2_INCLUDED

#include "bool.h"

#define T Bzip2_T
typedef struct T *T;

extern T
Bzip2_new (char *filename);

extern void
Bzip2_free (T *old);

extern int
bzgetc (T this);

extern bool
bzeof (T this);

extern char *
bzgets (T this, char *buffer, int maxlength);

#undef T
#endif
