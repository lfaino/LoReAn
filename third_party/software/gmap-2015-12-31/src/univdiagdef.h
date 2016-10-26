/* $Id: univdiagdef.h 166641 2015-05-29 21:13:04Z twu $ */
#ifndef UNIVDIAGDEF_INCLUDED
#define UNIVDIAGDEF_INCLUDED

#include "bool.h"

#define T Univdiag_T
struct T {
  Univcoord_T univdiagonal;	/* Used by sarray-read.c */
  int querystart;
  int queryend;
  int nconsecutive;
  bool nmismatches_known_p;

  int intscore;	/* Used for dynamic programming of diagonals in sarray-read.c */
  int nlinked; /* Used for dynamic programming of diagonals in sarray-read.c */
  struct T *prev; /* Used for dynamic programming of diagonals in sarray-read.c */
};

#undef T
#endif

