/* $Id: diagdef.h 157221 2015-01-22 18:38:57Z twu $ */
#ifndef DIAGDEF_INCLUDED
#define DIAGDEF_INCLUDED

#include "bool.h"

#define T Diag_T
struct T {
  Chrpos_T diagonal;
  int querystart;
  int queryend;
  int nconsecutive;
  bool dominatedp;
  double score;
};

#undef T
#endif

