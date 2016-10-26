/* $Id: splicestringpool.h 135405 2014-05-07 19:16:09Z twu $ */
#ifndef SPLICESTRINGPOOL_INCLUDED
#define SPLICESTRINGPOOL_INCLUDED

#include "list.h"
#include "types.h"

typedef struct Splicestring_T *Splicestring_T;
struct Splicestring_T {
  Genomecomp_T string;
  Genomecomp_T splicesite;
  Genomecomp_T splicesite_i;
};


#define T Splicestringpool_T
typedef struct T *T;

extern void
Splicestringpool_free (T *old);
extern void
Splicestringpool_free_memory (T this);
extern void
Splicestringpool_report_memory (T this);
extern T
Splicestringpool_new (void);
extern void
Splicestringpool_reset (T this);
extern List_T
Splicestringpool_push (List_T list, T this, Genomecomp_T string, Genomecomp_T splicesite,
		       Genomecomp_T splicesite_i);
extern List_T
Splicestringpool_pop (List_T list, Splicestring_T *x);

#undef T
#endif


