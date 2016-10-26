/* $Id: diag.h 166641 2015-05-29 21:13:04Z twu $ */
#ifndef DIAG_INCLUDED
#define DIAG_INCLUDED

#include "bool.h"
#include "list.h"
#include "genomicpos.h"
#include "types.h"

#define T Diag_T
typedef struct T *T;

extern T
Diag_new (int querystart, int queryend, Univcoord_T univdiagonal);
extern void
Diag_free (T *old);
extern void
Diag_gc (List_T *list);


extern Chrpos_T
Diag_diagonal (T this);
extern int
Diag_querystart (T this);
extern int
Diag_queryend (T this);
extern int
Diag_nconsecutive (T this);
extern bool
Diag_dominatedp (T this);
extern void
Diag_set_dominatedp (T this);
extern int
Diag_ascending_cmp (const void *a, const void *b);
extern int
Diag_descending_cmp (const void *a, const void *b);
extern double
Diag_update_coverage (bool *coveredp, int *ncovered, List_T diagonals, int querylength);
extern int
Diag_compare_querystart (const void *x, const void *y);
extern void
Diag_print_segments (List_T diagonals, char *queryseq_ptr, char *genomicseg_ptr);
extern void
Diag_range (int *start, int *end, List_T diagonals, int querylength);
extern int
Diag_compute_bounds (int *diag_querystart, int *diag_queryend,
		     Chrpos_T *minactive, Chrpos_T *maxactive, List_T diagonals,
		     int querylength, bool debug_graphic_p,
		     Chrpos_T chrstart, Chrpos_T chrend,
		     Univcoord_T chroffset, Univcoord_T chrhigh, bool plusp);
extern void
Diag_max_bounds (Chrpos_T *minactive, Chrpos_T *maxactive,
		 int querylength, Chrpos_T chrstart, Chrpos_T chrend,
		 Univcoord_T chroffset, Univcoord_T chrhigh, bool plusp);


#undef T
#endif

