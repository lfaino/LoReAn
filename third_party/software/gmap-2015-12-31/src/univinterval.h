/* $Id: univinterval.h 168395 2015-06-26 17:13:13Z twu $ */
#ifndef UNIVINTERVAL_INCLUDED
#define UNIVINTERVAL_INCLUDED

#include "bool.h"
#include "genomicpos.h"
#include "types.h"

#define T Univinterval_T
typedef struct T *T;
struct T {
  Univ_IIT_coord_T low;		/* low <= high */
  Univ_IIT_coord_T high;
  int sign;
  int type;
};

extern T
Univinterval_new (Univcoord_T low, Univcoord_T high, int type);
extern T
Univinterval_copy (T old);
extern void
Univinterval_free (T *old);
extern void
Univinterval_table_free (void **x);
extern void
Univinterval_print (T this);

extern Univcoord_T
Univinterval_low (T this);
extern Univcoord_T
Univinterval_high (T this);
extern void
Univinterval_store_length (T this, Chrpos_T length);
extern int
Univinterval_sign (T this);
extern Chrpos_T
Univinterval_length (T this);
extern int
Univinterval_type (T this);

extern Univcoord_T
Univinterval_array_low (struct T *intervals, int index);
extern Univcoord_T
Univinterval_array_high (struct T *intervals, int index);

extern bool
Univinterval_is_contained (Univcoord_T x, struct T *intervals, int index);
extern bool
Univinterval_overlap_p (Univcoord_T x, Univcoord_T y, struct T *intervals, int index);

extern void
Univinterval_qsort_by_sigma (int *table, int i, int j, struct T *intervals);
extern void
Univinterval_qsort_by_omega (int *table, int i, int j, struct T *intervals);

extern int
Univinterval_cmp (const void *a, const void *b);
extern int
Univinterval_cmp_low (const void *a, const void *b);
extern int
Univinterval_cmp_high (const void *a, const void *b);

#if 0
extern int
Univinterval_table_cmp (const void *a, const void *b);
extern unsigned int
Univinterval_table_hash (const void *a);
#endif

bool
Univinterval_equal (T x, T y);

#undef T
#endif


