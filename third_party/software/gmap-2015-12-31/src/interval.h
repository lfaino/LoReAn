/* $Id: interval.h 157221 2015-01-22 18:38:57Z twu $ */
#ifndef INTERVAL_INCLUDED
#define INTERVAL_INCLUDED

#include "bool.h"
#include "genomicpos.h"
#include "types.h"

#define T Interval_T
typedef struct T *T;
struct T {
  Chrpos_T low;		/* low <= high */
  Chrpos_T high;
  int sign;
  int type;
};

extern T
Interval_new (Chrpos_T low, Chrpos_T high, int type);
extern T
Interval_copy (T old);
extern void
Interval_copy_existing (T dest, T src);
extern void
Interval_free (T *old);
extern void
Interval_print (T this);

extern Chrpos_T
Interval_low (T this);
extern Chrpos_T
Interval_high (T this);
extern void
Interval_store_length (T this, Chrpos_T length);
extern int
Interval_sign (T this);
extern Chrpos_T
Interval_length (T this);
extern int
Interval_type (T this);

extern Chrpos_T
Interval_array_low (struct T *intervals, int index);
extern Chrpos_T
Interval_array_high (struct T *intervals, int index);

extern bool
Interval_is_contained (Chrpos_T x, struct T *intervals, int index);
extern bool
Interval_overlap_p (Chrpos_T x, Chrpos_T y, struct T *intervals, int index);

extern void
Interval_qsort_by_sigma (int *table, int i, int j, struct T *intervals);
extern void
Interval_qsort_by_omega (int *table, int i, int j, struct T *intervals);

extern int
Interval_cmp (const void *a, const void *b);
extern int
Interval_cmp_low (const void *a, const void *b);
extern int
Interval_cmp_high (const void *a, const void *b);
extern int
Interval_cmp_low_struct (const void *a, const void *b);


struct Interval_windex_T {
  int index;
  T interval;
};

extern int
Interval_windex_cmp (const void *a, const void *b);

#undef T
#endif


