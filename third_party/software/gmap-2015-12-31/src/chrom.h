/* $Id: chrom.h 157221 2015-01-22 18:38:57Z twu $ */
#ifndef CHROM_INCLUDED
#define CHROM_INCLUDED

#include "bool.h"
#include "genomicpos.h"
#include "types.h"

#define T Chrom_T
typedef struct T *T;

extern void
Chrom_free (T *old);
extern char *
Chrom_string (T this);
extern bool
Chrom_circularp (T this);
extern T
Chrom_from_string (char *string, char *mitochondrial_string, Univcoord_T order, bool circularp);

extern int
Chrom_cmp_alpha (T a, T b);
extern int
Chrom_cmp_numeric_alpha (T a, T b);
extern int
Chrom_cmp_chrom (T a, T b);
extern int
Chrom_cmp_order (T a, T b);

extern int
Chrom_compare_order (const void *x, const void *y);
extern int
Chrom_compare_alpha (const void *x, const void *y);
extern int
Chrom_compare_numeric_alpha (const void *x, const void *y);
extern int
Chrom_compare_chrom (const void *x, const void *y);
extern int
Chrom_compare_order (const void *x, const void *y);

extern int
Chrom_compare_table (const void *x, const void *y);
extern unsigned int
Chrom_hash_table (const void *key);

#undef T
#endif
