/* $Id: bitpack64-write.h 179363 2015-11-20 22:13:52Z twu $ */
#ifndef BITPACK64_WRITE_INCLUDED
#define BITPACK64_WRITE_INCLUDED
#include <stdio.h>
#include "types.h"

/* Stores the $(n+1)$ values [0..n] */
extern void
Bitpack64_write_differential (char *ptrsfile, char *compfile, UINT4 *ascending, Oligospace_T n);
extern void
Bitpack64_write_differential_paired (char *ptrsfile, char *compfile, UINT4 *ascending, Oligospace_T n);
extern void
Bitpack64_write_fixed10 (char *ptrsfile, char *compfile, UINT4 *ascending, Oligospace_T n);
extern void
Bitpack64_write_differential_huge (char *pagesfile, char *ptrsfile, char *compfile,
				   UINT8 *ascending, Oligospace_T n);
extern void
Bitpack64_write_fixed10_huge (char *pagesfile, char *ptrsfile, char *compfile,
			      UINT8 *ascending, Oligospace_T n);

/* Stores the $n$ values [0..(n-1)] */
extern void
Bitpack64_write_direct (char *ptrsfile, char *compfile, UINT4 *direct, Oligospace_T n);

#endif
