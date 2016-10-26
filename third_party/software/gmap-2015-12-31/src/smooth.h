/* $Id: smooth.h 106198 2013-08-28 23:07:34Z twu $ */
#ifndef SMOOTH_INCLUDED
#define SMOOTH_INCLUDED
#include "list.h"
#include "pairpool.h"

extern void
Smooth_reset (List_T pairs);
extern List_T
Smooth_pairs_by_netgap (bool *deletep, List_T pairs, Pairpool_T pairpool);
extern List_T
Smooth_pairs_by_size (bool *shortp, bool *deletep, List_T pairs, Pairpool_T pairpool, int stage2_indexsize);
extern List_T
Smooth_mark_short_exons (List_T pairs, Pairpool_T pairpool, int stage2_indexsize);
extern List_T
Smooth_pairs_by_intronprobs (bool *badp, List_T pairs, Pairpool_T pairpool);

#endif


