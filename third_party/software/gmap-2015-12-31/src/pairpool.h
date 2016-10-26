/* $Id: pairpool.h 147823 2014-09-15 23:13:11Z twu $ */
#ifndef PAIRPOOL_INCLUDED
#define PAIRPOOL_INCLUDED

typedef struct Pairpool_T *Pairpool_T;

#include "bool.h"
#include "pair.h"
#include "list.h"

#define T Pairpool_T

extern void
Pairpool_free (T *old);
extern void
Pairpool_free_memory (T this);
extern void
Pairpool_report_memory (T this);
extern T
Pairpool_new (void);
extern void
Pairpool_reset (T this);
extern List_T
Pairpool_push (List_T list, T this, int querypos, int genomepos, char cdna, char comp,
	       char genome, char genomealt, int dynprogindex);
extern List_T
Pairpool_push_copy (List_T list, T this, Pair_T orig);
extern List_T
Pairpool_push_gapalign (List_T list, T this, int querypos, int genomepos, char cdna, char comp,
			int introntype, char genome, char genomealt, bool extraexonp);
extern List_T
Pairpool_push_gapholder (List_T list, T this, int queryjump, int genomejump,
			 Pair_T leftpair, Pair_T rightpair, bool knownp);
extern List_T
Pairpool_push_existing (List_T list, T this, Pair_T pair);
extern List_T
Pairpool_pop (List_T list, Pair_T *x);
extern List_T
Pairpool_transfer (List_T dest, List_T source);
extern List_T
Pairpool_transfer_n (List_T dest, List_T source, int n);
extern int
Pairpool_count_bounded (int *nstart, List_T source, int minpos, int maxpos);
extern List_T
Pairpool_copy (List_T source, T this);
extern struct Pair_T *
Pairpool_copy_array (struct Pair_T *source, int npairs);
extern void
Pairpool_clean_join (List_T *left_path, List_T *right_pairs);
extern List_T
Pairpool_remove_gapholders (List_T pairs);
extern List_T
Pairpool_join_end3 (List_T path_orig, List_T end3_pairs_orig, Pairpool_T pairpool,
		    bool copy_end_p);
extern List_T
Pairpool_join_end5 (List_T pairs_orig, List_T end5_path_orig, Pairpool_T pairpool,
		    bool copy_end_p);

extern List_T
Pairpool_add_queryskip (List_T pairs, int r, int c, int dist, char *querysequence,
			int queryoffset, int genomeoffset, Pairpool_T pairpool, bool revp, int dynprogindex);
extern List_T
Pairpool_add_genomeskip (bool *add_dashes_p, List_T pairs, int r, int c, int dist,
			 char *genomesequence, char *genomesequenceuc,
			 int queryoffset, int genomeoffset, Pairpool_T pairpool, bool revp,
			 Univcoord_T chroffset, Univcoord_T chrhigh,
			 int cdna_direction, bool watsonp, int dynprogindex, bool use_genomicseg_p);

extern List_T
Pairpool_compact_copy (List_T list, T dest);

#undef T
#endif
