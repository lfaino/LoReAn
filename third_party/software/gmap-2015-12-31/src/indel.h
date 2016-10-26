/* $Id: indel.h 166641 2015-05-29 21:13:04Z twu $ */
#ifndef INDEL_INCLUDED
#define INDEL_INCLUDED

#include "bool.h"
#include "list.h"
#include "chrnum.h"
#include "genomicpos.h"
#include "compress.h"
#include "genome.h"

extern void
Indel_setup (int min_indel_end_matches_in, int indel_penalty_middle_in);

extern int
Indel_resolve_middle_insertion (int *best_nmismatches_i, int *best_nmismatches_j,
				Univcoord_T left, int indels, Compress_T query_compress,
				int querystart, int queryend, int querylength,
				int max_mismatches_allowed,
				bool plusp, int genestrand, bool first_read_p);

extern int
Indel_resolve_middle_deletion (int *best_nmismatches_i, int *best_nmismatches_j,
			       Univcoord_T left, int indels, Compress_T query_compress,
			       int querystart, int queryend, int querylength,
			       int max_mismatches_allowed,
			       bool plusp, int genestrand, bool first_read_p);


extern List_T
Indel_solve_middle_insertion (bool *foundp, int *found_score, int *nhits, List_T hits,
			      Univcoord_T left, Chrnum_T chrnum, Univcoord_T chroffset,
			      Univcoord_T chrhigh, Chrpos_T chrlength,
			      int indels, Compress_T query_compress,
			      int querylength, int max_mismatches_allowed,
			      bool plusp, int genestrand, bool first_read_p, bool sarrayp);

extern List_T
Indel_solve_middle_deletion (bool *foundp, int *found_score, int *nhits, List_T hits,
			     Univcoord_T left, Chrnum_T chrnum, Univcoord_T chroffset,
			     Univcoord_T chrhigh, Chrpos_T chrlength,
			     int indels, Compress_T query_compress, int querylength,
			     int max_mismatches_allowed,
			     bool plusp, int genestrand, bool first_read_p, bool sarrayp);

#endif

