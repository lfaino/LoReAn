/* $Id: dynprog_end.h 132731 2014-04-08 21:19:57Z twu $ */
#ifndef DYNPROG_END_INCLUDED
#define DYNPROG_END_INCLUDED

#include "bool.h"
#include "list.h"
#include "pairpool.h"
#include "chrnum.h"
#include "iit-read.h"
#include "splicetrie_build.h"	/* For splicetype */
#include "types.h"
#include "dynprog.h"

#define T Dynprog_T

extern void
Dynprog_end_setup (Univcoord_T *splicesites_in, Splicetype_T *splicetypes_in,
		   Chrpos_T *splicedists_in, int nsplicesites_in,
		   Trieoffset_T *trieoffsets_obs_in, Triecontent_T *triecontents_obs_in,
		   Trieoffset_T *trieoffsets_max_in, Triecontent_T *triecontents_max_in);


extern List_T
Dynprog_end5_gap (int *dynprogindex, int *finalscore, int *nmatches, int *nmismatches,
		  int *nopens, int *nindels, T dynprog,
		  char *revsequence1, char *revsequenceuc1,
		  int length1, int length2, int revoffset1, int revoffset2, 
		  Univcoord_T chroffset, Univcoord_T chrhigh,
		  int cdna_direction, bool watsonp, bool jump_late_p, Pairpool_T pairpool,
		  int extraband_end, double defect_rate, Endalign_T endalign);

extern List_T
Dynprog_end5_splicejunction (int *dynprogindex, int *finalscore, int *missscore,
			     int *nmatches, int *nmismatches, int *nopens, int *nindels, T dynprog, 
			     char *revsequence1, char *revsequenceuc1,
			     char *revsequence2, char *revsequenceuc2, char *revsequencealt2,
			     int length1, int length2, int revoffset1, int revoffset2_anchor, int revoffset2_far,
			     Univcoord_T chroffset, Univcoord_T chrhigh,
			     int cdna_direction, bool watsonp, bool jump_late_p, Pairpool_T pairpool,
			     int extraband_end, double defect_rate, int contlength);

extern List_T
Dynprog_end3_gap (int *dynprogindex, int *finalscore, int *nmatches, int *nmismatches,
		  int *nopens, int *nindels, T dynprog,
		  char *sequence1, char *sequenceuc1,
		  int length1, int length2, int offset1, int offset2, 
		  Univcoord_T chroffset, Univcoord_T chrhigh,
		  int cdna_direction, bool watsonp, bool jump_late_p, Pairpool_T pairpool,
		  int extraband_end, double defect_rate, Endalign_T endalign);

extern List_T
Dynprog_end3_splicejunction (int *dynprogindex, int *finalscore, int *missscore,
			     int *nmatches, int *nmismatches, int *nopens, int *nindels, T dynprog,
			     char *sequence1, char *sequenceuc1,
			     char *sequence2, char *sequenceuc2, char *sequencealt2,
			     int length1, int length2, int offset1, int offset2_anchor, int offset2_far,
			     Univcoord_T chroffset, Univcoord_T chrhigh,
			     int cdna_direction, bool watsonp, bool jump_late_p, Pairpool_T pairpool,
			     int extraband_end, double defect_rate, int contlength);

extern bool
Dynprog_make_splicejunction_5 (char *splicejunction, char *splicejunction_alt, Univcoord_T splicecoord,
			       int splicelength, int contlength, Splicetype_T far_splicetype,
			       bool watsonp);

extern bool
Dynprog_make_splicejunction_3 (char *splicejunction, char *splicejunction_alt, Univcoord_T splicecoord,
			       int splicelength, int contlength, Splicetype_T far_splicetype,
			       bool watsonp);

#if 0
extern List_T
Dynprog_add_known_splice_5 (int *nmatches_distal, List_T pairs, Univcoord_T anchor_splicesite, Univcoord_T far_splicesite,
			    Univcoord_T chroffset, Univcoord_T chrhigh, bool watsonp, Pairpool_T pairpool);

extern List_T
Dynprog_add_known_splice_3 (int *nmatches_distal, List_T pairs, Univcoord_T anchor_splicesite, Univcoord_T far_splicesite,
			    Univcoord_T chroffset, Univcoord_T chrhigh, bool watsonp, Pairpool_T pairpool);
#endif


extern List_T
Dynprog_end5_known (bool *knownsplicep, int *dynprogindex, int *finalscore,
		    int *ambig_end_length, Splicetype_T *ambig_splicetype,
		    int *nmatches, int *nmismatches, int *nopens, int *nindels, T dynprog, 
		    char *revsequence1, char *revsequenceuc1,
		    int length1, int length2, int revoffset1, int revoffset2, 
		    Univcoord_T chroffset, Univcoord_T chrhigh,
		    Univcoord_T knownsplice_limit_low, Univcoord_T knownsplice_limit_high,
#ifdef GSNAP
#ifdef END_KNOWNSPLICING_SHORTCUT
		    int cutoff_level, char *queryptr, int querylength, Compress_T query_compress,
#endif
#endif
		    int cdna_direction, bool watsonp, bool jump_late_p,
		    Pairpool_T pairpool, int extraband_end, double defect_rate);

extern List_T
Dynprog_end3_known (bool *knownsplicep, int *dynprogindex, int *finalscore,
		    int *ambig_end_length, Splicetype_T *ambig_splicetype,
		    int *nmatches, int *nmismatches, int *nopens, int *nindels, T dynprog, 
		    char *sequence1, char *sequenceuc1,
		    int length1, int length2, int offset1, int offset2, int querylength,
		    Univcoord_T chroffset, Univcoord_T chrhigh,
		    Univcoord_T knownsplice_limit_low, Univcoord_T knownsplice_limit_high,
#ifdef GSNAP
#ifdef END_KNOWNSPLICING_SHORTCUT
		    int cutoff_level, char *queryptr, int querylength, Compress_T query_compress,
#endif
#endif
		    int cdna_direction, bool watsonp, bool jump_late_p,
		    Pairpool_T pairpool, int extraband_end, double defect_rate);


#undef T
#endif

