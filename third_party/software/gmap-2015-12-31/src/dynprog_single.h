/* $Id: dynprog_single.h 132731 2014-04-08 21:19:57Z twu $ */
#ifndef DYNPROG_SINGLE_INCLUDED
#define DYNPROG_SINGLE_INCLUDED

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
Dynprog_single_setup (bool homopolymerp_in);


extern List_T
Dynprog_single_gap (int *dynprogindex, int *finalscore,
		    int *nmatches, int *nmismatches, int *nopens, int *nindels,
		    T dynprog, char *sequence1, char *sequenceuc1,
		    int length1, int length2, int offset1, int offset2,
		    Univcoord_T chroffset, Univcoord_T chrhigh,
		    bool watsonp, bool jump_late_p, Pairpool_T pairpool,
		    int extraband_single, double defect_rate, bool widebandp);


extern List_T
Dynprog_microexon_int (double *bestprob2, double *bestprob3, int *dynprogindex, int *microintrontype,
		       char *sequence1, char *sequenceuc1,
		       int length1, int length2L, int length2R,
		       int offset1, int offset2L, int revoffset2R, int cdna_direction,
		       char *queryseq, char *queryuc,
		       Univcoord_T chroffset, Univcoord_T chrhigh, bool watsonp,
		       Pairpool_T pairpool, double defect_rate);

extern List_T
Dynprog_microexon_5 (int *dynprogindex, int *microintrontype, int *microexonlength,
		     char *revsequence1, char *revsequenceuc1,
		     char *revsequence2, char *revsequenceuc2,
		     int length1, int length2, int revoffset1, int revoffset2, int cdna_direction,
		     char *queryseq, char *queryuc, char *genomicseg, char *genomicuc,
		     Pairpool_T pairpool, bool end_microexons_p);

extern List_T
Dynprog_microexon_3 (int *dynprogindex, int *microintrontype, int *microexonlength,
		     char *sequence1, char *sequenceuc1,
		     char *sequence2, char *sequenceuc2,
		     int length1, int length2, int offset1, int offset2, int cdna_direction,
		     char *queryseq, char *queryuc, char *genomicseg, char *genomicuc,
		     int genomiclength, Pairpool_T pairpool, bool end_microexons_p);

#undef T
#endif

