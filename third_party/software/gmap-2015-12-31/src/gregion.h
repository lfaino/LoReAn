/* $Id: gregion.h 157221 2015-01-22 18:38:57Z twu $ */
#ifndef GREGION_INCLUDED
#define GREGION_INCLUDED

#include "bool.h"
#include "genomicpos.h"
#include "types.h"
#include "chrnum.h"
#include "iit-read-univ.h"
#include "match.h"

#define T Gregion_T
typedef struct T *T;

extern void
Gregion_print (T this);

extern void
Gregion_free (T *old);

extern Univcoord_T
Gregion_genomicstart (T this);

extern Univcoord_T
Gregion_genomicend (T this);

extern Chrpos_T
Gregion_chrstart (T this);

extern Chrpos_T
Gregion_chrend (T this);

extern Chrpos_T
Gregion_genomiclength (T this);

extern bool
Gregion_plusp (T this);

extern bool
Gregion_revcompp (T this);

extern int
Gregion_genestrand (T this);

extern Chrnum_T
Gregion_chrnum (T this);

extern char *
Gregion_chr (T this, Univ_IIT_T chromosome_iit);

extern Univcoord_T
Gregion_chroffset (T this);

extern Univcoord_T
Gregion_chrhigh (T this);

extern Chrpos_T
Gregion_chrlength (T this);

extern int
Gregion_querystart (T this);

extern int
Gregion_queryend (T this);

extern int
Gregion_matchsize (T this);

extern double
Gregion_weight (T this);

extern int
Gregion_support (T this);

extern bool 
Gregion_extendedp (T this);

extern void
Gregion_set_ncovered (T this, int ncovered, int source);

extern int
Gregion_ncovered (T this);


extern T
Gregion_new (int nexons, Univcoord_T genomicstart, Univcoord_T genomicend,
	     bool plusp, int genestrand, Univ_IIT_T chromosome_iit, int querystart, int queryend, 
	     int querylength, int matchsize, int trimstart, int trimend, int circular_typeint);

extern T
Gregion_new_from_matches (Match_T match5, Match_T match3, int genestrand, Univ_IIT_T chromosome_iit,
			  int querylength, int matchsize, int trimstart, int trimend, int circular_typeint);

extern List_T
Gregion_filter_unique (List_T gregionlist);

extern List_T
Gregion_filter_support (List_T gregionlist, int boundary, double pct_max, int diff_max);

extern double
Gregion_best_weight (List_T gregionlist);

extern List_T
Gregion_filter_by_evidence (List_T gregionlist);

extern bool
Gregion_sufficient_support (T this);

extern void
Gregion_extend (T this, Chrpos_T extension5, Chrpos_T extension3, int querylength,
		int min_extra_end);

extern int
Gregion_cmp (const void *a, const void *b);

extern void
Gregion_filter_clean (List_T gregionlist, int nchrs);

#undef T
#endif


