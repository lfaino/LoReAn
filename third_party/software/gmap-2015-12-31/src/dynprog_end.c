static char rcsid[] = "$Id: dynprog_end.c 145990 2014-08-25 21:47:32Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "dynprog_end.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>		/* For ceil, log, pow */
#include <ctype.h>		/* For tolower */
#ifdef HAVE_SSE2
#include <emmintrin.h>
#endif
#ifdef HAVE_SSE4_1
#include <smmintrin.h>
#endif


#include "bool.h"
#include "except.h"
#include "assert.h"
#include "mem.h"
#include "comp.h"
#include "pair.h"
#include "pairdef.h"
#include "listdef.h"
#include "complement.h"
#include "splicetrie.h"
#include "dynprog_simd.h"


/* Tests whether get_genomic_nt == genomicseg in compute_scores procedures */
/* #define EXTRACT_GENOMICSEG 1 */


/* Prints parameters and results */
#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif

/* Ends */
#ifdef DEBUG6
#define debug6(x) x
#else
#define debug6(x)
#endif

/* Ends, known splicing.  May want to turn on DEBUG3 in splicetrie.c  */
#ifdef DEBUG7
#define debug7(x) x
#else
#define debug7(x)
#endif

/* Getting genomic nt */
#ifdef DEBUG8
#define debug8(x) x
#else
#define debug8(x)
#endif

/* Binary search */
#ifdef DEBUG10
#define debug10(x) x
#else
#define debug10(x)
#endif

/* traceback_nogaps */
#ifdef DEBUG11
#define debug11(x) x
#else
#define debug11(x)
#endif



#define FULLMATCH 3

/* Previously allowed lower mismatch scores on end to allow more
   complete alignments to the end, and because ends are typically of
   lower quality.  Previously made equal to FULLMATCH, because
   criterion is nmatches >= nmismatches.  However, extensions at ends
   appear to defeat purpose of trimming, so increase mismatch at end
   from -3 to -4. */
#define MISMATCH_ENDQ -5

/* Ends tend to be of lower quality, so we don't want to introduce gaps.
   Also, we make then indifferent to the quality of the rest of the
   sequence. */
/* was -10 open and -3 extend */
#define END_OPEN_HIGHQ -10
#define END_OPEN_MEDQ -8
#define END_OPEN_LOWQ -6

#define END_EXTEND_HIGHQ -2
#define END_EXTEND_MEDQ -2
#define END_EXTEND_LOWQ -2

#define LAZY_INDEL 1		/* Don't advance to next coordinate on final indel, since could go over chromosome bounds. */


#define T Dynprog_T


static Univcoord_T *splicesites;
static Splicetype_T *splicetypes;
static Chrpos_T *splicedists;
static int nsplicesites;

static Trieoffset_T *trieoffsets_obs;
static Triecontent_T *triecontents_obs;
static Trieoffset_T *trieoffsets_max;
static Triecontent_T *triecontents_max;

bool homopolymerp;

void
Dynprog_end_setup (Univcoord_T *splicesites_in, Splicetype_T *splicetypes_in,
		   Chrpos_T *splicedists_in, int nsplicesites_in,
		   Trieoffset_T *trieoffsets_obs_in, Triecontent_T *triecontents_obs_in,
		   Trieoffset_T *trieoffsets_max_in, Triecontent_T *triecontents_max_in) {

  splicesites = splicesites_in;
  splicetypes = splicetypes_in;
  splicedists = splicedists_in;
  nsplicesites = nsplicesites_in;
  trieoffsets_obs = trieoffsets_obs_in;
  triecontents_obs = triecontents_obs_in;
  trieoffsets_max = trieoffsets_max_in;
  triecontents_max = triecontents_max_in;

  return;
}


#if defined(HAVE_SSE4_1) || defined(HAVE_SSE2)
static void
find_best_endpoint_8 (int *finalscore, int *bestr, int *bestc, Score8_T **matrix_upper,
		      Score8_T **matrix_lower, int rlength, int glength,
		      int lband, int uband, bool jump_late_p) {
  Score8_T bestscore = 0;
  int r, c;
  int clo, chigh;
  /* No need for loffset or cmid because they apply only for cdnaend
     == FIVE, which doesn't require searching */

  *bestr = *bestc = 0;

#if 0
  /* Just go along main diagonal */
  uband = extraband_end_or_paired;
  lband = extraband_end_or_paired;
#endif

  if (jump_late_p == false) {
    /* use > bestscore */
    for (r = 1; r <= rlength; r++) {
      if ((clo = r - lband) < 1) {
	clo = 1;
      }
      if ((chigh = r + uband) > glength) {
	chigh = glength;
      }
      for (c = clo; c < /*to main diagonal*/r; c++) {
	if (matrix_lower[r][c] > bestscore) {
	  *bestr = r;
	  *bestc = c;
	  bestscore = matrix_lower[r][c];
	}
      }
      for (/*at main diagonal*/; c <= chigh; c++) {
	if (matrix_upper[c][r] > bestscore) {
	  *bestr = r;
	  *bestc = c;
	  bestscore = matrix_upper[c][r];
	}
      }

    }
  } else {
    /* use >= bestscore */
    for (r = 1; r <= rlength; r++) {
      if ((clo = r - lband) < 1) {
	clo = 1;
      }
      if ((chigh = r + uband) > glength) {
	chigh = glength;
      }
      for (c = clo; c < /*to main diagonal*/r; c++) {
	if (matrix_lower[r][c] >= bestscore) {
	  *bestr = r;
	  *bestc = c;
	  bestscore = matrix_lower[r][c];
	}
      }
      for (/*at main diagonal*/; c <= chigh; c++) {
	if (matrix_upper[c][r] >= bestscore) {
	  *bestr = r;
	  *bestc = c;
	  bestscore = matrix_upper[c][r];
	}
      }
    }
  }


  *finalscore = (int) bestscore;
  return;
}
#endif

static void
find_best_endpoint_16 (int *finalscore, int *bestr, int *bestc,
		       Score16_T **matrix_upper, Score16_T **matrix_lower,
		       int rlength, int glength, int lband, int uband,
		       bool jump_late_p) {
  Score16_T bestscore = 0;
  int r, c;
  int clo, chigh;
  /* No need for loffset or cmid because they apply only for cdnaend
     == FIVE, which doesn't require searching */

  *bestr = *bestc = 0;

#if 0
  /* Just go along main diagonal */
  uband = extraband_end_or_paired;
  lband = extraband_end_or_paired;
#endif

  if (jump_late_p == false) {
    /* use > bestscore */
    for (r = 1; r <= rlength; r++) {
      if ((clo = r - lband) < 1) {
	clo = 1;
      }
      if ((chigh = r + uband) > glength) {
	chigh = glength;
      }
      for (c = clo; c < /*to main diagonal*/r; c++) {
	if (matrix_lower[r][c] > bestscore) {
	  *bestr = r;
	  *bestc = c;
	  bestscore = matrix_lower[r][c];
	}
      }
      for (/*at main diagonal*/; c <= chigh; c++) {
	if (matrix_upper[c][r] > bestscore) {
	  *bestr = r;
	  *bestc = c;
	  bestscore = matrix_upper[c][r];
	}
      }
    }
  } else {
    /* use >= bestscore */
    for (r = 1; r <= rlength; r++) {
      if ((clo = r - lband) < 1) {
	clo = 1;
      }
      if ((chigh = r + uband) > glength) {
	chigh = glength;
      }
      for (c = clo; c < /*to main diagonal*/r; c++) {
	if (matrix_lower[r][c] >= bestscore) {
	  *bestr = r;
	  *bestc = c;
	  bestscore = matrix_lower[r][c];
	}
      }
      for (/*at main diagonal*/; c <= chigh; c++) {
	if (matrix_upper[c][r] >= bestscore) {
	  *bestr = r;
	  *bestc = c;
	  bestscore = matrix_upper[c][r];
	}
      }
    }
  }


  *finalscore = (int) bestscore;
  return;
}

static void
find_best_endpoint_std (int *finalscore, int *bestr, int *bestc, Score32_T **matrix,
			int rlength, int glength, int lband, int uband,
			bool jump_late_p) {
  Score32_T bestscore = 0;
  int r, c;
  int clo, chigh;
  /* No need for loffset or cmid because they apply only for cdnaend
     == FIVE, which doesn't require searching */

  *bestr = *bestc = 0;

#if 0
  /* Just go along main diagonal */
  uband = extraband_end_or_paired;
  lband = extraband_end_or_paired;
#endif

  if (jump_late_p == false) {
    /* use > bestscore */
    for (r = 1; r <= rlength; r++) {
      if ((clo = r - lband) < 1) {
	clo = 1;
      }
      if ((chigh = r + uband) > glength) {
	chigh = glength;
      }
      for (c = clo; c <= chigh; c++) {
	if (matrix[c][r] > bestscore) {
	  *bestr = r;
	  *bestc = c;
	  bestscore = matrix[c][r];
	}
      }
    }
  } else {
    /* use >= bestscore */
    for (r = 1; r <= rlength; r++) {
      if ((clo = r - lband) < 1) {
	clo = 1;
      }
      if ((chigh = r + uband) > glength) {
	chigh = glength;
      }
      for (c = clo; c <= chigh; c++) {
	if (matrix[c][r] >= bestscore) {
	  *bestr = r;
	  *bestc = c;
	  bestscore = matrix[c][r];
	}
      }
    }
  }


  *finalscore = (int) bestscore;
  return;
}

#if defined(HAVE_SSE4_1) || defined(HAVE_SSE2)
static void
find_best_endpoint_to_queryend_indels_8 (int *finalscore, int *bestr, int *bestc,
					 Score8_T **matrix_upper, Score8_T **matrix_lower,
					 int rlength, int glength, int lband, int uband,
					 bool jump_late_p) {
  Score8_T bestscore = NEG_INFINITY_8;
  int r, c;
  int clo, chigh;
  /* No need for loffset or cmid because they apply only for cdnaend
     == FIVE, which doesn't require searching */

#if 0
  if (glength >= rlength) {
    /* Widen band to right to reach destination */
    uband = glength - rlength + extraband_end_or_paired;
    lband = extraband_end_or_paired;
  } else {
    /* Widen band to left to reach destination */
    lband = rlength - glength + extraband_end_or_paired;
    uband = extraband_end_or_paired;
  }
#endif

  *bestr = r = rlength;
  *bestc = 0;

  if (jump_late_p == false) {
    if ((clo = r - lband) < 1) {
      clo = 1;
    }
    if ((chigh = r + uband) > glength) {
      chigh = glength;
    }
    for (c = clo; c < /*to main diagonal*/r; c++) {
      if (matrix_lower[r][c] > bestscore) {
	*bestr = r;
	*bestc = c;
	bestscore = matrix_lower[r][c];
      }
    }
    for (/*at main diagonal*/; c <= chigh; c++) {
      if (matrix_upper[c][r] > bestscore) {
	*bestr = r;
	*bestc = c;
	bestscore = matrix_upper[c][r];
      }
    }

  } else {
    if ((clo = r - lband) < 1) {
      clo = 1;
    }
    if ((chigh = r + uband) > glength) {
      chigh = glength;
    }
    for (c = clo; c < /*to main diagonal*/r; c++) {
      if (matrix_lower[r][c] >= bestscore) {
	*bestr = r;
	*bestc = c;
	bestscore = matrix_lower[r][c];
      }
    }
    for (/*at main diagonal*/; c <= chigh; c++) {
      if (matrix_upper[c][r] >= bestscore) {
	*bestr = r;
	*bestc = c;
	bestscore = matrix_upper[c][r];
      }
    }
  }

  *finalscore = (int) bestscore;
  return;
}
#endif


static void
find_best_endpoint_to_queryend_indels_16 (int *finalscore, int *bestr, int *bestc,
					  Score16_T **matrix_upper, Score16_T **matrix_lower,
					  int rlength, int glength, int lband, int uband,
					  bool jump_late_p) {
  Score16_T bestscore = NEG_INFINITY_16;
  int r, c;
  int clo, chigh;
  /* No need for loffset or cmid because they apply only for cdnaend
     == FIVE, which doesn't require searching */

#if 0
  if (glength >= rlength) {
    /* Widen band to right to reach destination */
    uband = glength - rlength + extraband_end_or_paired;
    lband = extraband_end_or_paired;
  } else {
    /* Widen band to left to reach destination */
    lband = rlength - glength + extraband_end_or_paired;
    uband = extraband_end_or_paired;
  }
#endif

  *bestr = r = rlength;
  *bestc = 0;

  if (jump_late_p == false) {
    if ((clo = r - lband) < 1) {
      clo = 1;
    }
    if ((chigh = r + uband) > glength) {
      chigh = glength;
    }
    for (c = clo; c < /*to main diagonal*/r; c++) {
      if (matrix_lower[r][c] > bestscore) {
	*bestr = r;
	*bestc = c;
	bestscore = matrix_lower[r][c];
      }
    }
    for (/*at main diagonal*/; c <= chigh; c++) {
      if (matrix_upper[c][r] > bestscore) {
	*bestr = r;
	*bestc = c;
	bestscore = matrix_upper[c][r];
      }
    }

  } else {
    if ((clo = r - lband) < 1) {
      clo = 1;
    }
    if ((chigh = r + uband) > glength) {
      chigh = glength;
    }
    for (c = clo; c < /*to main diagonal*/r; c++) {
      if (matrix_lower[r][c] >= bestscore) {
	*bestr = r;
	*bestc = c;
	bestscore = matrix_lower[r][c];
      }
    }
    for (/*at main diagonal*/; c <= chigh; c++) {
      if (matrix_upper[c][r] >= bestscore) {
	*bestr = r;
	*bestc = c;
	bestscore = matrix_upper[c][r];
      }
    }
  }

  *finalscore = (int) bestscore;
  return;
}


static void
find_best_endpoint_to_queryend_indels_std (int *finalscore, int *bestr, int *bestc, Score32_T **matrix, 
					   int rlength, int glength, int lband, int uband,
					   bool jump_late_p) {
  Score32_T bestscore = NEG_INFINITY_32;
  int r, c;
  int clo, chigh;
  /* No need for loffset or cmid because they apply only for cdnaend
     == FIVE, which doesn't require searching */

#if 0
  if (glength >= rlength) {
    /* Widen band to right to reach destination */
    uband = glength - rlength + extraband_end_or_paired;
    lband = extraband_end_or_paired;
  } else {
    /* Widen band to left to reach destination */
    lband = rlength - glength + extraband_end_or_paired;
    uband = extraband_end_or_paired;
  }
#endif

  *bestr = r = rlength;
  *bestc = 0;

  if (jump_late_p == false) {
    if ((clo = r - lband) < 1) {
      clo = 1;
    }
    if ((chigh = r + uband) > glength) {
      chigh = glength;
    }
    for (c = clo; c <= chigh; c++) {
      if (matrix[c][r] > bestscore) {
	*bestr = r;
	*bestc = c;
	bestscore = matrix[c][r];
      }
    }

  } else {
    if ((clo = r - lband) < 1) {
      clo = 1;
    }
    if ((chigh = r + uband) > glength) {
      chigh = glength;
    }
    for (c = clo; c <= chigh; c++) {
      if (matrix[c][r] >= bestscore) {
	*bestr = r;
	*bestc = c;
	bestscore = matrix[c][r];
      }
    }
  }

  *finalscore = (int) bestscore;
  return;
}


static void
find_best_endpoint_to_queryend_nogaps (int *bestr, int *bestc, int rlength, int glength) {
  if (glength < rlength) {
    *bestr = glength;
    *bestc = glength;
  } else {
    *bestr = rlength;
    *bestc = rlength;
  }

  return;
}


/* revp means both rev1p and rev2p, which must have equal values */
/* Iterative version */
static List_T
traceback_nogaps (List_T pairs, int *nmatches, int *nmismatches,
		  int r, int c, char *rsequence, char *rsequenceuc,
		  char *gsequence, char *gsequence_alt,
		  int queryoffset, int genomeoffset, Pairpool_T pairpool, 
#ifdef DEBUG14
		  Univcoord_T chroffset, Univcoord_T chrhigh, bool watsonp,
#endif
		  bool revp, int dynprogindex) {
  char c1, c1_uc, c2, c2_alt;
  int querycoord, genomecoord;
#ifdef DEBUG14
  char c2_single;
#endif

  debug11(printf("Starting traceback_nogaps at r=%d,c=%d (roffset=%d, goffset=%d), revp %d\n",
		 r,c,queryoffset,genomeoffset,revp));

  /* printf("genome sequence is %s\n",genomesequence); */
  while (r > 0 && c > 0) {
    querycoord = r-1;
    genomecoord = c-1;
    if (revp == true) {
      querycoord = -querycoord;
      genomecoord = -genomecoord;
    }

    c1 = rsequence[querycoord];
    c1_uc = rsequenceuc[querycoord];
    c2 = gsequence[genomecoord];
    c2_alt = gsequence_alt[genomecoord];

#ifdef DEBUG14
    c2_single = get_genomic_nt(&c2_alt,genomeoffset+genomecoord,chroffset,chrhigh,watsonp);
    if (c2 != c2_single) {
      abort();
    }
#endif

#ifdef EXTRACT_GENOMICSEG
    debug8(printf("genome sequence at %d is %c\n",genomecoord,genomesequence[genomecoord]));
    assert(c2 == genomesequence[genomecoord]);
#endif

    if (c2 == '*') {
      /* Don't push pairs past end of chromosome */

    } else if (/*querysequenceuc[querycoord]*/c1_uc == c2 || c1_uc == c2_alt) {
      debug11(printf("D: Pushing %d,%d [%d,%d] (%c,%c) - match\n",
		   r,c,queryoffset+querycoord,genomeoffset+genomecoord,c1_uc,c2));
      *nmatches += 1;
      pairs = Pairpool_push(pairs,pairpool,queryoffset+querycoord,genomeoffset+genomecoord,
			    c1,DYNPROG_MATCH_COMP,c2,c2_alt,dynprogindex);
      
    } else if (consistent_array[(int) c1_uc][(int) c2] == true || consistent_array[(int) c1_uc][(int) c2_alt] == true) {
      debug11(printf("D: Pushing %d,%d [%d,%d] (%c,%c) - ambiguous\n",
		   r,c,queryoffset+querycoord,genomeoffset+genomecoord,c1_uc,c2));
      *nmatches += 1;
      pairs = Pairpool_push(pairs,pairpool,queryoffset+querycoord,genomeoffset+genomecoord,
			    c1,AMBIGUOUS_COMP,c2,c2_alt,dynprogindex);

    } else {
      debug11(printf("D: Pushing %d,%d [%d,%d] (%c,%c) - mismatch\n",
		   r,c,queryoffset+querycoord,genomeoffset+genomecoord,c1_uc,c2));
      *nmismatches += 1;
      pairs = Pairpool_push(pairs,pairpool,queryoffset+querycoord,genomeoffset+genomecoord,
			    c1,MISMATCH_COMP,c2,c2_alt,dynprogindex);
    }
    r--; c--;
  }

  return pairs;
}


#if defined(HAVE_SSE4_1) || defined(HAVE_SSE2)
/* Want to keep pointers to r and c because traceback is interrupted */
static List_T
traceback_local_8_upper (List_T pairs, int *nmatches, int *nmismatches, int *nopens, int *nindels,
			 Direction8_T **directions_nogap, Direction8_T **directions_Egap,
			 int *r, int *c, int endc, char *rsequence, char *rsequenceuc,
			 char *genomesequence, char *genomesequenceuc, char *genomesequencealt,
			 int queryoffset, int genomeoffset, Pairpool_T pairpool, bool revp,
			 Univcoord_T chroffset, Univcoord_T chrhigh,
			 int cdna_direction, bool watsonp, int dynprogindex) {
  char c1, c1_uc, c2, c2_alt;
  int dist;
  bool add_dashes_p;
  int querycoord, genomecoord;
  Direction8_T dir;

  debug(printf("Starting traceback_local_8_upper at r=%d,c=%d (roffset=%d, goffset=%d)\n",
	       *r,*c,queryoffset,genomeoffset));

  /* We care only only about genomic coordinate c */

  while (*r > 0 && *c > endc) {
    if ((dir = directions_nogap[*c][*r]) != DIAG) {
      /* Must be HORIZ */
      dist = 1;
      while (*c > endc && directions_Egap[(*c)--][*r] != DIAG) {
	dist++;
      }
      /* assert(*c != endc); */

      debug(printf("H%d: ",dist));
      pairs = Pairpool_add_genomeskip(&add_dashes_p,pairs,*r,(*c)+dist,dist,
				      genomesequence,genomesequenceuc,
				      queryoffset,genomeoffset,pairpool,revp,chroffset,chrhigh,
				      cdna_direction,watsonp,dynprogindex,/*use_genomicseg_p*/true);
      if (add_dashes_p == true) {
	*nopens += 1;
	*nindels += dist;
      }
      debug(printf("\n"));

    } else {
      querycoord = (*r)-1;
      genomecoord = (*c)-1;
      if (revp == true) {
	querycoord = -querycoord;
	genomecoord = -genomecoord;
      }

      c1 = rsequence[querycoord];
      c1_uc = rsequenceuc[querycoord];
      c2 = genomesequence[genomecoord];
      c2_alt = genomesequencealt[genomecoord];

      if (/*querysequenceuc[querycoord]*/c1_uc == c2 || c1_uc == c2_alt) {
	debug(printf("Pushing %d,%d [%d,%d] (%c,%c) - match\n",
		     *r,*c,queryoffset+querycoord,genomeoffset+genomecoord,c1_uc,c2));
	*nmatches += 1;
	pairs = Pairpool_push(pairs,pairpool,queryoffset+querycoord,genomeoffset+genomecoord,
			      c1,DYNPROG_MATCH_COMP,c2,c2_alt,dynprogindex);

      } else if (consistent_array[(int) c1_uc][(int) c2] == true || consistent_array[(int) c1_uc][(int) c2_alt] == true) {
	debug(printf("Pushing %d,%d [%d,%d] (%c,%c) - ambiguous\n",
		     *r,*c,queryoffset+querycoord,genomeoffset+genomecoord,c1_uc,c2));
	*nmatches += 1;
	pairs = Pairpool_push(pairs,pairpool,queryoffset+querycoord,genomeoffset+genomecoord,
			      c1,AMBIGUOUS_COMP,c2,c2_alt,dynprogindex);

      } else {
	debug(printf("Pushing %d,%d [%d,%d] (%c,%c) - mismatch\n",
		     *r,*c,queryoffset+querycoord,genomeoffset+genomecoord,c1_uc,c2));
	*nmismatches += 1;
	pairs = Pairpool_push(pairs,pairpool,queryoffset+querycoord,genomeoffset+genomecoord,
			      c1,MISMATCH_COMP,c2,c2_alt,dynprogindex);
      }

      (*r)--; (*c)--;
    }
  }

  /* assert(*r == 0); */
  if (/* *r == 0 && */ *c == endc) {
    /* Finished with a diagonal step */

  } else {
    dist = (*c) - endc;
    debug(printf("H%d: ",dist));
    pairs = Pairpool_add_genomeskip(&add_dashes_p,pairs,/**r*/0+LAZY_INDEL,*c,dist,genomesequence,genomesequenceuc,
				    queryoffset,genomeoffset,pairpool,revp,chroffset,chrhigh,
				    cdna_direction,watsonp,dynprogindex,
				    /*use_genomicseg_p*/true);
    if (add_dashes_p == true) {
      *nopens += 1;
      *nindels += dist;
    }
    debug(printf("\n"));
  }

  return pairs;
}
#endif


#if defined(HAVE_SSE4_1) || defined(HAVE_SSE2)
static List_T
traceback_local_8_lower (List_T pairs, int *nmatches, int *nmismatches, int *nopens, int *nindels,
			 Direction8_T **directions_nogap, Direction8_T **directions_Egap,
			 int *r, int *c, int endc, char *rsequence, char *rsequenceuc,
			 char *genomesequence, char *genomesequenceuc, char *genomesequencealt,
			 int queryoffset, int genomeoffset, Pairpool_T pairpool, bool revp,
			 Univcoord_T chroffset, Univcoord_T chrhigh,
			 int cdna_direction, bool watsonp, int dynprogindex) {
  char c1, c1_uc, c2, c2_alt;
  int dist;
  int querycoord, genomecoord;
  Direction8_T dir;

  debug(printf("Starting traceback_local_8_lower at r=%d,c=%d (roffset=%d, goffset=%d)\n",*r,*c,queryoffset,genomeoffset));

  /* We care only only about genomic coordinate c */

  while (*r > 0 && *c > endc) {
    if ((dir = directions_nogap[*r][*c]) != DIAG) {
      /* Must be VERT */
      dist = 1;
      /* Should not need to check for r > 0 if the main diagonal is populated with DIAG */
      while (/* *r > 0 && */ directions_Egap[(*r)--][*c] != DIAG) {
	dist++;
      }
      /* assert(*r != 0); */
			     
      debug(printf("V%d: ",dist));
      pairs = Pairpool_add_queryskip(pairs,(*r)+dist,*c,dist,rsequence,
				     queryoffset,genomeoffset,pairpool,revp,
				     dynprogindex);
      *nopens += 1;
      *nindels += dist;
      debug(printf("\n"));
      
    } else {
      querycoord = (*r)-1;
      genomecoord = (*c)-1;
      if (revp == true) {
	querycoord = -querycoord;
	genomecoord = -genomecoord;
      }

      c1 = rsequence[querycoord];
      c1_uc = rsequenceuc[querycoord];
      c2 = genomesequence[genomecoord];
      c2_alt = genomesequencealt[genomecoord];

      if (/*querysequenceuc[querycoord]*/c1_uc == c2 || c1_uc == c2_alt) {
	debug(printf("Pushing %d,%d [%d,%d] (%c,%c) - match\n",
		     *r,*c,queryoffset+querycoord,genomeoffset+genomecoord,c1_uc,c2));
	*nmatches += 1;
	pairs = Pairpool_push(pairs,pairpool,queryoffset+querycoord,genomeoffset+genomecoord,
			      c1,DYNPROG_MATCH_COMP,c2,c2_alt,dynprogindex);

      } else if (consistent_array[(int) c1_uc][(int) c2] == true || consistent_array[(int) c1_uc][(int) c2_alt] == true) {
	debug(printf("Pushing %d,%d [%d,%d] (%c,%c) - ambiguous\n",
		     *r,*c,queryoffset+querycoord,genomeoffset+genomecoord,c1_uc,c2));
	*nmatches += 1;
	pairs = Pairpool_push(pairs,pairpool,queryoffset+querycoord,genomeoffset+genomecoord,
			      c1,AMBIGUOUS_COMP,c2,c2_alt,dynprogindex);
	
      } else {
	debug(printf("Pushing %d,%d [%d,%d] (%c,%c) - mismatch\n",
		     *r,*c,queryoffset+querycoord,genomeoffset+genomecoord,c1_uc,c2));
	*nmismatches += 1;
	pairs = Pairpool_push(pairs,pairpool,queryoffset+querycoord,genomeoffset+genomecoord,
			      c1,MISMATCH_COMP,c2,c2_alt,dynprogindex);
      }

      (*r)--; (*c)--;
    }
  }

  /* assert(*c == endc); */
  if (*r == 0 /* && *c == endc */) {
    /* Finished with a diagonal step */

  } else {
    /* Must be VERT */
    dist = *r;
    debug(printf("V%d: ",dist));
    pairs = Pairpool_add_queryskip(pairs,*r,/**c*/endc+LAZY_INDEL,dist,rsequence,
				   queryoffset,genomeoffset,pairpool,revp,
				   dynprogindex);
    *nopens += 1;
    *nindels += dist;
    debug(printf("\n"));
  }

  return pairs;
}
#endif


static List_T
traceback_local_16_upper (List_T pairs, int *nmatches, int *nmismatches, int *nopens, int *nindels,
			  Direction16_T **directions_nogap, Direction16_T **directions_Egap,
			  int *r, int *c, int endc, char *rsequence, char *rsequenceuc,
			  char *genomesequence, char *genomesequenceuc, char *genomesequencealt,
			  int queryoffset, int genomeoffset, Pairpool_T pairpool, bool revp,
			  Univcoord_T chroffset, Univcoord_T chrhigh,
			  int cdna_direction, bool watsonp, int dynprogindex) {
  char c1, c1_uc, c2, c2_alt;
  int dist;
  bool add_dashes_p;
  int querycoord, genomecoord;
  Direction16_T dir;

  debug(printf("Starting traceback_local_upper at r=%d,c=%d (roffset=%d, goffset=%d)\n",*r,*c,queryoffset,genomeoffset));

  /* We care only only about genomic coordinate c */

  while (*r > 0 && *c > endc) {
    if ((dir = directions_nogap[*c][*r]) != DIAG) {
      /* Must be HORIZ */
      dist = 1;
      while (*c > endc && directions_Egap[(*c)--][*r] != DIAG) {
	dist++;
      }

      debug(printf("H%d: ",dist));
      pairs = Pairpool_add_genomeskip(&add_dashes_p,pairs,*r,(*c)+dist,dist,
				      genomesequence,genomesequenceuc,
				      queryoffset,genomeoffset,pairpool,revp,chroffset,chrhigh,
				      cdna_direction,watsonp,dynprogindex,/*use_genomicseg_p*/true);
      if (add_dashes_p == true) {
	*nopens += 1;
	*nindels += dist;
      }
      debug(printf("\n"));

    } else {
      querycoord = (*r)-1;
      genomecoord = (*c)-1;
      if (revp == true) {
	querycoord = -querycoord;
	genomecoord = -genomecoord;
      }

      c1 = rsequence[querycoord];
      c1_uc = rsequenceuc[querycoord];
      c2 = genomesequence[genomecoord];
      c2_alt = genomesequencealt[genomecoord];

      if (/*querysequenceuc[querycoord]*/c1_uc == c2 || c1_uc == c2_alt) {
	debug(printf("Pushing %d,%d [%d,%d] (%c,%c) - match\n",
		     *r,*c,queryoffset+querycoord,genomeoffset+genomecoord,c1_uc,c2));
	*nmatches += 1;
	pairs = Pairpool_push(pairs,pairpool,queryoffset+querycoord,genomeoffset+genomecoord,
			      c1,DYNPROG_MATCH_COMP,c2,c2_alt,dynprogindex);

      } else if (consistent_array[(int) c1_uc][(int) c2] == true || consistent_array[(int) c1_uc][(int) c2_alt] == true) {
	debug(printf("Pushing %d,%d [%d,%d] (%c,%c) - ambiguous\n",
		     *r,*c,queryoffset+querycoord,genomeoffset+genomecoord,c1_uc,c2));
	*nmatches += 1;
	pairs = Pairpool_push(pairs,pairpool,queryoffset+querycoord,genomeoffset+genomecoord,
			      c1,AMBIGUOUS_COMP,c2,c2_alt,dynprogindex);

      } else {
	debug(printf("Pushing %d,%d [%d,%d] (%c,%c) - mismatch\n",
		     *r,*c,queryoffset+querycoord,genomeoffset+genomecoord,c1_uc,c2));
	*nmismatches += 1;
	pairs = Pairpool_push(pairs,pairpool,queryoffset+querycoord,genomeoffset+genomecoord,
			      c1,MISMATCH_COMP,c2,c2_alt,dynprogindex);
      }

      (*r)--; (*c)--;
    }
  }

  /* assert(*r == 0); */
  if (/* *r == 0 && */ *c == endc) {
    /* Finished with a diagonal step */

  } else {
    /* Must be HORIZ */
    dist = (*c) - endc;
    debug(printf("H%d: ",dist));
    pairs = Pairpool_add_genomeskip(&add_dashes_p,pairs,/**r*/0+LAZY_INDEL,*c,dist,genomesequence,genomesequenceuc,
				      queryoffset,genomeoffset,pairpool,revp,chroffset,chrhigh,
				      cdna_direction,watsonp,dynprogindex,
				      /*use_genomicseg_p*/true);
    if (add_dashes_p == true) {
      *nopens += 1;
      *nindels += dist;
    }
    debug(printf("\n"));
  }

  return pairs;
}


static List_T
traceback_local_16_lower (List_T pairs, int *nmatches, int *nmismatches, int *nopens, int *nindels,
			  Direction16_T **directions_nogap, Direction16_T **directions_Egap,
			  int *r, int *c, int endc, char *rsequence, char *rsequenceuc,
			  char *genomesequence, char *genomesequenceuc, char *genomesequencealt,
			  int queryoffset, int genomeoffset, Pairpool_T pairpool, bool revp,
			  Univcoord_T chroffset, Univcoord_T chrhigh,
			  int cdna_direction, bool watsonp, int dynprogindex) {
  char c1, c1_uc, c2, c2_alt;
  int dist;
  int querycoord, genomecoord;
  Direction16_T dir;

  debug(printf("Starting traceback_local at r=%d,c=%d (roffset=%d, goffset=%d)\n",*r,*c,queryoffset,genomeoffset));

  /* We care only only about genomic coordinate c */

  while (*r > 0 && *c > endc) {
    if ((dir = directions_nogap[*r][*c]) != DIAG) {
      /* Must be VERT */
      dist = 1;
      /* Should not need to check for r > 0 if the main diagonal is populated with DIAG */
      while (/* *r > 0 && */ directions_Egap[(*r)--][*c] != DIAG) {
	dist++;
      }
      /* assert(*r != 0); */

      debug(printf("V%d: ",dist));
      pairs = Pairpool_add_queryskip(pairs,(*r)+dist,*c,dist,rsequence,
				     queryoffset,genomeoffset,pairpool,revp,
				     dynprogindex);
      *nopens += 1;
      *nindels += dist;
      debug(printf("\n"));

    } else {
      querycoord = (*r)-1;
      genomecoord = (*c)-1;
      if (revp == true) {
	querycoord = -querycoord;
	genomecoord = -genomecoord;
      }

      c1 = rsequence[querycoord];
      c1_uc = rsequenceuc[querycoord];
      c2 = genomesequence[genomecoord];
      c2_alt = genomesequencealt[genomecoord];

      if (/*querysequenceuc[querycoord]*/c1_uc == c2 || c1_uc == c2_alt) {
	debug(printf("Pushing %d,%d [%d,%d] (%c,%c) - match\n",
		     *r,*c,queryoffset+querycoord,genomeoffset+genomecoord,c1_uc,c2));
	*nmatches += 1;
	pairs = Pairpool_push(pairs,pairpool,queryoffset+querycoord,genomeoffset+genomecoord,
			      c1,DYNPROG_MATCH_COMP,c2,c2_alt,dynprogindex);

      } else if (consistent_array[(int) c1_uc][(int) c2] == true || consistent_array[(int) c1_uc][(int) c2_alt] == true) {
	debug(printf("Pushing %d,%d [%d,%d] (%c,%c) - ambiguous\n",
		     *r,*c,queryoffset+querycoord,genomeoffset+genomecoord,c1_uc,c2));
	*nmatches += 1;
	pairs = Pairpool_push(pairs,pairpool,queryoffset+querycoord,genomeoffset+genomecoord,
			      c1,AMBIGUOUS_COMP,c2,c2_alt,dynprogindex);

      } else {
	debug(printf("Pushing %d,%d [%d,%d] (%c,%c) - mismatch\n",
		     *r,*c,queryoffset+querycoord,genomeoffset+genomecoord,c1_uc,c2));
	*nmismatches += 1;
	pairs = Pairpool_push(pairs,pairpool,queryoffset+querycoord,genomeoffset+genomecoord,
			      c1,MISMATCH_COMP,c2,c2_alt,dynprogindex);
      }

      (*r)--; (*c)--;
    }
  }

  /* assert(*c == endc); */
  if (*r == 0 /* && *c == endc */) {
    /* Finished with a diagonal step */

  } else {
    /* Must be VERT */
    dist = *r;
    debug(printf("V%d: ",dist));
    pairs = Pairpool_add_queryskip(pairs,*r,/**c*/endc+LAZY_INDEL,dist,rsequence,
				   queryoffset,genomeoffset,pairpool,revp,
				   dynprogindex);
    *nopens += 1;
    *nindels += dist;
    debug(printf("\n"));
  }

  return pairs;
}


static List_T
traceback_local_std (List_T pairs, int *nmatches, int *nmismatches, int *nopens, int *nindels,
		     Direction32_T **directions_nogap, Direction32_T **directions_Egap, Direction32_T **directions_Fgap,
		     int *r, int *c, int endc, char *rsequence, char *rsequenceuc,
		     char *genomesequence, char *genomesequenceuc, char *genomesequencealt,
		     int queryoffset, int genomeoffset, Pairpool_T pairpool, bool revp,
		     Univcoord_T chroffset, Univcoord_T chrhigh,
		     int cdna_direction, bool watsonp, int dynprogindex) {
  char c1, c1_uc, c2, c2_alt;
  int dist;
  bool add_dashes_p;
  int querycoord, genomecoord;
  Direction32_T dir;

  debug(printf("Starting traceback_local at r=%d,c=%d (roffset=%d, goffset=%d)\n",*r,*c,queryoffset,genomeoffset));

  /* We care only only about genomic coordinate c */

  if (*c <= endc) {
    /* Do nothing */

  } else if ((dir = directions_nogap[*c][*r]) == DIAG) {
    /* Not an indel.  Do nothing. */

  } else if (dir == HORIZ) {
    dist = 1;
    while (*c > 1 && directions_Egap[*c][*r] != DIAG) {
      dist++;
      (*c)--;
    }
    (*c)--;
    /* dir = directions_nogap[*c][*r]; */

    debug(printf("H%d: ",dist));
    pairs = Pairpool_add_genomeskip(&add_dashes_p,pairs,*r,(*c)+dist,dist,genomesequence,genomesequenceuc,
				    queryoffset,genomeoffset,pairpool,revp,chroffset,chrhigh,
				    cdna_direction,watsonp,dynprogindex,/*use_genomicseg_p*/true);
    if (add_dashes_p == true) {
      *nopens += 1;
      *nindels += dist;
    }
    debug(printf("\n"));

  } else {
    /* Must be VERT */
    dist = 1;
    while (*r > 1 && directions_Fgap[*c][*r] != DIAG) {
      dist++;
      (*r)--;
    }
    (*r)--;
    /* dir = directions_nogap[*c][*r]; */

    debug(printf("V%d: ",dist));
    pairs = Pairpool_add_queryskip(pairs,(*r)+dist,*c,dist,rsequence,
				   queryoffset,genomeoffset,pairpool,revp,
				   dynprogindex);
    *nopens += 1;
    *nindels += dist;
    debug(printf("\n"));
  }

  while (*r > 0 && *c > endc) {
    querycoord = (*r)-1;
    genomecoord = (*c)-1;
    if (revp == true) {
      querycoord = -querycoord;
      genomecoord = -genomecoord;
    }

    c1 = rsequence[querycoord];
    c1_uc = rsequenceuc[querycoord];
    c2 = genomesequence[genomecoord];
    c2_alt = genomesequencealt[genomecoord];

    if (/*querysequenceuc[querycoord]*/c1_uc == c2 || c1_uc == c2_alt) {
      debug(printf("Pushing %d,%d [%d,%d] (%c,%c) - match\n",
		   r,c,queryoffset+querycoord,genomeoffset+genomecoord,c1_uc,c2));
      *nmatches += 1;
      pairs = Pairpool_push(pairs,pairpool,queryoffset+querycoord,genomeoffset+genomecoord,
			    c1,DYNPROG_MATCH_COMP,c2,c2_alt,dynprogindex);

    } else if (consistent_array[(int) c1_uc][(int) c2] == true || consistent_array[(int) c1_uc][(int) c2_alt] == true) {
      debug(printf("Pushing %d,%d [%d,%d] (%c,%c) - ambiguous\n",
		   r,c,queryoffset+querycoord,genomeoffset+genomecoord,c1_uc,c2));
      *nmatches += 1;
      pairs = Pairpool_push(pairs,pairpool,queryoffset+querycoord,genomeoffset+genomecoord,
			    c1,AMBIGUOUS_COMP,c2,c2_alt,dynprogindex);

    } else {
      debug(printf("Pushing %d,%d [%d,%d] (%c,%c) - mismatch\n",
		   r,c,queryoffset+querycoord,genomeoffset+genomecoord,c1_uc,c2));
      *nmismatches += 1;
      pairs = Pairpool_push(pairs,pairpool,queryoffset+querycoord,genomeoffset+genomecoord,
			    c1,MISMATCH_COMP,c2,c2_alt,dynprogindex);
    }

    (*r)--; (*c)--;
    if (*r == 0 && *c == 0) {
      /* STOP condition.  Do nothing. */

    } else if ((dir = directions_nogap[*c][*r]) == DIAG) {
      /* Do nothing */

    } else if (dir == HORIZ) {
      dist = 1;
      while (*c > 1 && directions_Egap[*c][*r] != DIAG) {
	dist++;
	(*c)--;
      }
      (*c)--;
      /* dir = directions_nogap[*c][*r]; */

      debug(printf("H%d: ",dist));
      pairs = Pairpool_add_genomeskip(&add_dashes_p,pairs,*r,(*c)+dist,dist,genomesequence,genomesequenceuc,
				      queryoffset,genomeoffset,pairpool,revp,chroffset,chrhigh,
				      cdna_direction,watsonp,dynprogindex,
				      /*use_genomicseg_p*/true);
      if (add_dashes_p == true) {
	*nopens += 1;
	*nindels += dist;
      }
      debug(printf("\n"));

    } else {
      /* Must be VERT */
      dist = 1;
      while (*r > 1 && directions_Fgap[*c][*r] != DIAG) {
	dist++;
	(*r)--;
      }
      (*r)--;
      /* dir = directions_nogap[*c][*r]; */

      debug(printf("V%d: ",dist));
      pairs = Pairpool_add_queryskip(pairs,(*r)+dist,*c,dist,rsequence,
				     queryoffset,genomeoffset,pairpool,revp,
				     dynprogindex);
      *nopens += 1;
      *nindels += dist;
      debug(printf("\n"));

    }
  }

  return pairs;
}



List_T
Dynprog_end5_gap (int *dynprogindex, int *finalscore, int *nmatches, int *nmismatches, 
		  int *nopens, int *nindels, T dynprog, 
		  char *rev_rsequence, char *rev_rsequenceuc,
		  int rlength, int glength, int rev_roffset, int rev_goffset, 
		  Univcoord_T chroffset, Univcoord_T chrhigh,
		  int cdna_direction, bool watsonp, bool jump_late_p, Pairpool_T pairpool,
		  int extraband_end, double defect_rate, Endalign_T endalign) {
  List_T pairs = NULL;
  char *rev_gsequence, *rev_gsequence_alt;
  Pair_T pair;
  Mismatchtype_T mismatchtype; 
  int bestr, bestc, lband, uband;
  int open, extend;
#if defined(HAVE_SSE4_1) || defined(HAVE_SSE2)
  bool use8p = false;
  Score8_T **matrix8_upper, **matrix8_lower;
  Direction8_T **directions8_upper_nogap, **directions8_upper_Egap,
    **directions8_lower_nogap, **directions8_lower_Egap;

  Score16_T **matrix16_upper, **matrix16_lower;
  Direction16_T **directions16_upper_nogap, **directions16_upper_Egap,
    **directions16_lower_nogap, **directions16_lower_Egap;
#else
  Score32_T **matrix;
  Direction32_T **directions_nogap, **directions_Egap, **directions_Fgap;
#endif
#ifdef PMAP
  int initpos, initmod;
#endif

  debug6(
	printf("%c:  ",*dynprogindex > 0 ? (*dynprogindex-1)%26+'a' : (-(*dynprogindex)-1)%26+'A');
	printf("Aligning 5' end gap with endalign %d\n",endalign);
	);

  mismatchtype = ENDQ;
  if (defect_rate < DEFECT_HIGHQ) {
    open = END_OPEN_HIGHQ;
    extend = END_EXTEND_HIGHQ;
  } else if (defect_rate < DEFECT_MEDQ) {
    open = END_OPEN_MEDQ;
    extend = END_EXTEND_MEDQ;
  } else {
    open = END_OPEN_LOWQ;
    extend = END_EXTEND_LOWQ;
  }

  /* We can just chop lengths to work, since we're not constrained on 5' end */
  if (rlength <= 0) {
    /* Needed to avoid abort by Matrix16_alloc */
    debug6(printf("rlength %d <= 0, so returning NULL\n",rlength));
    *nmatches = *nmismatches = *nopens = *nindels = 0;
    *finalscore = 0;
    return (List_T) NULL;
  } else if (endalign == QUERYEND_NOGAPS) {
    /* Don't shorten rlength */
  } else if (rlength > dynprog->max_rlength) {
    debug6(printf("rlength %d is too long.  Chopping to %d\n",rlength,dynprog->max_rlength));
    rlength = dynprog->max_rlength;
  }
  if (glength <= 0) {
    /* Needed to avoid abort by Matrix16_alloc */
    debug6(printf("glength %d <= 0, so returning NULL\n",glength));
    *nmatches = *nmismatches = *nopens = *nindels = 0;
    *finalscore = 0;
    return (List_T) NULL;
  } else if (endalign == QUERYEND_NOGAPS) {
    /* Don't shorten glength */
  } else if (glength > dynprog->max_glength) {
    debug6(printf("glength %d is too long.  Chopping to %d\n",glength,dynprog->max_glength));
    glength = dynprog->max_glength;
  }

  debug6(printf("At query offset %d-%d, %.*s\n",rev_roffset-rlength+1,rev_roffset,rlength,&(rev_rsequence[-rlength+1])));

#ifdef EXTRACT_GENOMICSEG
  debug6(printf("At genomic offset %d-%d, %.*s\n",
		rev_goffset-glength+1,rev_goffset,glength,&(rev_gsequence[-glength+1])));
#endif


  rev_gsequence = (char *) MALLOCA((glength+1) * sizeof(char));
  rev_gsequence_alt = (char *) MALLOCA((glength+1) * sizeof(char));

  if (watsonp) {
    Genome_get_segment_blocks_left(rev_gsequence,rev_gsequence_alt,/*right*/chroffset+rev_goffset+1,
				   glength,chroffset,/*revcomp*/false);
  } else {
    Genome_get_segment_blocks_right(rev_gsequence,rev_gsequence_alt,/*left*/chrhigh-rev_goffset,
				    glength,chrhigh,/*revcomp*/true);
  }
  if (rev_gsequence[0] == '\0') {
    *nmatches = *nmismatches = *nopens = *nindels = 0;
    *finalscore = 0;
    FREEA(rev_gsequence_alt);
    FREEA(rev_gsequence);
    return (List_T) NULL;
  }

  if (endalign == QUERYEND_GAP || endalign == BEST_LOCAL) {
    Dynprog_compute_bands(&lband,&uband,rlength,glength,extraband_end,/*widebandp*/true);
#if defined(HAVE_SSE4_1) || defined(HAVE_SSE2)
    /* Use || because we want the minimum length (which determines the diagonal length) to achieve a score less than 128 */
    if (rlength <= SIMD_MAXLENGTH_EPI8 || glength <= SIMD_MAXLENGTH_EPI8) {
      use8p = true;
      matrix8_upper = Dynprog_simd_8_upper(&directions8_upper_nogap,&directions8_upper_Egap,dynprog,
					   rev_rsequence,&(rev_gsequence[glength-1]),&(rev_gsequence_alt[glength-1]),
					   rlength,glength,
#ifdef DEBUG14
					   rev_goffset,chroffset,chrhigh,watsonp,
#endif
					   mismatchtype,open,extend,
					   uband,/*for revp true*/!jump_late_p,/*revp*/true);
      matrix8_lower = Dynprog_simd_8_lower(&directions8_lower_nogap,&directions8_lower_Egap,dynprog,
					   rev_rsequence,&(rev_gsequence[glength-1]),&(rev_gsequence_alt[glength-1]),
					   rlength,glength,
#ifdef DEBUG14
					   rev_goffset,chroffset,chrhigh,watsonp,
#endif
					   mismatchtype,open,extend,
					   lband,/*for revp true*/!jump_late_p,/*revp*/true);

      find_best_endpoint_8(&(*finalscore),&bestr,&bestc,matrix8_upper,matrix8_lower,
			   rlength,glength,lband,uband,!jump_late_p);

    } else {
      matrix16_upper = Dynprog_simd_16_upper(&directions16_upper_nogap,&directions16_upper_Egap,dynprog,
					     rev_rsequence,&(rev_gsequence[glength-1]),&(rev_gsequence_alt[glength-1]),
					     rlength,glength,
#ifdef DEBUG14
					     rev_goffset,chroffset,chrhigh,watsonp,
#endif
					     mismatchtype,open,extend,
					     uband,/*for revp true*/!jump_late_p,/*revp*/true);
      matrix16_lower = Dynprog_simd_16_lower(&directions16_lower_nogap,&directions16_lower_Egap,dynprog,
					     rev_rsequence,&(rev_gsequence[glength-1]),&(rev_gsequence_alt[glength-1]),
					     rlength,glength,
#ifdef DEBUG14
					     rev_goffset,chroffset,chrhigh,watsonp,
#endif
					     mismatchtype,open,extend,
					     lband,/*for revp true*/!jump_late_p,/*revp*/true);

      find_best_endpoint_16(&(*finalscore),&bestr,&bestc,matrix16_upper,matrix16_lower,
			    rlength,glength,lband,uband,!jump_late_p);
    }

#else
    /* Non-SIMD methods */
    matrix = Dynprog_standard(&directions_nogap,&directions_Egap,&directions_Fgap,dynprog,
			      rev_rsequence,&(rev_gsequence[glength-1]),&(rev_gsequence_alt[glength-1]),
			      rlength,glength,rev_goffset,chroffset,chrhigh,watsonp,
			      mismatchtype,open,extend,lband,uband,
			      /*for revp true*/!jump_late_p,/*revp*/true,/*saturation*/NEG_INFINITY_INT);
    find_best_endpoint_std(&(*finalscore),&bestr,&bestc,matrix,rlength,glength,lband,uband,
			   !jump_late_p);
#endif

  } else if (endalign == QUERYEND_INDELS) {
    Dynprog_compute_bands(&lband,&uband,rlength,glength,extraband_end,/*widebandp*/true);
#if defined(HAVE_SSE4_1) || defined(HAVE_SSE2)
    /* Use || because we want the minimum length (which determines the diagonal length) to achive a score less than 128 */
    if (rlength <= SIMD_MAXLENGTH_EPI8 || glength <= SIMD_MAXLENGTH_EPI8) {


      use8p = true;
      matrix8_upper = Dynprog_simd_8_upper(&directions8_upper_nogap,&directions8_upper_Egap,dynprog,
					   rev_rsequence,&(rev_gsequence[glength-1]),&(rev_gsequence_alt[glength-1]),
					   rlength,glength,
#ifdef DEBUG14
					   rev_goffset,chroffset,chrhigh,watsonp,
#endif
					   mismatchtype,open,extend,
					   uband,/*for revp true*/!jump_late_p,/*revp*/true);
      matrix8_lower = Dynprog_simd_8_lower(&directions8_lower_nogap,&directions8_lower_Egap,dynprog,
					   rev_rsequence,&(rev_gsequence[glength-1]),&(rev_gsequence_alt[glength-1]),
					   rlength,glength,
#ifdef DEBUG14
					   rev_goffset,chroffset,chrhigh,watsonp,
#endif
					   mismatchtype,open,extend,
					   lband,/*for revp true*/!jump_late_p,/*revp*/true);
      find_best_endpoint_to_queryend_indels_8(&(*finalscore),&bestr,&bestc,matrix8_upper,matrix8_lower,
					      rlength,glength,lband,uband,!jump_late_p);
      /* *finalscore = 0 -- Splicetrie procedures need to know finalscore */

    } else {
      matrix16_upper = Dynprog_simd_16_upper(&directions16_upper_nogap,&directions16_upper_Egap,dynprog,
					     rev_rsequence,&(rev_gsequence[glength-1]),&(rev_gsequence_alt[glength-1]),
					     rlength,glength,
#ifdef DEBUG14
					     rev_goffset,chroffset,chrhigh,watsonp,
#endif
					     mismatchtype,open,extend,
					     uband,/*for revp true*/!jump_late_p,/*revp*/true);
      matrix16_lower = Dynprog_simd_16_lower(&directions16_lower_nogap,&directions16_lower_Egap,dynprog,
					     rev_rsequence,&(rev_gsequence[glength-1]),&(rev_gsequence_alt[glength-1]),
					     rlength,glength,
#ifdef DEBUG14
					     rev_goffset,chroffset,chrhigh,watsonp,
#endif
					     mismatchtype,open,extend,
					     lband,/*for revp true*/!jump_late_p,/*revp*/true);
      find_best_endpoint_to_queryend_indels_16(&(*finalscore),&bestr,&bestc,matrix16_upper,matrix16_lower,
					       rlength,glength,lband,uband,!jump_late_p);
      /* *finalscore = 0 -- Splicetrie procedures need to know finalscore */
    }

#else
    /* Non-SIMD methods */
    matrix = Dynprog_standard(&directions_nogap,&directions_Egap,&directions_Fgap,dynprog,
			      rev_rsequence,&(rev_gsequence[glength-1]),&(rev_gsequence_alt[glength-1]),
			      rlength,glength,rev_goffset,chroffset,chrhigh,watsonp,
			      mismatchtype,open,extend,lband,uband,
			      /*for revp true*/!jump_late_p,/*revp*/true,/*saturation*/NEG_INFINITY_INT);
    find_best_endpoint_to_queryend_indels_std(&(*finalscore),&bestr,&bestc,matrix,rlength,glength,lband,uband,
					      !jump_late_p);
    /* *finalscore = 0 -- Splicetrie procedures need to know finalscore */

#endif


  } else if (endalign == QUERYEND_NOGAPS) {
    find_best_endpoint_to_queryend_nogaps(&bestr,&bestc,rlength,glength);
    /* *finalscore = 0;	-- Splicetrie procedures need to know finalscore */

  } else {
    fprintf(stderr,"Unexpected endalign value %d\n",endalign);
    abort();
  }


#ifdef PMAP
  initpos = rev_roffset-(bestc-1);
  debug6(printf("Initial query pos is %d\n",initpos));
  if ((initmod = initpos % 3) > 0) {
    if (bestr + initmod < rlength && bestc + initmod < glength) {
      debug6(printf("Rounding down by %d\n",initmod));
      bestr += initmod;
      bestc += initmod;
    }
  }
#endif

  *nmatches = *nmismatches = *nopens = *nindels = 0;
  if (endalign == QUERYEND_NOGAPS) {
    pairs = traceback_nogaps(NULL,&(*nmatches),&(*nmismatches),bestr,bestc,
			     rev_rsequence,rev_rsequenceuc,
			     &(rev_gsequence[glength-1]),&(rev_gsequence_alt[glength-1]),
			     rev_roffset,rev_goffset,pairpool,
#ifdef DEBUG14
			     chroffset,chrhigh,watsonp,
#endif
			     /*revp*/true,*dynprogindex);
    *finalscore = (*nmatches)*FULLMATCH + (*nmismatches)*MISMATCH_ENDQ;

#if defined(HAVE_SSE4_1) || defined(HAVE_SSE2)
  } else if (use8p == true) {
    if (bestc >= bestr) {
      pairs = Dynprog_traceback_8_upper(NULL,&(*nmatches),&(*nmismatches),&(*nopens),&(*nindels),
					directions8_upper_nogap,directions8_upper_Egap,
					bestr,bestc,rev_rsequence,rev_rsequenceuc,
					&(rev_gsequence[glength-1]),&(rev_gsequence_alt[glength-1]),
					rev_roffset,rev_goffset,pairpool,/*revp*/true,
					chroffset,chrhigh,cdna_direction,watsonp,*dynprogindex);
    } else {
      pairs = Dynprog_traceback_8_lower(NULL,&(*nmatches),&(*nmismatches),&(*nopens),&(*nindels),
					directions8_lower_nogap,directions8_lower_Egap,
					bestr,bestc,rev_rsequence,rev_rsequenceuc,
					&(rev_gsequence[glength-1]),&(rev_gsequence_alt[glength-1]),
					rev_roffset,rev_goffset,pairpool,/*revp*/true,
					chroffset,chrhigh,cdna_direction,watsonp,*dynprogindex);
    }

  } else {
    if (bestc >= bestr) {
      pairs = Dynprog_traceback_16_upper(NULL,&(*nmatches),&(*nmismatches),&(*nopens),&(*nindels),
					 directions16_upper_nogap,directions16_upper_Egap,
					 bestr,bestc,rev_rsequence,rev_rsequenceuc,
					 &(rev_gsequence[glength-1]),&(rev_gsequence_alt[glength-1]),
					 rev_roffset,rev_goffset,pairpool,/*revp*/true,
					 chroffset,chrhigh,cdna_direction,watsonp,*dynprogindex);
    } else {
      pairs = Dynprog_traceback_16_lower(NULL,&(*nmatches),&(*nmismatches),&(*nopens),&(*nindels),
					 directions16_lower_nogap,directions16_lower_Egap,
					 bestr,bestc,rev_rsequence,rev_rsequenceuc,
					 &(rev_gsequence[glength-1]),&(rev_gsequence_alt[glength-1]),
					 rev_roffset,rev_goffset,pairpool,/*revp*/true,
					 chroffset,chrhigh,cdna_direction,watsonp,*dynprogindex);
    }

#else
  } else {
    pairs = Dynprog_traceback_std(NULL,&(*nmatches),&(*nmismatches),&(*nopens),&(*nindels),
				  directions_nogap,directions_Egap,directions_Fgap,bestr,bestc,
				  rev_rsequence,rev_rsequenceuc,
				  &(rev_gsequence[glength-1]),&(rev_gsequence_alt[glength-1]),
				  rev_roffset,rev_goffset,pairpool,/*revp*/true,
				  chroffset,chrhigh,cdna_direction,watsonp,*dynprogindex);
#endif
  }


  if ((endalign == QUERYEND_GAP || endalign == BEST_LOCAL) && (*nmatches + 1) < *nmismatches) {
    *finalscore = 0;
    /* No need to free pairs */
    pairs = NULL;
  } else {
    /* Add 1 to count the match already in the alignment */
    pairs = List_reverse(pairs); /* Look at 5' end to remove excess gaps */
    while (pairs != NULL && (pair = List_head(pairs)) && pair->comp == INDEL_COMP) {
      pairs = List_next(pairs);
    }
  }

  /*
    Directions_free(directions);
    Matrix_free(matrix);
  */
  
  FREEA(rev_gsequence_alt);
  FREEA(rev_gsequence);

  debug6(printf("End of dynprog end5 gap\n\n"));

  *dynprogindex += (*dynprogindex > 0 ? +1 : -1);
  return List_reverse(pairs);
}



/* rev_gsequence is the splicejunction */
List_T
Dynprog_end5_splicejunction (int *dynprogindex, int *finalscore, int *missscore,
			     int *nmatches, int *nmismatches, int *nopens, int *nindels, T dynprog, 
			     char *rev_rsequence, char *rev_rsequenceuc,
			     char *rev_gsequence, char *rev_gsequence_uc, char *rev_gsequence_alt,
			     int rlength, int glength, int rev_roffset, int rev_goffset_anchor, int rev_goffset_far,
			     Univcoord_T chroffset, Univcoord_T chrhigh,
			     int cdna_direction, bool watsonp, bool jump_late_p, Pairpool_T pairpool,
			     int extraband_end, double defect_rate, int contlength) {
  List_T pairs = NULL;
  Pair_T pair;
  Mismatchtype_T mismatchtype;
  int bestr, bestc, lband, uband;
  int open, extend;
#if defined(HAVE_SSE4_1) || defined(HAVE_SSE2)
  bool use8p = false;
  Score8_T **matrix8_upper, **matrix8_lower;
  Direction8_T **directions8_upper_nogap, **directions8_upper_Egap,
    **directions8_lower_nogap, **directions8_lower_Egap;

  Score16_T **matrix16_upper, **matrix16_lower;
  Direction16_T **directions16_upper_nogap, **directions16_upper_Egap,
    **directions16_lower_nogap, **directions16_lower_Egap;
#else
  Score32_T **matrix;
  Direction32_T **directions_nogap, **directions_Egap, **directions_Fgap;
#endif
#ifdef PMAP
  int initpos, initmod;
#endif

  debug6(
	printf("%c:  ",*dynprogindex > 0 ? (*dynprogindex-1)%26+'a' : (-(*dynprogindex)-1)%26+'A');
	printf("Aligning 5' end gap with endalign QUERYEND_NOGAPS\n");
	);

  mismatchtype = ENDQ;
  if (defect_rate < DEFECT_HIGHQ) {
    open = END_OPEN_HIGHQ;
    extend = END_EXTEND_HIGHQ;
  } else if (defect_rate < DEFECT_MEDQ) {
    open = END_OPEN_MEDQ;
    extend = END_EXTEND_MEDQ;
  } else {
    open = END_OPEN_LOWQ;
    extend = END_EXTEND_LOWQ;
  }

  /* We can just chop lengths to work, since we're not constrained on 5' end */
  if (rlength <= 0 || rlength > dynprog->max_rlength) {
    /* Needed to avoid abort by Matrix16_alloc */
    debug6(printf("rlength %d <= 0, so returning NULL\n",rlength));
    *nmatches = *nmismatches = *nopens = *nindels = 0;
    *finalscore = 0;
    *missscore = -100;
    return (List_T) NULL;
  }
  if (glength <= 0 || glength > dynprog->max_glength) {
    /* Needed to avoid abort by Matrix16_alloc */
    debug6(printf("glength %d <= 0, so returning NULL\n",glength));
    *nmatches = *nmismatches = *nopens = *nindels = 0;
    *finalscore = 0;
    *missscore = -100;
    return (List_T) NULL;
  }

  debug6(printf("At query offset %d-%d, %.*s\n",rev_roffset-rlength+1,rev_roffset,rlength,&(rev_rsequence[-rlength+1])));
  debug6(printf("At genomic offset %d-%d, %.*s\n",
		rev_goffset_anchor-glength+1,rev_goffset_anchor,glength,&(rev_gsequence[-glength+1])));
  
#ifdef PMAP
  initpos = rev_roffset-(bestc-1);
  debug6(printf("Initial query pos is %d\n",initpos));
  if ((initmod = initpos % 3) > 0) {
    if (bestr + initmod < rlength && bestc + initmod < glength) {
      debug6(printf("Rounding down by %d\n",initmod));
      bestr += initmod;
      bestc += initmod;
    }
  }
#endif

  Dynprog_compute_bands(&lband,&uband,rlength,glength,extraband_end,/*widebandp*/true);
#if defined(HAVE_SSE4_1) || defined(HAVE_SSE2)
  /* Use || because we want the minimum length (which determines the diagonal length) to achieve a score less than 128 */
  if (rlength <= SIMD_MAXLENGTH_EPI8 || glength <= SIMD_MAXLENGTH_EPI8) {
    use8p = true;
    matrix8_upper = Dynprog_simd_8_upper(&directions8_upper_nogap,&directions8_upper_Egap,dynprog,
					 rev_rsequence,rev_gsequence_uc,rev_gsequence_alt,
					 rlength,glength,
#ifdef DEBUG14
					 /*goffset*/0,chroffset,chrhigh,watsonp,
#endif
					 mismatchtype,open,extend,
					 uband,/*for revp true*/!jump_late_p,/*revp*/true);
    matrix8_lower = Dynprog_simd_8_lower(&directions8_lower_nogap,&directions8_lower_Egap,dynprog,
					 rev_rsequence,rev_gsequence_uc,rev_gsequence_alt,
					 rlength,glength,
#ifdef DEBUG14
					 /*goffset*/0,chroffset,chrhigh,watsonp,
#endif
					 mismatchtype,open,extend,
					 lband,/*for revp true*/!jump_late_p,/*revp*/true);

    find_best_endpoint_to_queryend_indels_8(&(*finalscore),&bestr,&bestc,matrix8_upper,matrix8_lower,
					    rlength,glength,lband,uband,!jump_late_p);

  } else {
    matrix16_upper = Dynprog_simd_16_upper(&directions16_upper_nogap,&directions16_upper_Egap,dynprog,
					   rev_rsequence,rev_gsequence_uc,rev_gsequence_alt,
					   rlength,glength,
#ifdef DEBUG14
					   /*goffset*/0,chroffset,chrhigh,watsonp,
#endif
					   mismatchtype,open,extend,
					   uband,/*for revp true*/!jump_late_p,/*revp*/true);
    matrix16_lower = Dynprog_simd_16_lower(&directions16_lower_nogap,&directions16_lower_Egap,dynprog,
					   rev_rsequence,rev_gsequence_uc,rev_gsequence_alt,
					   rlength,glength,
#ifdef DEBUG14
					   /*goffset*/0,chroffset,chrhigh,watsonp,
#endif
					   mismatchtype,open,extend,
					   lband,/*for revp true*/!jump_late_p,/*revp*/true);
    find_best_endpoint_to_queryend_indels_16(&(*finalscore),&bestr,&bestc,matrix16_upper,matrix16_lower,
					     rlength,glength,lband,uband,!jump_late_p);
  }

#else
  /* Non-SIMD methods */
  matrix = Dynprog_standard(&directions_nogap,&directions_Egap,&directions_Fgap,dynprog,
			    rev_rsequence,rev_gsequence,rev_gsequence_alt,
			    rlength,glength,/*goffset*/0,chroffset,chrhigh,watsonp,
			    mismatchtype,open,extend,lband,uband,
			    /*for revp true*/!jump_late_p,/*revp*/true,/*saturation*/NEG_INFINITY_INT);
  find_best_endpoint_to_queryend_indels_std(&(*finalscore),&bestr,&bestc,matrix,rlength,glength,lband,uband,
					    !jump_late_p);
#endif

  if (*finalscore < 0) {
    /* Need a reasonable alignment to call a splice */
    return (List_T) NULL;

  } else {
    *nmatches = *nmismatches = *nopens = *nindels = 0;
#if defined(HAVE_SSE4_1) || defined(HAVE_SSE2)
    if (use8p == true) {
      if (bestc >= bestr) {
	pairs = traceback_local_8_upper(NULL,&(*nmatches),&(*nmismatches),&(*nopens),&(*nindels),
					directions8_upper_nogap,directions8_upper_Egap,&bestr,&bestc,/*endc*/contlength,
					rev_rsequence,rev_rsequenceuc,
					rev_gsequence,rev_gsequence_uc,rev_gsequence_alt,
					rev_roffset,rev_goffset_far,pairpool,/*revp*/true,
					chroffset,chrhigh,cdna_direction,watsonp,*dynprogindex);
      } else {
	pairs = traceback_local_8_lower(NULL,&(*nmatches),&(*nmismatches),&(*nopens),&(*nindels),
					directions8_lower_nogap,directions8_lower_Egap,&bestr,&bestc,/*endc*/contlength,
					rev_rsequence,rev_rsequenceuc,
					rev_gsequence,rev_gsequence_uc,rev_gsequence_alt,
					rev_roffset,rev_goffset_far,pairpool,/*revp*/true,
					chroffset,chrhigh,cdna_direction,watsonp,*dynprogindex);
      }

      pairs = Pairpool_push_gapholder(pairs,pairpool,/*queryjump*/0,/*genomejump*/rev_goffset_anchor - rev_goffset_far,
				      /*leftpair*/NULL,/*rightpair*/NULL,/*knownp*/true);

      if (bestc >= bestr) {
	pairs = traceback_local_8_upper(pairs,&(*nmatches),&(*nmismatches),&(*nopens),&(*nindels),
					directions8_upper_nogap,directions8_upper_Egap,&bestr,&bestc,/*endc*/0,
					rev_rsequence,rev_rsequenceuc,
					rev_gsequence,rev_gsequence_uc,rev_gsequence_alt,
					rev_roffset,rev_goffset_anchor,pairpool,/*revp*/true,
					chroffset,chrhigh,cdna_direction,watsonp,*dynprogindex);
      } else {
	pairs = traceback_local_8_lower(pairs,&(*nmatches),&(*nmismatches),&(*nopens),&(*nindels),
					directions8_lower_nogap,directions8_lower_Egap,&bestr,&bestc,/*endc*/0,
					rev_rsequence,rev_rsequenceuc,
					rev_gsequence,rev_gsequence_uc,rev_gsequence_alt,
					rev_roffset,rev_goffset_anchor,pairpool,/*revp*/true,
					chroffset,chrhigh,cdna_direction,watsonp,*dynprogindex);
      }

    } else {
      if (bestc >= bestr) {
	pairs = traceback_local_16_upper(NULL,&(*nmatches),&(*nmismatches),&(*nopens),&(*nindels),
					 directions16_upper_nogap,directions16_upper_Egap,&bestr,&bestc,/*endc*/contlength,
					 rev_rsequence,rev_rsequenceuc,
					 rev_gsequence,rev_gsequence_uc,rev_gsequence_alt,
					 rev_roffset,rev_goffset_far,pairpool,/*revp*/true,
					 chroffset,chrhigh,cdna_direction,watsonp,*dynprogindex);
      } else {
	pairs = traceback_local_16_lower(NULL,&(*nmatches),&(*nmismatches),&(*nopens),&(*nindels),
					 directions16_lower_nogap,directions16_lower_Egap,&bestr,&bestc,/*endc*/contlength,
					 rev_rsequence,rev_rsequenceuc,
					 rev_gsequence,rev_gsequence_uc,rev_gsequence_alt,
					 rev_roffset,rev_goffset_far,pairpool,/*revp*/true,
					 chroffset,chrhigh,cdna_direction,watsonp,*dynprogindex);
      }

      pairs = Pairpool_push_gapholder(pairs,pairpool,/*queryjump*/0,/*genomejump*/rev_goffset_anchor - rev_goffset_far,
				      /*leftpair*/NULL,/*rightpair*/NULL,/*knownp*/true);

      if (bestc >= bestr) {
	pairs = traceback_local_16_upper(pairs,&(*nmatches),&(*nmismatches),&(*nopens),&(*nindels),
					 directions16_upper_nogap,directions16_upper_Egap,&bestr,&bestc,/*endc*/0,
					 rev_rsequence,rev_rsequenceuc,
					 rev_gsequence,rev_gsequence_uc,rev_gsequence_alt,
					 rev_roffset,rev_goffset_anchor,pairpool,/*revp*/true,
					 chroffset,chrhigh,cdna_direction,watsonp,*dynprogindex);
      } else {
	pairs = traceback_local_16_lower(pairs,&(*nmatches),&(*nmismatches),&(*nopens),&(*nindels),
					 directions16_lower_nogap,directions16_lower_Egap,&bestr,&bestc,/*endc*/0,
					 rev_rsequence,rev_rsequenceuc,
					 rev_gsequence,rev_gsequence_uc,rev_gsequence_alt,
					 rev_roffset,rev_goffset_anchor,pairpool,/*revp*/true,
					 chroffset,chrhigh,cdna_direction,watsonp,*dynprogindex);
      }
    }

#else
    pairs = traceback_local_std(NULL,&(*nmatches),&(*nmismatches),&(*nopens),&(*nindels),
				directions_nogap,directions_Egap,directions_Fgap,&bestr,&bestc,/*endc*/contlength,
				rev_rsequence,rev_rsequenceuc,
				rev_gsequence,rev_gsequence_uc,rev_gsequence_alt,
				rev_roffset,rev_goffset_far,pairpool,/*revp*/true,
				chroffset,chrhigh,cdna_direction,watsonp,*dynprogindex);

    pairs = Pairpool_push_gapholder(pairs,pairpool,/*queryjump*/0,/*genomejump*/rev_goffset_anchor - rev_goffset_far,
				    /*leftpair*/NULL,/*rightpair*/NULL,/*knownp*/true);

    pairs = traceback_local_std(pairs,&(*nmatches),&(*nmismatches),&(*nopens),&(*nindels),
				directions_nogap,directions_Egap,directions_Fgap,&bestr,&bestc,/*endc*/0,
				rev_rsequence,rev_rsequenceuc,
				rev_gsequence,rev_gsequence_uc,rev_gsequence_alt,
				rev_roffset,rev_goffset_anchor,pairpool,/*revp*/true,
				chroffset,chrhigh,cdna_direction,watsonp,*dynprogindex);
#endif


    /* Score compared with perfect score, so heavy weight on mismatches may not be necessary */
    *finalscore = (*nmatches)*FULLMATCH + (*nmismatches)*MISMATCH_ENDQ + (*nopens)*open + (*nindels)*extend;
    *missscore = (*finalscore) - rlength*FULLMATCH;
    debug6(printf("finalscore %d = %d*%d matches + %d*%d mismatches + %d*%d opens + %d*%d extends\n",
		  *finalscore,FULLMATCH,*nmatches,MISMATCH_ENDQ,*nmismatches,open,*nopens,extend,*nindels));
    debug6(printf("missscore = %d\n",*missscore));

    /* Add 1 to count the match already in the alignment */
    pairs = List_reverse(pairs); /* Look at 5' end to remove excess gaps */
    while (pairs != NULL && (pair = List_head(pairs)) && pair->comp == INDEL_COMP) {
      pairs = List_next(pairs);
    }

    debug6(Pair_dump_list(pairs,true));
    debug6(printf("End of dynprog end5 gap splicejunction\n\n"));

    *dynprogindex += (*dynprogindex > 0 ? +1 : -1);
    return List_reverse(pairs);
  }
}



List_T
Dynprog_end3_gap (int *dynprogindex, int *finalscore, int *nmatches, int *nmismatches, 
		  int *nopens, int *nindels, T dynprog, 
		  char *rsequence, char *rsequenceuc,
		  int rlength, int glength, int roffset, int goffset, 
		  Univcoord_T chroffset, Univcoord_T chrhigh,
		  int cdna_direction, bool watsonp, bool jump_late_p, Pairpool_T pairpool,
		  int extraband_end, double defect_rate, Endalign_T endalign) {
  List_T pairs = NULL;
  char *gsequence, *gsequence_alt;
  Pair_T pair;
  Mismatchtype_T mismatchtype;
  int bestr, bestc, lband, uband;
  int open, extend;
#if defined(HAVE_SSE4_1) || defined(HAVE_SSE2)
  bool use8p = false;
  Score8_T **matrix8_upper, **matrix8_lower;
  Direction8_T **directions8_upper_nogap, **directions8_upper_Egap,
    **directions8_lower_nogap, **directions8_lower_Egap;

  Score16_T **matrix16_upper, **matrix16_lower;
  Direction16_T **directions16_upper_nogap, **directions16_upper_Egap,
    **directions16_lower_nogap, **directions16_lower_Egap;
#else
  Score32_T **matrix;
  Direction32_T **directions_nogap, **directions_Egap, **directions_Fgap;
#endif
#ifdef PMAP
  int termpos, termmod;
#endif

  debug6(
	printf("%c:  ",*dynprogindex > 0 ? (*dynprogindex-1)%26+'a' : (-(*dynprogindex)-1)%26+'A');
	printf("Aligning 3' end gap with endalign = %d\n",endalign);
	);

  mismatchtype = ENDQ;
  if (defect_rate < DEFECT_HIGHQ) {
    open = END_OPEN_HIGHQ;
    extend = END_EXTEND_HIGHQ;
  } else if (defect_rate < DEFECT_MEDQ) {
    open = END_OPEN_MEDQ;
    extend = END_EXTEND_MEDQ;
  } else {
    open = END_OPEN_LOWQ;
    extend = END_EXTEND_LOWQ;
  }

  /* We can just chop lengths to work, since we're not constrained on 3' end */
  if (rlength <= 0) {
    /* Needed to avoid abort by Matrix16_alloc */
    *nmatches = *nmismatches = *nopens = *nindels = 0;
    *finalscore = 0;
    return (List_T) NULL;
  } else if (endalign == QUERYEND_NOGAPS) {
    /* Don't shorten rlength */
  } else if (rlength > dynprog->max_rlength) {
    debug6(printf("rlength %d is too long.  Chopping to %d\n",rlength,dynprog->max_rlength));
    rlength = dynprog->max_rlength;
  }
  if (glength <= 0) {
    /* Needed to avoid abort by Matrix16_alloc */
    *nmatches = *nmismatches = *nopens = *nindels = 0;
    *finalscore = 0;
    return (List_T) NULL;
  } else if (endalign == QUERYEND_NOGAPS) {
    /* Don't shorten glength */
  } else if (glength > dynprog->max_glength) {
    debug6(printf("glength %d is too long.  Chopping to %d\n",glength,dynprog->max_glength));
    glength = dynprog->max_glength;
  }

  debug6(printf("At query offset %d-%d, %.*s\n",roffset,roffset+rlength-1,rlength,rsequence));
#ifdef EXTRACT_GENOMICSEG
  debug6(printf("At genomic offset %d-%d, %.*s\n",goffset,goffset+glength-1,glength,gsequence));
#endif


  gsequence = (char *) MALLOCA((glength+1) * sizeof(char));
  gsequence_alt = (char *) MALLOCA((glength+1) * sizeof(char));

  if (watsonp) {
    Genome_get_segment_blocks_right(gsequence,gsequence_alt,/*left*/chroffset+goffset,
				    glength,chrhigh,/*revcomp*/false);
  } else {
    Genome_get_segment_blocks_left(gsequence,gsequence_alt,/*right*/chrhigh-goffset+1,
				   glength,chroffset,/*revcomp*/true);
  }
  if (gsequence[0] == '\0') {
    *nmatches = *nmismatches = *nopens = *nindels = 0;
    *finalscore = 0;
    FREEA(gsequence_alt);
    FREEA(gsequence);
    return (List_T) NULL;
  }

  if (endalign == QUERYEND_GAP || endalign == BEST_LOCAL) {
    Dynprog_compute_bands(&lband,&uband,rlength,glength,extraband_end,/*widebandp*/true);
#if defined(HAVE_SSE4_1) || defined(HAVE_SSE2)
    /* Use || because we want the minimum length (which determines the diagonal length) to achieve a score less than 128 */
    if (rlength <= SIMD_MAXLENGTH_EPI8 || glength <= SIMD_MAXLENGTH_EPI8) {
      use8p = true;
      matrix8_upper = Dynprog_simd_8_upper(&directions8_upper_nogap,&directions8_upper_Egap,dynprog,
					   rsequenceuc,gsequence,gsequence_alt,rlength,glength,
#ifdef DEBUG14
					   goffset,chroffset,chrhigh,watsonp,
#endif
					   mismatchtype,open,extend,
					   uband,jump_late_p,/*revp*/false);
      matrix8_lower = Dynprog_simd_8_lower(&directions8_lower_nogap,&directions8_lower_Egap,dynprog,
					   rsequenceuc,gsequence,gsequence_alt,rlength,glength,
#ifdef DEBUG14
					   goffset,chroffset,chrhigh,watsonp,
#endif
					   mismatchtype,open,extend,
					   lband,jump_late_p,/*revp*/false);
      find_best_endpoint_8(&(*finalscore),&bestr,&bestc,matrix8_upper,matrix8_lower,
			   rlength,glength,lband,uband,jump_late_p);

    } else {
      matrix16_upper = Dynprog_simd_16_upper(&directions16_upper_nogap,&directions16_upper_Egap,dynprog,
					     rsequenceuc,gsequence,gsequence_alt,rlength,glength,
#ifdef DEBUG14
					     goffset,chroffset,chrhigh,watsonp,
#endif
					     mismatchtype,open,extend,
					     uband,jump_late_p,/*revp*/false);
      matrix16_lower = Dynprog_simd_16_lower(&directions16_lower_nogap,&directions16_lower_Egap,dynprog,
					     rsequenceuc,gsequence,gsequence_alt,rlength,glength,
#ifdef DEBUG14
					     goffset,chroffset,chrhigh,watsonp,
#endif
					     mismatchtype,open,extend,
					     lband,jump_late_p,/*revp*/false);
      find_best_endpoint_16(&(*finalscore),&bestr,&bestc,matrix16_upper,matrix16_lower,
			    rlength,glength,lband,uband,jump_late_p);
    }

#else
    /* Non-SIMD methods */
    matrix = Dynprog_standard(&directions_nogap,&directions_Egap,&directions_Fgap,dynprog,
			      rsequenceuc,gsequence,gsequence_alt,rlength,glength,
			      goffset,chroffset,chrhigh,watsonp,mismatchtype,open,extend,
			      lband,uband,jump_late_p,/*revp*/false,/*saturation*/NEG_INFINITY_INT);
    find_best_endpoint_std(&(*finalscore),&bestr,&bestc,matrix,rlength,glength,lband,uband,
			   jump_late_p);
#endif

  } else if (endalign == QUERYEND_INDELS) {
    Dynprog_compute_bands(&lband,&uband,rlength,glength,extraband_end,/*widebandp*/true);
#if defined(HAVE_SSE4_1) || defined(HAVE_SSE2)
    /* Use || because we want the minimum length (which determines the diagonal length) to achieve a score less than 128 */
    if (rlength <= SIMD_MAXLENGTH_EPI8 || glength <= SIMD_MAXLENGTH_EPI8) {
      use8p = true;
      matrix8_upper = Dynprog_simd_8_upper(&directions8_upper_nogap,&directions8_upper_Egap,dynprog,
					   rsequenceuc,gsequence,gsequence_alt,rlength,glength,
#ifdef DEBUG14
					   goffset,chroffset,chrhigh,watsonp,
#endif
					   mismatchtype,open,extend,
					   uband,jump_late_p,/*revp*/false);
      matrix8_lower = Dynprog_simd_8_lower(&directions8_lower_nogap,&directions8_lower_Egap,dynprog,
					   rsequenceuc,gsequence,gsequence_alt,rlength,glength,
#ifdef DEBUG14
					   goffset,chroffset,chrhigh,watsonp,
#endif
					   mismatchtype,open,extend,
					   lband,jump_late_p,/*revp*/false);
      find_best_endpoint_to_queryend_indels_8(&(*finalscore),&bestr,&bestc,matrix8_upper,matrix8_lower,
					      rlength,glength,lband,uband,jump_late_p);
      /* *finalscore = 0; -- Splicetrie procedures need to know finalscore */

    } else {
      matrix16_upper = Dynprog_simd_16_upper(&directions16_upper_nogap,&directions16_upper_Egap,dynprog,
					     rsequenceuc,gsequence,gsequence_alt,rlength,glength,
#ifdef DEBUG14
					     goffset,chroffset,chrhigh,watsonp,
#endif
					     mismatchtype,open,extend,
					     uband,jump_late_p,/*revp*/false);
      matrix16_lower = Dynprog_simd_16_lower(&directions16_lower_nogap,&directions16_lower_Egap,dynprog,
					     rsequenceuc,gsequence,gsequence_alt,rlength,glength,
#ifdef DEBUG14
					     goffset,chroffset,chrhigh,watsonp,
#endif
					     mismatchtype,open,extend,
					     lband,jump_late_p,/*revp*/false);
      find_best_endpoint_to_queryend_indels_16(&(*finalscore),&bestr,&bestc,matrix16_upper,matrix16_lower,
					       rlength,glength,lband,uband,jump_late_p);
      /* *finalscore = 0; -- Splicetrie procedures need to know finalscore */
    }

#else
    /* Non-SIMD methods */
    matrix = Dynprog_standard(&directions_nogap,&directions_Egap,&directions_Fgap,dynprog,
			      rsequenceuc,gsequence,gsequence_alt,rlength,glength,
			      goffset,chroffset,chrhigh,watsonp,mismatchtype,open,extend,
			      lband,uband,jump_late_p,/*revp*/false,/*saturation*/NEG_INFINITY_INT);
    find_best_endpoint_to_queryend_indels_std(&(*finalscore),&bestr,&bestc,matrix,rlength,glength,lband,uband,
					      jump_late_p);
    /* *finalscore = 0; -- Splicetrie procedures need to know finalscore */
#endif

  } else if (endalign == QUERYEND_NOGAPS) {
    find_best_endpoint_to_queryend_nogaps(&bestr,&bestc,rlength,glength);
    /* *finalscore = 0; -- Splicetrie procedures need to know finalscore */

  } else {
    fprintf(stderr,"Unexpected endalign value %d\n",endalign);
    abort();
  }

#ifdef PMAP
  termpos = roffset+(bestc-1);
  debug6(printf("Final query pos is %d\n",termpos));
  if ((termmod = termpos % 3) < 2) {
    if (bestr + (2 - termmod) < rlength && bestc + (2 - termmod) < glength) {
      debug6(printf("Rounding up by %d\n",2 - termmod));
      bestr += 2 - termmod;
      bestc += 2 - termmod;
    }
  }
#endif

  *nmatches = *nmismatches = *nopens = *nindels = 0;
  if (endalign == QUERYEND_NOGAPS) {
    pairs = traceback_nogaps(NULL,&(*nmatches),&(*nmismatches),bestr,bestc,
			     rsequence,rsequenceuc,
			     gsequence,gsequence_alt,roffset,goffset,pairpool,
#ifdef DEBUG14
			     chroffset,chrhigh,watsonp,
#endif
			     /*revp*/false,*dynprogindex);
    *finalscore = (*nmatches)*FULLMATCH + (*nmismatches)*MISMATCH_ENDQ;

#if defined(HAVE_SSE4_1) || defined(HAVE_SSE2)
  } else if (use8p == true) {
    if (bestc >= bestr) {
      pairs = Dynprog_traceback_8_upper(NULL,&(*nmatches),&(*nmismatches),&(*nopens),&(*nindels),
					directions8_upper_nogap,directions8_upper_Egap,bestr,bestc,
					rsequence,rsequenceuc,
					gsequence,gsequence_alt,roffset,goffset,pairpool,/*revp*/false,
					chroffset,chrhigh,cdna_direction,watsonp,*dynprogindex);
    } else {
      pairs = Dynprog_traceback_8_lower(NULL,&(*nmatches),&(*nmismatches),&(*nopens),&(*nindels),
					directions8_lower_nogap,directions8_lower_Egap,bestr,bestc,
					rsequence,rsequenceuc,
					gsequence,gsequence_alt,roffset,goffset,pairpool,/*revp*/false,
					chroffset,chrhigh,cdna_direction,watsonp,*dynprogindex);
    }

  } else {
    if (bestc >= bestr) {
      pairs = Dynprog_traceback_16_upper(NULL,&(*nmatches),&(*nmismatches),&(*nopens),&(*nindels),
					 directions16_upper_nogap,directions16_upper_Egap,bestr,bestc,
					 rsequence,rsequenceuc,
					 gsequence,gsequence_alt,roffset,goffset,pairpool,/*revp*/false,
					 chroffset,chrhigh,cdna_direction,watsonp,*dynprogindex);
    } else {
      pairs = Dynprog_traceback_16_lower(NULL,&(*nmatches),&(*nmismatches),&(*nopens),&(*nindels),
					 directions16_lower_nogap,directions16_lower_Egap,bestr,bestc,
					 rsequence,rsequenceuc,
					 gsequence,gsequence_alt,roffset,goffset,pairpool,/*revp*/false,
					 chroffset,chrhigh,cdna_direction,watsonp,*dynprogindex);
    }

#else
  } else {
    pairs = Dynprog_traceback_std(NULL,&(*nmatches),&(*nmismatches),&(*nopens),&(*nindels),
				  directions_nogap,directions_Egap,directions_Fgap,bestr,bestc,
				  rsequence,rsequenceuc,
				  gsequence,gsequence_alt,roffset,goffset,pairpool,/*revp*/false,
				  chroffset,chrhigh,cdna_direction,watsonp,*dynprogindex);
#endif
  }

  if ((endalign == QUERYEND_GAP || endalign == BEST_LOCAL) && (*nmatches + 1) < *nmismatches) {
    *finalscore = 0;
    /* No need to free pairs */
    pairs = NULL;

  } else {
    /* Add 1 to count the match already in the alignment */
    pairs = List_reverse(pairs); /* Look at 3' end to remove excess gaps */
    while (pairs != NULL && (pair = List_head(pairs)) && pair->comp == INDEL_COMP) {
      pairs = List_next(pairs);
    }
  }

  /*
    Directions_free(directions);
    Matrix_free(matrix);
  */

  FREEA(gsequence_alt);
  FREEA(gsequence);

  debug6(printf("End of dynprog end3 gap\n\n"));

  *dynprogindex += (*dynprogindex > 0 ? +1 : -1);
  return pairs;			/* not List_reverse(pairs) */
}


/* gsequence is the splicejunction */
List_T
Dynprog_end3_splicejunction (int *dynprogindex, int *finalscore, int *missscore,
			     int *nmatches, int *nmismatches, int *nopens, int *nindels, T dynprog, 
			     char *rsequence, char *rsequenceuc,
			     char *gsequence, char *gsequence_uc, char *gsequence_alt,
			     int rlength, int glength, int roffset, int goffset_anchor, int goffset_far,
			     Univcoord_T chroffset, Univcoord_T chrhigh,
			     int cdna_direction, bool watsonp, bool jump_late_p, Pairpool_T pairpool,
			     int extraband_end, double defect_rate, int contlength) {
  List_T pairs = NULL;
  Pair_T pair;
  Mismatchtype_T mismatchtype;
  int bestr, bestc, lband, uband;
  int open, extend;
#if defined(HAVE_SSE4_1) || defined(HAVE_SSE2)
  bool use8p = false;
  Score8_T **matrix8_upper, **matrix8_lower;
  Direction8_T **directions8_upper_nogap, **directions8_upper_Egap,
    **directions8_lower_nogap, **directions8_lower_Egap;

  Score16_T **matrix16_upper, **matrix16_lower;
  Direction16_T **directions16_upper_nogap, **directions16_upper_Egap,
    **directions16_lower_nogap, **directions16_lower_Egap;
#else
  Score32_T **matrix;
  Direction32_T **directions_nogap, **directions_Egap, **directions_Fgap;
#endif
#ifdef PMAP
  int termpos, termmod;
#endif

  debug6(
	printf("%c:  ",*dynprogindex > 0 ? (*dynprogindex-1)%26+'a' : (-(*dynprogindex)-1)%26+'A');
	printf("Aligning 3' end gap splicejunction\n");
	);

  mismatchtype = ENDQ;
  if (defect_rate < DEFECT_HIGHQ) {
    open = END_OPEN_HIGHQ;
    extend = END_EXTEND_HIGHQ;
  } else if (defect_rate < DEFECT_MEDQ) {
    open = END_OPEN_MEDQ;
    extend = END_EXTEND_MEDQ;
  } else {
    open = END_OPEN_LOWQ;
    extend = END_EXTEND_LOWQ;
  }


  /* We can just chop lengths to work, since we're not constrained on 3' end */
  if (rlength <= 0 || rlength > dynprog->max_rlength) {
    /* Needed to avoid abort by Matrix16_alloc */
    *nmatches = *nmismatches = *nopens = *nindels = 0;
    *finalscore = 0;
    *missscore = -100;
    return (List_T) NULL;
  }
  if (glength <= 0 || glength > dynprog->max_glength) {
    /* Needed to avoid abort by Matrix16_alloc */
    *nmatches = *nmismatches = *nopens = *nindels = 0;
    *finalscore = 0;
    *missscore = -100;
    return (List_T) NULL;
  }

  debug6(printf("At query offset %d-%d, %.*s\n",roffset,roffset+rlength-1,rlength,rsequence));
  debug6(printf("At genomic offset %d-%d, %.*s\n",
		goffset_anchor,goffset_anchor+glength-1,glength,gsequence));


  /* find_best_endpoint_to_queryend_nogaps(bestr,bestc,rlength,glength); */
  /* bestr = bestc = rlength; */
  /* *finalscore = 0; -- Splicetrie procedures need to know finalscore */

#ifdef PMAP
  termpos = roffset+(bestc-1);
  debug6(printf("Final query pos is %d\n",termpos));
  if ((termmod = termpos % 3) < 2) {
    if (bestr + (2 - termmod) < rlength && bestc + (2 - termmod) < glength) {
      debug6(printf("Rounding up by %d\n",2 - termmod));
      bestr += 2 - termmod;
      bestc += 2 - termmod;
    }
  }
#endif

  Dynprog_compute_bands(&lband,&uband,rlength,glength,extraband_end,/*widebandp*/true);
#if defined(HAVE_SSE4_1) || defined(HAVE_SSE2)
  /* Use || because we want the minimum length (which determines the diagonal length) to achieve a score less than 128 */
  if (rlength <= SIMD_MAXLENGTH_EPI8 || glength <= SIMD_MAXLENGTH_EPI8) {
    use8p = true;
    matrix8_upper = Dynprog_simd_8_upper(&directions8_upper_nogap,&directions8_upper_Egap,dynprog,
					 rsequenceuc,gsequence_uc,gsequence_alt,rlength,glength,
#ifdef DEBUG14
					 /*goffset*/0,chroffset,chrhigh,watsonp,
#endif
					 mismatchtype,open,extend,
					 uband,jump_late_p,/*revp*/false);
    matrix8_lower = Dynprog_simd_8_lower(&directions8_lower_nogap,&directions8_lower_Egap,dynprog,
					 rsequenceuc,gsequence_uc,gsequence_alt,rlength,glength,
#ifdef DEBUG14
					 /*goffset*/0,chroffset,chrhigh,watsonp,
#endif
					 mismatchtype,open,extend,
					 lband,jump_late_p,/*revp*/false);


    find_best_endpoint_to_queryend_indels_8(&(*finalscore),&bestr,&bestc,matrix8_upper,matrix8_lower,
					    rlength,glength,lband,uband,jump_late_p);

  } else {
    matrix16_upper = Dynprog_simd_16_upper(&directions16_upper_nogap,&directions16_upper_Egap,dynprog,
					   rsequenceuc,gsequence_uc,gsequence_alt,rlength,glength,
#ifdef DEBUG14
					   /*goffset*/0,chroffset,chrhigh,watsonp,
#endif
					   mismatchtype,open,extend,
					   uband,jump_late_p,/*revp*/false);
    matrix16_lower = Dynprog_simd_16_lower(&directions16_lower_nogap,&directions16_lower_Egap,dynprog,
					   rsequenceuc,gsequence_uc,gsequence_alt,rlength,glength,
#ifdef DEBUG14
					   /*goffset*/0,chroffset,chrhigh,watsonp,
#endif
					   mismatchtype,open,extend,
					   lband,jump_late_p,/*revp*/false);

    find_best_endpoint_to_queryend_indels_16(&(*finalscore),&bestr,&bestc,matrix16_upper,matrix16_lower,
					     rlength,glength,lband,uband,jump_late_p);
  }

#else
  /* Non-SIMD methods */
  matrix = Dynprog_standard(&directions_nogap,&directions_Egap,&directions_Fgap,dynprog,
			    rsequenceuc,gsequence_uc,gsequence_alt,rlength,glength,
			    /*goffset*/0,chroffset,chrhigh,watsonp,mismatchtype,open,extend,
			    lband,uband,jump_late_p,/*revp*/false,/*saturation*/NEG_INFINITY_INT);
  find_best_endpoint_to_queryend_indels_std(&(*finalscore),&bestr,&bestc,matrix,rlength,glength,lband,uband,
					    jump_late_p);
#endif

  if (*finalscore < 0) {
    /* Need a reasonable alignment to call a splice */
    return (List_T) NULL;

  } else {
    *nmatches = *nmismatches = *nopens = *nindels = 0;
#if defined(HAVE_SSE4_1) || defined(HAVE_SSE2)
    if (use8p == true) {
      if (bestc >= bestr) {
	pairs = traceback_local_8_upper(NULL,&(*nmatches),&(*nmismatches),&(*nopens),&(*nindels),
					directions8_upper_nogap,directions8_upper_Egap,&bestr,&bestc,/*endc*/contlength,
					rsequence,rsequenceuc,gsequence,gsequence_uc,gsequence_alt,
					roffset,goffset_far,pairpool,/*revp*/false,
					chroffset,chrhigh,cdna_direction,watsonp,*dynprogindex);
      } else {
	pairs = traceback_local_8_lower(NULL,&(*nmatches),&(*nmismatches),&(*nopens),&(*nindels),
					directions8_lower_nogap,directions8_lower_Egap,&bestr,&bestc,/*endc*/contlength,
					rsequence,rsequenceuc,gsequence,gsequence_uc,gsequence_alt,
					roffset,goffset_far,pairpool,/*revp*/false,
					chroffset,chrhigh,cdna_direction,watsonp,*dynprogindex);
      }
      pairs = Pairpool_push_gapholder(pairs,pairpool,/*queryjump*/0,/*genomejump*/goffset_far - goffset_anchor,
				      /*leftpair*/NULL,/*rightpair*/NULL,/*knownp*/true);
      if (bestc >= bestr) {
	pairs = traceback_local_8_upper(pairs,&(*nmatches),&(*nmismatches),&(*nopens),&(*nindels),
					directions8_upper_nogap,directions8_upper_Egap,&bestr,&bestc,/*endc*/0,
					rsequence,rsequenceuc,gsequence,gsequence_uc,gsequence_alt,
					roffset,goffset_anchor,pairpool,/*revp*/false,
					chroffset,chrhigh,cdna_direction,watsonp,*dynprogindex);
      } else {
	pairs = traceback_local_8_lower(pairs,&(*nmatches),&(*nmismatches),&(*nopens),&(*nindels),
					directions8_lower_nogap,directions8_lower_Egap,&bestr,&bestc,/*endc*/0,
					rsequence,rsequenceuc,gsequence,gsequence_uc,gsequence_alt,
					roffset,goffset_anchor,pairpool,/*revp*/false,
					chroffset,chrhigh,cdna_direction,watsonp,*dynprogindex);
      }

    } else {
      if (bestc >= bestr) {
	pairs = traceback_local_16_upper(NULL,&(*nmatches),&(*nmismatches),&(*nopens),&(*nindels),
					 directions16_upper_nogap,directions16_upper_Egap,&bestr,&bestc,/*endc*/contlength,
					 rsequence,rsequenceuc,gsequence,gsequence_uc,gsequence_alt,
					 roffset,goffset_far,pairpool,/*revp*/false,
					 chroffset,chrhigh,cdna_direction,watsonp,*dynprogindex);
      } else {
	pairs = traceback_local_16_lower(NULL,&(*nmatches),&(*nmismatches),&(*nopens),&(*nindels),
					 directions16_lower_nogap,directions16_lower_Egap,&bestr,&bestc,/*endc*/contlength,
					 rsequence,rsequenceuc,gsequence,gsequence_uc,gsequence_alt,
					 roffset,goffset_far,pairpool,/*revp*/false,
					 chroffset,chrhigh,cdna_direction,watsonp,*dynprogindex);
      }
      pairs = Pairpool_push_gapholder(pairs,pairpool,/*queryjump*/0,/*genomejump*/goffset_far - goffset_anchor,
				      /*leftpair*/NULL,/*rightpair*/NULL,/*knownp*/true);
      if (bestc >= bestr) {
	pairs = traceback_local_16_upper(pairs,&(*nmatches),&(*nmismatches),&(*nopens),&(*nindels),
					 directions16_upper_nogap,directions16_upper_Egap,&bestr,&bestc,/*endc*/0,
					 rsequence,rsequenceuc,gsequence,gsequence_uc,gsequence_alt,
					 roffset,goffset_anchor,pairpool,/*revp*/false,
					 chroffset,chrhigh,cdna_direction,watsonp,*dynprogindex);
      } else {
	pairs = traceback_local_16_lower(pairs,&(*nmatches),&(*nmismatches),&(*nopens),&(*nindels),
					 directions16_lower_nogap,directions16_lower_Egap,&bestr,&bestc,/*endc*/0,
					 rsequence,rsequenceuc,gsequence,gsequence_uc,gsequence_alt,
					 roffset,goffset_anchor,pairpool,/*revp*/false,
					 chroffset,chrhigh,cdna_direction,watsonp,*dynprogindex);
      }
    }

#else
    /* Non-SIMD methods */
    pairs = traceback_local_std(NULL,&(*nmatches),&(*nmismatches),&(*nopens),&(*nindels),
				directions_nogap,directions_Egap,directions_Fgap,&bestr,&bestc,/*endc*/contlength,
				rsequence,rsequenceuc,gsequence,gsequence_uc,gsequence_alt,
				roffset,goffset_far,pairpool,/*revp*/false,
				chroffset,chrhigh,cdna_direction,watsonp,*dynprogindex);
  
    pairs = Pairpool_push_gapholder(pairs,pairpool,/*queryjump*/0,/*genomejump*/goffset_far - goffset_anchor,
				    /*leftpair*/NULL,/*rightpair*/NULL,/*knownp*/true);

    pairs = traceback_local_std(pairs,&(*nmatches),&(*nmismatches),&(*nopens),&(*nindels),
				directions_nogap,directions_Egap,directions_Fgap,&bestr,&bestc,/*endc*/0,
				rsequence,rsequenceuc,gsequence,gsequence_uc,gsequence_alt,
				roffset,goffset_anchor,pairpool,/*revp*/false,
				chroffset,chrhigh,cdna_direction,watsonp,*dynprogindex);
#endif

    /* Score compared with perfect score, so heavy weight on mismatches may not be necessary */
    *finalscore = (*nmatches)*FULLMATCH + (*nmismatches)*MISMATCH_ENDQ + (*nopens)*open + (*nindels)*extend;
    *missscore = (*finalscore) - rlength*FULLMATCH;
    debug6(printf("finalscore %d = %d*%d matches + %d*%d mismatches + %d*%d opens + %d*%d extends\n",
		  *finalscore,FULLMATCH,*nmatches,MISMATCH_ENDQ,*nmismatches,open,*nopens,extend,*nindels));
    debug6(printf("missscore = %d\n",*missscore));

    /* Add 1 to count the match already in the alignment */
    pairs = List_reverse(pairs); /* Look at 3' end to remove excess gaps */
    while (pairs != NULL && (pair = List_head(pairs)) && pair->comp == INDEL_COMP) {
      pairs = List_next(pairs);
    }

    debug6(Pair_dump_list(pairs,true));
    debug6(printf("End of dynprog end3 gap splicejunction\n\n"));

    *dynprogindex += (*dynprogindex > 0 ? +1 : -1);
    return pairs;			/* not List_reverse(pairs) */
  }
}


static char complCode[128] = COMPLEMENT_LC;

static void
make_complement_inplace (char *sequence, unsigned int length) {
  char temp;
  unsigned int i, j;

  for (i = 0, j = length-1; i < length/2; i++, j--) {
    temp = complCode[(int) sequence[i]];
    sequence[i] = complCode[(int) sequence[j]];
    sequence[j] = temp;
  }
  if (i == j) {
    sequence[i] = complCode[(int) sequence[i]];
  }

  return;
}


static void
make_contjunction_5 (char *splicejunction, char *splicejunction_alt, Univcoord_T splicecoord,
		     int splicelength, int contlength, Splicetype_T anchor_splicetype,
		     bool watsonp) {
  char *proximal, *proximal_alt;

  debug7(printf("make_contjunction_5 at %u, splice (%s), contlength %d, splicelength %d:",
		splicecoord,Splicetype_string(anchor_splicetype),contlength, splicelength));

  proximal = &(splicejunction[splicelength]);
  proximal_alt = &(splicejunction_alt[splicelength]);

  if (anchor_splicetype == ACCEPTOR) {
    Genome_fill_buffer_blocks_noterm(splicecoord,contlength,proximal,proximal_alt);

  } else if (anchor_splicetype == ANTIDONOR) {
    Genome_fill_buffer_blocks_noterm(splicecoord,contlength,proximal,proximal_alt);

  } else if (anchor_splicetype == ANTIACCEPTOR) {
    Genome_fill_buffer_blocks_noterm(splicecoord-contlength,contlength,proximal,proximal_alt);

  } else if (anchor_splicetype == DONOR) {
    Genome_fill_buffer_blocks_noterm(splicecoord-contlength,contlength,proximal,proximal_alt);
    
  } else {
    fprintf(stderr,"Unexpected anchor_splicetype value %d\n",anchor_splicetype);
    abort();
  }

  if (watsonp == false) {
    make_complement_inplace(proximal,contlength);
    make_complement_inplace(proximal_alt,contlength);
  }

#ifdef DEBUG7
  if (watsonp == true) {
    printf(" (fwd)  contjunction    : %.*s\n",contlength,proximal);
    printf(" (fwd)  contjunction_alt: %.*s\n",contlength,proximal_alt);
  } else {
    printf(" (rev)  contjunction    : %.*s\n",contlength,proximal);
    printf(" (rev)  contjunction_alt: %.*s\n",contlength,proximal_alt);
  }
#endif

  return;
}



/* Fills in just the distal part, keeping the proximal part same for contlength */
bool
Dynprog_make_splicejunction_5 (char *splicejunction, char *splicejunction_alt, Univcoord_T splicecoord,
			       int splicelength, int contlength, Splicetype_T far_splicetype,
			       bool watsonp) {
  char *distal, *distal_alt;

  debug7(printf("make_splicejunction_5 at %u, splice (%s), contlength %d, splicelength %d:\n",
		splicecoord,Splicetype_string(far_splicetype),contlength, splicelength));

  distal = &(splicejunction[0]);
  distal_alt = &(splicejunction_alt[0]);

  if (far_splicetype == ACCEPTOR) {
    Genome_fill_buffer_blocks_noterm(splicecoord,splicelength,distal,distal_alt);

  } else if (far_splicetype == ANTIDONOR) {
    Genome_fill_buffer_blocks_noterm(splicecoord,splicelength,distal,distal_alt);

  } else if (splicecoord <= splicelength) {
    return false;

  } else if (far_splicetype == ANTIACCEPTOR) {
    Genome_fill_buffer_blocks_noterm(splicecoord-splicelength,splicelength,distal,distal_alt);

  } else if (far_splicetype == DONOR) {
    Genome_fill_buffer_blocks_noterm(splicecoord-splicelength,splicelength,distal,distal_alt);
    
  } else {
    fprintf(stderr,"Unexpected far_splicetype value %d\n",far_splicetype);
    abort();
  }

  if (watsonp == false) {
    make_complement_inplace(distal,splicelength);
    make_complement_inplace(distal_alt,splicelength);
  }

#ifdef DEBUG7
  if (watsonp == true) {
    printf(" (fwd)  splicejunction    : %s\n",splicejunction);
    printf(" (fwd)  splicejunction_alt: %s\n",splicejunction_alt);
  } else {
    printf(" (rev)  splicejunction    : %s\n",splicejunction);
    printf(" (rev)  splicejunction_alt: %s\n",splicejunction_alt);
  }
#endif

  return true;
}


static void
make_contjunction_3 (char *splicejunction, char *splicejunction_alt, Univcoord_T splicecoord,
		     int splicelength, int contlength, Splicetype_T anchor_splicetype,
		     bool watsonp) {
  char *proximal, *proximal_alt;

  debug7(printf("make_contjunction_3 at %u, splice (%s), contlength %d, splicelength %d:\n",
		splicecoord,Splicetype_string(anchor_splicetype),contlength,splicelength));

  proximal = &(splicejunction[0]);
  proximal_alt = &(splicejunction_alt[0]);

  if (anchor_splicetype == DONOR) {
    Genome_fill_buffer_blocks_noterm(splicecoord-contlength,contlength,proximal,proximal_alt);

  } else if (anchor_splicetype == ANTIACCEPTOR) {
    Genome_fill_buffer_blocks_noterm(splicecoord-contlength,contlength,proximal,proximal_alt);

  } else if (anchor_splicetype == ANTIDONOR) {
    Genome_fill_buffer_blocks_noterm(splicecoord,contlength,proximal,proximal_alt);

  } else if (anchor_splicetype == ACCEPTOR) {
    Genome_fill_buffer_blocks_noterm(splicecoord,contlength,proximal,proximal_alt);
    
  } else {
    fprintf(stderr,"Unexpected anchor_splicetype value %d\n",anchor_splicetype);
    abort();
  }

  if (watsonp == false) {
    make_complement_inplace(proximal,contlength);
    make_complement_inplace(proximal_alt,contlength);
  }

#ifdef DEBUG7
  if (watsonp == true) {
    printf(" (fwd)  contjunction    : %.*s\n",contlength,proximal);
    printf(" (fwd)  contjunction_alt: %.*s\n",contlength,proximal_alt);
  } else {
    printf(" (rev)  contjunction    : %.*s\n",contlength,proximal);
    printf(" (rev)  contjunction_alt: %.*s\n",contlength,proximal_alt);
  }
#endif

  return;
}



/* Fills in just the distal part, keeping the proximal part same for contlength */
bool
Dynprog_make_splicejunction_3 (char *splicejunction, char *splicejunction_alt, Univcoord_T splicecoord,
			       int splicelength, int contlength, Splicetype_T far_splicetype,
			       bool watsonp) {
  char *distal, *distal_alt;

  debug7(printf("make_splicejunction_3 at %u, splice (%s), contlength %d, splicelength %d:\n",
		splicecoord,Splicetype_string(far_splicetype),contlength,splicelength));

  distal = &(splicejunction[contlength]);
  distal_alt = &(splicejunction_alt[contlength]);

  if (far_splicetype == ANTIDONOR) {
    Genome_fill_buffer_blocks_noterm(splicecoord,splicelength,distal,distal_alt);

  } else if (far_splicetype == ACCEPTOR) {
    Genome_fill_buffer_blocks_noterm(splicecoord,splicelength,distal,distal_alt);
    
  } else if (splicecoord <= splicelength) {
    return false;

  } else if (far_splicetype == DONOR) {
    Genome_fill_buffer_blocks_noterm(splicecoord-splicelength,splicelength,distal,distal_alt);

  } else if (far_splicetype == ANTIACCEPTOR) {
    Genome_fill_buffer_blocks_noterm(splicecoord-splicelength,splicelength,distal,distal_alt);


  } else {
    fprintf(stderr,"Unexpected far_splicetype value %d\n",far_splicetype);
    abort();
  }

  if (watsonp == false) {
    make_complement_inplace(distal,splicelength);
    make_complement_inplace(distal_alt,splicelength);
  }

#ifdef DEBUG7
  if (watsonp == true) {
    printf(" (fwd)  splicejunction    : %s\n",splicejunction);
    printf(" (fwd)  splicejunction_alt: %s\n",splicejunction_alt);
  } else {
    printf(" (rev)  splicejunction    : %s\n",splicejunction);
    printf(" (rev)  splicejunction_alt: %s\n",splicejunction_alt);
  }
#endif

  return true;
}


static int
binary_search (int lowi, int highi, Univcoord_T *positions, Univcoord_T goal) {
  int middlei;

  debug10(printf("entered binary search with lowi=%d, highi=%d, goal=%u\n",lowi,highi,goal));

  while (lowi < highi) {
    middlei = lowi + ((highi - lowi) / 2);
    debug10(printf("  binary: %d:%u %d:%u %d:%u   vs. %u\n",
		  lowi,positions[lowi],middlei,positions[middlei],
		  highi,positions[highi],goal));
    if (goal < positions[middlei]) {
      highi = middlei;
    } else if (goal > positions[middlei]) {
      lowi = middlei + 1;
    } else {
      debug10(printf("binary search returns %d\n",middlei));
      return middlei;
    }
  }

  debug10(printf("binary search returns %d\n",highi));
  return highi;
}


List_T
Dynprog_end5_known (bool *knownsplicep, int *dynprogindex, int *finalscore,
		    int *ambig_end_length, Splicetype_T *ambig_splicetype,
		    int *nmatches, int *nmismatches, int *nopens, int *nindels, T dynprog, 
		    char *rev_rsequence, char *rev_rsequenceuc,
		    int rlength, int glength, int rev_roffset, int rev_goffset, 
		    Univcoord_T chroffset, Univcoord_T chrhigh,
		    Univcoord_T knownsplice_limit_low, Univcoord_T knownsplice_limit_high,
		    int cdna_direction, bool watsonp, bool jump_late_p, Pairpool_T pairpool,
		    int extraband_end, double defect_rate) {
  List_T best_pairs = NULL, orig_pairs;
  Pair_T pair;
  Univcoord_T low, high, far_limit_low, far_limit_high;
  Splicetype_T anchor_splicetype, far_splicetype;
  int contlength, splicelength, endlength;
  char *splicejunction, *splicejunction_alt;
#ifdef EXTRACT_GENOMICSEG
  char *splicejunction_test;
#endif
  int jstart, j;

  int orig_score, threshold_miss_score, perfect_score;
  int obsmax_penalty;


  assert(glength >= rlength);

  debug7(
	printf("%c:  ",*dynprogindex > 0 ? (*dynprogindex-1)%26+'a' : (-(*dynprogindex)-1)%26+'A');
	printf("Aligning 5' end gap, known\n")
	);

  *ambig_end_length = 0;

  /* We can just chop lengths to work, since we're not constrained on 5' end */
  if (rlength <= 0) {
    /* Needed to avoid abort by Matrix16_alloc */
    *finalscore = 0;
    *knownsplicep = false;
    return (List_T) NULL;
  }
  if (glength <= 0) {
    /* Needed to avoid abort by Matrix16_alloc */
    *finalscore = 0;
    *knownsplicep = false;
    return (List_T) NULL;
  }

  debug7(printf("At query offset %d-%d, %.*s\n",rev_roffset-rlength+1,rev_roffset,rlength,&(rev_rsequence[-rlength+1])));
#ifdef EXTRACT_GENOMICSEG
  debug7(printf("At genomic offset %d-%d, %.*s\n",rev_goffset-glength+1,rev_goffset,glength,&(rev_gsequence[-glength+1])));
#endif

  perfect_score = rlength*FULLMATCH;

  /* Try without splicing, all the way to query end */
  best_pairs = Dynprog_end5_gap(&(*dynprogindex),&(*finalscore),&(*nmatches),&(*nmismatches),
				&(*nopens),&(*nindels),dynprog,rev_rsequence,rev_rsequenceuc,
				rlength,glength,rev_roffset,rev_goffset,chroffset,chrhigh,
				cdna_direction,watsonp,jump_late_p,pairpool,
				extraband_end,defect_rate,/*endalign*/QUERYEND_NOGAPS);
  if (*finalscore < 0) {
    orig_score = 0;
    orig_pairs = best_pairs = (List_T) NULL;
  } else {
    orig_score = *finalscore;
    orig_pairs = best_pairs;
  }
  threshold_miss_score = orig_score - perfect_score;
  debug7(printf("score %d - perfect score %d = threshold %d",
		orig_score,perfect_score,threshold_miss_score));
  if (threshold_miss_score < -2*FULLMATCH) {
    /* Don't allow more than 2 mismatches in a distant splice */
    threshold_miss_score = -2*FULLMATCH;
    debug7(printf(", but revising to %d\n",threshold_miss_score));
  }
  debug7(printf("\n"));
  *knownsplicep = false;


  if (threshold_miss_score < 0 && glength > 0) {
    /* Try known splicing */
    splicejunction = (char *) MALLOCA((glength+1) * sizeof(char));
    splicejunction_alt = (char *) MALLOCA((glength+1) * sizeof(char));
#ifdef EXTRACT_GENOMICSEG
    splicejunction_test = (char *) MALLOCA((glength+1) * sizeof(char));
#endif

    endlength = rlength;
    if (watsonp == true) {
      low = chroffset + rev_goffset-endlength + 2;
      high = chroffset + rev_goffset + 1;
      debug7(printf("5' watson\n"));
      debug7(printf("Calculating low %u (%u) = %u + %d-%d + 2\n",
		    low,low-chroffset,chroffset,rev_goffset,endlength));
      debug7(printf("Calculating high %u (%u) = %u + %d + 1\n",
		    high,high-chroffset,chroffset,rev_goffset));
      if (cdna_direction > 0) {
	anchor_splicetype = ACCEPTOR;
	far_splicetype = DONOR;
      } else {
	anchor_splicetype = ANTIDONOR;
	far_splicetype = ANTIACCEPTOR;
      }
    } else {
      low = chrhigh - rev_goffset;
      high = chrhigh - (rev_goffset-endlength) - 1;
      debug7(printf("5' crick\n"));
      debug7(printf("Calculating low %u (%u) = %u - %d\n",
		    low,low-chroffset,chrhigh,rev_goffset));
      debug7(printf("Calculating high %u (%u) = %u - (%d-%d) - 1\n",
		    high,high-chroffset,chrhigh,rev_goffset,endlength));
      if (cdna_direction > 0) {
	anchor_splicetype = ANTIACCEPTOR;
	far_splicetype = ANTIDONOR;
      } else {
	anchor_splicetype = DONOR;
	far_splicetype = ACCEPTOR;
      }
    }

    far_limit_low = knownsplice_limit_low;
    far_limit_high = knownsplice_limit_high;
    debug7(printf("Genomic positions: %u..%u (%u..%u), looking for anchor splicetype %s\n",
		  low,high,low-chroffset,high-chroffset,Splicetype_string(anchor_splicetype)));
    j = jstart = binary_search(0,nsplicesites,splicesites,low);
    while (j < nsplicesites && splicesites[j] <= high) {
      if (splicetypes[j] == anchor_splicetype) {
	debug7(printf("Found one at %u (%u)\n",splicesites[j],splicesites[j]-chroffset));
	if (watsonp == true) {
	  contlength = high - splicesites[j];
	} else {
	  contlength = splicesites[j] - low;
	}
	debug7(printf("contlength %d, splicelength %d, rlength %d, glength %d\n",
		      contlength,glength-contlength,rlength,glength));
	assert(contlength >= 0 && contlength < rlength);

#ifdef EXTRACT_GENOMICSEG
	debug7(printf("cont: %.*s\n",contlength,&(rev_gsequence[-contlength+1])));
#endif
	splicelength = glength - contlength;
	assert(splicelength > 0);
	debug7(printf("  Saw %u (%u) of type %s (cont length %d, splice length %d)\n",
		      splicesites[j],splicesites[j]-chroffset,Splicetype_string(splicetypes[j]),contlength,splicelength));

	make_contjunction_5(splicejunction,splicejunction_alt,splicesites[j],splicelength,contlength,anchor_splicetype,watsonp);
#ifdef EXTRACT_GENOMICSEG
	strncpy(&(splicejunction_test[splicelength]),&(rev_gsequence[-contlength+1]),contlength);
	debug7(printf("contjunction_gen:  %s\n",&(splicejunction[splicelength])));
	debug7(printf("contjunction_test: %s\n",&(splicejunction_test[splicelength])));
	assert(!strncmp(&(splicejunction[splicelength]),&(splicejunction_test[splicelength]),contlength));
#endif

	if (watsonp) {
	  far_limit_high = splicesites[j];
	} else {
	  far_limit_low = splicesites[j];
	}

	obsmax_penalty = 0;
	if (trieoffsets_obs != NULL) {
	  debug7(printf("  Running Splicetrie_solve_end5 on observed splice sites with rev_goffset %d\n",rev_goffset));
	  best_pairs = Splicetrie_solve_end5(best_pairs,triecontents_obs,trieoffsets_obs,j,
					     far_limit_low,far_limit_high,
					     &(*finalscore),&(*nmatches),&(*nmismatches),
					     &(*nopens),&(*nindels),&(*knownsplicep),&(*ambig_end_length),
					     &threshold_miss_score,/*obsmax_penalty*/0,perfect_score,
					     /*anchor_splicesite*/splicesites[j],splicejunction,splicejunction_alt,
					     splicelength,contlength,far_splicetype,
					     chroffset,chrhigh,&(*dynprogindex),dynprog,
					     rev_rsequence,rev_rsequenceuc,rlength,glength,rev_roffset,rev_goffset,
					     cdna_direction,watsonp,jump_late_p,pairpool,extraband_end,defect_rate);
	  debug7(printf("  Result on obs with ambig_end_length_5 %d\n",*ambig_end_length));
	  debug7(Pair_dump_list(best_pairs,/*zerobasedp*/true));
	  obsmax_penalty += FULLMATCH;
	}

	if (threshold_miss_score + obsmax_penalty < 0 && trieoffsets_max != NULL) {
	  debug7(printf("  Running Splicetrie_solve_end5 on maxdistance splice sites with rev_goffset %d\n",rev_goffset));
	  best_pairs = Splicetrie_solve_end5(best_pairs,triecontents_max,trieoffsets_max,j,
					     far_limit_low,far_limit_high,
					     &(*finalscore),&(*nmatches),&(*nmismatches),
					     &(*nopens),&(*nindels),&(*knownsplicep),&(*ambig_end_length),
					     &threshold_miss_score,obsmax_penalty,perfect_score,
					     /*anchor_splicesite*/splicesites[j],splicejunction,splicejunction_alt,
					     splicelength,contlength,far_splicetype,
					     chroffset,chrhigh,&(*dynprogindex),dynprog,
					     rev_rsequence,rev_rsequenceuc,rlength,glength,rev_roffset,rev_goffset,
					     cdna_direction,watsonp,jump_late_p,pairpool,extraband_end,defect_rate);
	  debug7(printf("  Result on max with ambig_end_length_5 %d\n",*ambig_end_length));
	  debug7(Pair_dump_list(best_pairs,/*zerobasedp*/true));
	}
      }
      j++;
    }

#ifdef EXTRACT_GENOMICSEG
    FREEA(splicejunction_test);
#endif
    FREEA(splicejunction_alt);
    FREEA(splicejunction);
  }


  if (best_pairs == NULL) {
    if (*ambig_end_length == 0) {
      /* Don't go to query end this time */
      if (rlength > dynprog->max_rlength) {
	debug7(printf("rlength %d is too long.  Chopping to %d\n",rlength,dynprog->max_rlength));
	rlength = dynprog->max_rlength;
      }
      if (glength > dynprog->max_glength) {
	debug7(printf("glength %d is too long.  Chopping to %d\n",glength,dynprog->max_glength));
	glength = dynprog->max_glength;
      }
      orig_pairs = Dynprog_end5_gap(&(*dynprogindex),&(*finalscore),&(*nmatches),&(*nmismatches),
				    &(*nopens),&(*nindels),dynprog,rev_rsequence,rev_rsequenceuc,
				    rlength,glength,rev_roffset,rev_goffset,chroffset,chrhigh,
				    cdna_direction,watsonp,jump_late_p,pairpool,
				    extraband_end,defect_rate,/*endalign*/BEST_LOCAL);
      debug7(Pair_dump_list(orig_pairs,/*zerobasedp*/true));
      debug7(printf("End of dynprog end5 known\n"));
      *knownsplicep = false;
      return orig_pairs;

    } else {
      *ambig_splicetype = anchor_splicetype;
      debug7(printf("Final result: best_pairs is NULL.  ambig_end_length is %d.  ambig_splicetype is %s.  Result after truncate:\n",
		    *ambig_end_length,Splicetype_string(*ambig_splicetype)));
      /* Truncate ambiguous part.  querypos is decreasing. */
      orig_pairs = List_reverse(orig_pairs);
      while (orig_pairs != NULL && ((Pair_T) orig_pairs->first)->querypos < *ambig_end_length) {
	orig_pairs = Pairpool_pop(orig_pairs,&pair);
      }
      orig_pairs = List_reverse(orig_pairs);
      debug7(Pair_dump_list(orig_pairs,/*zerobasedp*/true));
      *knownsplicep = false;
      *finalscore = orig_score;
      debug7(printf("End of dynprog end5 known\n"));
      return orig_pairs;
    }

  } else {
    debug7(printf("Found a best splice\n"));
    *ambig_end_length = 0;
    debug7(printf("End of dynprog end5 known\n"));
    if (*knownsplicep == true) {
      return Pair_protect_end5(best_pairs,pairpool);
    } else {
      return best_pairs;
    }
  }
}


List_T
Dynprog_end3_known (bool *knownsplicep, int *dynprogindex, int *finalscore,
		    int *ambig_end_length, Splicetype_T *ambig_splicetype,
		    int *nmatches, int *nmismatches, int *nopens, int *nindels, T dynprog, 
		    char *rsequence, char *rsequenceuc,
		    int rlength, int glength, int roffset, int goffset, int querylength,
		    Univcoord_T chroffset, Univcoord_T chrhigh,
		    Univcoord_T knownsplice_limit_low, Univcoord_T knownsplice_limit_high,
		    int cdna_direction, bool watsonp, bool jump_late_p, Pairpool_T pairpool,
		    int extraband_end, double defect_rate) {
  List_T best_pairs = NULL, orig_pairs;
  Pair_T pair;
  Univcoord_T low, high, far_limit_low, far_limit_high;
  Splicetype_T anchor_splicetype, far_splicetype;
  int contlength, splicelength, endlength;
  char *splicejunction, *splicejunction_alt;
#ifdef EXTRACT_GENOMICSEG
  char *splicejunction_test;
#endif
  int jstart, j;

  int orig_score, threshold_miss_score, perfect_score;
  int obsmax_penalty;


  assert(glength >= rlength);

  debug7(
	printf("%c:  ",*dynprogindex > 0 ? (*dynprogindex-1)%26+'a' : (-(*dynprogindex)-1)%26+'A');
	printf("Aligning 3' end gap, known\n")
	);

  *ambig_end_length = 0;

  /* We can just chop lengths to work, since we're not constrained on 3' end */
  if (rlength <= 0) {
    /* Needed to avoid abort by Matrix16_alloc */
    *finalscore = 0;
    *knownsplicep = false;
    return (List_T) NULL;
  }
  if (glength <= 0) {
    /* Needed to avoid abort by Matrix16_alloc */
    *finalscore = 0;
    *knownsplicep = false;
    return (List_T) NULL;
  }


  debug7(printf("At query offset %d-%d, %.*s\n",roffset,roffset+rlength-1,rlength,rsequence));
#ifdef EXTRACT_GENOMICSEG
  debug7(printf("At genomic offset %d-%d, %.*s\n",goffset,goffset+glength-1,glength,gsequence));
#endif

  perfect_score = rlength*FULLMATCH;

  /* Try without splicing, all the way to query end */
  best_pairs = Dynprog_end3_gap(&(*dynprogindex),&(*finalscore),&(*nmatches),&(*nmismatches),
				&(*nopens),&(*nindels),dynprog,rsequence,rsequenceuc,
				rlength,glength,roffset,goffset,chroffset,chrhigh,
				cdna_direction,watsonp,jump_late_p,pairpool,
				extraband_end,defect_rate,/*endalign*/QUERYEND_NOGAPS);
  if (*finalscore < 0) {
    orig_score = 0;
    orig_pairs = best_pairs = (List_T) NULL;
  } else {
    orig_score = *finalscore;
    orig_pairs = best_pairs;
  }
  threshold_miss_score = orig_score - perfect_score;
  debug7(printf("score %d - perfect score %d = threshold %d",
		orig_score,perfect_score,threshold_miss_score));
  if (threshold_miss_score < -2*FULLMATCH) {
    /* Don't allow more than 2 mismatches in a distant splice */
    threshold_miss_score = -2*FULLMATCH;
    debug7(printf(", but revising to %d\n",threshold_miss_score));
  }
  debug7(printf("\n"));
  *knownsplicep = false;


  if (threshold_miss_score < 0 && glength > 0) {
    /* Try known splicing */
    splicejunction = (char *) MALLOCA((glength+1) * sizeof(char));
    splicejunction_alt = (char *) MALLOCA((glength+1) * sizeof(char));
#ifdef EXTRACT_GENOMICSEG
    splicejunction_test = (char *) MALLOCA((glength+1) * sizeof(char));
#endif

    endlength = rlength;
    if (watsonp == true) {
      low = chroffset + goffset;
      high = chroffset + goffset+endlength - 1;
      debug7(printf("3' watson\n"));
      debug7(printf("Calculating low %u (%u) = %u + %d\n",
		    low,low-chroffset,chroffset,goffset));
      debug7(printf("Calculating high %u (%u) = %u + %d+%d - 1\n",
		    high,high-chroffset,chroffset,goffset,endlength));
      if (cdna_direction > 0) {
	anchor_splicetype = DONOR;
	far_splicetype = ACCEPTOR;
      } else {
	anchor_splicetype = ANTIACCEPTOR;
	far_splicetype = ANTIDONOR;
      }
    } else {
      low = chrhigh - (goffset+endlength) + 2;
      high = chrhigh - goffset + 1;
      debug7(printf("3' crick\n"));
      debug7(printf("Calculating low %u (%u) = %u - (%d+%d)\n",
		    low,low-chroffset,chrhigh,goffset,endlength));
      debug7(printf("Calculating high %u (%u) = %u - %d + 1\n",
		    high,high-chroffset,chrhigh,goffset));
      if (cdna_direction > 0) {
	anchor_splicetype = ANTIDONOR;
	far_splicetype = ANTIACCEPTOR;
      } else {
	anchor_splicetype = ACCEPTOR;
	far_splicetype = DONOR;
      }
    }

    far_limit_low = knownsplice_limit_low;
    far_limit_high = knownsplice_limit_high;
    debug7(printf("Genomic positions: %u..%u (%u..%u), looking for anchor splicetype %s\n",
		  low,high,low-chroffset,high-chroffset,Splicetype_string(anchor_splicetype)));
    j = jstart = binary_search(0,nsplicesites,splicesites,low);
    while (j < nsplicesites && splicesites[j] <= high) {
      if (splicetypes[j] == anchor_splicetype) {
	debug7(printf("Found one at %u (%u)\n",splicesites[j],splicesites[j]-chroffset));
	if (watsonp == true) {
	  contlength = splicesites[j] - low;
	} else {
	  contlength = high - splicesites[j];
	}
	debug7(printf("contlength %d, splicelength %d, rlength %d, glength %d\n",
		      contlength,glength-contlength,rlength,glength));
	assert(contlength >= 0 && contlength < rlength);

#ifdef EXTRACT_GENOMICSEG
	debug7(printf("cont: %.*s\n",contlength,gsequence));
#endif
	splicelength = glength - contlength;
	assert(splicelength > 0);
	debug7(printf("  Saw %u (%u) of type %s (cont length %d, splice length %d)\n",
		      splicesites[j],splicesites[j]-chroffset,Splicetype_string(splicetypes[j]),contlength,splicelength));

	make_contjunction_3(splicejunction,splicejunction_alt,splicesites[j],splicelength,contlength,anchor_splicetype,watsonp);
#ifdef EXTRACT_GENOMICSEG
	strncpy(splicejunction_test,gsequence,contlength);
	debug7(printf("contjunction_gen:  %s\n",splicejunction));
	debug7(printf("contjunction_test: %s\n",splicejunction_test));
	assert(!strncmp(splicejunction,splicejunction_test,contlength));
#endif

	if (watsonp) {
	  far_limit_low = splicesites[j];
	} else {
	  far_limit_high = splicesites[j];
	}

	obsmax_penalty = 0;
	if (trieoffsets_obs != NULL) {
	  debug7(printf("  Running Splicetrie_solve_end3 on observed splice sites with goffset %d\n",goffset));
	  best_pairs = Splicetrie_solve_end3(best_pairs,triecontents_obs,trieoffsets_obs,j,
					     far_limit_low,far_limit_high,
					     &(*finalscore),&(*nmatches),&(*nmismatches),
					     &(*nopens),&(*nindels),&(*knownsplicep),&(*ambig_end_length),
					     &threshold_miss_score,/*obsmax_penalty*/0,perfect_score,
					     /*anchor_splicesite*/splicesites[j],splicejunction,splicejunction_alt,
					     splicelength,contlength,far_splicetype,
					     chroffset,chrhigh,&(*dynprogindex),dynprog,
					     rsequence,rsequenceuc,rlength,glength,roffset,goffset,
					     cdna_direction,watsonp,jump_late_p,pairpool,extraband_end,defect_rate);
	  debug7(printf("  Result on obs with ambig_end_length_3 %d\n",*ambig_end_length));
	  debug7(Pair_dump_list(best_pairs,/*zerobasedp*/true));
	  obsmax_penalty += FULLMATCH;
	}

	if (threshold_miss_score + obsmax_penalty < 0 && trieoffsets_max != NULL) {
	  debug7(printf("  Running Splicetrie_solve_end3 on maxdistance splice sites with goffset %d\n",goffset));
	  best_pairs = Splicetrie_solve_end3(best_pairs,triecontents_max,trieoffsets_max,j,
					     far_limit_low,far_limit_high,
					     &(*finalscore),&(*nmatches),&(*nmismatches),
					     &(*nopens),&(*nindels),&(*knownsplicep),&(*ambig_end_length),
					     &threshold_miss_score,obsmax_penalty,perfect_score,
					     /*anchor_splicesite*/splicesites[j],splicejunction,splicejunction_alt,
					     splicelength,contlength,far_splicetype,
					     chroffset,chrhigh,&(*dynprogindex),dynprog,
					     rsequence,rsequenceuc,rlength,glength,roffset,goffset,
					     cdna_direction,watsonp,jump_late_p,pairpool,extraband_end,defect_rate);
	  debug7(printf("  Result on max with ambig_end_length_3 %d\n",*ambig_end_length));
	  debug7(Pair_dump_list(best_pairs,/*zerobasedp*/true));
	}
      }
      j++;
    }

#ifdef EXTRACT_GENOMICSEG
    FREEA(splicejunction_test);
#endif
    FREEA(splicejunction_alt);
    FREEA(splicejunction);
  }


  if (best_pairs == NULL) {
    if (*ambig_end_length == 0) {
      /* Don't go to query end this time */
      if (rlength > dynprog->max_rlength) {
	debug7(printf("rlength %d is too long.  Chopping to %d\n",rlength,dynprog->max_rlength));
	rlength = dynprog->max_rlength;
      }
      if (glength > dynprog->max_glength) {
	debug7(printf("glength %d is too long.  Chopping to %d\n",glength,dynprog->max_glength));
	glength = dynprog->max_glength;
      }
      orig_pairs = Dynprog_end3_gap(&(*dynprogindex),&(*finalscore),&(*nmatches),&(*nmismatches),
				    &(*nopens),&(*nindels),dynprog,rsequence,rsequenceuc,
				    rlength,glength,roffset,goffset,chroffset,chrhigh,
				    cdna_direction,watsonp,jump_late_p,pairpool,
				    extraband_end,defect_rate,/*endalign*/BEST_LOCAL);
      debug7(Pair_dump_list(orig_pairs,/*zerobasedp*/true));
      *knownsplicep = false;
      debug7(printf("End of dynprog end5 known\n"));
      return orig_pairs;

    } else {
      *ambig_splicetype = anchor_splicetype;
      debug7(printf("Final result: best_pairs is NULL.  ambig_end_length is %d.  ambig_splicetype is %s.  Result after truncate:\n",
		    *ambig_end_length,Splicetype_string(*ambig_splicetype)));
      /* Truncate ambiguous part.  querypos is decreasing */
      while (orig_pairs != NULL && ((Pair_T) orig_pairs->first)->querypos >= querylength - *ambig_end_length) {
	orig_pairs = Pairpool_pop(orig_pairs,&pair);
      }
      debug7(Pair_dump_list(orig_pairs,/*zerobasedp*/true));
      *knownsplicep = false;
      *finalscore = orig_score;
      debug7(printf("End of dynprog end5 known\n"));
      return orig_pairs;
    }

  } else {
    debug7(printf("Found a best splice\n"));
    *ambig_end_length = 0;
    debug7(printf("End of dynprog end3 known\n"));
    if (*knownsplicep == true) {
      return Pair_protect_end3(best_pairs,pairpool);
    } else {
      return best_pairs;
    }
  }
}


