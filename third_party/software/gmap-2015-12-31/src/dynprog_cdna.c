static char rcsid[] = "$Id: dynprog_cdna.c 145990 2014-08-25 21:47:32Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "dynprog_genome.h"
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
#include "dynprog_simd.h"

/* These values were set to -5, -4, -3, but this led to chopped ends
   in GMAP alignments, and failure to find chimeras */
#define MISMATCH_HIGHQ -3
#define MISMATCH_MEDQ -2
#define MISMATCH_LOWQ -1

/* cDNA insertions are biologically not meaningful, so look for a good
   gap opening somewhere */
#define CDNA_OPEN_HIGHQ -10
#define CDNA_OPEN_MEDQ -10
#define CDNA_OPEN_LOWQ -10

#define CDNA_EXTEND_HIGHQ -7
#define CDNA_EXTEND_MEDQ -7
#define CDNA_EXTEND_LOWQ -7

#define INSERT_PAIRS 9


#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif

#ifdef DEBUG3
#define debug3(x) x
#else
#define debug3(x)
#endif

#ifdef DEBUG8
#define debug8(x) x
#else
#define debug8(x)
#endif

#define T Dynprog_T


/************************************************************************
 * get_genomic_nt
 ************************************************************************/

static char complCode[128] = COMPLEMENT_LC;

static char
get_genomic_nt (char *g_alt, int genomicpos, Univcoord_T chroffset, Univcoord_T chrhigh,
		bool watsonp) {
  char c2, c2_alt;
  Univcoord_T pos;

#if 0
  /* If the read has a deletion, then we will extend beyond 0 or genomiclength, so do not restrict. */
  if (genomicpos < 0) {
    return '*';

  } else if (genomicpos >= genomiclength) {
    return '*';

  }
#endif

  if (watsonp) {
    if ((pos = chroffset + genomicpos) < chroffset) { /* Must be <, and not <=, or dynamic programming will fail */
      *g_alt = '*';
      return '*';

    } else if (pos >= chrhigh) {
      *g_alt = '*';
      return '*';

#if 0
    } else if (genome) {
      /* Not necessary, because Genome_get_char_blocks should work */
      debug8(printf("At %u, genomicnt is %c\n",
		    genomicpos,Genome_get_char(genome,pos)));
      return Genome_get_char(genome,pos);
#endif

    } else {
      /* GMAP with user-supplied genomic segment */
      debug8(printf("At %u, genomicnt is %c\n",
		    genomicpos,Genome_get_char_blocks(&(g_alt),pos)));
      return Genome_get_char_blocks(&(*g_alt),pos);
    }

  } else {
    if ((pos = chrhigh - genomicpos) < chroffset) { /* Must be <, and not <=, or dynamic programming will fail */
      *g_alt = '*';
      return '*';

    } else if (pos >= chrhigh) {
      *g_alt = '*';
      return '*';

#if 0
    } else if (genome) {
      /* Not necessary, because Genome_get_char_blocks should work */
      c2 = Genome_get_char(genome,pos);
#endif

    } else {
      /* GMAP with user-supplied genomic segment */
      c2 = Genome_get_char_blocks(&c2_alt,pos);
    }
    debug8(printf("At %u, genomicnt is %c\n",genomicpos,complCode[(int) c2]));
    *g_alt = complCode[(int) c2_alt];
    return complCode[(int) c2];
  }
}


#if defined(HAVE_SSE4_1) || defined(HAVE_SSE2)
/* Columns are always genomic.  Rows are always query.  Bridging across common columns */
static void
bridge_cdna_gap_8_ud (int *finalscore, int *bestcL, int *bestcR, int *bestrL, int *bestrR,
		      Score8_T **matrixL_upper, Score8_T **matrixL_lower,
		      Score8_T **matrixR_upper, Score8_T **matrixR_lower,
		      int glength, int rlengthL, int rlengthR,
		      int lbandL, int ubandL, int lbandR, int ubandR,
		      int open, int extend, int leftoffset, int rightoffset, bool jump_late_p) {
  int bestscore = NEG_INFINITY_16, score, scoreL, scoreR, pen, end_reward = 0;
  int rL, rR, cL, cR;
  int rloL, rhighL;
  int rloR, rhighR;

#if 0
  /* Perform computations */
  lbandL = rlengthL - glength + extraband_paired;
  ubandL = extraband_paired;

  lbandR = rlengthR - glength + extraband_paired;
  ubandR = extraband_paired;
#endif

  /* Need a double loop on rows here, in contrast with a single loop
     for introns, because we allow a genomic insertion that doesn't
     match the cDNA.  So we need to add a penalty for a genomic
     insertion */

  if (jump_late_p) {
    for (cL = 1; cL < glength; cL++) {

      /* Note: opening penalty is added at the bottom of the loop */
      for (cR = glength-cL, pen = 0; cR >= 0; cR--, pen += extend) {
	/* debug3(printf("\nAt row %d on left and %d on right\n",cL,cR)); */
	if ((rloL = cL - ubandL) < 1) {
	  rloL = 1;
	}
	if ((rhighL = cL + lbandL) > rlengthL-1) {
	  rhighL = rlengthL-1;
	}

	if ((rloR = cR - ubandR) < 1) {
	  rloR = 1;
	}
	if ((rhighR = cR + lbandR) > rlengthR-1) {
	  rhighR = rlengthR-1;
	}

	for (rL = rloL; rL < /*to main diagonal*/cL; rL++) {
	  scoreL = (int) matrixL_upper[cL][rL];
	
	  /* Disallow leftoffset + rL >= rightoffset - rR, or rR >= rightoffset - leftoffset - rL */
	  debug3(printf("  Disallowing rR to be >= %d\n",rightoffset-leftoffset-rL));
	  for (rR = rloR; rR < /*to main diagonal*/cR && rR < rightoffset-leftoffset-rL; rR++) {
	    scoreR = (int) matrixR_upper[cR][rR];

	    if ((score = scoreL + scoreR + pen + end_reward) >= bestscore) {  /* Use >= for jump late */
	      debug3(printf("At %d left to %d right, score is (%d)+(%d)+(%d)+(%d) = %d (bestscore)\n",
			    rL,rR,scoreL,scoreR,pen,end_reward,score));

	      bestscore = score;
	      *bestcL = cL;
	      *bestcR = cR;
	      *bestrL = rL;
	      *bestrR = rR;

	    } else {
	      debug3(printf("At %d left to %d right, score is (%d)+(%d)+(%d) = %d\n",
			    rL,rR,scoreL,scoreR,pen,scoreL+scoreR+pen));
	    }
	  }

	  for (/*at main diagonal*/; rR <= rhighR && rR < rightoffset-leftoffset-rL; rR++) {
	    scoreR = (int) matrixR_lower[rR][cR];

	    if ((score = scoreL + scoreR + pen + end_reward) >= bestscore) {  /* Use >= for jump late */
	      debug3(printf("At %d left to %d right, score is (%d)+(%d)+(%d)+(%d) = %d (bestscore)\n",
			    rL,rR,scoreL,scoreR,pen,end_reward,score));

	      bestscore = score;
	      *bestcL = cL;
	      *bestcR = cR;
	      *bestrL = rL;
	      *bestrR = rR;

	    } else {
	      debug3(printf("At %d left to %d right, score is (%d)+(%d)+(%d) = %d\n",
			    rL,rR,scoreL,scoreR,pen,scoreL+scoreR+pen));
	    }
	  }
	}

	for (/*at main diagonal*/; rL <= rhighL; rL++) {
	  scoreL = (int) matrixL_lower[rL][cL];
	
	  /* Disallow leftoffset + rL >= rightoffset - rR, or rR >= rightoffset - leftoffset - rL */
	  debug3(printf("  Disallowing rR to be >= %d\n",rightoffset-leftoffset-rL));
	  for (rR = rloR; rR < /*to main diagonal*/cR && rR < rightoffset-leftoffset-rL; rR++) {
	    scoreR = (int) matrixR_upper[cR][rR];

	    if ((score = scoreL + scoreR + pen + end_reward) >= bestscore) {  /* Use >= for jump late */
	      debug3(printf("At %d left to %d right, score is (%d)+(%d)+(%d)+(%d) = %d (bestscore)\n",
			    rL,rR,scoreL,scoreR,pen,end_reward,score));

	      bestscore = score;
	      *bestcL = cL;
	      *bestcR = cR;
	      *bestrL = rL;
	      *bestrR = rR;

	    } else {
	      debug3(printf("At %d left to %d right, score is (%d)+(%d)+(%d) = %d\n",
			    rL,rR,scoreL,scoreR,pen,scoreL+scoreR+pen));
	    }
	  }

	  for (/*at main diagonal*/; rR <= rhighR && rR < rightoffset-leftoffset-rL; rR++) {
	    scoreR = (int) matrixR_lower[rR][cR];

	    if ((score = scoreL + scoreR + pen + end_reward) >= bestscore) {  /* Use >= for jump late */
	      debug3(printf("At %d left to %d right, score is (%d)+(%d)+(%d)+(%d) = %d (bestscore)\n",
			    rL,rR,scoreL,scoreR,pen,end_reward,score));

	      bestscore = score;
	      *bestcL = cL;
	      *bestcR = cR;
	      *bestrL = rL;
	      *bestrR = rR;

	    } else {
	      debug3(printf("At %d left to %d right, score is (%d)+(%d)+(%d) = %d\n",
			    rL,rR,scoreL,scoreR,pen,scoreL+scoreR+pen));
	    }
	  }
	}

	pen = open - extend;	/* Subtract extend to compensate for
                                   its addition in the for loop */
      }
    }

  } else {
    /* Do not jump late */
    for (cL = 1; cL < glength; cL++) {

      /* Note: opening penalty is added at the bottom of the loop */
      for (cR = glength-cL, pen = 0; cR >= 0; cR--, pen += extend) {
	/* debug3(printf("\nAt row %d on left and %d on right\n",cL,cR)); */
	if ((rloL = cL - ubandL) < 1) {
	  rloL = 1;
	}
	if ((rhighL = cL + lbandL) > rlengthL-1) {
	  rhighL = rlengthL-1;
	}

	if ((rloR = cR - ubandR) < 1) {
	  rloR = 1;
	}
	if ((rhighR = cR + lbandR) > rlengthR-1) {
	  rhighR = rlengthR-1;
	}

	for (rL = rloL; rL < /*to main diagonal*/cL; rL++) {
	  scoreL = (int) matrixL_upper[cL][rL];
	
	  /* Disallow leftoffset + rL >= rightoffset - rR, or rR >= rightoffset - leftoffset - rL */
	  debug3(printf("  Disallowing rR to be >= %d\n",rightoffset-leftoffset-rL));
	  for (rR = rloR; rR < /*to main diagonal*/cR && rR < rightoffset-leftoffset-rL; rR++) {
	    scoreR = (int) matrixR_upper[cR][rR];

	    if ((score = scoreL + scoreR + pen + end_reward) > bestscore) {  /* Use > for jump early */
	      debug3(printf("At %d left to %d right, score is (%d)+(%d)+(%d)+(%d) = %d (bestscore)\n",
			    rL,rR,scoreL,scoreR,pen,end_reward,score));

	      bestscore = score;
	      *bestcL = cL;
	      *bestcR = cR;
	      *bestrL = rL;
	      *bestrR = rR;

	    } else {
	      debug3(printf("At %d left to %d right, score is (%d)+(%d)+(%d) = %d\n",
			    rL,rR,scoreL,scoreR,pen,scoreL+scoreR+pen));
	    }
	  }

	  for (/*at main diagonal*/; rR <= rhighR && rR < rightoffset-leftoffset-rL; rR++) {
	    scoreR = (int) matrixR_lower[rR][cR];

	    if ((score = scoreL + scoreR + pen + end_reward) > bestscore) {  /* Use > for jump early */
	      debug3(printf("At %d left to %d right, score is (%d)+(%d)+(%d)+(%d) = %d (bestscore)\n",
			    rL,rR,scoreL,scoreR,pen,end_reward,score));

	      bestscore = score;
	      *bestcL = cL;
	      *bestcR = cR;
	      *bestrL = rL;
	      *bestrR = rR;

	    } else {
	      debug3(printf("At %d left to %d right, score is (%d)+(%d)+(%d) = %d\n",
			    rL,rR,scoreL,scoreR,pen,scoreL+scoreR+pen));
	    }
	  }
	}

	for (/*at main diagonal*/; rL <= rhighL; rL++) {
	  scoreL = (int) matrixL_lower[rL][cL];
	
	  /* Disallow leftoffset + rL >= rightoffset - rR, or rR >= rightoffset - leftoffset - rL */
	  debug3(printf("  Disallowing rR to be >= %d\n",rightoffset-leftoffset-rL));
	  for (rR = rloR; rR < /*to main diagonal*/cR && rR < rightoffset-leftoffset-rL; rR++) {
	    scoreR = (int) matrixR_upper[cR][rR];

	    if ((score = scoreL + scoreR + pen + end_reward) > bestscore) {  /* Use > for jump early */
	      debug3(printf("At %d left to %d right, score is (%d)+(%d)+(%d)+(%d) = %d (bestscore)\n",
			    rL,rR,scoreL,scoreR,pen,end_reward,score));

	      bestscore = score;
	      *bestcL = cL;
	      *bestcR = cR;
	      *bestrL = rL;
	      *bestrR = rR;

	    } else {
	      debug3(printf("At %d left to %d right, score is (%d)+(%d)+(%d) = %d\n",
			    rL,rR,scoreL,scoreR,pen,scoreL+scoreR+pen));
	    }
	  }

	  for (/*at main diagonal*/; rR <= rhighR && rR < rightoffset-leftoffset-rL; rR++) {
	    scoreR = (int) matrixR_lower[rR][cR];

	    if ((score = scoreL + scoreR + pen + end_reward) > bestscore) {  /* Use > for jump early */
	      debug3(printf("At %d left to %d right, score is (%d)+(%d)+(%d)+(%d) = %d (bestscore)\n",
			    rL,rR,scoreL,scoreR,pen,end_reward,score));

	      bestscore = score;
	      *bestcL = cL;
	      *bestcR = cR;
	      *bestrL = rL;
	      *bestrR = rR;

	    } else {
	      debug3(printf("At %d left to %d right, score is (%d)+(%d)+(%d) = %d\n",
			    rL,rR,scoreL,scoreR,pen,scoreL+scoreR+pen));
	    }
	  }
	}

	pen = open - extend;	/* Subtract extend to compensate for
                                   its addition in the for loop */
      }
    }
  }
      
  *finalscore = (int) bestscore;
  debug3(printf("Returning final score of %d at (%d,%d) left to (%d,%d) right\n",
		*finalscore,*bestcL,*bestrL,*bestcR,*bestrR));

  return;
}
#endif

#if defined(HAVE_SSE4_1) || defined(HAVE_SSE2)
/* Columns are always genomic.  Rows are always query.  Bridging across common columns */
static void
bridge_cdna_gap_16_ud (int *finalscore, int *bestcL, int *bestcR, int *bestrL, int *bestrR,
		      Score16_T **matrixL_upper, Score16_T **matrixL_lower,
		      Score16_T **matrixR_upper, Score16_T **matrixR_lower,
		      int glength, int rlengthL, int rlengthR,
		      int lbandL, int ubandL, int lbandR, int ubandR,
		      int open, int extend, int leftoffset, int rightoffset, bool jump_late_p) {
  int bestscore = NEG_INFINITY_16, score, scoreL, scoreR, pen, end_reward = 0;
  int rL, rR, cL, cR;
  int rloL, rhighL;
  int rloR, rhighR;

#if 0
  /* Perform computations */
  lbandL = rlengthL - glength + extraband_paired;
  ubandL = extraband_paired;

  lbandR = rlengthR - glength + extraband_paired;
  ubandR = extraband_paired;
#endif

  /* Need a double loop on rows here, in contrast with a single loop
     for introns, because we allow a genomic insertion that doesn't
     match the cDNA.  So we need to add a penalty for a genomic
     insertion */

  if (jump_late_p) {
    for (cL = 1; cL < glength; cL++) {

      /* Note: opening penalty is added at the bottom of the loop */
      for (cR = glength-cL, pen = 0; cR >= 0; cR--, pen += extend) {
	/* debug3(printf("\nAt row %d on left and %d on right\n",cL,cR)); */
	if ((rloL = cL - ubandL) < 1) {
	  rloL = 1;
	}
	if ((rhighL = cL + lbandL) > rlengthL-1) {
	  rhighL = rlengthL-1;
	}

	if ((rloR = cR - ubandR) < 1) {
	  rloR = 1;
	}
	if ((rhighR = cR + lbandR) > rlengthR-1) {
	  rhighR = rlengthR-1;
	}

	for (rL = rloL; rL < /*to main diagonal*/cL; rL++) {
	  scoreL = (int) matrixL_upper[cL][rL];
	
	  /* Disallow leftoffset + rL >= rightoffset - rR, or rR >= rightoffset - leftoffset - rL */
	  debug3(printf("  Disallowing rR to be >= %d\n",rightoffset-leftoffset-rL));
	  for (rR = rloR; rR < /*to main diagonal*/cR && rR < rightoffset-leftoffset-rL; rR++) {
	    scoreR = (int) matrixR_upper[cR][rR];

	    if ((score = scoreL + scoreR + pen + end_reward) >= bestscore) {  /* Use >= for jump late */
	      debug3(printf("At %d left to %d right, score is (%d)+(%d)+(%d)+(%d) = %d (bestscore)\n",
			    rL,rR,scoreL,scoreR,pen,end_reward,score));

	      bestscore = score;
	      *bestcL = cL;
	      *bestcR = cR;
	      *bestrL = rL;
	      *bestrR = rR;

	    } else {
	      debug3(printf("At %d left to %d right, score is (%d)+(%d)+(%d) = %d\n",
			    rL,rR,scoreL,scoreR,pen,scoreL+scoreR+pen));
	    }
	  }

	  for (/*at main diagonal*/; rR <= rhighR && rR < rightoffset-leftoffset-rL; rR++) {
	    scoreR = (int) matrixR_lower[rR][cR];

	    if ((score = scoreL + scoreR + pen + end_reward) >= bestscore) {  /* Use >= for jump late */
	      debug3(printf("At %d left to %d right, score is (%d)+(%d)+(%d)+(%d) = %d (bestscore)\n",
			    rL,rR,scoreL,scoreR,pen,end_reward,score));

	      bestscore = score;
	      *bestcL = cL;
	      *bestcR = cR;
	      *bestrL = rL;
	      *bestrR = rR;

	    } else {
	      debug3(printf("At %d left to %d right, score is (%d)+(%d)+(%d) = %d\n",
			    rL,rR,scoreL,scoreR,pen,scoreL+scoreR+pen));
	    }
	  }
	}

	for (/*at main diagonal*/; rL <= rhighL; rL++) {
	  scoreL = (int) matrixL_lower[rL][cL];
	
	  /* Disallow leftoffset + rL >= rightoffset - rR, or rR >= rightoffset - leftoffset - rL */
	  debug3(printf("  Disallowing rR to be >= %d\n",rightoffset-leftoffset-rL));
	  for (rR = rloR; rR < /*to main diagonal*/cR && rR < rightoffset-leftoffset-rL; rR++) {
	    scoreR = (int) matrixR_upper[cR][rR];

	    if ((score = scoreL + scoreR + pen + end_reward) >= bestscore) {  /* Use >= for jump late */
	      debug3(printf("At %d left to %d right, score is (%d)+(%d)+(%d)+(%d) = %d (bestscore)\n",
			    rL,rR,scoreL,scoreR,pen,end_reward,score));

	      bestscore = score;
	      *bestcL = cL;
	      *bestcR = cR;
	      *bestrL = rL;
	      *bestrR = rR;

	    } else {
	      debug3(printf("At %d left to %d right, score is (%d)+(%d)+(%d) = %d\n",
			    rL,rR,scoreL,scoreR,pen,scoreL+scoreR+pen));
	    }
	  }

	  for (/*at main diagonal*/; rR <= rhighR && rR < rightoffset-leftoffset-rL; rR++) {
	    scoreR = (int) matrixR_lower[rR][cR];

	    if ((score = scoreL + scoreR + pen + end_reward) >= bestscore) {  /* Use >= for jump late */
	      debug3(printf("At %d left to %d right, score is (%d)+(%d)+(%d)+(%d) = %d (bestscore)\n",
			    rL,rR,scoreL,scoreR,pen,end_reward,score));

	      bestscore = score;
	      *bestcL = cL;
	      *bestcR = cR;
	      *bestrL = rL;
	      *bestrR = rR;

	    } else {
	      debug3(printf("At %d left to %d right, score is (%d)+(%d)+(%d) = %d\n",
			    rL,rR,scoreL,scoreR,pen,scoreL+scoreR+pen));
	    }
	  }
	}

	pen = open - extend;	/* Subtract extend to compensate for
                                   its addition in the for loop */
      }
    }

  } else {
    /* Do not jump late */
    for (cL = 1; cL < glength; cL++) {

      /* Note: opening penalty is added at the bottom of the loop */
      for (cR = glength-cL, pen = 0; cR >= 0; cR--, pen += extend) {
	/* debug3(printf("\nAt row %d on left and %d on right\n",cL,cR)); */
	if ((rloL = cL - ubandL) < 1) {
	  rloL = 1;
	}
	if ((rhighL = cL + lbandL) > rlengthL-1) {
	  rhighL = rlengthL-1;
	}

	if ((rloR = cR - ubandR) < 1) {
	  rloR = 1;
	}
	if ((rhighR = cR + lbandR) > rlengthR-1) {
	  rhighR = rlengthR-1;
	}

	for (rL = rloL; rL < /*to main diagonal*/cL; rL++) {
	  scoreL = (int) matrixL_upper[cL][rL];
	
	  /* Disallow leftoffset + rL >= rightoffset - rR, or rR >= rightoffset - leftoffset - rL */
	  debug3(printf("  Disallowing rR to be >= %d\n",rightoffset-leftoffset-rL));
	  for (rR = rloR; rR < /*to main diagonal*/cR && rR < rightoffset-leftoffset-rL; rR++) {
	    scoreR = (int) matrixR_upper[cR][rR];

	    if ((score = scoreL + scoreR + pen + end_reward) > bestscore) {  /* Use > for jump early */
	      debug3(printf("At %d left to %d right, score is (%d)+(%d)+(%d)+(%d) = %d (bestscore)\n",
			    rL,rR,scoreL,scoreR,pen,end_reward,score));

	      bestscore = score;
	      *bestcL = cL;
	      *bestcR = cR;
	      *bestrL = rL;
	      *bestrR = rR;

	    } else {
	      debug3(printf("At %d left to %d right, score is (%d)+(%d)+(%d) = %d\n",
			    rL,rR,scoreL,scoreR,pen,scoreL+scoreR+pen));
	    }
	  }

	  for (/*at main diagonal*/; rR <= rhighR && rR < rightoffset-leftoffset-rL; rR++) {
	    scoreR = (int) matrixR_lower[rR][cR];

	    if ((score = scoreL + scoreR + pen + end_reward) > bestscore) {  /* Use > for jump early */
	      debug3(printf("At %d left to %d right, score is (%d)+(%d)+(%d)+(%d) = %d (bestscore)\n",
			    rL,rR,scoreL,scoreR,pen,end_reward,score));

	      bestscore = score;
	      *bestcL = cL;
	      *bestcR = cR;
	      *bestrL = rL;
	      *bestrR = rR;

	    } else {
	      debug3(printf("At %d left to %d right, score is (%d)+(%d)+(%d) = %d\n",
			    rL,rR,scoreL,scoreR,pen,scoreL+scoreR+pen));
	    }
	  }
	}

	for (/*at main diagonal*/; rL <= rhighL; rL++) {
	  scoreL = (int) matrixL_lower[rL][cL];
	
	  /* Disallow leftoffset + rL >= rightoffset - rR, or rR >= rightoffset - leftoffset - rL */
	  debug3(printf("  Disallowing rR to be >= %d\n",rightoffset-leftoffset-rL));
	  for (rR = rloR; rR < /*to main diagonal*/cR && rR < rightoffset-leftoffset-rL; rR++) {
	    scoreR = (int) matrixR_upper[cR][rR];

	    if ((score = scoreL + scoreR + pen + end_reward) > bestscore) {  /* Use > for jump early */
	      debug3(printf("At %d left to %d right, score is (%d)+(%d)+(%d)+(%d) = %d (bestscore)\n",
			    rL,rR,scoreL,scoreR,pen,end_reward,score));

	      bestscore = score;
	      *bestcL = cL;
	      *bestcR = cR;
	      *bestrL = rL;
	      *bestrR = rR;

	    } else {
	      debug3(printf("At %d left to %d right, score is (%d)+(%d)+(%d) = %d\n",
			    rL,rR,scoreL,scoreR,pen,scoreL+scoreR+pen));
	    }
	  }

	  for (/*at main diagonal*/; rR <= rhighR && rR < rightoffset-leftoffset-rL; rR++) {
	    scoreR = (int) matrixR_lower[rR][cR];

	    if ((score = scoreL + scoreR + pen + end_reward) > bestscore) {  /* Use > for jump early */
	      debug3(printf("At %d left to %d right, score is (%d)+(%d)+(%d)+(%d) = %d (bestscore)\n",
			    rL,rR,scoreL,scoreR,pen,end_reward,score));

	      bestscore = score;
	      *bestcL = cL;
	      *bestcR = cR;
	      *bestrL = rL;
	      *bestrR = rR;

	    } else {
	      debug3(printf("At %d left to %d right, score is (%d)+(%d)+(%d) = %d\n",
			    rL,rR,scoreL,scoreR,pen,scoreL+scoreR+pen));
	    }
	  }
	}

	pen = open - extend;	/* Subtract extend to compensate for
                                   its addition in the for loop */
      }
    }
  }
      
  *finalscore = (int) bestscore;
  debug3(printf("Returning final score of %d at (%d,%d) left to (%d,%d) right\n",
		*finalscore,*bestcL,*bestrL,*bestcR,*bestrR));

  return;
}
#endif



#ifndef HAVE_SSE2
/* Columns are always genomic.  Rows are always query.  Bridging across common columns */
static void
bridge_cdna_gap (int *finalscore, int *bestcL, int *bestcR, int *bestrL, int *bestrR,
#ifdef HAVE_SSE2
		 Score16_T **matrixL, Score16_T **matrixR,
#else
		 Score32_T **matrixL, Score32_T **matrixR,
#endif
		 int glength, int rlengthL, int rlengthR, int extraband_paired,
		 int open, int extend, int leftoffset, int rightoffset, bool jump_late_p) {
  int bestscore = NEG_INFINITY_32, score, scoreL, scoreR, pen, end_reward = 0;
  int rL, rR, cL, cR;
  int lbandL, ubandL, rloL, rhighL;
  int lbandR, ubandR, rloR, rhighR;

  /* Perform computations */
  lbandL = rlengthL - glength + extraband_paired;
  ubandL = extraband_paired;

  lbandR = rlengthR - glength + extraband_paired;
  ubandR = extraband_paired;

  /* Need a double loop on rows here, in contrast with a single loop
     for introns, because we allow a genomic insertion that doesn't
     match the cDNA.  So we need to add a penalty for a genomic
     insertion */

  if (jump_late_p) {
    for (cL = 1; cL < glength; cL++) {

      /* Note: opening penalty is added at the bottom of the loop */
      for (cR = glength-cL, pen = 0; cR >= 0; cR--, pen += extend) {
	/* debug3(printf("\nAt row %d on left and %d on right\n",cL,cR)); */
	if ((rloL = cL - ubandL) < 1) {
	  rloL = 1;
	}
	if ((rhighL = cL + lbandL) > rlengthL-1) {
	  rhighL = rlengthL-1;
	}

	if ((rloR = cR - ubandR) < 1) {
	  rloR = 1;
	}
	if ((rhighR = cR + lbandR) > rlengthR-1) {
	  rhighR = rlengthR-1;
	}

	for (rL = rloL; rL <= rhighL; rL++) {
	  scoreL = (int) matrixL[cL][rL];
	
	  /* Disallow leftoffset + rL >= rightoffset - rR, or rR >= rightoffset - leftoffset - rL */
	  debug3(printf("  Disallowing rR to be >= %d\n",rightoffset-leftoffset-rL));
	  for (rR = rloR; rR <= rhighR && rR < rightoffset-leftoffset-rL; rR++) {
	    scoreR = (int) matrixR[cR][rR];

	    if ((score = scoreL + scoreR + pen + end_reward) >= bestscore) {  /* Use >= for jump late */
	      debug3(printf("At %d left to %d right, score is (%d)+(%d)+(%d)+(%d) = %d (bestscore)\n",
			    rL,rR,scoreL,scoreR,pen,end_reward,score));

	      bestscore = score;
	      *bestcL = cL;
	      *bestcR = cR;
	      *bestrL = rL;
	      *bestrR = rR;

	    } else {
	      debug3(printf("At %d left to %d right, score is (%d)+(%d)+(%d) = %d\n",
			    rL,rR,scoreL,scoreR,pen,scoreL+scoreR+pen));
	    }
	  }
	}
	pen = open - extend;	/* Subtract extend to compensate for
                                   its addition in the for loop */
      }
    }

  } else {
    /* Do not jump late */
    for (cL = 1; cL < glength; cL++) {

      /* Note: opening penalty is added at the bottom of the loop */
      for (cR = glength-cL, pen = 0; cR >= 0; cR--, pen += extend) {
	/* debug3(printf("\nAt row %d on left and %d on right\n",cL,cR)); */
	if ((rloL = cL - ubandL) < 1) {
	  rloL = 1;
	}
	if ((rhighL = cL + lbandL) > rlengthL-1) {
	  rhighL = rlengthL-1;
	}

	if ((rloR = cR - ubandR) < 1) {
	  rloR = 1;
	}
	if ((rhighR = cR + lbandR) > rlengthR-1) {
	  rhighR = rlengthR-1;
	}

	for (rL = rloL; rL <= rhighL; rL++) {
	  scoreL = (int) matrixL[cL][rL];
	
	  /* Disallow leftoffset + rL >= rightoffset - rR, or rR >= rightoffset - leftoffset - rL */
	  debug3(printf("  Disallowing rR to be >= %d\n",rightoffset-leftoffset-rL));
	  for (rR = rloR; rR <= rhighR && rR < rightoffset-leftoffset-rL; rR++) {
	    scoreR = (int) matrixR[cR][rR];

	    if ((score = scoreL + scoreR + pen + end_reward) > bestscore) {  /* Use > for jump early */
	      debug3(printf("At %d left to %d right, score is (%d)+(%d)+(%d)+(%d) = %d (bestscore)\n",
			    rL,rR,scoreL,scoreR,pen,end_reward,score));

	      bestscore = score;
	      *bestcL = cL;
	      *bestcR = cR;
	      *bestrL = rL;
	      *bestrR = rR;

	    } else {
	      debug3(printf("At %d left to %d right, score is (%d)+(%d)+(%d) = %d\n",
			    rL,rR,scoreL,scoreR,pen,scoreL+scoreR+pen));
	    }
	  }
	}
	pen = open - extend;	/* Subtract extend to compensate for
                                   its addition in the for loop */
      }
    }
  }
      
  *finalscore = (int) bestscore;
  debug3(printf("Returning final score of %d at (%d,%d) left to (%d,%d) right\n",
		*finalscore,*bestcL,*bestrL,*bestcR,*bestrR));

  return;
}
#endif


/* Sequences rsequenceL and rsequenceR represent the two ends of the cDNA insertion */
List_T
Dynprog_cdna_gap (int *dynprogindex, int *finalscore, bool *incompletep,
		  T dynprogL, T dynprogR, char *rsequenceL, char *rsequence_ucL, 
		  char *rev_rsequenceR, char *rev_rsequence_ucR,
#if 0
		  char *gsequence, char *gsequence_uc,
#endif
		  int rlengthL, int rlengthR, int glength,
		  int roffsetL, int rev_roffsetR, int goffset,
		  Univcoord_T chroffset, Univcoord_T chrhigh,
		  int cdna_direction, bool watsonp, bool jump_late_p, Pairpool_T pairpool,
		  int extraband_paired, double defect_rate) {
  List_T pairs = NULL;
  char *gsequence, *gsequence_alt, *rev_gsequence, *rev_gsequence_alt;
  Mismatchtype_T mismatchtype;
  int lbandL, ubandL, lbandR, ubandR;
  int open, extend;
#if defined(HAVE_SSE4_1) || defined(HAVE_SSE2)
  Score8_T **matrix8L_upper, **matrix8L_lower, **matrix8R_upper, **matrix8R_lower;
  Direction8_T **directions8L_upper_nogap, **directions8L_upper_Egap,
    **directions8L_lower_nogap, **directions8L_lower_Egap,
    **directions8R_upper_nogap, **directions8R_upper_Egap,
    **directions8R_lower_nogap, **directions8R_lower_Egap;
  bool use8p;

  Score16_T **matrix16L_upper, **matrix16L_lower, **matrix16R_upper, **matrix16R_lower;
  Direction16_T **directions16L_upper_nogap, **directions16L_upper_Egap,
    **directions16L_lower_nogap, **directions16L_lower_Egap,
    **directions16R_upper_nogap, **directions16R_upper_Egap,
    **directions16R_lower_nogap, **directions16R_lower_Egap;
#else
  Score32_T **matrixL, **matrixR;
  Direction32_T **directionsL_nogap, **directionsL_Egap, **directionsL_Fgap,
    **directionsR_nogap, **directionsR_Egap, **directionsR_Fgap;
#endif
  int rev_goffset, bestrL, bestrR, bestcL, bestcR, k;
  int nmatches, nmismatches, nopens, nindels;
  int queryjump, genomejump;
  char c2, c2_alt;


  if (glength <= 1) {
    return NULL;
  }

  debug(printf("\n"));
  debug(printf("%c:  ",*dynprogindex > 0 ? (*dynprogindex-1)%26+'a' : (-(*dynprogindex)-1)%26+'A'));
  debug(printf("Aligning cdna gap\n"));
#ifdef EXTRACT_GENOMICSEG
  debug(printf("At genomic offset %d-%d, %.*s\n",goffset,goffset+glength-1,glength,gsequence));
#endif
  debug(printf("\n"));

  /* ?check if offsets are too close.  But this eliminates a segment
     of the cDNA.  Should check in stage 3, and do single gap instead. */
  /*
  if (roffsetL+rlengthL-1 >= rev_roffsetR-rlengthR+1) {
    debug(printf("Bounds don't make sense\n"));
    *finalscore = NEG_INFINITY_16;
    return NULL;
  }
  */

  if (defect_rate < DEFECT_HIGHQ) {
    mismatchtype = HIGHQ;
    /* mismatch = MISMATCH_HIGHQ; */
    open = CDNA_OPEN_HIGHQ;
    extend = CDNA_EXTEND_HIGHQ;
  } else if (defect_rate < DEFECT_MEDQ) {
    mismatchtype = MEDQ;
    /* mismatch = MISMATCH_MEDQ; */
    open = CDNA_OPEN_MEDQ;
    extend = CDNA_EXTEND_MEDQ;
  } else {
    mismatchtype = LOWQ;
    /* mismatch = MISMATCH_LOWQ; */
    open = CDNA_OPEN_LOWQ;
    extend = CDNA_EXTEND_LOWQ;
  }

  if (glength > dynprogR->max_glength || rlengthR > dynprogR->max_rlength) {
    debug(printf("glength %d or rlengthR %d is too long.  Returning NULL\n",glength,rlengthR));
#if 0
    rev_goffset = goffset + glength - 1;
    queryjump = rev_roffsetR - roffsetL + 1;
    genomejump = rev_goffset - goffset + 1;
    pairs = Pairpool_push_gapholder(NULL,pairpool,queryjump,genomejump,
				    /*leftpair*/NULL,/*rightpair*/NULL,/*knownp*/false);
#endif
    *dynprogindex += (*dynprogindex > 0 ? +1 : -1);
    return (List_T) NULL;
  }

  if (glength > dynprogL->max_glength || rlengthL > dynprogL->max_rlength) {
    debug(printf("glength %d or rlengthL %d is too long.  Returning NULL\n",glength,rlengthL));
#if 0
    rev_goffset = goffset + glength - 1;
    queryjump = rev_roffsetR - roffsetL + 1;
    genomejump = rev_goffset - goffset + 1;
    pairs = Pairpool_push_gapholder(NULL,pairpool,queryjump,genomejump,
				    /*leftpair*/NULL,/*rightpair*/NULL,/*knownp*/false);
#endif
    *dynprogindex += (*dynprogindex > 0 ? +1 : -1);
    return (List_T) NULL;
  }

#if 0
  /* Right side looks like 5' end */
  /* Note: sequence 1 and 2 flipped, because 1 has extramaterial */
  rev_gsequence = &(gsequence[glength-1]);
  rev_gsequence_uc = &(gsequence_uc[glength-1]);
#endif
  rev_goffset = goffset+glength-1;

  debug(printf("At query offset %d-%d, %.*s\n",roffsetL,roffsetL+rlengthL-1,rlengthL,rsequenceL));
  debug(printf("At query offset %d-%d, %.*s\n",rev_roffsetR-rlengthR+1,rev_roffsetR,rlengthR,&(rev_rsequenceR[-rlengthR+1])));
  debug(printf("Whole piece at query offset %d-%d, %.*s\n",roffsetL,rev_roffsetR,rev_roffsetR-roffsetL+1,rsequenceL));

#if defined(HAVE_SSE4_1) || defined(HAVE_SSE2)
  /* Use || because we want the minimum length (which determines the diagonal length) to achieve a score less than 128 */
  if (glength <= SIMD_MAXLENGTH_EPI8 || (rlengthL <= SIMD_MAXLENGTH_EPI8 && rlengthR <= SIMD_MAXLENGTH_EPI8)) {
    use8p = true;
  } else {
    use8p = false;
  }
#endif


  rev_gsequence = (char *) MALLOCA((glength+1) * sizeof(char));
  rev_gsequence_alt = (char *) MALLOCA((glength+1) * sizeof(char));
  gsequence = (char *) MALLOCA((glength+1) * sizeof(char));
  gsequence_alt = (char *) MALLOCA((glength+1) * sizeof(char));

  if (watsonp) {
    Genome_get_segment_blocks_left(rev_gsequence,rev_gsequence_alt,/*right*/chroffset+rev_goffset+1,
				   glength,chroffset,/*revcomp*/false);
    Genome_get_segment_blocks_right(gsequence,gsequence_alt,/*left*/chroffset+goffset,
				    glength,chrhigh,/*revcomp*/false);
  } else {
    Genome_get_segment_blocks_right(rev_gsequence,rev_gsequence_alt,/*left*/chrhigh-rev_goffset,
				    glength,chrhigh,/*revcomp*/true);
    Genome_get_segment_blocks_left(gsequence,gsequence_alt,/*right*/chrhigh-goffset+1,
				   glength,chroffset,/*revcomp*/true);
  }
  if (gsequence[0] == '\0' || rev_gsequence[0] == '\0') {
    FREEA(gsequence_alt);
    FREEA(gsequence);
    FREEA(rev_gsequence_alt);
    FREEA(rev_gsequence);
    return (List_T) NULL;
  }


#if defined(HAVE_SSE4_1) || defined(HAVE_SSE2)
  if (use8p == true) {
    Dynprog_compute_bands(&lbandL,&ubandL,rlengthL,glength,extraband_paired,/*widebandp*/true);
    matrix8L_upper = Dynprog_simd_8_upper(&directions8L_upper_nogap,&directions8L_upper_Egap,dynprogL,
					  rsequenceL,gsequence,gsequence_alt,rlengthL,glength,
#ifdef DEBUG14
					  goffset,chroffset,chrhigh,watsonp,
#endif
					  mismatchtype,open,extend,ubandL,jump_late_p,/*revp*/false);

    matrix8L_lower = Dynprog_simd_8_lower(&directions8L_lower_nogap,&directions8L_lower_Egap,dynprogL,
					  rsequenceL,gsequence,gsequence_alt,rlengthL,glength,
#ifdef DEBUG14
					  goffset,chroffset,chrhigh,watsonp,
#endif
					  mismatchtype,open,extend,lbandL,jump_late_p,/*revp*/false);
    

    Dynprog_compute_bands(&lbandR,&ubandR,rlengthR,glength,extraband_paired,/*widebandp*/true);
    matrix8R_upper = Dynprog_simd_8_upper(&directions8R_upper_nogap,&directions8R_upper_Egap,dynprogR,
					  rev_rsequenceR,&(rev_gsequence[glength-1]),&(rev_gsequence_alt[glength-1]),
					  rlengthR,glength,
#ifdef DEBUG14
					  rev_goffset,chroffset,chrhigh,watsonp,
#endif
					  mismatchtype,open,extend,ubandR,/*for revp true*/!jump_late_p,/*revp*/true);

    matrix8R_lower = Dynprog_simd_8_lower(&directions8R_lower_nogap,&directions8R_lower_Egap,dynprogR,
					  rev_rsequenceR,&(rev_gsequence[glength-1]),&(rev_gsequence_alt[glength-1]),
					  rlengthR,glength,
#ifdef DEBUG14
					  rev_goffset,chroffset,chrhigh,watsonp,
#endif
					  mismatchtype,open,extend,lbandR,/*for revp true*/!jump_late_p,/*revp*/true);


    nmatches = nmismatches = nopens = nindels = 0;
    bridge_cdna_gap_8_ud(&(*finalscore),&bestcL,&bestcR,&bestrL,&bestrR,
			 matrix8L_upper,matrix8L_lower,matrix8R_upper,matrix8R_lower,
			 glength,rlengthL,rlengthR,lbandL,ubandL,lbandR,ubandR,
			 open,extend,roffsetL,rev_roffsetR,jump_late_p);

    if (bestcR >= bestrR) {
      pairs = Dynprog_traceback_8_upper(NULL,&nmatches,&nmismatches,&nopens,&nindels,
					directions8R_upper_nogap,directions8R_upper_Egap,
					bestrR,bestcR,rev_rsequenceR,rev_rsequence_ucR,
					&(rev_gsequence[glength-1]),&(rev_gsequence_alt[glength-1]),
					rev_roffsetR,rev_goffset,pairpool,/*revp*/true,
					chroffset,chrhigh,cdna_direction,watsonp,*dynprogindex);
    } else {
      pairs = Dynprog_traceback_8_lower(NULL,&nmatches,&nmismatches,&nopens,&nindels,
					directions8R_lower_nogap,directions8R_lower_Egap,
					bestrR,bestcR,rev_rsequenceR,rev_rsequence_ucR,
					&(rev_gsequence[glength-1]),&(rev_gsequence_alt[glength-1]),
					rev_roffsetR,rev_goffset,pairpool,/*revp*/true,
					chroffset,chrhigh,cdna_direction,watsonp,*dynprogindex);
    }

    pairs = List_reverse(pairs);

    queryjump = (rev_roffsetR-bestrR) - (roffsetL+bestrL) + 1;
    genomejump = (rev_goffset-bestcR) - (goffset+bestcL) + 1;
    /* No need to revise queryjump or genomejump, because the above
       coordinates are internal to the gap. */

    if (queryjump == INSERT_PAIRS && genomejump == INSERT_PAIRS) {
      /* Add cDNA insertion, if any */
      for (k = rev_roffsetR-bestrR; k >= roffsetL+bestrL; k--) {
	debug(printf("cDNA insertion, Pushing [%d,%d] (%c,-)\n",k,rev_goffset-bestcR+1,rsequenceL[k-roffsetL]));
	pairs = Pairpool_push(pairs,pairpool,k,rev_goffset-bestcR+1,rsequenceL[k-roffsetL],SHORTGAP_COMP,
			      /*genome*/' ',/*genomealt*/' ',*dynprogindex);
      }
      debug(printf("\n"));


      /* This loop not yet checked for get_genomic_nt giving correct answer */
      for (k = rev_goffset-bestcR; k >= goffset+bestcL; k--) {
	c2 = get_genomic_nt(&c2_alt,k,chroffset,chrhigh,watsonp);
	debug(printf("genome insertion, Pushing [%d,%d] (-,%c)\n",roffsetL+bestrL,k,c2));
#if 0
	assert(c2 == gsequence[k-goffset]);
	pairs = Pairpool_push(pairs,pairpool,roffsetL+bestrL,k,' ',SHORTGAP_COMP,
			      gsequence[k-goffset],/*genomealt*/GENOMEALT_DEFERRED,*dynprogindex);
#else
	pairs = Pairpool_push(pairs,pairpool,roffsetL+bestrL,k,' ',SHORTGAP_COMP,c2,c2_alt,*dynprogindex);
#endif
      }
      debug(printf("\n"));

    } else {

      /* Add gapholder to be solved in the future */
#ifndef NOGAPHOLDER
      debug(printf("Pushing a gap with queryjump = %d, genomejump = %d\n",queryjump,genomejump));
      pairs = Pairpool_push_gapholder(pairs,pairpool,queryjump,genomejump,
				      /*leftpair*/NULL,/*rightpair*/NULL,/*knownp*/false);
#endif
      *incompletep = true;
    }

    if (bestcL >= bestrL) {
      pairs = Dynprog_traceback_8_upper(pairs,&nmatches,&nmismatches,&nopens,&nindels,
					directions8L_upper_nogap,directions8L_upper_Egap,
					bestrL,bestcL,rsequenceL,rsequence_ucL,
					gsequence,gsequence_alt,roffsetL,goffset,pairpool,/*revp*/false,
					chroffset,chrhigh,cdna_direction,watsonp,*dynprogindex);
    } else {
      pairs = Dynprog_traceback_8_lower(pairs,&nmatches,&nmismatches,&nopens,&nindels,
					directions8L_lower_nogap,directions8L_lower_Egap,
					bestrL,bestcL,rsequenceL,rsequence_ucL,
					gsequence,gsequence_alt,roffsetL,goffset,pairpool,/*revp*/false,
					chroffset,chrhigh,cdna_direction,watsonp,*dynprogindex);
    }

    if (List_length(pairs) == 1) {
      /* Only a gap added */
      pairs = (List_T) NULL;
    }

    FREEA(gsequence_alt);
    FREEA(gsequence);
    FREEA(rev_gsequence_alt);
    FREEA(rev_gsequence);

    debug(printf("End of dynprog cDNA gap\n"));

    *dynprogindex += (*dynprogindex > 0 ? +1 : -1);
    return List_reverse(pairs);

  } else {
    /* Use 16-mers */
    Dynprog_compute_bands(&lbandL,&ubandL,rlengthL,glength,extraband_paired,/*widebandp*/true);
    matrix16L_upper = Dynprog_simd_16_upper(&directions16L_upper_nogap,&directions16L_upper_Egap,dynprogL,
					    rsequenceL,gsequence,gsequence_alt,rlengthL,glength,
#ifdef DEBUG14
					    goffset,chroffset,chrhigh,watsonp,
#endif
					    mismatchtype,open,extend,ubandL,jump_late_p,/*revp*/false);

    matrix16L_lower = Dynprog_simd_16_lower(&directions16L_lower_nogap,&directions16L_lower_Egap,dynprogL,
					    rsequenceL,gsequence,gsequence_alt,rlengthL,glength,
#ifdef DEBUG14
					    goffset,chroffset,chrhigh,watsonp,
#endif
					    mismatchtype,open,extend,lbandL,jump_late_p,/*revp*/false);

    Dynprog_compute_bands(&lbandR,&ubandR,rlengthR,glength,extraband_paired,/*widebandp*/true);
    matrix16R_upper = Dynprog_simd_16_upper(&directions16R_upper_nogap,&directions16R_upper_Egap,dynprogR,
					    rev_rsequenceR,&(rev_gsequence[glength-1]),&(rev_gsequence_alt[glength-1]),
					    rlengthR,glength,
#ifdef DEBUG14
					    rev_goffset,chroffset,chrhigh,watsonp,
#endif
					    mismatchtype,open,extend,ubandR,/*for revp true*/!jump_late_p,/*revp*/true);

    matrix16R_lower = Dynprog_simd_16_lower(&directions16R_lower_nogap,&directions16R_lower_Egap,dynprogR,
					    rev_rsequenceR,&(rev_gsequence[glength-1]),&(rev_gsequence_alt[glength-1]),
					    rlengthR,glength,
#ifdef DEBUG14
					    rev_goffset,chroffset,chrhigh,watsonp,
#endif
					    mismatchtype,open,extend,lbandR,/*for revp true*/!jump_late_p,/*revp*/true);

    nmatches = nmismatches = nopens = nindels = 0;

    bridge_cdna_gap_16_ud(&(*finalscore),&bestcL,&bestcR,&bestrL,&bestrR,
			 matrix16L_upper,matrix16L_lower,matrix16R_upper,matrix16R_lower,
			 glength,rlengthL,rlengthR,lbandL,ubandL,lbandR,ubandR,
			 open,extend,roffsetL,rev_roffsetR,jump_late_p);

    if (bestcR >= bestrR) {
      pairs = Dynprog_traceback_16_upper(NULL,&nmatches,&nmismatches,&nopens,&nindels,
					 directions16R_upper_nogap,directions16R_upper_Egap,
					 bestrR,bestcR,rev_rsequenceR,rev_rsequence_ucR,
					 &(rev_gsequence[glength-1]),&(rev_gsequence_alt[glength-1]),
					 rev_roffsetR,rev_goffset,pairpool,/*revp*/true,
					 chroffset,chrhigh,cdna_direction,watsonp,*dynprogindex);
    } else {
      pairs = Dynprog_traceback_16_lower(NULL,&nmatches,&nmismatches,&nopens,&nindels,
					 directions16R_lower_nogap,directions16R_lower_Egap,
					 bestrR,bestcR,rev_rsequenceR,rev_rsequence_ucR,
					 &(rev_gsequence[glength-1]),&(rev_gsequence_alt[glength-1]),
					 rev_roffsetR,rev_goffset,pairpool,/*revp*/true,
					 chroffset,chrhigh,cdna_direction,watsonp,*dynprogindex);
    }

    pairs = List_reverse(pairs);

    queryjump = (rev_roffsetR-bestrR) - (roffsetL+bestrL) + 1;
    genomejump = (rev_goffset-bestcR) - (goffset+bestcL) + 1;
    /* No need to revise queryjump or genomejump, because the above
       coordinates are internal to the gap. */

    if (queryjump == INSERT_PAIRS && genomejump == INSERT_PAIRS) {
      /* Add cDNA insertion, if any */
      for (k = rev_roffsetR-bestrR; k >= roffsetL+bestrL; k--) {
	debug(printf("cDNA insertion, Pushing [%d,%d] (%c,-)\n",k,rev_goffset-bestcR+1,rsequenceL[k-roffsetL]));
	pairs = Pairpool_push(pairs,pairpool,k,rev_goffset-bestcR+1,rsequenceL[k-roffsetL],SHORTGAP_COMP,
			      /*genome*/' ',/*genomealt*/' ',*dynprogindex);
      }
      debug(printf("\n"));


      /* This loop not yet checked for get_genomic_nt giving correct answer */
      for (k = rev_goffset-bestcR; k >= goffset+bestcL; k--) {
	c2 = get_genomic_nt(&c2_alt,k,chroffset,chrhigh,watsonp);
	debug(printf("genome insertion, Pushing [%d,%d] (-,%c)\n",roffsetL+bestrL,k,c2));
#if 0
	assert(c2 == gsequence[k-goffset]);
	pairs = Pairpool_push(pairs,pairpool,roffsetL+bestrL,k,' ',SHORTGAP_COMP,
			      gsequence[k-goffset],/*genomealt*/GENOMEALT_DEFERRED,*dynprogindex);
#else
	pairs = Pairpool_push(pairs,pairpool,roffsetL+bestrL,k,' ',SHORTGAP_COMP,c2,c2_alt,*dynprogindex);
#endif
      }
      debug(printf("\n"));

    } else {

      /* Add gapholder to be solved in the future */
#ifndef NOGAPHOLDER
      debug(printf("Pushing a gap with queryjump = %d, genomejump = %d\n",queryjump,genomejump));
      pairs = Pairpool_push_gapholder(pairs,pairpool,queryjump,genomejump,
				      /*leftpair*/NULL,/*rightpair*/NULL,/*knownp*/false);
#endif
      *incompletep = true;
    }

    if (bestcL >= bestrL) {
      pairs = Dynprog_traceback_16_upper(pairs,&nmatches,&nmismatches,&nopens,&nindels,
					 directions16L_upper_nogap,directions16L_upper_Egap,
					 bestrL,bestcL,rsequenceL,rsequence_ucL,
					 gsequence,gsequence_alt,roffsetL,goffset,pairpool,/*revp*/false,
					 chroffset,chrhigh,cdna_direction,watsonp,*dynprogindex);
    } else {
      pairs = Dynprog_traceback_16_lower(pairs,&nmatches,&nmismatches,&nopens,&nindels,
					 directions16L_lower_nogap,directions16L_lower_Egap,
					 bestrL,bestcL,rsequenceL,rsequence_ucL,
					 gsequence,gsequence_alt,roffsetL,goffset,pairpool,/*revp*/false,
					 chroffset,chrhigh,cdna_direction,watsonp,*dynprogindex);
    }

    if (List_length(pairs) == 1) {
      /* Only a gap added */
      pairs = (List_T) NULL;
    }

    FREEA(gsequence_alt);
    FREEA(gsequence);
    FREEA(rev_gsequence_alt);
    FREEA(rev_gsequence);

    debug(printf("End of dynprog cDNA gap\n"));

    *dynprogindex += (*dynprogindex > 0 ? +1 : -1);
    return List_reverse(pairs);
  }

#else
  /* Non-SIMD methods */
  Dynprog_compute_bands(&lbandL,&ubandL,rlengthL,glength,extraband_paired,/*widebandp*/true);
  matrixL = Dynprog_standard(&directionsL_nogap,&directionsL_Egap,&directionsL_Fgap,dynprogL,
			     rsequenceL,gsequence,gsequence_alt,rlengthL,glength,
			     goffset,chroffset,chrhigh,watsonp,
			     mismatchtype,open,extend,lbandL,ubandL,
			     jump_late_p,/*revp*/false,/*saturation*/NEG_INFINITY_INT);
  
  Dynprog_compute_bands(&lbandR,&ubandR,rlengthR,glength,extraband_paired,/*widebandp*/true);
  matrixR = Dynprog_standard(&directionsR_nogap,&directionsR_Egap,&directionsR_Fgap,dynprogR,
			     rev_rsequenceR,&(rev_gsequence[glength-1]),&(rev_gsequence_alt[glength-1]),
			     rlengthR,glength,rev_goffset,chroffset,chrhigh,watsonp,
			     mismatchtype,open,extend,lbandL,ubandR,
			     /*for revp true*/!jump_late_p,/*revp*/true,/*saturation*/NEG_INFINITY_INT);

  nmatches = nmismatches = nopens = nindels = 0;

  bridge_cdna_gap(&(*finalscore),&bestcL,&bestcR,&bestrL,&bestrR,matrixL,matrixR,
		  glength,rlengthL,rlengthR,extraband_paired,
		  open,extend,roffsetL,rev_roffsetR,jump_late_p);

  pairs = Dynprog_traceback_std(NULL,&nmatches,&nmismatches,&nopens,&nindels,
				directionsR_nogap,directionsR_Egap,directionsR_Fgap,bestrR,bestcR,
				rev_rsequenceR,rev_rsequence_ucR,
				&(rev_gsequence[glength-1]),&(rev_gsequence_alt[glength-1]),
				rev_roffsetR,rev_goffset,pairpool,/*revp*/true,
				chroffset,chrhigh,cdna_direction,watsonp,*dynprogindex);

  pairs = List_reverse(pairs);

  queryjump = (rev_roffsetR-bestrR) - (roffsetL+bestrL) + 1;
  genomejump = (rev_goffset-bestcR) - (goffset+bestcL) + 1;
  /* No need to revise queryjump or genomejump, because the above
     coordinates are internal to the gap. */

  if (queryjump == INSERT_PAIRS && genomejump == INSERT_PAIRS) {
    /* Add cDNA insertion, if any */
    for (k = rev_roffsetR-bestrR; k >= roffsetL+bestrL; k--) {
      debug(printf("cDNA insertion, Pushing [%d,%d] (%c,-)\n",k,rev_goffset-bestcR+1,rsequenceL[k-roffsetL]));
      pairs = Pairpool_push(pairs,pairpool,k,rev_goffset-bestcR+1,rsequenceL[k-roffsetL],SHORTGAP_COMP,
			    /*genome*/' ',/*genomealt*/' ',*dynprogindex);
    }
    debug(printf("\n"));


    /* This loop not yet checked for get_genomic_nt giving correct answer */
    for (k = rev_goffset-bestcR; k >= goffset+bestcL; k--) {
      c2 = get_genomic_nt(&c2_alt,k,chroffset,chrhigh,watsonp);
      debug(printf("genome insertion, Pushing [%d,%d] (-,%c)\n",roffsetL+bestrL,k,c2));
#if 0
      assert(c2 == gsequence[k-goffset]);
      pairs = Pairpool_push(pairs,pairpool,roffsetL+bestrL,k,' ',SHORTGAP_COMP,
			    gsequence[k-goffset],/*genomealt*/GENOMEALT_DEFERRED,*dynprogindex);
#else
      pairs = Pairpool_push(pairs,pairpool,roffsetL+bestrL,k,' ',SHORTGAP_COMP,c2,c2_alt,*dynprogindex);
#endif
    }
    debug(printf("\n"));

  } else {

    /* Add gapholder to be solved in the future */
#ifndef NOGAPHOLDER
    debug(printf("Pushing a gap with queryjump = %d, genomejump = %d\n",queryjump,genomejump));
    pairs = Pairpool_push_gapholder(pairs,pairpool,queryjump,genomejump,
				    /*leftpair*/NULL,/*rightpair*/NULL,/*knownp*/false);
#endif
    *incompletep = true;
  }

  pairs = Dynprog_traceback_std(pairs,&nmatches,&nmismatches,&nopens,&nindels,
				directionsL_nogap,directionsL_Egap,directionsL_Fgap,bestrL,bestcL,
				rsequenceL,rsequence_ucL,
				gsequence,gsequence_alt,roffsetL,goffset,pairpool,/*revp*/false,
				chroffset,chrhigh,cdna_direction,watsonp,*dynprogindex);

  if (List_length(pairs) == 1) {
    /* Only a gap added */
    pairs = (List_T) NULL;
  }

  FREEA(gsequence_alt);
  FREEA(gsequence);
  FREEA(rev_gsequence_alt);
  FREEA(rev_gsequence);

  debug(printf("End of dynprog cDNA gap\n"));

  *dynprogindex += (*dynprogindex > 0 ? +1 : -1);
  return List_reverse(pairs);

#endif

}

