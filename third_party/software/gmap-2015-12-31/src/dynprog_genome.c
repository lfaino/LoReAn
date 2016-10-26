static char rcsid[] = "$Id: dynprog_genome.c 182795 2016-01-21 19:51:16Z twu $";
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
#include "intron.h"
#include "complement.h"
#include "maxent.h"
#include "maxent_hr.h"
#include "dynprog.h"		/* For parameters */
#include "dynprog_simd.h"


#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif

/* Prints all winning bridge scores */
#ifdef DEBUG3
#define debug3(x) x
#else
#define debug3(x)
#endif

/* Prints all losing bridge scores */
#ifdef DEBUG3A
#define debug3a(x) x
#else
#define debug3a(x)
#endif

/* Known splicing */
#ifdef DEBUG5
#define debug5(x) x
#else
#define debug5(x)
#endif

/* Getting genomic nt */
#ifdef DEBUG8
#define debug8(x) x
#else
#define debug8(x)
#endif

/* Splice site probabilities */
#ifdef DEBUG9
#define debug9(x) x
#else
#define debug9(x)
#endif

/* print_vector */
#ifdef DEBUG15
#define debug15(x) x
#else
#define debug15(x)
#endif


#define USE_SCOREI 1
#define USE_WEAK_SCOREI 1

#define PROB_CEILING 0.85
#define PROB_FLOOR 0.75
#define PROB_BAD 0.50

/* Prefer alternate intron to other non-canonicals, but don't
   introduce mismatches or gaps to identify */
#ifdef USE_WEAK_SCOREI
#define CANONICAL_INTRON 6
#define GCAG_INTRON 4
#define ATAC_INTRON 2
#define FINAL_GCAG_INTRON 4    /* Amount above regular should approximately
				   match FINAL_CANONICAL_INTRON - CANONICAL_INTRON */
#define FINAL_ATAC_INTRON 2
#else
#define GCAG_INTRON 15
#define ATAC_INTRON 12
#define FINAL_GCAG_INTRON 20    /* Amount above regular should approximately
				   match FINAL_CANONICAL_INTRON - CANONICAL_INTRON */
#define FINAL_ATAC_INTRON 12
#endif


/* Don't want to make too high, otherwise we will harm evaluation of
   dual introns vs. single intron */
#define CANONICAL_INTRON_HIGHQ 10 /* GT-AG */
#define CANONICAL_INTRON_MEDQ  16
#define CANONICAL_INTRON_LOWQ  22

#define FINAL_CANONICAL_INTRON_HIGHQ 30 /* GT-AG */
#define FINAL_CANONICAL_INTRON_MEDQ  36
#define FINAL_CANONICAL_INTRON_LOWQ  42

#define KNOWN_SPLICESITE_REWARD 20



static bool novelsplicingp;

static IIT_T splicing_iit;
static int *splicing_divint_crosstable;
static int donor_typeint;
static int acceptor_typeint;

#define T Dynprog_T

void
Dynprog_genome_setup (bool novelsplicingp_in,
		      IIT_T splicing_iit_in, int *splicing_divint_crosstable_in,
		      int donor_typeint_in, int acceptor_typeint_in) {
  novelsplicingp = novelsplicingp_in;

  splicing_iit = splicing_iit_in;
  splicing_divint_crosstable = splicing_divint_crosstable_in;
  donor_typeint = donor_typeint_in;
  acceptor_typeint = acceptor_typeint_in;

  return;
}


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
		    genomicpos,Genome_get_char_blocks(pos)));
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



static int
intron_score (int *introntype, int leftdi, int rightdi, int cdna_direction, int canonical_reward, 
	      bool finalp) {
  int scoreI;

#ifdef USE_WEAK_SCOREI
  canonical_reward = CANONICAL_INTRON;
#endif

#ifdef PMAP
  if ((*introntype = leftdi & rightdi) == NONINTRON) {
    scoreI = 0.0;
  } else {
    switch (*introntype) {
    case GTAG_FWD: scoreI = canonical_reward; break;
    case GCAG_FWD: scoreI = finalp == true ? FINAL_GCAG_INTRON : GCAG_INTRON; break;
    case ATAC_FWD: scoreI = finalp == true ? FINAL_ATAC_INTRON : ATAC_INTRON; break;
    default: *introntype = NONINTRON; scoreI = 0.0;
    }
  }
#else
  if ((*introntype = leftdi & rightdi) == NONINTRON) {
    scoreI = 0.0;
  } else if (cdna_direction > 0) {
    switch (*introntype) {
    case GTAG_FWD: scoreI = canonical_reward; break;
    case GCAG_FWD: scoreI = finalp == true ? FINAL_GCAG_INTRON : GCAG_INTRON; break;
    case ATAC_FWD: scoreI = finalp == true ? FINAL_ATAC_INTRON : ATAC_INTRON; break;
    default: *introntype = NONINTRON; scoreI = 0.0;
    }
  } else if (cdna_direction < 0) {
    switch (*introntype) {
    case GTAG_REV: scoreI = canonical_reward; break;
    case GCAG_REV: scoreI = finalp == true ? FINAL_GCAG_INTRON : GCAG_INTRON; break;
    case ATAC_REV: scoreI = finalp == true ? FINAL_ATAC_INTRON : ATAC_INTRON; break;
    default: *introntype = NONINTRON; scoreI = 0.0;
    }
  } else {
    switch (*introntype) {
    case GTAG_FWD: case GTAG_REV: scoreI = canonical_reward; break;
    case GCAG_FWD: case GCAG_REV: scoreI = finalp == true ? FINAL_GCAG_INTRON : GCAG_INTRON; break;
    case ATAC_FWD: case ATAC_REV: scoreI = finalp == true ? FINAL_ATAC_INTRON : ATAC_INTRON; break;
    default: *introntype = NONINTRON; scoreI = 0.0;
    }
  }
#endif

  return scoreI;
}


static void
get_splicesite_probs (double *left_prob, double *right_prob, int cL, int cR,
		      int *left_known, int *right_known, Univcoord_T leftoffset, Univcoord_T rightoffset,
		      Univcoord_T chroffset, Univcoord_T chrhigh, int cdna_direction, bool watsonp) {
  Univcoord_T splicesitepos;
  
  if (left_known[cL] > 0) {
    debug9(printf("left position is known, so prob is 1.0\n"));
    *left_prob = 1.0;

  } else if (watsonp == true) {
    splicesitepos = chroffset + leftoffset + cL;
    if (cdna_direction > 0) {
      *left_prob = Maxent_hr_donor_prob(splicesitepos,chroffset);
      debug9(printf("1. donor splicesitepos is %u, prob %f, known %d\n",
		    splicesitepos,*left_prob,left_known[cL]));

    } else {
      *left_prob = Maxent_hr_antiacceptor_prob(splicesitepos,chroffset);
      debug9(printf("2. antiacceptor splicesitepos is %u, prob %f, known %d\n",
		    splicesitepos,*left_prob,left_known[cL]));

    }
  } else {
    splicesitepos = chrhigh - leftoffset - cL + 1;
    if (cdna_direction > 0) {
      *left_prob = Maxent_hr_antidonor_prob(splicesitepos,chroffset);
      debug9(printf("3. antidonor splicesitepos is %u, prob %f, known %d\n",
		    splicesitepos,*left_prob,left_known[cL]));

    } else {
      *left_prob = Maxent_hr_acceptor_prob(splicesitepos,chroffset);
      debug9(printf("4. acceptor splicesitepos is %u, prob %f, known %d\n",
		    splicesitepos,*left_prob,left_known[cL]));
    }
  }

  if (right_known[cR] > 0) {
    debug9(printf("right position is known, so prob is 1.0\n"));
    *right_prob = 1.0;

  } else if (watsonp == true) {
    splicesitepos = chroffset + rightoffset - cR + 1;
    if (cdna_direction > 0) {
      *right_prob = Maxent_hr_acceptor_prob(splicesitepos,chroffset);
      debug9(printf("5. acceptor splicesitepos is %u, prob %f, known %d\n",
		    splicesitepos,*right_prob,right_known[cR]));
    } else {
      *right_prob = Maxent_hr_antidonor_prob(splicesitepos,chroffset);
      debug9(printf("6. antidonor splicesitepos is %u, prob %f, known %d\n",
		    splicesitepos,*right_prob,right_known[cR]));

    }
  } else {
    splicesitepos = chrhigh - rightoffset + cR;
    if (cdna_direction > 0) {
      *right_prob = Maxent_hr_antiacceptor_prob(splicesitepos,chroffset);
      debug9(printf("7. antiacceptor splicesitepos is %u, prob %f, known %d\n",
		    splicesitepos,*right_prob,right_known[cR]));

    } else {
      *right_prob = Maxent_hr_donor_prob(splicesitepos,chroffset);
      debug9(printf("8. donor splicesitepos is %u, prob %f, known %d\n",
		    splicesitepos,*right_prob,right_known[cR]));
    }
  }

  return;
}


static void
get_known_splicesites (int *left_known, int *right_known, int glengthL, int glengthR,
		       int leftoffset, int rightoffset, int cdna_direction, bool watsonp,
		       Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh) {

  int *matches, nmatches, i;
  Univcoord_T splicesitepos;

  if (splicing_iit != NULL && donor_typeint >= 0 && acceptor_typeint >= 0) {
    /* Handle known splicing, splice site level */
    if (watsonp == true) {
      if (cdna_direction > 0) {
	/* splicesitepos = leftoffset + cL;  cL = 0 to < glengthL - 1 */
	matches = IIT_get_typed_signed_with_divno(&nmatches,splicing_iit,splicing_divint_crosstable[chrnum],
						  leftoffset+1,leftoffset+glengthL-2,donor_typeint,/*sign*/+1,/*sortp*/false);
	for (i = 0; i < nmatches; i++) {
	  splicesitepos = IIT_interval_low(splicing_iit,matches[i]);
	  debug5(printf("1. Found known donor at %u\n",splicesitepos));
	  left_known[splicesitepos - leftoffset] = KNOWN_SPLICESITE_REWARD;
	}
	FREE(matches);

	/* splicesitepos = rightoffset - cR + 1; cR = 0 to < glengthR - 1 */
	matches = IIT_get_typed_signed_with_divno(&nmatches,splicing_iit,splicing_divint_crosstable[chrnum],
						  rightoffset-glengthR+4,rightoffset+1,acceptor_typeint,/*sign*/+1,/*sortp*/false);
	for (i = 0; i < nmatches; i++) {
	  splicesitepos = IIT_interval_low(splicing_iit,matches[i]);
	  debug5(printf("2. Found known acceptor at %u\n",splicesitepos));
	  right_known[rightoffset - splicesitepos + 1] = KNOWN_SPLICESITE_REWARD;
	}
	FREE(matches);

      } else {
	/* splicesitepos = leftoffset + cL;  cL = 0 to < glengthL - 1 */
	matches = IIT_get_typed_signed_with_divno(&nmatches,splicing_iit,splicing_divint_crosstable[chrnum],
						  leftoffset+1,leftoffset+glengthL-2,acceptor_typeint,/*sign*/-1,/*sortp*/false);
	for (i = 0; i < nmatches; i++) {
	  splicesitepos = IIT_interval_low(splicing_iit,matches[i]);
	  debug5(printf("3. Found known antiacceptor at %u\n",splicesitepos));
	  left_known[splicesitepos - leftoffset] = KNOWN_SPLICESITE_REWARD;
	}
	FREE(matches);

	/* splicesitepos = rightoffset - cR + 1; cR = 0 to < glengthR - 1 */
	matches = IIT_get_typed_signed_with_divno(&nmatches,splicing_iit,splicing_divint_crosstable[chrnum],
						  rightoffset-glengthR+4,rightoffset+1,donor_typeint,/*sign*/-1,/*sortp*/false);
	for (i = 0; i < nmatches; i++) {
	  splicesitepos = IIT_interval_low(splicing_iit,matches[i]);
	  debug5(printf("4. Found known antidonor at %u\n",splicesitepos));
	  right_known[rightoffset - splicesitepos + 1] = KNOWN_SPLICESITE_REWARD;
	}
	FREE(matches);

      }

    } else {
      if (cdna_direction > 0) {
	/* splicesitepos = (chrhigh - chroffset) - leftoffset - cL + 1; cL = 0 to < glengthL - 1 */
	matches = IIT_get_typed_signed_with_divno(&nmatches,splicing_iit,splicing_divint_crosstable[chrnum],
						  (chrhigh - chroffset) - leftoffset - glengthL + 4,
						  (chrhigh - chroffset) - leftoffset + 1,
						  donor_typeint,/*sign*/-1,/*sortp*/false);
	for (i = 0; i < nmatches; i++) {
	  splicesitepos = IIT_interval_low(splicing_iit,matches[i]);
	  debug5(printf("5. Found known antidonor at %u\n",splicesitepos));
	  left_known[(chrhigh - chroffset) - leftoffset - splicesitepos + 1] = KNOWN_SPLICESITE_REWARD;
	}
	FREE(matches);

	/* splicesitepos = (chrhigh - chroffset) - rightoffset + cR; cR = 0 to < glengthR - 1 */
	matches = IIT_get_typed_signed_with_divno(&nmatches,splicing_iit,splicing_divint_crosstable[chrnum],
						  (chrhigh - chroffset) - rightoffset + 1,
						  (chrhigh - chroffset) - rightoffset + glengthR - 2,
						  acceptor_typeint,/*sign*/-1,/*sortp*/false);
	for (i = 0; i < nmatches; i++) {
	  splicesitepos = IIT_interval_low(splicing_iit,matches[i]);
	  debug5(printf("6. Found known antiacceptor at %u\n",splicesitepos));
	  right_known[splicesitepos - (chrhigh - chroffset) + rightoffset] = KNOWN_SPLICESITE_REWARD;
	}
	FREE(matches);

      } else {
	/* splicesitepos = (chrhigh - chroffset) - leftoffset - cL + 1; cL = 0 to < glengthL - 1 */
	matches = IIT_get_typed_signed_with_divno(&nmatches,splicing_iit,splicing_divint_crosstable[chrnum],
						  (chrhigh - chroffset) - leftoffset - glengthL + 4,
						  (chrhigh - chroffset) - leftoffset + 1,
						  acceptor_typeint,/*sign*/+1,/*sortp*/false);
	for (i = 0; i < nmatches; i++) {
	  splicesitepos = IIT_interval_low(splicing_iit,matches[i]);
	  debug5(printf("7. Found known acceptor at %u\n",splicesitepos));
	  left_known[(chrhigh - chroffset) - leftoffset - splicesitepos + 1] = KNOWN_SPLICESITE_REWARD;
	}
	FREE(matches);

	/* splicesitepos = (chrhigh - chroffset) - rightoffset + cR; cR = 0 to < glengthR - 1 */
	matches = IIT_get_typed_signed_with_divno(&nmatches,splicing_iit,splicing_divint_crosstable[chrnum],
						  (chrhigh - chroffset) - rightoffset + 1,
						  (chrhigh - chroffset) - rightoffset + glengthR - 2,
						  donor_typeint,/*sign*/+1,/*sortp*/false);
	for (i = 0; i < nmatches; i++) {
	  splicesitepos = IIT_interval_low(splicing_iit,matches[i]);
	  debug5(printf("8. Found known donor at %u\n",splicesitepos));
	  right_known[splicesitepos - (chrhigh - chroffset) + rightoffset] = KNOWN_SPLICESITE_REWARD;
	}
	FREE(matches);
	  
      }
    }

  } else if (splicing_iit != NULL) {
    /* Handle known splicing, intron level */
    if (watsonp == true) {
      if (cdna_direction > 0) {
	/* splicesitepos = leftoffset + cL; cL = 0 to < glengthL - 1 */
	matches = IIT_get_lows_signed(&nmatches,splicing_iit,splicing_divint_crosstable[chrnum],
				      leftoffset,leftoffset+glengthL-2,/*sign*/+1);
	for (i = 0; i < nmatches; i++) {
	  splicesitepos = IIT_interval_low(splicing_iit,matches[i]);
	  debug5(printf("1. Found known donor at %u\n",splicesitepos));
	  left_known[splicesitepos - leftoffset] = KNOWN_SPLICESITE_REWARD;
	}
	FREE(matches);

	/* splicesitepos+1U = rightoffset - cR + 2; cR = 0 to < glengthR - 1 */
	matches = IIT_get_highs_signed(&nmatches,splicing_iit,splicing_divint_crosstable[chrnum],
				       rightoffset-glengthR+4,rightoffset+2,/*sign*/+1);
	for (i = 0; i < nmatches; i++) {
	  splicesitepos = IIT_interval_high(splicing_iit,matches[i]);
	  debug5(printf("2. Found known acceptor at %u\n",splicesitepos));
	  right_known[rightoffset - splicesitepos + 2] = KNOWN_SPLICESITE_REWARD;
	}
	FREE(matches);

      } else {
	/* splicesitepos = leftoffset + cL; cL = 0 to < glengthL - 1 */
	matches = IIT_get_lows_signed(&nmatches,splicing_iit,splicing_divint_crosstable[chrnum],
				      leftoffset,leftoffset+glengthL-2,/*sign*/-1);
	for (i = 0; i < nmatches; i++) {
	  splicesitepos = IIT_interval_low(splicing_iit,matches[i]);
	  debug5(printf("3. Found known antiacceptor at %u\n",splicesitepos));
	  left_known[splicesitepos - leftoffset] = KNOWN_SPLICESITE_REWARD;
	}
	FREE(matches);

	/* splicesitepos+1U = rightoffset - cR + 2; cR = 0 to < glengthR - 1 */
	matches = IIT_get_highs_signed(&nmatches,splicing_iit,splicing_divint_crosstable[chrnum],
				       rightoffset-glengthR+4,rightoffset+2,/*sign*/-1);
	for (i = 0; i < nmatches; i++) {
	  splicesitepos = IIT_interval_high(splicing_iit,matches[i]);
	  debug5(printf("4. Found known antidonor at %u\n",splicesitepos));
	  right_known[rightoffset - splicesitepos + 2] = KNOWN_SPLICESITE_REWARD;
	}
	FREE(matches);

      }
    } else {
      if (cdna_direction > 0) {
	/* splicesitepos+1U = (chrhigh - chroffset) - leftoffset - cL + 2; cL = 0 to < glengthL - 1 */
	matches = IIT_get_highs_signed(&nmatches,splicing_iit,splicing_divint_crosstable[chrnum],
				       (chrhigh - chroffset) - leftoffset - glengthL + 4,
				       (chrhigh - chroffset) - leftoffset + 2,/*sign*/-1);
	for (i = 0; i < nmatches; i++) {
	  splicesitepos = IIT_interval_high(splicing_iit,matches[i]);
	  debug5(printf("5. Found known antidonor at %u\n",splicesitepos));
	  left_known[(chrhigh - chroffset) - leftoffset - splicesitepos + 2] = KNOWN_SPLICESITE_REWARD;
	}
	FREE(matches);

	/* splicesitepos = (chrhigh - chroffset) - rightoffset + cR; cR = 0 to < glengthR - 1 */
	matches = IIT_get_lows_signed(&nmatches,splicing_iit,splicing_divint_crosstable[chrnum],
				      (chrhigh - chroffset) - rightoffset,
				      (chrhigh - chroffset) - rightoffset + glengthR - 2,/*sign*/-1);
	for (i = 0; i < nmatches; i++) {
	  splicesitepos = IIT_interval_low(splicing_iit,matches[i]);
	  debug5(printf("6. Found known antiacceptor at %u\n",splicesitepos));
	  right_known[splicesitepos - (chrhigh - chroffset) + rightoffset] = KNOWN_SPLICESITE_REWARD;
	}
	FREE(matches);

      } else {
	/* splicesitepos+1U = (chrhigh - chroffset) - leftoffset - cL + 2; cL = 0 to < glengthL - 1 */
	matches = IIT_get_highs_signed(&nmatches,splicing_iit,splicing_divint_crosstable[chrnum],
				       (chrhigh - chroffset) - leftoffset - glengthL + 4,
				       (chrhigh - chroffset) - leftoffset + 2,/*sign*/+1);
	for (i = 0; i < nmatches; i++) {
	  splicesitepos = IIT_interval_high(splicing_iit,matches[i]);
	  debug5(printf("7. Found known acceptor at %u\n",splicesitepos));
	  left_known[(chrhigh - chroffset) - leftoffset - splicesitepos + 2] = KNOWN_SPLICESITE_REWARD;
	}
	FREE(matches);

	/* splicesitepos = (chrhigh - chroffset) - rightoffset + cR; cR = 0 to < glengthR - 1 */
	matches = IIT_get_lows_signed(&nmatches,splicing_iit,splicing_divint_crosstable[chrnum],
				      (chrhigh - chroffset) - rightoffset,
				      (chrhigh - chroffset) - rightoffset + glengthR - 2,/*sign*/+1);
	for (i = 0; i < nmatches; i++) {
	  splicesitepos = IIT_interval_low(splicing_iit,matches[i]);
	  debug5(printf("8. Found known donor at %u\n",splicesitepos));
	  right_known[splicesitepos - (chrhigh - chroffset) + rightoffset] = KNOWN_SPLICESITE_REWARD;
	}
	FREE(matches);
      }
    }
  }

  return;
}


#if defined(HAVE_SSE4_1) || defined(HAVE_SSE2)

/* Returns finalscore */
static int
bridge_intron_gap_8_intron_level (int *bestrL, int *bestrR, int *bestcL, int *bestcR,
				  int *best_introntype,
				  Score8_T **matrixL_upper, Score8_T **matrixL_lower,
				  Score8_T **matrixR_upper, Score8_T **matrixR_lower,
				  Direction8_T **directionsL_upper_nogap, Direction8_T **directionsL_lower_nogap, 
				  Direction8_T **directionsR_upper_nogap, Direction8_T **directionsR_lower_nogap,
				  int *left_known, int *right_known,
				  int rlength, int glengthL, int glengthR,
				  int cdna_direction, bool watsonp, int lbandL, int ubandL, int lbandR, int ubandR,
				  int leftoffset, int rightoffset,
				  Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh,
				  bool jump_late_p) {
  int rL, rR, cL, cR;
  int cloL, chighL;
  int cloR, chighR;
  int bestscore = NEG_INFINITY_8, score, scoreL, scoreR;
  Univcoord_T splicesitepos1, splicesitepos2;
  bool bestp;


  for (rL = 1, rR = rlength-1; rL < rlength; rL++, rR--) {
    debug3(printf("\nGenomic insert: At row %d on left and %d on right\n",rL,rR));
    if ((cloL = rL - lbandL) < 1) {
      cloL = 1;
    }
    if ((chighL = rL + ubandL) > glengthL-1) {
      chighL = glengthL-1;
    }
    
    if ((cloR = rR - lbandR) < 1) {
      cloR = 1;
    }
    if ((chighR = rR + ubandR) > glengthR-1) {
      chighR = glengthR-1;
    }

    /* Test indels on left and right */
    for (cL = cloL; cL < /* left of main diagonal*/rL; cL++) {
      /* The following check limits genomic inserts (horizontal) and
	 multiple cDNA inserts (vertical). */
      if (left_known[cL] > 0) {
	scoreL = (int) matrixL_lower[rL][cL];
	if (directionsL_lower_nogap[rL][cL] != DIAG) {
	  /* Favor gaps away from intron if possible */
	  scoreL -= 1;
	}

	/* Disallow leftoffset + cL >= rightoffset - cR, or cR >= rightoffset - leftoffset - cL */
	for (cR = cloR; cR < /* left of main diagonal*/rR && cR < rightoffset-leftoffset-cL; cR++) {
	  if (right_known[cR] > 0) {
	    scoreR = (int) matrixR_lower[rR][cR];
	    if (directionsR_lower_nogap[rR][cR] != DIAG) {
	      /* Favor gaps away from intron if possible */
	      scoreR -= 1;
	    }

	    if ((score = scoreL + scoreR) > bestscore ||
		(score >= bestscore && jump_late_p)) { /* Use >= for jump late */
	      bestp = false;
	      if (watsonp == true) {
		splicesitepos1 = leftoffset + cL;
		splicesitepos2 = rightoffset - cR + 1;
		if (IIT_exists_with_divno_signed(splicing_iit,splicing_divint_crosstable[chrnum],
						 splicesitepos1,splicesitepos2+1U,/*sign*/cdna_direction) == true) {
		  bestp = true;
		}
	      } else {
		splicesitepos1 = (chrhigh - chroffset) - leftoffset - cL + 1;
		splicesitepos2 = (chrhigh - chroffset) - rightoffset + cR;
		if (IIT_exists_with_divno_signed(splicing_iit,splicing_divint_crosstable[chrnum],
						 splicesitepos2,splicesitepos1+1U,/*sign*/-cdna_direction) == true) {
		  bestp = true;
		}
	      }
	      if (bestp == true) {
		debug3(printf("At %d left to %d right, score is (%d)+(%d) = %d (bestscore)\n",
			      cL,cR,scoreL,scoreR,score));
		bestscore = score;
		*bestrL = rL;
		*bestrR = rR;
		*bestcL = cL;
		*bestcR = cR;
	      } else {
		debug3a(printf("At %d left to %d right, score is (%d)+(%d) = %d\n",
			       cL,cR,scoreL,scoreR,score));
	      }
	    }
	  }
	}

	for (/* at main diagonal*/; cR < chighR && cR < rightoffset-leftoffset-cL; cR++) {
	  if (right_known[cR] > 0) {
	    scoreR = (int) matrixR_upper[cR][rR];
	    if (directionsR_upper_nogap[cR][rR] != DIAG) {
	      /* Favor gaps away from intron if possible */
	      scoreR -= 1;
	    }

	    if ((score = scoreL + scoreR) > bestscore ||
		(score >= bestscore && jump_late_p)) {  /* Use >= for jump late */
	      bestp = false;
	      if (watsonp == true) {
		splicesitepos1 = leftoffset + cL;
		splicesitepos2 = rightoffset - cR + 1;
		if (IIT_exists_with_divno_signed(splicing_iit,splicing_divint_crosstable[chrnum],
						 splicesitepos1,splicesitepos2+1U,/*sign*/cdna_direction) == true) {
		  bestp = true;
		}
	      } else {
		splicesitepos1 = (chrhigh - chroffset) - leftoffset - cL + 1;
		splicesitepos2 = (chrhigh - chroffset) - rightoffset + cR;
		if (IIT_exists_with_divno_signed(splicing_iit,splicing_divint_crosstable[chrnum],
						 splicesitepos2,splicesitepos1+1U,/*sign*/-cdna_direction) == true) {
		  bestp = true;
		}
	      }
	      if (bestp == true) {
		debug3(printf("At %d left to %d right, score is (%d)+(%d) = %d (bestscore)\n",
			      cL,cR,scoreL,scoreR,score));
		bestscore = score;
		*bestrL = rL;
		*bestrR = rR;
		*bestcL = cL;
		*bestcR = cR;
	      } else {
		debug3a(printf("At %d left to %d right, score is (%d)+(%d) = %d\n",
			       cL,cR,scoreL,scoreR,score));
	      }
	    }
	  }
	}
      }
    }

    for (/* at main diagonal*/; cL < chighL; cL++) {
      /* The following check limits genomic inserts (horizontal) and
	 multiple cDNA inserts (vertical). */
      if (left_known[cL] > 0) {
	scoreL = (int) matrixL_upper[cL][rL];
	if (directionsL_upper_nogap[cL][rL] != DIAG) {
	  /* Favor gaps away from intron if possible */
	  scoreL -= 1;
	}

	/* Disallow leftoffset + cL >= rightoffset - cR, or cR >= rightoffset - leftoffset - cL */
	for (cR = cloR; cR < /* left of main diagonal*/rR && cR < rightoffset-leftoffset-cL; cR++) {
	  if (right_known[cR] > 0) {
	    scoreR = (int) matrixR_lower[rR][cR];
	    if (directionsR_lower_nogap[rR][cR] != DIAG) {
	      /* Favor gaps away from intron if possible */
	      scoreR -= 1;
	    }

	    if ((score = scoreL + scoreR) > bestscore ||
		(score >= bestscore && jump_late_p)) {  /* Use >= for jump late */
	      bestp = false;
	      if (watsonp == true) {
		splicesitepos1 = leftoffset + cL;
		splicesitepos2 = rightoffset - cR + 1;
		if (IIT_exists_with_divno_signed(splicing_iit,splicing_divint_crosstable[chrnum],
						 splicesitepos1,splicesitepos2+1U,/*sign*/cdna_direction) == true) {
		  bestp = true;
		}
	      } else {
		splicesitepos1 = (chrhigh - chroffset) - leftoffset - cL + 1;
		splicesitepos2 = (chrhigh - chroffset) - rightoffset + cR;
		if (IIT_exists_with_divno_signed(splicing_iit,splicing_divint_crosstable[chrnum],
						 splicesitepos2,splicesitepos1+1U,/*sign*/-cdna_direction) == true) {
		  bestp = true;
		}
	      }
	      if (bestp == true) {
		debug3(printf("At %d left to %d right, score is (%d)+(%d) = %d (bestscore)\n",
			      cL,cR,scoreL,scoreR,score));
		bestscore = score;
		*bestrL = rL;
		*bestrR = rR;
		*bestcL = cL;
		*bestcR = cR;
	      } else {
		debug3a(printf("At %d left to %d right, score is (%d)+(%d) = %d\n",
			       cL,cR,scoreL,scoreR,score));
	      }
	    }
	  }
	}

	for (/* at main diagonal*/; cR < chighR && cR < rightoffset-leftoffset-cL; cR++) {
	  if (right_known[cR] > 0) {
	    scoreR = (int) matrixR_upper[cR][rR];
	    if (directionsR_upper_nogap[cR][rR] != DIAG) {
	      /* Favor gaps away from intron if possible */
	      scoreR -= 1;
	    }

	    if ((score = scoreL + scoreR) > bestscore ||
		(score >= bestscore && jump_late_p)) {  /* Use >= for jump late */
	      bestp = false;
	      if (watsonp == true) {
		splicesitepos1 = leftoffset + cL;
		splicesitepos2 = rightoffset - cR + 1;
		if (IIT_exists_with_divno_signed(splicing_iit,splicing_divint_crosstable[chrnum],
						 splicesitepos1,splicesitepos2+1U,/*sign*/cdna_direction) == true) {
		  bestp = true;
		}
	      } else {
		splicesitepos1 = (chrhigh - chroffset) - leftoffset - cL + 1;
		splicesitepos2 = (chrhigh - chroffset) - rightoffset + cR;
		if (IIT_exists_with_divno_signed(splicing_iit,splicing_divint_crosstable[chrnum],
						 splicesitepos2,splicesitepos1+1U,/*sign*/-cdna_direction) == true) {
		  bestp = true;
		}
	      }
	      if (bestp == true) {
		debug3(printf("At %d left to %d right, score is (%d)+(%d) = %d (bestscore)\n",
			      cL,cR,scoreL,scoreR,score));
		bestscore = score;
		*bestrL = rL;
		*bestrR = rR;
		*bestcL = cL;
		*bestcR = cR;
	      } else {
		debug3a(printf("At %d left to %d right, score is (%d)+(%d) = %d\n",
			       cL,cR,scoreL,scoreR,score));
	      }
	    }
	  }
	}
      }
    }
  }

  *best_introntype = NONINTRON;
  return bestscore;
}


/* Returns finalscore */
static int
bridge_intron_gap_8_site_level (int *bestrL, int *bestrR, int *bestcL, int *bestcR,
				Score8_T **matrixL_upper, Score8_T **matrixL_lower,
				Score8_T **matrixR_upper, Score8_T **matrixR_lower,
				Direction8_T **directionsL_upper_nogap, Direction8_T **directionsL_lower_nogap, 
				Direction8_T **directionsR_upper_nogap, Direction8_T **directionsR_lower_nogap,
				char *gsequenceL, char *gsequenceL_alt, char *rev_gsequenceR, char *rev_gsequenceR_alt,
				int goffsetL, int rev_goffsetR, int *left_known, int *right_known,
				int rlength, int glengthL, int glengthR,
				int cdna_direction, bool watsonp, int lbandL, int ubandL, int lbandR, int ubandR,
				int canonical_reward, int leftoffset, int rightoffset,
				Univcoord_T chroffset, Univcoord_T chrhigh,
				bool halfp, bool finalp) {
  int rL, rR, cL, cR;
  int bestrL_with_prob, bestrR_with_prob, bestcL_with_prob, bestcR_with_prob;
  int cloL, chighL;
  int cloR, chighR;
  int introntype;
  int bestscore = NEG_INFINITY_8, score, scoreL, scoreR, scoreI;
  int bestscore_with_prob = NEG_INFINITY_8;
  double *left_probabilities, *right_probabilities, probL, probR, probL_trunc, probR_trunc, bestprob, bestprob_trunc;
  Univcoord_T splicesitepos;
  char left1, left2, right2, right1, left1_alt, left2_alt, right2_alt, right1_alt;
  int *leftdi, *rightdi;
  bool use_prob_p;


  /* Read dinucleotides */
  leftdi = (int *) MALLOCA((glengthL+1) * sizeof(int));
  rightdi = (int *) MALLOCA((glengthR+1) * sizeof(int));

  for (cL = 0; cL < glengthL - 1; cL++) {
    left1 = gsequenceL[cL];
    left1_alt = gsequenceL_alt[cL];
    left2 = gsequenceL[cL+1];
    left2_alt = gsequenceL_alt[cL+1];
    assert(left1 == get_genomic_nt(&left1_alt,goffsetL+cL,chroffset,chrhigh,watsonp));
    assert(left2 == get_genomic_nt(&left2_alt,goffsetL+cL+1,chroffset,chrhigh,watsonp));

    if ((left1 == 'G' || left1_alt == 'G') && (left2 == 'T' || left2_alt == 'T')) {
      leftdi[cL] = LEFT_GT;
    } else if ((left1 == 'G' || left1_alt == 'G') && (left2 == 'C' || left2_alt == 'C')) {
      leftdi[cL] = LEFT_GC;
    } else if ((left1 == 'A' || left1_alt == 'A') && (left2 == 'T' || left2_alt == 'T')) {
      leftdi[cL] = LEFT_AT;
#ifndef PMAP
    } else if ((left1 == 'C' || left1_alt == 'C') && (left2 == 'T' || left2_alt == 'T')) {
      leftdi[cL] = LEFT_CT;
#endif
    } else {
      leftdi[cL] = 0x00;
    }
  }
  leftdi[glengthL-1] = leftdi[glengthL] = 0x00;

  for (cR = 0; cR < glengthR - 1; cR++) {
    right2 = rev_gsequenceR[-cR-1];
    right2_alt = rev_gsequenceR_alt[-cR-1];
    right1 = rev_gsequenceR[-cR];
    right1_alt = rev_gsequenceR_alt[-cR];
    assert(right2 == get_genomic_nt(&right2_alt,rev_goffsetR-cR-1,chroffset,chrhigh,watsonp));
    assert(right1 == get_genomic_nt(&right1_alt,rev_goffsetR-cR,chroffset,chrhigh,watsonp));

    if ((right2 == 'A' || right2_alt == 'A') && (right1 == 'G' || right1_alt == 'G')) {
      rightdi[cR] = RIGHT_AG;
    } else if ((right2 == 'A' || right2_alt == 'A') && (right1 == 'C' || right1_alt == 'C')) {
      rightdi[cR] = RIGHT_AC;
#ifndef PMAP
    } else if ((right2 == 'G' || right2_alt == 'G') && (right1 == 'C' || right1_alt == 'C')) {
      rightdi[cR] = RIGHT_GC;
    } else if ((right2 == 'A' || right2_alt == 'A') && (right1 == 'T' || right1_alt == 'T')) {
      rightdi[cR] = RIGHT_AT;
#endif
    } else {
      rightdi[cR] = 0x00;
    }
  }
  rightdi[glengthR-1] = rightdi[glengthR] = 0x00;


  left_probabilities = (double *) MALLOCA(glengthL * sizeof(double));
  right_probabilities = (double *) MALLOCA(glengthR * sizeof(double));

  debug3(printf("watsonp is %d.  cdna_direction is %d\n",watsonp,cdna_direction));
  if (watsonp == true) {
    if (cdna_direction > 0) {
      for (cL = 0; cL < glengthL - 1; cL++) {
	splicesitepos = chroffset + leftoffset + cL;
	if (left_known[cL]) {
	  left_probabilities[cL] = 1.0;
	} else {
	  left_probabilities[cL] = Maxent_hr_donor_prob(splicesitepos,chroffset);
	  debug3(printf("left donor probability at cL %d is %f\n",cL,left_probabilities[cL]));
	}
      }

      for (cR = 0; cR < glengthR - 1; cR++) {
	splicesitepos = chroffset + rightoffset - cR + 1;
	if (right_known[cR]) {
	  right_probabilities[cR] = 1.0;
	} else {
	  right_probabilities[cR] = Maxent_hr_acceptor_prob(splicesitepos,chroffset);
	  debug3(printf("right acceptor probability at cR %d is %f\n",cR,right_probabilities[cR]));
	}
      }

    } else {
      for (cL = 0; cL < glengthL - 1; cL++) {
	splicesitepos = chroffset + leftoffset + cL;
	if (left_known[cL]) {
	  left_probabilities[cL] = 1.0;
	} else {
	  left_probabilities[cL] = Maxent_hr_antiacceptor_prob(splicesitepos,chroffset);
	  debug3(printf("left antiacceptor probability at cL %d is %f\n",cL,left_probabilities[cL]));
	}
      }

      for (cR = 0; cR < glengthR - 1; cR++) {
	splicesitepos = chroffset + rightoffset - cR + 1;
	if (right_known[cR]) {
	  right_probabilities[cR] = 1.0;
	} else {
	  right_probabilities[cR] = Maxent_hr_antidonor_prob(splicesitepos,chroffset);
	  debug3(printf("right antidonor probability at cR %d is %f\n",cR,right_probabilities[cR]));
	}
      }
    }

  } else {
    if (cdna_direction > 0) {
      for (cL = 0; cL < glengthL - 1; cL++) {
	splicesitepos = chrhigh - leftoffset - cL + 1;
	if (left_known[cL]) {
	  left_probabilities[cL] = 1.0;
	} else {
	  left_probabilities[cL] = Maxent_hr_antidonor_prob(splicesitepos,chroffset);
	  debug3(printf("left antidonor probability at cL %d is %f\n",cL,left_probabilities[cL]));
	}
      }

      for (cR = 0; cR < glengthR - 1; cR++) {
	splicesitepos = chrhigh - rightoffset + cR;
	if (right_known[cR]) {
	  right_probabilities[cR] = 1.0;
	} else {
	  right_probabilities[cR] = Maxent_hr_antiacceptor_prob(splicesitepos,chroffset);
	  debug3(printf("right antiacceptor probability at cR %d is %f\n",cR,right_probabilities[cR]));
	}
      }

    } else {
      for (cL = 0; cL < glengthL - 1; cL++) {
	splicesitepos = chrhigh - leftoffset - cL + 1;
	if (left_known[cL]) {
	  left_probabilities[cL] = 1.0;
	} else {
	  left_probabilities[cL] = Maxent_hr_acceptor_prob(splicesitepos,chroffset);
	  debug3(printf("left acceptor probability at cL %d is %f\n",cL,left_probabilities[cL]));
	}
      }

      for (cR = 0; cR < glengthR - 1; cR++) {
	splicesitepos = chrhigh - rightoffset + cR;
	if (right_known[cR]) {
	  right_probabilities[cR] = 1.0;
	} else {
	  right_probabilities[cR] = Maxent_hr_donor_prob(splicesitepos,chroffset);
	  debug3(printf("right donor probability at cR %d is %f\n",cR,right_probabilities[cR]));
	}
      }
    }
  }

  /* Search using probs and without simultaneously */
  bestscore = NEG_INFINITY_8;
  bestprob = bestprob_trunc = 0.0;
  for (rL = 1, rR = rlength-1; rL < rlength; rL++, rR--) {
    debug3(printf("\nAt row %d on left and %d on right\n",rL,rR));
    if ((cloL = rL - lbandL) < 1) {
      cloL = 1;
    }
    if ((chighL = rL + ubandL) > glengthL-1) {
      chighL = glengthL-1;
    }

    if ((cloR = rR - lbandR) < 1) {
      cloR = 1;
    }
    if ((chighR = rR + ubandR) > glengthR-1) {
      chighR = glengthR-1;
    }

    debug3(printf("A. Test no indels\n"));
    cL = rL;
    probL = left_probabilities[cL];
    if (probL > PROB_CEILING) {
      probL_trunc = PROB_CEILING;
    } else if (probL < PROB_FLOOR) {
      probL_trunc = PROB_FLOOR;
    } else {
      probL_trunc = probL;
    }
    scoreL = (int) matrixL_upper[cL][rL];

    cR = rR;
    probR = right_probabilities[cR];
    if (probR > PROB_CEILING) {
      probR_trunc = PROB_CEILING;
    } else if (probR < PROB_FLOOR) {
      probR_trunc = PROB_FLOOR;
    } else {
      probR_trunc = probR;
    }
    scoreR = (int) matrixR_upper[cR][rR];
    

#ifdef USE_SCOREI
    scoreI = intron_score(&introntype,leftdi[cL],rightdi[cR],cdna_direction,canonical_reward,finalp);
#else
    scoreI = 0;
#endif

    if ((score = scoreL + scoreI + scoreR) > bestscore) {
      debug3(printf("Best score: At %d left to %d right, score is (%d)+(%d)+(%d) = %d (bestscore, prob %f + %f)\n",
		    cL,cR,scoreL,scoreI,scoreR,scoreL+scoreI+scoreR,probL,probR));
      debug3(printf("probL %f, probR %f\n",left_probabilities[cL],right_probabilities[cR]));
      bestscore = score;
      *bestrL = rL;
      *bestrR = rR;
      *bestcL = cL;
      *bestcR = cR;
      bestprob = probL + probR;
    } else if (score == bestscore && probL + probR > bestprob) {
      debug3(printf("Improved prob: At %d left to %d right, score is (%d)+(%d)+(%d) = %d (bestscore, prob %f + %f)\n",
		    cL,cR,scoreL,scoreI,scoreR,scoreL+scoreI+scoreR,probL,probR));
      debug3(printf("probL %f, probR %f\n",left_probabilities[cL],right_probabilities[cR]));
      *bestrL = rL;
      *bestrR = rR;
      *bestcL = cL;
      *bestcR = cR;
      bestprob = probL + probR;
    } else {
      debug3a(printf("Not best score: At %d left to %d right, score is (%d)+(%d)+(%d) = %d (bestscore, prob %f + %f)\n",
		     cL,cR,scoreL,scoreI,scoreR,scoreL+scoreI+scoreR,probL,probR));
    }
	
    if (probL_trunc + probR_trunc < bestprob_trunc) {
      debug3a(printf("At %d left to %d right, prob is %f + %f = %f\n",
		     cL,cR,probL_trunc,probR_trunc,probL_trunc+probR_trunc));
      
    } else if (probL_trunc + probR_trunc == bestprob_trunc) {
      debug3(printf("At %d left to %d right, prob is %f + %f = %f\n",
		    cL,cR,probL_trunc,probR_trunc,probL_trunc+probR_trunc));
      
      if (scoreL + scoreI + scoreR > bestscore_with_prob) {
	debug3(printf(" (bestscore %d)\n",scoreL+scoreR));
	bestprob_trunc = probL_trunc + probR_trunc;
	bestcL_with_prob = cL;
	bestcR_with_prob = cR;
	bestrL_with_prob = rL;
	bestrR_with_prob = rR;
	bestscore_with_prob = scoreL + scoreI + scoreR;
      }
	  
    } else {
      /* probL_trunc + probR_trunc > bestprob_trunc */
      debug3(printf("At %d left to %d right, prob is %f + %f = %f\n",
		    cL,cR,probL_trunc,probR_trunc,probL_trunc+probR_trunc));
      
      debug3(printf(" (bestscore %d)\n",scoreL+scoreR));
      bestprob_trunc = probL_trunc + probR_trunc;
      bestcL_with_prob = cL;
      bestcR_with_prob = cR;
      bestrL_with_prob = rL;
      bestrR_with_prob = rR;
      bestscore_with_prob = scoreL + scoreI + scoreR;
    }


    debug3(printf("B. Test indel on right\n"));
    /* Test indel on right */
    cL = rL;
    probL = left_probabilities[cL];
    if (probL > PROB_CEILING) {
      probL_trunc = PROB_CEILING;
    } else if (probL < PROB_FLOOR) {
      probL_trunc = PROB_FLOOR;
    } else {
      probL_trunc = probL;
    }
    scoreL = (int) matrixL_upper[cL][rL];
    if (directionsL_upper_nogap[cL][rL] != DIAG) {
      /* Favor gaps away from intron if possible */
      scoreL -= 1;
    }

    /* Disallow leftoffset + cL >= rightoffset - cR, or cR >= rightoffset - leftoffset - cL */
    for (cR = cloR; cR < /*to main diagonal*/rR && cR < rightoffset-leftoffset-cL; cR++) {
      probR = right_probabilities[cR];
      if (probR > PROB_CEILING) {
	probR_trunc = PROB_CEILING;
      } else if (probR < PROB_FLOOR) {
	probR_trunc = PROB_FLOOR;
      } else {
	probR_trunc = probR;
      }
      scoreR = (int) matrixR_lower[rR][cR];
      if (directionsR_lower_nogap[rR][cR] != DIAG) {
	/* Favor gaps away from intron if possible */
	scoreR -= 1;
      }
	
#ifdef USE_SCOREI
      scoreI = intron_score(&introntype,leftdi[cL],rightdi[cR],cdna_direction,canonical_reward,finalp);
#else
      scoreI = 0;
#endif
	
      if ((score = scoreL + scoreI + scoreR) > bestscore) {
	debug3(printf("Best score: At %d left to %d right, score is (%d)+(%d)+(%d) = %d (bestscore, prob %f + %f)\n",
		      cL,cR,scoreL,scoreI,scoreR,scoreL+scoreI+scoreR,probL,probR));
	debug3(printf("probL %f, probR %f\n",left_probabilities[cL],right_probabilities[cR]));
	bestscore = score;
	*bestrL = rL;
	*bestrR = rR;
	*bestcL = cL;
	*bestcR = cR;
	bestprob = probL + probR;
      } else if (score == bestscore && probL + probR > bestprob) {
	debug3(printf("Improved prob: At %d left to %d right, score is (%d)+(%d)+(%d) = %d (bestscore, prob %f + %f)\n",
		      cL,cR,scoreL,scoreI,scoreR,scoreL+scoreI+scoreR,probL,probR));
	debug3(printf("probL %f, probR %f\n",left_probabilities[cL],right_probabilities[cR]));
	*bestrL = rL;
	*bestrR = rR;
	*bestcL = cL;
	*bestcR = cR;
	bestprob = probL + probR;
      } else {
	debug3a(printf("Not best score: At %d left to %d right, score is (%d)+(%d)+(%d) = %d (bestscore, prob %f + %f)\n",
		       cL,cR,scoreL,scoreI,scoreR,scoreL+scoreI+scoreR,probL,probR));
      }
	
      if (probL_trunc + probR_trunc < bestprob_trunc) {
	debug3a(printf("At %d left to %d right, prob is %f + %f = %f\n",
		       cL,cR,probL_trunc,probR_trunc,probL_trunc+probR_trunc));
	  
      } else if (probL_trunc + probR_trunc == bestprob_trunc) {
	debug3(printf("At %d left to %d right, prob is %f + %f = %f\n",
		      cL,cR,probL_trunc,probR_trunc,probL_trunc+probR_trunc));
	  
	if (scoreL + scoreI + scoreR > bestscore_with_prob) {
	  debug3(printf(" (bestscore %d)\n",scoreL+scoreR));
	  bestprob_trunc = probL_trunc + probR_trunc;
	  bestcL_with_prob = cL;
	  bestcR_with_prob = cR;
	  bestrL_with_prob = rL;
	  bestrR_with_prob = rR;
	  bestscore_with_prob = scoreL + scoreI + scoreR;
	}
		  
      } else {
	/* probL_trunc + probR_trunc > bestprob_trunc */
	debug3(printf("At %d left to %d right, prob is %f + %f = %f\n",
		      cL,cR,probL_trunc,probR_trunc,probL_trunc+probR_trunc));
	  
	debug3(printf(" (bestscore %d)\n",scoreL+scoreR));
	bestprob_trunc = probL_trunc + probR_trunc;
	bestcL_with_prob = cL;
	bestcR_with_prob = cR;
	bestrL_with_prob = rL;
	bestrR_with_prob = rR;
	bestscore_with_prob = scoreL + scoreI + scoreR;
      }
    }

    debug3(printf("Skip main diagonal\n"));
    for (/*skip main diagonal*/cR++; cR < chighR && cR < rightoffset-leftoffset-cL; cR++) {
      probR = right_probabilities[cR];
      if (probR > PROB_CEILING) {
	probR_trunc = PROB_CEILING;
      } else if (probR < PROB_FLOOR) {
	probR_trunc = PROB_FLOOR;
      } else {
	probR_trunc = probR;
      }
      scoreR = (int) matrixR_upper[cR][rR];
      if (directionsR_upper_nogap[cR][rR] != DIAG) {
	/* Favor gaps away from intron if possible */
	scoreR -= 1;
      }
	
#ifdef USE_SCOREI
      scoreI = intron_score(&introntype,leftdi[cL],rightdi[cR],cdna_direction,canonical_reward,finalp);
#else
      scoreI = 0;
#endif
	
      if ((score = scoreL + scoreI + scoreR) > bestscore) {
	debug3(printf("Best score: At %d left to %d right, score is (%d)+(%d)+(%d) = %d (bestscore, prob %f + %f)\n",
		      cL,cR,scoreL,scoreI,scoreR,scoreL+scoreI+scoreR,probL,probR));
	debug3(printf("probL %f, probR %f\n",left_probabilities[cL],right_probabilities[cR]));
	bestscore = score;
	*bestrL = rL;
	*bestrR = rR;
	*bestcL = cL;
	*bestcR = cR;
	bestprob = probL + probR;
      } else if (score == bestscore && probL + probR > bestprob) {
	debug3(printf("Improved prob: At %d left to %d right, score is (%d)+(%d)+(%d) = %d (bestscore, prob %f + %f)\n",
		      cL,cR,scoreL,scoreI,scoreR,scoreL+scoreI+scoreR,probL,probR));
	debug3(printf("probL %f, probR %f\n",left_probabilities[cL],right_probabilities[cR]));
	*bestrL = rL;
	*bestrR = rR;
	*bestcL = cL;
	*bestcR = cR;
	bestprob = probL + probR;
      } else {
	debug3a(printf("Not best score: At %d left to %d right, score is (%d)+(%d)+(%d) = %d (bestscore, prob %f + %f)\n",
		       cL,cR,scoreL,scoreI,scoreR,scoreL+scoreI+scoreR,probL,probR));
      }
	
      if (probL_trunc + probR_trunc < bestprob_trunc) {
	debug3a(printf("At %d left to %d right, prob is %f + %f = %f\n",
		       cL,cR,probL_trunc,probR_trunc,probL_trunc+probR_trunc));
	  
      } else if (probL_trunc + probR_trunc == bestprob_trunc) {
	debug3(printf("At %d left to %d right, prob is %f + %f = %f\n",
		      cL,cR,probL_trunc,probR_trunc,probL_trunc+probR_trunc));
	  
	if (scoreL + scoreI + scoreR > bestscore_with_prob) {
	  debug3(printf(" (bestscore %d)\n",scoreL+scoreR));
	  bestprob_trunc = probL_trunc + probR_trunc;
	  bestcL_with_prob = cL;
	  bestcR_with_prob = cR;
	  bestrL_with_prob = rL;
	  bestrR_with_prob = rR;
	  bestscore_with_prob = scoreL + scoreI + scoreR;
	}
		  
      } else {
	/* probL_trunc + probR_trunc > bestprob_trunc */
	debug3(printf("At %d left to %d right, prob is %f + %f = %f\n",
		      cL,cR,probL_trunc,probR_trunc,probL_trunc+probR_trunc));
	  
	debug3(printf(" (bestscore %d)\n",scoreL+scoreR));
	bestprob_trunc = probL_trunc + probR_trunc;
	bestcL_with_prob = cL;
	bestcR_with_prob = cR;
	bestrL_with_prob = rL;
	bestrR_with_prob = rR;
	bestscore_with_prob = scoreL + scoreI + scoreR;
      }
    }

    debug3(printf("C. Test indel on left\n"));
    /* Test indel on left */
    cR = rR;
    probR = right_probabilities[cR];
    if (probR > PROB_CEILING) {
      probR_trunc = PROB_CEILING;
    } else if (probR < PROB_FLOOR) {
      probR_trunc = PROB_FLOOR;
    } else {
      probR_trunc = probR;
    }
    scoreR = (int) matrixR_upper[cR][rR];
    if (directionsR_upper_nogap[cR][rR] != DIAG) {
      /* Favor gaps away from intron if possible */
      scoreR -= 1;
    }

    /* Disallow leftoffset + cL >= rightoffset - cR, or cR >= rightoffset - leftoffset - cL */
    for (cL = cloL; cL < /*to main diagonal*/rL && cL < rightoffset-leftoffset-cR; cL++) {
      probL = left_probabilities[cL];
      if (probL > PROB_CEILING) {
	probL_trunc = PROB_CEILING;
      } else if (probL < PROB_FLOOR) {
	probL_trunc = PROB_FLOOR;
      } else {
	probL_trunc = probL;
      }
      scoreL = (int) matrixL_lower[rL][cL];
      if (directionsL_lower_nogap[rL][cL] != DIAG) {
	/* Favor gaps away from intron if possible */
	scoreL -= 1;
      }

#ifdef USE_SCOREI
      scoreI = intron_score(&introntype,leftdi[cL],rightdi[cR],cdna_direction,canonical_reward,finalp);
#else
      scoreI = 0;
#endif

      if ((score = scoreL + scoreI + scoreR) > bestscore) {
	debug3(printf("Best score: At %d left to %d right, score is (%d)+(%d)+(%d) = %d (bestscore, prob %f + %f)\n",
		      cL,cR,scoreL,scoreI,scoreR,scoreL+scoreI+scoreR,probL,probR));
	debug3(printf("probL %f, probR %f\n",left_probabilities[cL],right_probabilities[cR]));
	bestscore = score;
	*bestrL = rL;
	*bestrR = rR;
	*bestcL = cL;
	*bestcR = cR;
	bestprob = probL + probR;
      } else if (score == bestscore && probL + probR > bestprob) {
	debug3(printf("Improved prob: At %d left to %d right, score is (%d)+(%d)+(%d) = %d (bestscore, prob %f + %f)\n",
		      cL,cR,scoreL,scoreI,scoreR,scoreL+scoreI+scoreR,probL,probR));
	debug3(printf("probL %f, probR %f\n",left_probabilities[cL],right_probabilities[cR]));
	*bestrL = rL;
	*bestrR = rR;
	*bestcL = cL;
	*bestcR = cR;
	bestprob = probL + probR;
      } else {
	debug3a(printf("Not best score: At %d left to %d right, score is (%d)+(%d)+(%d) = %d (bestscore, prob %f + %f)\n",
		       cL,cR,scoreL,scoreI,scoreR,scoreL+scoreI+scoreR,probL,probR));
      }
	
      if (probL_trunc + probR_trunc < bestprob_trunc) {
	debug3a(printf("At %d left to %d right, prob is %f + %f = %f\n",
		       cL,cR,probL_trunc,probR_trunc,probL_trunc+probR_trunc));
	  
      } else if (probL_trunc + probR_trunc == bestprob_trunc) {
	debug3(printf("At %d left to %d right, prob is %f + %f = %f\n",
		      cL,cR,probL_trunc,probR_trunc,probL_trunc+probR_trunc));
	  
	if (scoreL + scoreI + scoreR > bestscore_with_prob) {
	  debug3(printf(" (bestscore %d)\n",scoreL+scoreR));
	  bestprob_trunc = probL_trunc + probR_trunc;
	  bestcL_with_prob = cL;
	  bestcR_with_prob = cR;
	  bestrL_with_prob = rL;
	  bestrR_with_prob = rR;
	  bestscore_with_prob = scoreL + scoreI + scoreR;
	}
	  
      } else {
	/* probL_trunc + probR_trunc > bestprob_trunc */
	debug3(printf("At %d left to %d right, prob is %f + %f = %f\n",
		      cL,cR,probL_trunc,probR_trunc,probL_trunc+probR_trunc));
	  
	debug3(printf(" (bestscore %d)\n",scoreL+scoreR));
	bestprob_trunc = probL_trunc + probR_trunc;
	bestcL_with_prob = cL;
	bestcR_with_prob = cR;
	bestrL_with_prob = rL;
	bestrR_with_prob = rR;
	bestscore_with_prob = scoreL + scoreI + scoreR;
      }
    }

    debug3(printf("Skip main diagonal\n"));
    for (/*Skip main diagonal*/cL++; cL < chighL && cL < rightoffset-leftoffset-cR; cL++) {
      probL = left_probabilities[cL];
      if (probL > PROB_CEILING) {
	probL_trunc = PROB_CEILING;
      } else if (probL < PROB_FLOOR) {
	probL_trunc = PROB_FLOOR;
      } else {
	probL_trunc = probL;
      }
      scoreL = (int) matrixL_upper[cL][rL];
      if (directionsL_upper_nogap[cL][rL] != DIAG) {
	/* Favor gaps away from intron if possible */
	scoreL -= 1;
      }

#ifdef USE_SCOREI
      scoreI = intron_score(&introntype,leftdi[cL],rightdi[cR],cdna_direction,canonical_reward,finalp);
#else
      scoreI = 0;
#endif

      if ((score = scoreL + scoreI + scoreR) > bestscore) {
	debug3(printf("Best score: At %d left to %d right, score is (%d)+(%d)+(%d) = %d (bestscore, prob %f + %f)\n",
		      cL,cR,scoreL,scoreI,scoreR,scoreL+scoreI+scoreR,probL,probR));
	debug3(printf("probL %f, probR %f\n",left_probabilities[cL],right_probabilities[cR]));
	bestscore = score;
	*bestrL = rL;
	*bestrR = rR;
	*bestcL = cL;
	*bestcR = cR;
	bestprob = probL + probR;
      } else if (score == bestscore && probL + probR > bestprob) {
	debug3(printf("Improved prob: At %d left to %d right, score is (%d)+(%d)+(%d) = %d (bestscore, prob %f + %f)\n",
		      cL,cR,scoreL,scoreI,scoreR,scoreL+scoreI+scoreR,probL,probR));
	debug3(printf("probL %f, probR %f\n",left_probabilities[cL],right_probabilities[cR]));
	*bestrL = rL;
	*bestrR = rR;
	*bestcL = cL;
	*bestcR = cR;
	bestprob = probL + probR;
      } else {
	debug3a(printf("Not best score: At %d left to %d right, score is (%d)+(%d)+(%d) = %d (bestscore, prob %f + %f)\n",
		       cL,cR,scoreL,scoreI,scoreR,scoreL+scoreI+scoreR,probL,probR));
      }
	
      if (probL_trunc + probR_trunc < bestprob_trunc) {
	debug3a(printf("At %d left to %d right, prob is %f + %f = %f\n",
		       cL,cR,probL_trunc,probR_trunc,probL_trunc+probR_trunc));
	  
      } else if (probL_trunc + probR_trunc == bestprob_trunc) {
	debug3(printf("At %d left to %d right, prob is %f + %f = %f\n",
		      cL,cR,probL_trunc,probR_trunc,probL_trunc+probR_trunc));
	  
	if (scoreL + scoreI + scoreR > bestscore_with_prob) {
	  debug3(printf(" (bestscore %d)\n",scoreL+scoreR));
	  bestprob_trunc = probL_trunc + probR_trunc;
	  bestcL_with_prob = cL;
	  bestcR_with_prob = cR;
	  bestrL_with_prob = rL;
	  bestrR_with_prob = rR;
	  bestscore_with_prob = scoreL + scoreI + scoreR;
	}
	  
      } else {
	/* probL_trunc + probR_trunc > bestprob_trunc */
	debug3(printf("At %d left to %d right, prob is %f + %f = %f\n",
		      cL,cR,probL_trunc,probR_trunc,probL_trunc+probR_trunc));
	  
	debug3(printf(" (bestscore %d)\n",scoreL+scoreR));
	bestprob_trunc = probL_trunc + probR_trunc;
	bestcL_with_prob = cL;
	bestcR_with_prob = cR;
	bestrL_with_prob = rL;
	bestrR_with_prob = rR;
	bestscore_with_prob = scoreL + scoreI + scoreR;
      }
    }
  }

  if (bestprob > 2*PROB_CEILING) {
    /* Probability is good with best alignment, so take that */
    debug3(printf("Best alignment has good probability\n"));
    use_prob_p = true;
  } else if (left_probabilities[bestcL_with_prob] < PROB_CEILING && right_probabilities[bestcR_with_prob] < PROB_CEILING) {
    /* Probability-based solution is bad, so use alignment */
    debug3(printf("Probability-based solution is bad\n"));
    use_prob_p = false;
  } else if (bestscore_with_prob < bestscore - 9) {
    debug3(printf("Probability-based solution requires very bad alignment, because bestscore_with_prob %d < bestscore %d - 9\n",
		  bestscore_with_prob,bestscore));
    use_prob_p = false;
  } else {
    use_prob_p = true;
  }

  if (use_prob_p == true) {
    /* Best alignment yields bad probability, and probability-based alignment yields good probability, so switch */
    debug3(printf("Switch to probability-based solution\n"));
    debug3(printf("SIMD 8. bestscore %d (bestprob %f) vs bestscore_with_prob %d (bestprob_trunc %f, actually %f and %f)\n",
		 bestscore,bestprob,bestscore_with_prob,bestprob_trunc,left_probabilities[bestcL_with_prob],right_probabilities[bestcR_with_prob]));
    *bestcL = bestcL_with_prob;
    *bestcR = bestcR_with_prob;
    *bestrL = bestrL_with_prob;
    *bestrR = bestrR_with_prob;
    bestscore = bestscore_with_prob;
  }
    
  FREEA(rightdi);
  FREEA(leftdi);
  FREEA(left_probabilities);
  FREEA(right_probabilities);

  if (bestscore < 0) {
    return bestscore;
  } else if (halfp == true) {
    scoreI = intron_score(&introntype,leftdi[*bestcL],rightdi[*bestcR],cdna_direction,canonical_reward,finalp);
    return (bestscore - scoreI/2);
  } else {
    return bestscore;
  }
}


static int
bridge_intron_gap_8_ud (int *bestrL, int *bestrR, int *bestcL, int *bestcR,
			int *best_introntype, double *left_prob, double *right_prob,
			Score8_T **matrixL_upper, Score8_T **matrixL_lower,
			Score8_T **matrixR_upper, Score8_T **matrixR_lower,
			Direction8_T **directionsL_upper_nogap, Direction8_T **directionsL_lower_nogap, 
			Direction8_T **directionsR_upper_nogap, Direction8_T **directionsR_lower_nogap,
			char *gsequenceL, char *gsequenceL_alt, char *rev_gsequenceR, char *rev_gsequenceR_alt,
			int goffsetL, int rev_goffsetR, int rlength, int glengthL, int glengthR,
			int cdna_direction, bool watsonp, int lbandL, int ubandL, int lbandR, int ubandR,
			double defect_rate, int canonical_reward, int leftoffset, int rightoffset,
			Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh,
			bool halfp, bool finalp, bool jump_late_p) {
  int finalscore;
  int *left_known, *right_known;

  debug(printf("Running bridge_intron_gap_8_ud\n"));

  if (glengthL+1 <= 0) {
    fprintf(stderr,"Problem with glengthL = %d\n",glengthL);
    abort();
  }

  if (glengthR+1 <= 0) {
    fprintf(stderr,"Problem with glengthR = %d\n",glengthR);
    abort();
  }

  left_known = (int *) CALLOCA(glengthL+1,sizeof(int));
  right_known = (int *) CALLOCA(glengthR+1,sizeof(int));
  get_known_splicesites(left_known,right_known,glengthL,glengthR,
			/*leftoffset*/goffsetL,/*rightoffset*/rev_goffsetR,
			cdna_direction,watsonp,chrnum,chroffset,chrhigh);

  if (novelsplicingp == false && splicing_iit != NULL && (donor_typeint < 0 || acceptor_typeint < 0)) {
    /* Constrain to given introns */
    finalscore = bridge_intron_gap_8_intron_level(&(*bestrL),&(*bestrR),&(*bestcL),&(*bestcR),&(*best_introntype),
						  matrixL_upper,matrixL_lower,matrixR_upper,matrixR_lower,
						  directionsL_upper_nogap,directionsL_lower_nogap, 
						  directionsR_upper_nogap,directionsR_lower_nogap,
						  left_known,right_known,rlength,glengthL,glengthR,
						  cdna_direction,watsonp,lbandL,ubandL,lbandR,ubandR,
						  leftoffset,rightoffset,chrnum,chroffset,chrhigh,jump_late_p);
  } else {
    finalscore = bridge_intron_gap_8_site_level(&(*bestrL),&(*bestrR),&(*bestcL),&(*bestcR),
						matrixL_upper,matrixL_lower,matrixR_upper,matrixR_lower,
						directionsL_upper_nogap,directionsL_lower_nogap, 
						directionsR_upper_nogap,directionsR_lower_nogap,
						gsequenceL,gsequenceL_alt,rev_gsequenceR,rev_gsequenceR_alt,goffsetL,rev_goffsetR,
						left_known,right_known,rlength,glengthL,glengthR,
						cdna_direction,watsonp,lbandL,ubandL,lbandR,ubandR,
						canonical_reward,leftoffset,rightoffset,
						chroffset,chrhigh,halfp,finalp);
  }


#if 0
  /* Determine if result meets given constraints */
  if (*finalscore < 0) {
    result = false;
  } else if (novelsplicingp == true) {
    result = true;
  } else if (splicing_iit == NULL) {
    result = true;
  } else if (donor_typeint >= 0 && acceptor_typeint >= 0) {
    /* If novelsplicingp is false and using splicing at splice site level, require sites to be known */
    if (left_known[*bestcL] == 0 || right_known[*bestcR] == 0) {
      debug(printf("Novel splicing not allowed, so bridge_intron_gap returning false\n"));
      result = false;
    } else {
      result = true;
    }
  } else {
    /* If novelsplicingp is false and using splicing at splice site level, result was already constrained */
    result = true;
  }
#endif

  if (finalscore >= 0) {
    get_splicesite_probs(&(*left_prob),&(*right_prob),*bestcL,*bestcR,
			 left_known,right_known,leftoffset,rightoffset,chroffset,chrhigh,
			 cdna_direction,watsonp);
  }

  FREEA(right_known);
  FREEA(left_known);

#if defined(DEBUG) || defined(DEBUG3)
  printf("Returning final score of %d at (%d,%d) left to (%d,%d) right, with probs %f and %f\n",
	 finalscore,*bestrL,*bestcL,*bestrR,*bestcR,*left_prob,*right_prob);
#endif

  return finalscore;
}
#endif


#if defined(HAVE_SSE2)
static int
bridge_intron_gap_16_intron_level (int *bestrL, int *bestrR, int *bestcL, int *bestcR,
				  int *best_introntype,
				  Score16_T **matrixL_upper, Score16_T **matrixL_lower,
				  Score16_T **matrixR_upper, Score16_T **matrixR_lower,
				  Direction16_T **directionsL_upper_nogap, Direction16_T **directionsL_lower_nogap, 
				  Direction16_T **directionsR_upper_nogap, Direction16_T **directionsR_lower_nogap,
				  int *left_known, int *right_known,
				  int rlength, int glengthL, int glengthR,
				  int cdna_direction, bool watsonp, int lbandL, int ubandL, int lbandR, int ubandR,
				  int leftoffset, int rightoffset,
				  Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh,
				  bool jump_late_p) {
  int rL, rR, cL, cR;
  int cloL, chighL;
  int cloR, chighR;
  int bestscore = NEG_INFINITY_16, score, scoreL, scoreR;
  Univcoord_T splicesitepos1, splicesitepos2;
  bool bestp;


  for (rL = 1, rR = rlength-1; rL < rlength; rL++, rR--) {
    debug3(printf("\nGenomic insert: At row %d on left and %d on right\n",rL,rR));
    if ((cloL = rL - lbandL) < 1) {
      cloL = 1;
    }
    if ((chighL = rL + ubandL) > glengthL-1) {
      chighL = glengthL-1;
    }

    if ((cloR = rR - lbandR) < 1) {
      cloR = 1;
    }
    if ((chighR = rR + ubandR) > glengthR-1) {
      chighR = glengthR-1;
    }

    /* Test indels on left and right */
    for (cL = cloL; cL < /* left of main diagonal*/rL; cL++) {
      /* The following check limits genomic inserts (horizontal) and
	 multiple cDNA inserts (vertical). */
      if (left_known[cL] > 0) {
	scoreL = (int) matrixL_lower[rL][cL];
	if (directionsL_lower_nogap[rL][cL] != DIAG) {
	  /* Favor gaps away from intron if possible */
	  scoreL -= 1;
	}

	/* Disallow leftoffset + cL >= rightoffset - cR, or cR >= rightoffset - leftoffset - cL */
	for (cR = cloR; cR < /* left of main diagonal*/rR && cR < rightoffset-leftoffset-cL; cR++) {
	  if (right_known[cR] > 0) {
	    scoreR = (int) matrixR_lower[rR][cR];
	    if (directionsR_lower_nogap[rR][cR] != DIAG) {
	      /* Favor gaps away from intron if possible */
	      scoreR -= 1;
	    }

	    if ((score = scoreL + scoreR) > bestscore ||
		(score >= bestscore && jump_late_p)) {  /* Use >= for jump late */
	      bestp = false;
	      if (watsonp == true) {
		splicesitepos1 = leftoffset + cL;
		splicesitepos2 = rightoffset - cR + 1;
		if (IIT_exists_with_divno_signed(splicing_iit,splicing_divint_crosstable[chrnum],
						 splicesitepos1,splicesitepos2+1U,/*sign*/cdna_direction) == true) {
		  bestp = true;
		}
	      } else {
		splicesitepos1 = (chrhigh - chroffset) - leftoffset - cL + 1;
		splicesitepos2 = (chrhigh - chroffset) - rightoffset + cR;
		if (IIT_exists_with_divno_signed(splicing_iit,splicing_divint_crosstable[chrnum],
						 splicesitepos2,splicesitepos1+1U,/*sign*/-cdna_direction) == true) {
		  bestp = true;
		}
	      }
	      if (bestp == true) {
		debug3(printf("At %d left to %d right, score is (%d)+(%d) = %d (bestscore)\n",
			      cL,cR,scoreL,scoreR,score));
		bestscore = score;
		*bestrL = rL;
		*bestrR = rR;
		*bestcL = cL;
		*bestcR = cR;
	      } else {
		debug3a(printf("At %d left to %d right, score is (%d)+(%d) = %d\n",
			       cL,cR,scoreL,scoreR,score));
	      }
	    }
	  }
	}

	for (/* at main diagonal*/; cR < chighR && cR < rightoffset-leftoffset-cL; cR++) {
	  if (right_known[cR] > 0) {
	    scoreR = (int) matrixR_upper[cR][rR];
	    if (directionsR_upper_nogap[cR][rR] != DIAG) {
	      /* Favor gaps away from intron if possible */
	      scoreR -= 1;
	    }

	    if ((score = scoreL + scoreR) > bestscore ||
		(score >= bestscore && jump_late_p)) {  /* Use >= for jump late */
	      bestp = false;
	      if (watsonp == true) {
		splicesitepos1 = leftoffset + cL;
		splicesitepos2 = rightoffset - cR + 1;
		if (IIT_exists_with_divno_signed(splicing_iit,splicing_divint_crosstable[chrnum],
						 splicesitepos1,splicesitepos2+1U,/*sign*/cdna_direction) == true) {
		  bestp = true;
		}
	      } else {
		splicesitepos1 = (chrhigh - chroffset) - leftoffset - cL + 1;
		splicesitepos2 = (chrhigh - chroffset) - rightoffset + cR;
		if (IIT_exists_with_divno_signed(splicing_iit,splicing_divint_crosstable[chrnum],
						 splicesitepos2,splicesitepos1+1U,/*sign*/-cdna_direction) == true) {
		  bestp = true;
		}
	      }
	      if (bestp == true) {
		debug3(printf("At %d left to %d right, score is (%d)+(%d) = %d (bestscore)\n",
			      cL,cR,scoreL,scoreR,score));
		bestscore = score;
		*bestrL = rL;
		*bestrR = rR;
		*bestcL = cL;
		*bestcR = cR;
	      } else {
		debug3a(printf("At %d left to %d right, score is (%d)+(%d) = %d\n",
			       cL,cR,scoreL,scoreR,score));
	      }
	    }
	  }
	}
      }
    }

    for (/* at main diagonal*/; cL < chighL; cL++) {
      /* The following check limits genomic inserts (horizontal) and
	 multiple cDNA inserts (vertical). */
      if (left_known[cL] > 0) {
	scoreL = (int) matrixL_upper[cL][rL];
	if (directionsL_upper_nogap[cL][rL] != DIAG) {
	  /* Favor gaps away from intron if possible */
	  scoreL -= 1;
	}

	/* Disallow leftoffset + cL >= rightoffset - cR, or cR >= rightoffset - leftoffset - cL */
	for (cR = cloR; cR < /* left of main diagonal*/rR && cR < rightoffset-leftoffset-cL; cR++) {
	  if (right_known[cR] > 0) {
	    scoreR = (int) matrixR_lower[rR][cR];
	    if (directionsR_lower_nogap[rR][cR] != DIAG) {
	      /* Favor gaps away from intron if possible */
	      scoreR -= 1;
	    }

	    if ((score = scoreL + scoreR) > bestscore ||
		(score >= bestscore && jump_late_p)) {  /* Use >= for jump late */
	      bestp = false;
	      if (watsonp == true) {
		splicesitepos1 = leftoffset + cL;
		splicesitepos2 = rightoffset - cR + 1;
		if (IIT_exists_with_divno_signed(splicing_iit,splicing_divint_crosstable[chrnum],
						 splicesitepos1,splicesitepos2+1U,/*sign*/cdna_direction) == true) {
		  bestp = true;
		}
	      } else {
		splicesitepos1 = (chrhigh - chroffset) - leftoffset - cL + 1;
		splicesitepos2 = (chrhigh - chroffset) - rightoffset + cR;
		if (IIT_exists_with_divno_signed(splicing_iit,splicing_divint_crosstable[chrnum],
						 splicesitepos2,splicesitepos1+1U,/*sign*/-cdna_direction) == true) {
		  bestp = true;
		}
	      }
	      if (bestp == true) {
		debug3(printf("At %d left to %d right, score is (%d)+(%d) = %d (bestscore)\n",
			      cL,cR,scoreL,scoreR,score));
		bestscore = score;
		*bestrL = rL;
		*bestrR = rR;
		*bestcL = cL;
		*bestcR = cR;
	      } else {
		debug3a(printf("At %d left to %d right, score is (%d)+(%d) = %d\n",
			       cL,cR,scoreL,scoreR,score));
	      }
	    }
	  }
	}

	for (/* at main diagonal*/; cR < chighR && cR < rightoffset-leftoffset-cL; cR++) {
	  if (right_known[cR] > 0) {
	    scoreR = (int) matrixR_upper[cR][rR];
	    if (directionsR_upper_nogap[cR][rR] != DIAG) {
	      /* Favor gaps away from intron if possible */
	      scoreR -= 1;
	    }

	    if ((score = scoreL + scoreR) > bestscore ||
		(score >= bestscore && jump_late_p)) {  /* Use >= for jump late */
	      bestp = false;
	      if (watsonp == true) {
		splicesitepos1 = leftoffset + cL;
		splicesitepos2 = rightoffset - cR + 1;
		if (IIT_exists_with_divno_signed(splicing_iit,splicing_divint_crosstable[chrnum],
						 splicesitepos1,splicesitepos2+1U,/*sign*/cdna_direction) == true) {
		  bestp = true;
		}
	      } else {
		splicesitepos1 = (chrhigh - chroffset) - leftoffset - cL + 1;
		splicesitepos2 = (chrhigh - chroffset) - rightoffset + cR;
		if (IIT_exists_with_divno_signed(splicing_iit,splicing_divint_crosstable[chrnum],
						 splicesitepos2,splicesitepos1+1U,/*sign*/-cdna_direction) == true) {
		  bestp = true;
		}
	      }
	      if (bestp == true) {
		debug3(printf("At %d left to %d right, score is (%d)+(%d) = %d (bestscore)\n",
			      cL,cR,scoreL,scoreR,score));
		bestscore = score;
		*bestrL = rL;
		*bestrR = rR;
		*bestcL = cL;
		*bestcR = cR;
	      } else {
		debug3a(printf("At %d left to %d right, score is (%d)+(%d) = %d\n",
			       cL,cR,scoreL,scoreR,score));
	      }
	    }
	  }
	}
      }
    }
  }

  *best_introntype = NONINTRON;
  return bestscore;
}


/* Returns finalscore */
static int
bridge_intron_gap_16_site_level (int *bestrL, int *bestrR, int *bestcL, int *bestcR,
				Score16_T **matrixL_upper, Score16_T **matrixL_lower,
				Score16_T **matrixR_upper, Score16_T **matrixR_lower,
				Direction16_T **directionsL_upper_nogap, Direction16_T **directionsL_lower_nogap, 
				Direction16_T **directionsR_upper_nogap, Direction16_T **directionsR_lower_nogap,
				char *gsequenceL, char *gsequenceL_alt, char *rev_gsequenceR, char *rev_gsequenceR_alt,
				int goffsetL, int rev_goffsetR, int *left_known, int *right_known,
				int rlength, int glengthL, int glengthR,
				int cdna_direction, bool watsonp, int lbandL, int ubandL, int lbandR, int ubandR,
				int canonical_reward, int leftoffset, int rightoffset,
				Univcoord_T chroffset, Univcoord_T chrhigh,
				bool halfp, bool finalp) {
  int rL, rR, cL, cR;
  int bestrL_with_prob, bestrR_with_prob, bestcL_with_prob, bestcR_with_prob;
  int cloL, chighL;
  int cloR, chighR;
  int introntype;
  int bestscore = NEG_INFINITY_16, score, scoreL, scoreR, scoreI;
  int bestscore_with_prob = NEG_INFINITY_16;
  double *left_probabilities, *right_probabilities, probL, probR, probL_trunc, probR_trunc, bestprob, bestprob_trunc;
  Univcoord_T splicesitepos;
  char left1, left2, right2, right1, left1_alt, left2_alt, right2_alt, right1_alt;
  int *leftdi, *rightdi;
  bool use_prob_p;


  /* Read dinucleotides */
  leftdi = (int *) MALLOCA((glengthL+1) * sizeof(int));
  rightdi = (int *) MALLOCA((glengthR+1) * sizeof(int));

  for (cL = 0; cL < glengthL - 1; cL++) {
    left1 = gsequenceL[cL];
    left1_alt = gsequenceL_alt[cL];
    left2 = gsequenceL[cL+1];
    left2_alt = gsequenceL_alt[cL+1];
    assert(left1 == get_genomic_nt(&left1_alt,goffsetL+cL,chroffset,chrhigh,watsonp));
    assert(left2 == get_genomic_nt(&left2_alt,goffsetL+cL+1,chroffset,chrhigh,watsonp));

    if ((left1 == 'G' || left1_alt == 'G') && (left2 == 'T' || left2_alt == 'T')) {
      leftdi[cL] = LEFT_GT;
    } else if ((left1 == 'G' || left1_alt == 'G') && (left2 == 'C' || left2_alt == 'C')) {
      leftdi[cL] = LEFT_GC;
    } else if ((left1 == 'A' || left1_alt == 'A') && (left2 == 'T' || left2_alt == 'T')) {
      leftdi[cL] = LEFT_AT;
#ifndef PMAP
    } else if ((left1 == 'C' || left1_alt == 'C') && (left2 == 'T' || left2_alt == 'T')) {
      leftdi[cL] = LEFT_CT;
#endif
    } else {
      leftdi[cL] = 0x00;
    }
  }
  leftdi[glengthL-1] = leftdi[glengthL] = 0x00;

  for (cR = 0; cR < glengthR - 1; cR++) {
    right2 = rev_gsequenceR[-cR-1];
    right2_alt = rev_gsequenceR_alt[-cR-1];
    right1 = rev_gsequenceR[-cR];
    right1_alt = rev_gsequenceR_alt[-cR];
    assert(right2 == get_genomic_nt(&right2_alt,rev_goffsetR-cR-1,chroffset,chrhigh,watsonp));
    assert(right1 == get_genomic_nt(&right1_alt,rev_goffsetR-cR,chroffset,chrhigh,watsonp));

    if ((right2 == 'A' || right2_alt == 'A') && (right1 == 'G' || right1_alt == 'G')) {
      rightdi[cR] = RIGHT_AG;
    } else if ((right2 == 'A' || right2_alt == 'A') && (right1 == 'C' || right1_alt == 'C')) {
      rightdi[cR] = RIGHT_AC;
#ifndef PMAP
    } else if ((right2 == 'G' || right2_alt == 'G') && (right1 == 'C' || right1_alt == 'C')) {
      rightdi[cR] = RIGHT_GC;
    } else if ((right2 == 'A' || right2_alt == 'A') && (right1 == 'T' || right1_alt == 'T')) {
      rightdi[cR] = RIGHT_AT;
#endif
    } else {
      rightdi[cR] = 0x00;
    }
  }
  rightdi[glengthR-1] = rightdi[glengthR] = 0x00;


  left_probabilities = (double *) MALLOCA(glengthL * sizeof(double));
  right_probabilities = (double *) MALLOCA(glengthR * sizeof(double));

  debug3(printf("watsonp is %d.  cdna_direction is %d\n",watsonp,cdna_direction));
  if (watsonp == true) {
    if (cdna_direction > 0) {
      for (cL = 0; cL < glengthL - 1; cL++) {
	splicesitepos = chroffset + leftoffset + cL;
	if (left_known[cL]) {
	  left_probabilities[cL] = 1.0;
	} else {
	  left_probabilities[cL] = Maxent_hr_donor_prob(splicesitepos,chroffset);
	  debug3(printf("left donor probability at cL %d is %f\n",cL,left_probabilities[cL]));
	}
      }

      for (cR = 0; cR < glengthR - 1; cR++) {
	splicesitepos = chroffset + rightoffset - cR + 1;
	if (right_known[cR]) {
	  right_probabilities[cR] = 1.0;
	} else {
	  right_probabilities[cR] = Maxent_hr_acceptor_prob(splicesitepos,chroffset);
	  debug3(printf("right acceptor probability at cR %d is %f\n",cR,right_probabilities[cR]));
	}
      }

    } else {
      for (cL = 0; cL < glengthL - 1; cL++) {
	splicesitepos = chroffset + leftoffset + cL;
	if (left_known[cL]) {
	  left_probabilities[cL] = 1.0;
	} else {
	  left_probabilities[cL] = Maxent_hr_antiacceptor_prob(splicesitepos,chroffset);
	  debug3(printf("left antiacceptor probability at cL %d is %f\n",cL,left_probabilities[cL]));
	}
      }

      for (cR = 0; cR < glengthR - 1; cR++) {
	splicesitepos = chroffset + rightoffset - cR + 1;
	if (right_known[cR]) {
	  right_probabilities[cR] = 1.0;
	} else {
	  right_probabilities[cR] = Maxent_hr_antidonor_prob(splicesitepos,chroffset);
	  debug3(printf("right antidonor probability at cR %d is %f\n",cR,right_probabilities[cR]));
	}
      }
    }

  } else {
    if (cdna_direction > 0) {
      for (cL = 0; cL < glengthL - 1; cL++) {
	splicesitepos = chrhigh - leftoffset - cL + 1;
	if (left_known[cL]) {
	  left_probabilities[cL] = 1.0;
	} else {
	  left_probabilities[cL] = Maxent_hr_antidonor_prob(splicesitepos,chroffset);
	  debug3(printf("left antidonor probability at cL %d is %f\n",cL,left_probabilities[cL]));
	}
      }

      for (cR = 0; cR < glengthR - 1; cR++) {
	splicesitepos = chrhigh - rightoffset + cR;
	if (right_known[cR]) {
	  right_probabilities[cR] = 1.0;
	} else {
	  right_probabilities[cR] = Maxent_hr_antiacceptor_prob(splicesitepos,chroffset);
	  debug3(printf("right antiacceptor probability at cR %d is %f\n",cR,right_probabilities[cR]));
	}
      }

    } else {
      for (cL = 0; cL < glengthL - 1; cL++) {
	splicesitepos = chrhigh - leftoffset - cL + 1;
	if (left_known[cL]) {
	  left_probabilities[cL] = 1.0;
	} else {
	  left_probabilities[cL] = Maxent_hr_acceptor_prob(splicesitepos,chroffset);
	  debug3(printf("left acceptor probability at cL %d is %f\n",cL,left_probabilities[cL]));
	}
      }

      for (cR = 0; cR < glengthR - 1; cR++) {
	splicesitepos = chrhigh - rightoffset + cR;
	if (right_known[cR]) {
	  right_probabilities[cR] = 1.0;
	} else {
	  right_probabilities[cR] = Maxent_hr_donor_prob(splicesitepos,chroffset);
	  debug3(printf("right donor probability at cR %d is %f\n",cR,right_probabilities[cR]));
	}
      }
    }
  }

  /* Search using probs and without simultaneously */
  bestscore = NEG_INFINITY_16;
  bestprob = bestprob_trunc = 0.0;
  for (rL = 1, rR = rlength-1; rL < rlength; rL++, rR--) {
    debug3(printf("\nAt row %d on left and %d on right\n",rL,rR));
    if ((cloL = rL - lbandL) < 1) {
      cloL = 1;
    }
    if ((chighL = rL + ubandL) > glengthL-1) {
      chighL = glengthL-1;
    }

    if ((cloR = rR - lbandR) < 1) {
      cloR = 1;
    }
    if ((chighR = rR + ubandR) > glengthR-1) {
      chighR = glengthR-1;
    }

    debug3(printf("A. Test no indels\n"));
    cL = rL;
    probL = left_probabilities[cL];
    if (probL > PROB_CEILING) {
      probL_trunc = PROB_CEILING;
    } else if (probL < PROB_FLOOR) {
      probL_trunc = PROB_FLOOR;
    } else {
      probL_trunc = probL;
    }
    scoreL = (int) matrixL_upper[cL][rL];

    cR = rR;
    probR = right_probabilities[cR];
    if (probR > PROB_CEILING) {
      probR_trunc = PROB_CEILING;
    } else if (probR < PROB_FLOOR) {
      probR_trunc = PROB_FLOOR;
    } else {
      probR_trunc = probR;
    }
    scoreR = (int) matrixR_upper[cR][rR];

#ifdef USE_SCOREI
    scoreI = intron_score(&introntype,leftdi[cL],rightdi[cR],cdna_direction,canonical_reward,finalp);
#else
    scoreI = 0;
#endif
	
    if ((score = scoreL + scoreI + scoreR) > bestscore) {
      debug3(printf("Best score: At %d left to %d right, score is (%d)+(%d)+(%d) = %d (bestscore, prob %f + %f)\n",
		    cL,cR,scoreL,scoreI,scoreR,scoreL+scoreI+scoreR,probL,probR));
      debug3(printf("probL %f, probR %f\n",left_probabilities[cL],right_probabilities[cR]));
      bestscore = score;
      *bestrL = rL;
      *bestrR = rR;
      *bestcL = cL;
      *bestcR = cR;
      bestprob = probL + probR;
    } else if (score == bestscore && probL + probR > bestprob) {
      debug3(printf("Improved prob: At %d left to %d right, score is (%d)+(%d)+(%d) = %d (bestscore, prob %f + %f)\n",
		    cL,cR,scoreL,scoreI,scoreR,scoreL+scoreI+scoreR,probL,probR));
      debug3(printf("probL %f, probR %f\n",left_probabilities[cL],right_probabilities[cR]));
      *bestrL = rL;
      *bestrR = rR;
      *bestcL = cL;
      *bestcR = cR;
      bestprob = probL + probR;
    } else {
      debug3a(printf("Not best score: At %d left to %d right, score is (%d)+(%d)+(%d) = %d (bestscore, prob %f + %f)\n",
		     cL,cR,scoreL,scoreI,scoreR,scoreL+scoreI+scoreR,probL,probR));
    }
	
    if (probL_trunc + probR_trunc < bestprob_trunc) {
      debug3a(printf("At %d left to %d right, prob is %f + %f = %f\n",
		     cL,cR,probL_trunc,probR_trunc,probL_trunc+probR_trunc));
	  
    } else if (probL_trunc + probR_trunc == bestprob_trunc) {
      debug3(printf("At %d left to %d right, prob is %f + %f = %f\n",
		    cL,cR,probL_trunc,probR_trunc,probL_trunc+probR_trunc));
	  
      if (scoreL + scoreI + scoreR > bestscore_with_prob) {
	debug3(printf(" (bestscore %d)\n",scoreL+scoreR));
	bestprob_trunc = probL_trunc + probR_trunc;
	bestcL_with_prob = cL;
	bestcR_with_prob = cR;
	bestrL_with_prob = rL;
	bestrR_with_prob = rR;
	bestscore_with_prob = scoreL + scoreI + scoreR;
      }
		  
    } else {
      /* probL_trunc + probR_trunc > bestprob_trunc */
      debug3(printf("At %d left to %d right, prob is %f + %f = %f\n",
		    cL,cR,probL_trunc,probR_trunc,probL_trunc+probR_trunc));
	  
      debug3(printf(" (bestscore %d)\n",scoreL+scoreR));
      bestprob_trunc = probL_trunc + probR_trunc;
      bestcL_with_prob = cL;
      bestcR_with_prob = cR;
      bestrL_with_prob = rL;
      bestrR_with_prob = rR;
      bestscore_with_prob = scoreL + scoreI + scoreR;
    }


    debug3(printf("B. Test indel on right\n"));
    /* Test indel on right */
    cL = rL;
    probL = left_probabilities[cL];
    if (probL > PROB_CEILING) {
      probL_trunc = PROB_CEILING;
    } else if (probL < PROB_FLOOR) {
      probL_trunc = PROB_FLOOR;
    } else {
      probL_trunc = probL;
    }
    scoreL = (int) matrixL_upper[cL][rL];
    if (directionsL_upper_nogap[cL][rL] != DIAG) {
      /* Favor gaps away from intron if possible */
      scoreL -= 1;
    }

    /* Disallow leftoffset + cL >= rightoffset - cR, or cR >= rightoffset - leftoffset - cL */
    for (cR = cloR; cR < /*to main diagonal*/rR && cR < rightoffset-leftoffset-cL; cR++) {
      probR = right_probabilities[cR];
      if (probR > PROB_CEILING) {
	probR_trunc = PROB_CEILING;
      } else if (probR < PROB_FLOOR) {
	probR_trunc = PROB_FLOOR;
      } else {
	probR_trunc = probR;
      }
      scoreR = (int) matrixR_lower[rR][cR];
      if (directionsR_lower_nogap[rR][cR] != DIAG) {
	/* Favor gaps away from intron if possible */
	scoreR -= 1;
      }
	
#ifdef USE_SCOREI
      scoreI = intron_score(&introntype,leftdi[cL],rightdi[cR],cdna_direction,canonical_reward,finalp);
#else
      scoreI = 0;
#endif
	
      if ((score = scoreL + scoreI + scoreR) > bestscore) {
	debug3(printf("Best score: At %d left to %d right, score is (%d)+(%d)+(%d) = %d (bestscore, prob %f + %f)\n",
		      cL,cR,scoreL,scoreI,scoreR,scoreL+scoreI+scoreR,probL,probR));
	debug3(printf("probL %f, probR %f\n",left_probabilities[cL],right_probabilities[cR]));
	bestscore = score;
	*bestrL = rL;
	*bestrR = rR;
	*bestcL = cL;
	*bestcR = cR;
	bestprob = probL + probR;
      } else if (score == bestscore && probL + probR > bestprob) {
	debug3(printf("Improved prob: At %d left to %d right, score is (%d)+(%d)+(%d) = %d (bestscore, prob %f + %f)\n",
		      cL,cR,scoreL,scoreI,scoreR,scoreL+scoreI+scoreR,probL,probR));
	debug3(printf("probL %f, probR %f\n",left_probabilities[cL],right_probabilities[cR]));
	*bestrL = rL;
	*bestrR = rR;
	*bestcL = cL;
	*bestcR = cR;
	bestprob = probL + probR;
      } else {
	debug3a(printf("Not best score: At %d left to %d right, score is (%d)+(%d)+(%d) = %d (bestscore, prob %f + %f)\n",
		       cL,cR,scoreL,scoreI,scoreR,scoreL+scoreI+scoreR,probL,probR));
      }
	
      if (probL_trunc + probR_trunc < bestprob_trunc) {
	debug3a(printf("At %d left to %d right, prob is %f + %f = %f\n",
		       cL,cR,probL_trunc,probR_trunc,probL_trunc+probR_trunc));
	  
      } else if (probL_trunc + probR_trunc == bestprob_trunc) {
	debug3(printf("At %d left to %d right, prob is %f + %f = %f\n",
		      cL,cR,probL_trunc,probR_trunc,probL_trunc+probR_trunc));
	  
	if (scoreL + scoreI + scoreR > bestscore_with_prob) {
	  debug3(printf(" (bestscore %d)\n",scoreL+scoreR));
	  bestprob_trunc = probL_trunc + probR_trunc;
	  bestcL_with_prob = cL;
	  bestcR_with_prob = cR;
	  bestrL_with_prob = rL;
	  bestrR_with_prob = rR;
	  bestscore_with_prob = scoreL + scoreI + scoreR;
	}
		  
      } else {
	/* probL_trunc + probR_trunc > bestprob_trunc */
	debug3(printf("At %d left to %d right, prob is %f + %f = %f\n",
		      cL,cR,probL_trunc,probR_trunc,probL_trunc+probR_trunc));
	  
	debug3(printf(" (bestscore %d)\n",scoreL+scoreR));
	bestprob_trunc = probL_trunc + probR_trunc;
	bestcL_with_prob = cL;
	bestcR_with_prob = cR;
	bestrL_with_prob = rL;
	bestrR_with_prob = rR;
	bestscore_with_prob = scoreL + scoreI + scoreR;
      }
    }

    debug3(printf("Skip main diagonal\n"));
    for (/*Skip main diagonal*/cR++; cR < chighR && cR < rightoffset-leftoffset-cL; cR++) {
      probR = right_probabilities[cR];
      if (probR > PROB_CEILING) {
	probR_trunc = PROB_CEILING;
      } else if (probR < PROB_FLOOR) {
	probR_trunc = PROB_FLOOR;
      } else {
	probR_trunc = probR;
      }
      scoreR = (int) matrixR_upper[cR][rR];
      if (directionsR_upper_nogap[cR][rR] != DIAG) {
	/* Favor gaps away from intron if possible */
	scoreR -= 1;
      }
	
#ifdef USE_SCOREI
      scoreI = intron_score(&introntype,leftdi[cL],rightdi[cR],cdna_direction,canonical_reward,finalp);
#else
      scoreI = 0;
#endif
	
      if ((score = scoreL + scoreI + scoreR) > bestscore) {
	debug3(printf("Best score: At %d left to %d right, score is (%d)+(%d)+(%d) = %d (bestscore, prob %f + %f)\n",
		      cL,cR,scoreL,scoreI,scoreR,scoreL+scoreI+scoreR,probL,probR));
	debug3(printf("probL %f, probR %f\n",left_probabilities[cL],right_probabilities[cR]));
	bestscore = score;
	*bestrL = rL;
	*bestrR = rR;
	*bestcL = cL;
	*bestcR = cR;
	bestprob = probL + probR;
      } else if (score == bestscore && probL + probR > bestprob) {
	debug3(printf("Improved prob: At %d left to %d right, score is (%d)+(%d)+(%d) = %d (bestscore, prob %f + %f)\n",
		      cL,cR,scoreL,scoreI,scoreR,scoreL+scoreI+scoreR,probL,probR));
	debug3(printf("probL %f, probR %f\n",left_probabilities[cL],right_probabilities[cR]));
	*bestrL = rL;
	*bestrR = rR;
	*bestcL = cL;
	*bestcR = cR;
	bestprob = probL + probR;
      } else {
	debug3a(printf("Not best score: At %d left to %d right, score is (%d)+(%d)+(%d) = %d (bestscore, prob %f + %f)\n",
		       cL,cR,scoreL,scoreI,scoreR,scoreL+scoreI+scoreR,probL,probR));
      }
	
      if (probL_trunc + probR_trunc < bestprob_trunc) {
	debug3a(printf("At %d left to %d right, prob is %f + %f = %f\n",
		       cL,cR,probL_trunc,probR_trunc,probL_trunc+probR_trunc));
	  
      } else if (probL_trunc + probR_trunc == bestprob_trunc) {
	debug3(printf("At %d left to %d right, prob is %f + %f = %f\n",
		      cL,cR,probL_trunc,probR_trunc,probL_trunc+probR_trunc));
	  
	if (scoreL + scoreI + scoreR > bestscore_with_prob) {
	  debug3(printf(" (bestscore %d)\n",scoreL+scoreR));
	  bestprob_trunc = probL_trunc + probR_trunc;
	  bestcL_with_prob = cL;
	  bestcR_with_prob = cR;
	  bestrL_with_prob = rL;
	  bestrR_with_prob = rR;
	  bestscore_with_prob = scoreL + scoreI + scoreR;
	}
		  
      } else {
	/* probL_trunc + probR_trunc > bestprob_trunc */
	debug3(printf("At %d left to %d right, prob is %f + %f = %f\n",
		      cL,cR,probL_trunc,probR_trunc,probL_trunc+probR_trunc));
	  
	debug3(printf(" (bestscore %d)\n",scoreL+scoreR));
	bestprob_trunc = probL_trunc + probR_trunc;
	bestcL_with_prob = cL;
	bestcR_with_prob = cR;
	bestrL_with_prob = rL;
	bestrR_with_prob = rR;
	bestscore_with_prob = scoreL + scoreI + scoreR;
      }
    }


    debug3(printf("C. Test indel on left\n"));
    /* Test indel on left */
    cR = rR;
    probR = right_probabilities[cR];
    if (probR > PROB_CEILING) {
      probR_trunc = PROB_CEILING;
    } else if (probR < PROB_FLOOR) {
      probR_trunc = PROB_FLOOR;
    } else {
      probR_trunc = probR;
    }
    scoreR = (int) matrixR_upper[cR][rR];
    if (directionsR_upper_nogap[cR][rR] != DIAG) {
      /* Favor gaps away from intron if possible */
      scoreR -= 1;
    }

    /* Disallow leftoffset + cL >= rightoffset - cR, or cR >= rightoffset - leftoffset - cL */
    for (cL = cloL; cL < /*to main diagonal*/rL && cL < rightoffset-leftoffset-cR; cL++) {
      probL = left_probabilities[cL];
      if (probL > PROB_CEILING) {
	probL_trunc = PROB_CEILING;
      } else if (probL < PROB_FLOOR) {
	probL_trunc = PROB_FLOOR;
      } else {
	probL_trunc = probL;
      }
      scoreL = (int) matrixL_lower[rL][cL];
      if (directionsL_lower_nogap[rL][cL] != DIAG) {
	/* Favor gaps away from intron if possible */
	scoreL -= 1;
      }

#ifdef USE_SCOREI
      scoreI = intron_score(&introntype,leftdi[cL],rightdi[cR],cdna_direction,canonical_reward,finalp);
#else
      scoreI = 0;
#endif

      if ((score = scoreL + scoreI + scoreR) > bestscore) {
	debug3(printf("Best score: At %d left to %d right, score is (%d)+(%d)+(%d) = %d (bestscore, prob %f + %f)\n",
		      cL,cR,scoreL,scoreI,scoreR,scoreL+scoreI+scoreR,probL,probR));
	debug3(printf("probL %f, probR %f\n",left_probabilities[cL],right_probabilities[cR]));
	bestscore = score;
	*bestrL = rL;
	*bestrR = rR;
	*bestcL = cL;
	*bestcR = cR;
	bestprob = probL + probR;
      } else if (score == bestscore && probL + probR > bestprob) {
	debug3(printf("Improved prob: At %d left to %d right, score is (%d)+(%d)+(%d) = %d (bestscore, prob %f + %f)\n",
		      cL,cR,scoreL,scoreI,scoreR,scoreL+scoreI+scoreR,probL,probR));
	debug3(printf("probL %f, probR %f\n",left_probabilities[cL],right_probabilities[cR]));
	*bestrL = rL;
	*bestrR = rR;
	*bestcL = cL;
	*bestcR = cR;
	bestprob = probL + probR;
      } else {
	debug3a(printf("Not best score: At %d left to %d right, score is (%d)+(%d)+(%d) = %d (bestscore, prob %f + %f)\n",
		       cL,cR,scoreL,scoreI,scoreR,scoreL+scoreI+scoreR,probL,probR));
      }
	
      if (probL_trunc + probR_trunc < bestprob_trunc) {
	debug3a(printf("At %d left to %d right, prob is %f + %f = %f\n",
		       cL,cR,probL_trunc,probR_trunc,probL_trunc+probR_trunc));
	  
      } else if (probL_trunc + probR_trunc == bestprob_trunc) {
	debug3(printf("At %d left to %d right, prob is %f + %f = %f\n",
		      cL,cR,probL_trunc,probR_trunc,probL_trunc+probR_trunc));
	  
	if (scoreL + scoreI + scoreR > bestscore_with_prob) {
	  debug3(printf(" (bestscore %d)\n",scoreL+scoreR));
	  bestprob_trunc = probL_trunc + probR_trunc;
	  bestcL_with_prob = cL;
	  bestcR_with_prob = cR;
	  bestrL_with_prob = rL;
	  bestrR_with_prob = rR;
	  bestscore_with_prob = scoreL + scoreI + scoreR;
	}
	  
      } else {
	/* probL_trunc + probR_trunc > bestprob_trunc */
	debug3(printf("At %d left to %d right, prob is %f + %f = %f\n",
		      cL,cR,probL_trunc,probR_trunc,probL_trunc+probR_trunc));
	  
	debug3(printf(" (bestscore %d)\n",scoreL+scoreR));
	bestprob_trunc = probL_trunc + probR_trunc;
	bestcL_with_prob = cL;
	bestcR_with_prob = cR;
	bestrL_with_prob = rL;
	bestrR_with_prob = rR;
	bestscore_with_prob = scoreL + scoreI + scoreR;
      }
    }

    debug3(printf("Skip main diagonal\n"));
    for (/*Skip main diagonal*/cL++; cL < chighL && cL < rightoffset-leftoffset-cR; cL++) {
      probL = left_probabilities[cL];
      if (probL > PROB_CEILING) {
	probL_trunc = PROB_CEILING;
      } else if (probL < PROB_FLOOR) {
	probL_trunc = PROB_FLOOR;
      } else {
	probL_trunc = probL;
      }
      scoreL = (int) matrixL_upper[cL][rL];
      if (directionsL_upper_nogap[cL][rL] != DIAG) {
	/* Favor gaps away from intron if possible */
	scoreL -= 1;
      }

#ifdef USE_SCOREI
      scoreI = intron_score(&introntype,leftdi[cL],rightdi[cR],cdna_direction,canonical_reward,finalp);
#else
      scoreI = 0;
#endif

      if ((score = scoreL + scoreI + scoreR) > bestscore) {
	debug3(printf("Best score: At %d left to %d right, score is (%d)+(%d)+(%d) = %d (bestscore, prob %f + %f)\n",
		      cL,cR,scoreL,scoreI,scoreR,scoreL+scoreI+scoreR,probL,probR));
	debug3(printf("probL %f, probR %f\n",left_probabilities[cL],right_probabilities[cR]));
	bestscore = score;
	*bestrL = rL;
	*bestrR = rR;
	*bestcL = cL;
	*bestcR = cR;
	bestprob = probL + probR;
      } else if (score == bestscore && probL + probR > bestprob) {
	debug3(printf("Improved prob: At %d left to %d right, score is (%d)+(%d)+(%d) = %d (bestscore, prob %f + %f)\n",
		      cL,cR,scoreL,scoreI,scoreR,scoreL+scoreI+scoreR,probL,probR));
	debug3(printf("probL %f, probR %f\n",left_probabilities[cL],right_probabilities[cR]));
	*bestrL = rL;
	*bestrR = rR;
	*bestcL = cL;
	*bestcR = cR;
	bestprob = probL + probR;
      } else {
	debug3a(printf("Not best score: At %d left to %d right, score is (%d)+(%d)+(%d) = %d (bestscore, prob %f + %f)\n",
		       cL,cR,scoreL,scoreI,scoreR,scoreL+scoreI+scoreR,probL,probR));
      }
	
      if (probL_trunc + probR_trunc < bestprob_trunc) {
	debug3a(printf("At %d left to %d right, prob is %f + %f = %f\n",
		       cL,cR,probL_trunc,probR_trunc,probL_trunc+probR_trunc));
	  
      } else if (probL_trunc + probR_trunc == bestprob_trunc) {
	debug3(printf("At %d left to %d right, prob is %f + %f = %f\n",
		      cL,cR,probL_trunc,probR_trunc,probL_trunc+probR_trunc));
	  
	if (scoreL + scoreI + scoreR > bestscore_with_prob) {
	  debug3(printf(" (bestscore %d)\n",scoreL+scoreR));
	  bestprob_trunc = probL_trunc + probR_trunc;
	  bestcL_with_prob = cL;
	  bestcR_with_prob = cR;
	  bestrL_with_prob = rL;
	  bestrR_with_prob = rR;
	  bestscore_with_prob = scoreL + scoreI + scoreR;
	}
	  
      } else {
	/* probL_trunc + probR_trunc > bestprob_trunc */
	debug3(printf("At %d left to %d right, prob is %f + %f = %f\n",
		      cL,cR,probL_trunc,probR_trunc,probL_trunc+probR_trunc));
	  
	debug3(printf(" (bestscore %d)\n",scoreL+scoreR));
	bestprob_trunc = probL_trunc + probR_trunc;
	bestcL_with_prob = cL;
	bestcR_with_prob = cR;
	bestrL_with_prob = rL;
	bestrR_with_prob = rR;
	bestscore_with_prob = scoreL + scoreI + scoreR;
      }
    }
  }

  if (bestprob > 2*PROB_CEILING) {
    /* Probability is good with best alignment, so take that */
    debug(printf("Best alignment has good probability\n"));
    use_prob_p = true;
  } else if (left_probabilities[bestcL_with_prob] < PROB_CEILING && right_probabilities[bestcR_with_prob] < PROB_CEILING) {
    /* Probability-based solution is bad, so use alignment */
    debug(printf("Probability-based solution is bad\n"));
    use_prob_p = false;
  } else if (bestscore_with_prob < bestscore - 9) {
    debug(printf("Probability-based solution requires very bad alignment\n"));
    use_prob_p = false;
  } else {
    use_prob_p = true;
  }

  if (use_prob_p == true) {
    /* Best alignment yields bad probability, and probability-based alignment yields good probability, so switch */
    debug(printf("Switch to probability-based solution\n"));
    debug(printf("SIMD 16. bestscore %d (bestprob %f) vs bestscore_with_prob %d (bestprob_trunc %f, actually %f and %f)\n",
		 bestscore,bestprob,bestscore_with_prob,bestprob_trunc,left_probabilities[bestcL_with_prob],right_probabilities[bestcR_with_prob]));
    *bestcL = bestcL_with_prob;
    *bestcR = bestcR_with_prob;
    *bestrL = bestrL_with_prob;
    *bestrR = bestrR_with_prob;
    bestscore = bestscore_with_prob;
  }
    
  FREEA(rightdi);
  FREEA(leftdi);
  FREEA(left_probabilities);
  FREEA(right_probabilities);

  if (bestscore < 0) {
    return bestscore;
  } else if (halfp == true) {
    scoreI = intron_score(&introntype,leftdi[*bestcL],rightdi[*bestcR],cdna_direction,canonical_reward,finalp);
    return (bestscore - scoreI/2);
  } else {
    return bestscore;
  }
}



static int
bridge_intron_gap_16_ud (int *bestrL, int *bestrR, int *bestcL, int *bestcR,
			 int *best_introntype, double *left_prob, double *right_prob,
			 Score16_T **matrixL_upper, Score16_T **matrixL_lower,
			 Score16_T **matrixR_upper, Score16_T **matrixR_lower,
			 Direction16_T **directionsL_upper_nogap, Direction16_T **directionsL_lower_nogap, 
			 Direction16_T **directionsR_upper_nogap, Direction16_T **directionsR_lower_nogap,
			 char *gsequenceL, char *gsequenceL_alt, char *rev_gsequenceR, char *rev_gsequenceR_alt,
			 int goffsetL, int rev_goffsetR, int rlength, int glengthL, int glengthR,
			 int cdna_direction, bool watsonp, int lbandL, int ubandL, int lbandR, int ubandR,
			 double defect_rate, int canonical_reward, int leftoffset, int rightoffset,
			 Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh,
			 bool halfp, bool finalp, bool jump_late_p) {
  int finalscore;
  int *left_known, *right_known;

  debug(printf("Running bridge_intron_gap_16_ud\n"));

  if (glengthL+1 <= 0) {
    fprintf(stderr,"Problem with glengthL = %d\n",glengthL);
    abort();
  }

  if (glengthR+1 <= 0) {
    fprintf(stderr,"Problem with glengthR = %d\n",glengthR);
    abort();
  }

  left_known = (int *) CALLOCA(glengthL+1,sizeof(int));
  right_known = (int *) CALLOCA(glengthR+1,sizeof(int));
  get_known_splicesites(left_known,right_known,glengthL,glengthR,
			/*leftoffset*/goffsetL,/*rightoffset*/rev_goffsetR,
			cdna_direction,watsonp,chrnum,chroffset,chrhigh);

  if (novelsplicingp == false && splicing_iit != NULL && (donor_typeint < 0 || acceptor_typeint < 0)) {
    /* Constrain to given introns */
    finalscore = bridge_intron_gap_16_intron_level(&(*bestrL),&(*bestrR),&(*bestcL),&(*bestcR),&(*best_introntype),
						  matrixL_upper,matrixL_lower,matrixR_upper,matrixR_lower,
						  directionsL_upper_nogap,directionsL_lower_nogap, 
						  directionsR_upper_nogap,directionsR_lower_nogap,
						  left_known,right_known,rlength,glengthL,glengthR,
						  cdna_direction,watsonp,lbandL,ubandL,lbandR,ubandR,
						  leftoffset,rightoffset,chrnum,chroffset,chrhigh,jump_late_p);

  } else {
    finalscore = bridge_intron_gap_16_site_level(&(*bestrL),&(*bestrR),&(*bestcL),&(*bestcR),
						matrixL_upper,matrixL_lower,matrixR_upper,matrixR_lower,
						directionsL_upper_nogap,directionsL_lower_nogap, 
						directionsR_upper_nogap,directionsR_lower_nogap,
						gsequenceL,gsequenceL_alt,rev_gsequenceR,rev_gsequenceR_alt,goffsetL,rev_goffsetR,
						left_known,right_known,rlength,glengthL,glengthR,
						cdna_direction,watsonp,lbandL,ubandL,lbandR,ubandR,
						canonical_reward,leftoffset,rightoffset,
						chroffset,chrhigh,halfp,finalp);
  }


#if 0
  /* Determine if result meets given constraints */
  if (*finalscore < 0) {
    result = false;
  } else if (novelsplicingp == true) {
    result = true;
  } else if (splicing_iit == NULL) {
    result = true;
  } else if (donor_typeint >= 0 && acceptor_typeint >= 0) {
    /* If novelsplicingp is false and using splicing at splice site level, require sites to be known */
    if (left_known[*bestcL] == 0 || right_known[*bestcR] == 0) {
      debug(printf("Novel splicing not allowed, so bridge_intron_gap returning false\n"));
      result = false;
    } else {
      result = true;
    }
  } else {
    /* If novelsplicingp is false and using splicing at splice site level, result was already constrained */
    result = true;
  }
#endif

  if (finalscore >= 0) {
    get_splicesite_probs(&(*left_prob),&(*right_prob),*bestcL,*bestcR,
			 left_known,right_known,leftoffset,rightoffset,chroffset,chrhigh,
			 cdna_direction,watsonp);
  }

  FREEA(right_known);
  FREEA(left_known);

#if defined(DEBUG) || defined(DEBUG3)
  printf("Returning final score of %d at (%d,%d) left to (%d,%d) right, with probs %f and %f\n",
	 finalscore,*bestrL,*bestcL,*bestrR,*bestcR,*left_prob,*right_prob);
#endif

  return finalscore;
}
#endif


#ifndef HAVE_SSE2
static int
bridge_intron_gap_intron_level (int *bestrL, int *bestrR, int *bestcL, int *bestcR,
				int *best_introntype,
				Score32_T **matrixL, Score32_T **matrixR,
				Direction32_T **directionsL_nogap, Direction32_T **directionsR_nogap,
				int *left_known, int *right_known,
				int rlength, int glengthL, int glengthR,
				int cdna_direction, bool watsonp, int lbandL, int ubandL, int lbandR, int ubandR,
				int leftoffset, int rightoffset,
				Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh,
				bool jump_late_p) {
  int rL, rR, cL, cR;
  int cloL, chighL;
  int cloR, chighR;
  int bestscore = NEG_INFINITY_32, score, scoreL, scoreR;
  Univcoord_T splicesitepos1, splicesitepos2;
  bool bestp;


  for (rL = 1, rR = rlength-1; rL < rlength; rL++, rR--) {
    debug3(printf("\nGenomic insert: At row %d on left and %d on right\n",rL,rR));
    if ((cloL = rL - lbandL) < 1) {
      cloL = 1;
    }
    if ((chighL = rL + ubandL) > glengthL-1) {
      chighL = glengthL-1;
    }

    if ((cloR = rR - lbandR) < 1) {
      cloR = 1;
    }
    if ((chighR = rR + ubandR) > glengthR-1) {
      chighR = glengthR-1;
    }

    /* Test indels on left and right */
    for (cL = cloL; cL < chighL; cL++) {
      /* The following check limits genomic inserts (horizontal) and
	 multiple cDNA inserts (vertical). */
      if (left_known[cL] > 0) {
	scoreL = (int) matrixL[cL][rL];
	if (directionsL_nogap[cL][rL] != DIAG) {
	  /* Favor gaps away from intron if possible */
	  scoreL -= 1;
	}

	/* Disallow leftoffset + cL >= rightoffset - cR, or cR >= rightoffset - leftoffset - cL */
	for (cR = cloR; cR < chighR && cR < rightoffset-leftoffset-cL; cR++) {
	  if (right_known[cR] > 0) {
	    scoreR = (int) matrixR[cR][rR];
	    if (directionsR_nogap[cR][rR] != DIAG) {
	      /* Favor gaps away from intron if possible */
	      scoreR -= 1;
	    }

	    if ((score = scoreL + scoreR) > bestscore ||
		(score >= bestscore && jump_late_p)) {  /* Use >= for jump late */
	      bestp = false;
	      if (watsonp == true) {
		splicesitepos1 = leftoffset + cL;
		splicesitepos2 = rightoffset - cR + 1;
		if (IIT_exists_with_divno_signed(splicing_iit,splicing_divint_crosstable[chrnum],
						 splicesitepos1,splicesitepos2+1U,/*sign*/cdna_direction) == true) {
		  bestp = true;
		}
	      } else {
		splicesitepos1 = (chrhigh - chroffset) - leftoffset - cL + 1;
		splicesitepos2 = (chrhigh - chroffset) - rightoffset + cR;
		if (IIT_exists_with_divno_signed(splicing_iit,splicing_divint_crosstable[chrnum],
						 splicesitepos2,splicesitepos1+1U,/*sign*/-cdna_direction) == true) {
		  bestp = true;
		}
	      }
	      if (bestp == true) {
		debug3(printf("At %d left to %d right, score is (%d)+(%d) = %d (bestscore)\n",
			      cL,cR,scoreL,scoreR,score));
		bestscore = score;
		*bestrL = rL;
		*bestrR = rR;
		*bestcL = cL;
		*bestcR = cR;
	      } else {
		debug3a(printf("At %d left to %d right, score is (%d)+(%d) = %d\n",
			       cL,cR,scoreL,scoreR,score));
	      }
	    }
	  }
	}
      }
    }
  }

  *best_introntype = NONINTRON;
  return bestscore;
}


/* Returns finalscore */
static int
bridge_intron_gap_site_level (int *bestrL, int *bestrR, int *bestcL, int *bestcR,
			      Score32_T **matrixL, Score32_T **matrixR,
			      Direction32_T **directionsL_nogap, Direction32_T **directionsR_nogap,
			      char *gsequenceL, char *gsequenceL_alt, char *rev_gsequenceR, char *rev_gsequenceR_alt,
			      int goffsetL, int rev_goffsetR, int *left_known, int *right_known,
			      int rlength, int glengthL, int glengthR,
			      int cdna_direction, bool watsonp, int lbandL, int ubandL, int lbandR, int ubandR,
			      int canonical_reward, int leftoffset, int rightoffset,
			      Univcoord_T chroffset, Univcoord_T chrhigh,
			      bool halfp, bool finalp) {
  int rL, rR, cL, cR;
  int bestrL_with_prob, bestrR_with_prob, bestcL_with_prob, bestcR_with_prob;
  int cloL, chighL;
  int cloR, chighR;
  int introntype;
  int bestscore = NEG_INFINITY_32, score, scoreL, scoreR, scoreI;
  int bestscore_with_prob = NEG_INFINITY_32;
  double *left_probabilities, *right_probabilities, probL, probR, probL_trunc, probR_trunc, bestprob, bestprob_trunc;
  Univcoord_T splicesitepos;
  char left1, left2, right2, right1, left1_alt, left2_alt, right2_alt, right1_alt;
  int *leftdi, *rightdi;
  bool use_prob_p;


  /* Read dinucleotides */
  leftdi = (int *) MALLOCA((glengthL+1) * sizeof(int));
  rightdi = (int *) MALLOCA((glengthR+1) * sizeof(int));

  for (cL = 0; cL < glengthL - 1; cL++) {
    left1 = gsequenceL[cL];
    left1_alt = gsequenceL_alt[cL];
    left2 = gsequenceL[cL+1];
    left2_alt = gsequenceL_alt[cL+1];
    assert(left1 == get_genomic_nt(&left1_alt,goffsetL+cL,chroffset,chrhigh,watsonp));
    assert(left2 == get_genomic_nt(&left2_alt,goffsetL+cL+1,chroffset,chrhigh,watsonp));

    if ((left1 == 'G' || left1_alt == 'G') && (left2 == 'T' || left2_alt == 'T')) {
      leftdi[cL] = LEFT_GT;
    } else if ((left1 == 'G' || left1_alt == 'G') && (left2 == 'C' || left2_alt == 'C')) {
      leftdi[cL] = LEFT_GC;
    } else if ((left1 == 'A' || left1_alt == 'A') && (left2 == 'T' || left2_alt == 'T')) {
      leftdi[cL] = LEFT_AT;
#ifndef PMAP
    } else if ((left1 == 'C' || left1_alt == 'C') && (left2 == 'T' || left2_alt == 'T')) {
      leftdi[cL] = LEFT_CT;
#endif
    } else {
      leftdi[cL] = 0x00;
    }
  }
  leftdi[glengthL-1] = leftdi[glengthL] = 0x00;

  for (cR = 0; cR < glengthR - 1; cR++) {
    right2 = rev_gsequenceR[-cR-1];
    right2_alt = rev_gsequenceR_alt[-cR-1];
    right1 = rev_gsequenceR[-cR];
    right1_alt = rev_gsequenceR_alt[-cR];
    assert(right2 == get_genomic_nt(&right2_alt,rev_goffsetR-cR-1,chroffset,chrhigh,watsonp));
    assert(right1 == get_genomic_nt(&right1_alt,rev_goffsetR-cR,chroffset,chrhigh,watsonp));

    if ((right2 == 'A' || right2_alt == 'A') && (right1 == 'G' || right1_alt == 'G')) {
      rightdi[cR] = RIGHT_AG;
    } else if ((right2 == 'A' || right2_alt == 'A') && (right1 == 'C' || right1_alt == 'C')) {
      rightdi[cR] = RIGHT_AC;
#ifndef PMAP
    } else if ((right2 == 'G' || right2_alt == 'G') && (right1 == 'C' || right1_alt == 'C')) {
      rightdi[cR] = RIGHT_GC;
    } else if ((right2 == 'A' || right2_alt == 'A') && (right1 == 'T' || right1_alt == 'T')) {
      rightdi[cR] = RIGHT_AT;
#endif
    } else {
      rightdi[cR] = 0x00;
    }
  }
  rightdi[glengthR-1] = rightdi[glengthR] = 0x00;


  left_probabilities = (double *) MALLOCA(glengthL * sizeof(double));
  right_probabilities = (double *) MALLOCA(glengthR * sizeof(double));

  debug3(printf("watsonp is %d.  cdna_direction is %d\n",watsonp,cdna_direction));
  if (watsonp == true) {
    if (cdna_direction > 0) {
      for (cL = 0; cL < glengthL - 1; cL++) {
	splicesitepos = chroffset + leftoffset + cL;
	if (left_known[cL]) {
	  left_probabilities[cL] = 1.0;
	} else {
	  left_probabilities[cL] = Maxent_hr_donor_prob(splicesitepos,chroffset);
	  debug3(printf("left donor probability at cL %d is %f\n",cL,left_probabilities[cL]));
	}
      }

      for (cR = 0; cR < glengthR - 1; cR++) {
	splicesitepos = chroffset + rightoffset - cR + 1;
	if (right_known[cR]) {
	  right_probabilities[cR] = 1.0;
	} else {
	  right_probabilities[cR] = Maxent_hr_acceptor_prob(splicesitepos,chroffset);
	  debug3(printf("right acceptor probability at cR %d is %f\n",cR,right_probabilities[cR]));
	}
      }

    } else {
      for (cL = 0; cL < glengthL - 1; cL++) {
	splicesitepos = chroffset + leftoffset + cL;
	if (left_known[cL]) {
	  left_probabilities[cL] = 1.0;
	} else {
	  left_probabilities[cL] = Maxent_hr_antiacceptor_prob(splicesitepos,chroffset);
	  debug3(printf("left antiacceptor probability at cL %d is %f\n",cL,left_probabilities[cL]));
	}
      }

      for (cR = 0; cR < glengthR - 1; cR++) {
	splicesitepos = chroffset + rightoffset - cR + 1;
	if (right_known[cR]) {
	  right_probabilities[cR] = 1.0;
	} else {
	  right_probabilities[cR] = Maxent_hr_antidonor_prob(splicesitepos,chroffset);
	  debug3(printf("right antidonor probability at cR %d is %f\n",cR,right_probabilities[cR]));
	}
      }
    }

  } else {
    if (cdna_direction > 0) {
      for (cL = 0; cL < glengthL - 1; cL++) {
	splicesitepos = chrhigh - leftoffset - cL + 1;
	if (left_known[cL]) {
	  left_probabilities[cL] = 1.0;
	} else {
	  left_probabilities[cL] = Maxent_hr_antidonor_prob(splicesitepos,chroffset);
	  debug3(printf("left antidonor probability at cL %d is %f\n",cL,left_probabilities[cL]));
	}
      }

      for (cR = 0; cR < glengthR - 1; cR++) {
	splicesitepos = chrhigh - rightoffset + cR;
	if (right_known[cR]) {
	  right_probabilities[cR] = 1.0;
	} else {
	  right_probabilities[cR] = Maxent_hr_antiacceptor_prob(splicesitepos,chroffset);
	  debug3(printf("right antiacceptor probability at cR %d is %f\n",cR,right_probabilities[cR]));
	}
      }

    } else {
      for (cL = 0; cL < glengthL - 1; cL++) {
	splicesitepos = chrhigh - leftoffset - cL + 1;
	if (left_known[cL]) {
	  left_probabilities[cL] = 1.0;
	} else {
	  left_probabilities[cL] = Maxent_hr_acceptor_prob(splicesitepos,chroffset);
	  debug3(printf("left acceptor probability at cL %d is %f\n",cL,left_probabilities[cL]));
	}
      }

      for (cR = 0; cR < glengthR - 1; cR++) {
	splicesitepos = chrhigh - rightoffset + cR;
	if (right_known[cR]) {
	  right_probabilities[cR] = 1.0;
	} else {
	  right_probabilities[cR] = Maxent_hr_donor_prob(splicesitepos,chroffset);
	  debug3(printf("right donor probability at cR %d is %f\n",cR,right_probabilities[cR]));
	}
      }
    }
  }

  /* Search using probs and without simultaneously */
  bestscore = NEG_INFINITY_32;
  bestprob = bestprob_trunc = 0.0;
  for (rL = 1, rR = rlength-1; rL < rlength; rL++, rR--) {
    debug3(printf("\nAt row %d on left and %d on right\n",rL,rR));
    if ((cloL = rL - lbandL) < 1) {
      cloL = 1;
    }
    if ((chighL = rL + ubandL) > glengthL-1) {
      chighL = glengthL-1;
    }

    if ((cloR = rR - lbandR) < 1) {
      cloR = 1;
    }
    if ((chighR = rR + ubandR) > glengthR-1) {
      chighR = glengthR-1;
    }

    debug3(printf("A. Test no indels\n"));
    cL = rL;
    probL = left_probabilities[cL];
    if (probL > PROB_CEILING) {
      probL_trunc = PROB_CEILING;
    } else if (probL < PROB_FLOOR) {
      probL_trunc = PROB_FLOOR;
    } else {
      probL_trunc = probL;
    }
    scoreL = (int) matrixL[cL][rL];

    cR = rR;
    probR = right_probabilities[cR];
    if (probR > PROB_CEILING) {
      probR_trunc = PROB_CEILING;
    } else if (probR < PROB_FLOOR) {
      probR_trunc = PROB_FLOOR;
    } else {
      probR_trunc = probR;
    }
    scoreR = (int) matrixR[cR][rR];


    debug3(printf("B. Test indel on right\n"));
    /* Test indel on right */
    cL = rL;
    probL = left_probabilities[cL];
    if (probL > PROB_CEILING) {
      probL_trunc = PROB_CEILING;
    } else if (probL < PROB_FLOOR) {
      probL_trunc = PROB_FLOOR;
    } else {
      probL_trunc = probL;
    }
    scoreL = (int) matrixL[cL][rL];
    if (directionsL_nogap[cL][rL] != DIAG) {
      /* Favor gaps away from intron if possible */
      scoreL -= 1;
    }

    /* Disallow leftoffset + cL >= rightoffset - cR, or cR >= rightoffset - leftoffset - cL */
    for (cR = cloR; cR < chighR && cR < rightoffset-leftoffset-cL; cR++) {
      probR = right_probabilities[cR];
      if (probR > PROB_CEILING) {
	probR_trunc = PROB_CEILING;
      } else if (probR < PROB_FLOOR) {
	probR_trunc = PROB_FLOOR;
      } else {
	probR_trunc = probR;
      }
      scoreR = (int) matrixR[cR][rR];
      if (directionsR_nogap[cR][rR] != DIAG) {
	/* Favor gaps away from intron if possible */
	scoreR -= 1;
      }
	      
#ifdef USE_SCOREI
      scoreI = intron_score(&introntype,leftdi[cL],rightdi[cR],cdna_direction,canonical_reward,finalp);
#else
      scoreI = 0;
#endif
	
      if ((score = scoreL + scoreI + scoreR) > bestscore) {
	debug3(printf("Best score: At %d left to %d right, score is (%d)+(%d)+(%d) = %d (bestscore, prob %f + %f)\n",
		      cL,cR,scoreL,scoreI,scoreR,scoreL+scoreI+scoreR,probL,probR));
	debug3(printf("probL %f, probR %f\n",left_probabilities[cL],right_probabilities[cR]));
	bestscore = score;
	*bestrL = rL;
	*bestrR = rR;
	*bestcL = cL;
	*bestcR = cR;
	bestprob = probL + probR;
      } else if (score == bestscore && probL + probR > bestprob) {
	debug3(printf("Improved prob: At %d left to %d right, score is (%d)+(%d)+(%d) = %d (bestscore, prob %f + %f)\n",
		      cL,cR,scoreL,scoreI,scoreR,scoreL+scoreI+scoreR,probL,probR));
	debug3(printf("probL %f, probR %f\n",left_probabilities[cL],right_probabilities[cR]));
	*bestrL = rL;
	*bestrR = rR;
	*bestcL = cL;
	*bestcR = cR;
	bestprob = probL + probR;
      } else {
	debug3a(printf("Not best score: At %d left to %d right, score is (%d)+(%d)+(%d) = %d (bestscore, prob %f + %f)\n",
		       cL,cR,scoreL,scoreI,scoreR,scoreL+scoreI+scoreR,probL,probR));
      }

      if (probL_trunc + probR_trunc < bestprob_trunc) {
	debug3a(printf("At %d left to %d right, prob is %f + %f = %f\n",
		       cL,cR,probL_trunc,probR_trunc,probL_trunc+probR_trunc));
	  
      } else if (probL_trunc + probR_trunc == bestprob_trunc) {
	debug3(printf("At %d left to %d right, prob is %f + %f = %f\n",
		      cL,cR,probL_trunc,probR_trunc,probL_trunc+probR_trunc));
	  
	if (scoreL + scoreI + scoreR > bestscore_with_prob) {
	  debug3(printf(" (bestscore %d)\n",scoreL+scoreR));
	  bestprob_trunc = probL_trunc + probR_trunc;
	  bestcL_with_prob = cL;
	  bestcR_with_prob = cR;
	  bestrL_with_prob = rL;
	  bestrR_with_prob = rR;
	  bestscore_with_prob = scoreL + scoreI + scoreR;
	}

      } else {
	/* probL_trunc + probR_trunc > bestprob_trunc */
	debug3(printf("At %d left to %d right, prob is %f + %f = %f\n",
		      cL,cR,probL_trunc,probR_trunc,probL_trunc+probR_trunc));
	  
	debug3(printf(" (bestscore %d)\n",scoreL+scoreR));
	bestprob_trunc = probL_trunc + probR_trunc;
	bestcL_with_prob = cL;
	bestcR_with_prob = cR;
	bestrL_with_prob = rL;
	bestrR_with_prob = rR;
	bestscore_with_prob = scoreL + scoreI + scoreR;
      }
    }

    debug3(printf("C. Test indel on left\n"));
    /* Test indel on left */
    cR = rR;
    probR = right_probabilities[cR];
    if (probR > PROB_CEILING) {
      probR_trunc = PROB_CEILING;
    } else if (probR < PROB_FLOOR) {
      probR_trunc = PROB_FLOOR;
    } else {
      probR_trunc = probR;
    }
    scoreR = (int) matrixR[cR][rR];
    if (directionsR_nogap[cR][rR] != DIAG) {
      /* Favor gaps away from intron if possible */
      scoreR -= 1;
    }

    /* Disallow leftoffset + cL >= rightoffset - cR, or cR >= rightoffset - leftoffset - cL */
    for (cL = cloL; cL < chighL && cL < rightoffset-leftoffset-cR; cL++) {
      probL = left_probabilities[cL];
      if (probL > PROB_CEILING) {
	probL_trunc = PROB_CEILING;
      } else if (probL < PROB_FLOOR) {
	probL_trunc = PROB_FLOOR;
      } else {
	probL_trunc = probL;
      }
      scoreL = (int) matrixL[cL][rL];
      if (directionsL_nogap[cL][rL] != DIAG) {
	/* Favor gaps away from intron if possible */
	scoreL -= 1;
      }

#ifdef USE_SCOREI
      scoreI = intron_score(&introntype,leftdi[cL],rightdi[cR],cdna_direction,canonical_reward,finalp);
#else
      scoreI = 0;
#endif
	
      if ((score = scoreL + scoreI + scoreR) > bestscore) {
	debug3(printf("Best score: At %d left to %d right, score is (%d)+(%d)+(%d) = %d (bestscore, prob %f + %f)\n",
		      cL,cR,scoreL,scoreI,scoreR,scoreL+scoreI+scoreR,probL,probR));
	debug3(printf("probL %f, probR %f\n",left_probabilities[cL],right_probabilities[cR]));
	bestscore = score;
	*bestrL = rL;
	*bestrR = rR;
	*bestcL = cL;
	*bestcR = cR;
	bestprob = probL + probR;
      } else if (score == bestscore && probL + probR > bestprob) {
	debug3(printf("Improved prob: At %d left to %d right, score is (%d)+(%d)+(%d) = %d (bestscore, prob %f + %f)\n",
		      cL,cR,scoreL,scoreI,scoreR,scoreL+scoreI+scoreR,probL,probR));
	debug3(printf("probL %f, probR %f\n",left_probabilities[cL],right_probabilities[cR]));
	*bestrL = rL;
	*bestrR = rR;
	*bestcL = cL;
	*bestcR = cR;
	bestprob = probL + probR;
      } else {
	debug3a(printf("Not best score: At %d left to %d right, score is (%d)+(%d)+(%d) = %d (bestscore, prob %f + %f)\n",
		       cL,cR,scoreL,scoreI,scoreR,scoreL+scoreI+scoreR,probL,probR));
      }

      if (probL_trunc + probR_trunc < bestprob_trunc) {
	debug3a(printf("At %d left to %d right, prob is %f + %f = %f\n",
		       cL,cR,probL_trunc,probR_trunc,probL_trunc+probR_trunc));
	  
      } else if (probL_trunc + probR_trunc == bestprob_trunc) {
	debug3(printf("At %d left to %d right, prob is %f + %f = %f\n",
		      cL,cR,probL_trunc,probR_trunc,probL_trunc+probR_trunc));
	  
	if (scoreL + scoreI + scoreR > bestscore_with_prob) {
	  debug3(printf(" (bestscore %d)\n",scoreL+scoreR));
	  bestprob_trunc = probL_trunc + probR_trunc;
	  bestcL_with_prob = cL;
	  bestcR_with_prob = cR;
	  bestrL_with_prob = rL;
	  bestrR_with_prob = rR;
	  bestscore_with_prob = scoreL + scoreI + scoreR;
	}

      } else {
	/* probL_trunc + probR_trunc > bestprob_trunc */
	debug3(printf("At %d left to %d right, prob is %f + %f = %f\n",
		      cL,cR,probL_trunc,probR_trunc,probL_trunc+probR_trunc));
	  
	debug3(printf(" (bestscore %d)\n",scoreL+scoreR));
	bestprob_trunc = probL_trunc + probR_trunc;
	bestcL_with_prob = cL;
	bestcR_with_prob = cR;
	bestrL_with_prob = rL;
	bestrR_with_prob = rR;
	bestscore_with_prob = scoreL + scoreI + scoreR;
      }
    }
  }

  if (bestprob > 2*PROB_CEILING) {
    /* Probability is good with best alignment, so take that */
    debug(printf("Best alignment has good probability\n"));
    use_prob_p = true;
  } else if (left_probabilities[bestcL_with_prob] < PROB_CEILING && right_probabilities[bestcR_with_prob] < PROB_CEILING) {
    /* Probability-based solution is bad, so use alignment */
    debug(printf("Probability-based solution is bad\n"));
    use_prob_p = false;
  } else if (bestscore_with_prob < bestscore - 9) {
    debug(printf("Probability-based solution requires very bad alignment\n"));
    use_prob_p = false;
  } else {
    use_prob_p = true;
  }

  if (use_prob_p == true) {
    /* Best alignment yields bad probability, and probability-based alignment yields good probability, so switch */
    debug(printf("Switch to probability-based solution\n"));
    debug(printf("Non-SIMD. bestscore %d (bestprob %f) vs bestscore_with_prob %d (bestprob_trunc %f, actually %f and %f)\n",
		 bestscore,bestprob,bestscore_with_prob,bestprob_trunc,left_probabilities[bestcL_with_prob],right_probabilities[bestcR_with_prob]));
    *bestcL = bestcL_with_prob;
    *bestcR = bestcR_with_prob;
    *bestrL = bestrL_with_prob;
    *bestrR = bestrR_with_prob;
    bestscore = bestscore_with_prob;
  }


  FREEA(rightdi);
  FREEA(leftdi);
  FREEA(left_probabilities);
  FREEA(right_probabilities);

  if (bestscore < 0) {
    return bestscore;
  } else if (halfp == true) {
    scoreI = intron_score(&introntype,leftdi[*bestcL],rightdi[*bestcR],cdna_direction,canonical_reward,finalp);
    return (bestscore - scoreI/2);
  } else {
    return bestscore;
  }
}


static int
bridge_intron_gap (int *bestrL, int *bestrR, int *bestcL, int *bestcR,
		   int *best_introntype, double *left_prob, double *right_prob,
		   Score32_T **matrixL, Score32_T **matrixR,
		   Direction32_T **directionsL_nogap, Direction32_T **directionsR_nogap, 
		   char *gsequenceL, char *gsequenceL_alt, char *rev_gsequenceR, char *rev_gsequenceR_alt,
		   int goffsetL, int rev_goffsetR, int rlength, int glengthL, int glengthR,
		   int cdna_direction, bool watsonp, int extraband_paired, double defect_rate, int canonical_reward,
		   int leftoffset, int rightoffset,
		   Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh,
		   bool halfp, bool finalp, bool jump_late_p) {
  int finalscore;
  int *left_known, *right_known;
  int ubandL, lbandL, ubandR, lbandR;


  if (glengthL+1 <= 0) {
    fprintf(stderr,"Problem with glengthL = %d\n",glengthL);
    abort();
  }

  if (glengthR+1 <= 0) {
    fprintf(stderr,"Problem with glengthR = %d\n",glengthR);
    abort();
  }

#if 1
  /* Allows unlimited indel lengths */
  ubandL = glengthL - rlength + extraband_paired;
  lbandL = extraband_paired;

  ubandR = glengthR - rlength + extraband_paired;
  lbandR = extraband_paired;
#else
  /* Limit indels to 3 bp around splice sites.  Doesn't work on PacBio reads. */
  ubandL = 3;
  lbandL = 3;

  ubandR = 3;
  lbandR = 3;
#endif

  left_known = (int *) CALLOCA(glengthL+1,sizeof(int));
  right_known = (int *) CALLOCA(glengthR+1,sizeof(int));
  get_known_splicesites(left_known,right_known,glengthL,glengthR,
			/*leftoffset*/goffsetL,/*rightoffset*/rev_goffsetR,
			cdna_direction,watsonp,chrnum,chroffset,chrhigh);

  if (novelsplicingp == false && splicing_iit != NULL && (donor_typeint < 0 || acceptor_typeint < 0)) {
    /* Constrain to given introns */
    finalscore = bridge_intron_gap_intron_level(&(*bestrL),&(*bestrR),&(*bestcL),&(*bestcR),&(*best_introntype),
						matrixL,matrixR,directionsL_nogap,directionsR_nogap,
						left_known,right_known,rlength,glengthL,glengthR,
						cdna_direction,watsonp,lbandL,ubandL,lbandR,ubandR,
						leftoffset,rightoffset,chrnum,chroffset,chrhigh,jump_late_p);

  } else {
    finalscore = bridge_intron_gap_site_level(&(*bestrL),&(*bestrR),&(*bestcL),&(*bestcR),
					      matrixL,matrixR,directionsL_nogap,directionsR_nogap,
					      gsequenceL,gsequenceL_alt,rev_gsequenceR,rev_gsequenceR_alt,goffsetL,rev_goffsetR,
					      left_known,right_known,rlength,glengthL,glengthR,
					      cdna_direction,watsonp,lbandL,ubandL,lbandR,ubandR,
					      canonical_reward,leftoffset,rightoffset,
					      chroffset,chrhigh,halfp,finalp);
  }


#if 0
  /* Determine if result meets given constraints */
  if (*finalscore < 0) {
    result = false;
  } else if (novelsplicingp == true) {
    result = true;
  } else if (splicing_iit == NULL) {
    result = true;
  } else if (donor_typeint >= 0 && acceptor_typeint >= 0) {
    /* If novelsplicingp is false and using splicing at splice site level, require sites to be known */
    if (left_known[*bestcL] == 0 || right_known[*bestcR] == 0) {
      debug(printf("Novel splicing not allowed, so bridge_intron_gap returning false\n"));
      result = false;
    } else {
      result = true;
    }
  } else {
    /* If novelsplicingp is false and using splicing at splice site level, result was already constrained */
    result = true;
  }
#endif

  if (finalscore >= 0) {
    get_splicesite_probs(&(*left_prob),&(*right_prob),*bestcL,*bestcR,
			 left_known,right_known,leftoffset,rightoffset,chroffset,chrhigh,
			 cdna_direction,watsonp);
  }

  FREEA(right_known);
  FREEA(left_known);

#if defined(DEBUG) || defined(DEBUG3)
  printf("Returning final score of %d at (%d,%d) left to (%d,%d) right, with probs %f and %f\n",
	 finalscore,*bestrL,*bestcL,*bestrR,*bestcR,*left_prob,*right_prob);
#endif

  return finalscore;
}
#endif


static List_T
genome_gap_simple (int *finalscore, int *best_introntype, int *new_leftgenomepos, int *new_rightgenomepos,
		   double *left_prob, double *right_prob, int *exonhead, int *nmatches, int *nmismatches,
		   char *rsequence, char *rsequenceuc, char *rev_rsequence, char *rev_rsequenceuc, int rlength, 
		   char *gsequenceL, char *gsequenceL_alt, char *rev_gsequenceR, char *rev_gsequenceR_alt,
		   int roffset, int rev_roffset, int leftoffset, int rightoffset,
		   Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh, Pairpool_T pairpool,
		   int mismatchtype, int canonical_reward,
		   int cdna_direction, bool watsonp, int dynprogindex, bool halfp) {
  List_T pairs = NULL;
  Pair_T gappair;
  bool result;
  int bestrL, bestrR, rL, rR, r;
  int querycoord, genomecoord;
  int na1, na2, na2_alt, c1, c1_uc, c2, c2_alt;
  char left1, left1_alt, left2, left2_alt, right2, right2_alt, right1, right1_alt;
  int bestscore = 0, bestscoreI = 0, scoreL, scoreI, scoreR, pairscore, score;
  int leftdi, rightdi;
  int *left_known, *right_known;
  int introntype;
  Pairdistance_T **pairdistance_array_type;

  debug(printf("Starting genome_gap_simple with cdna_direction %d and watsonp %d\n",cdna_direction,watsonp));
  pairdistance_array_type = pairdistance_array[mismatchtype];

  left_known = (int *) CALLOCA(rlength+1,sizeof(int));
  right_known = (int *) CALLOCA(rlength+1,sizeof(int));
  get_known_splicesites(left_known,right_known,/*glengthL*/rlength,/*glengthR*/rlength,
			leftoffset,rightoffset,
			cdna_direction,watsonp,chrnum,chroffset,chrhigh);

  scoreR = 0;
  for (rR = 1; rR < rlength; rR++) {
    na1 = rev_rsequenceuc[1-rR];
    na2 = rev_gsequenceR[1-rR];
    na2_alt = rev_gsequenceR_alt[1-rR];
    pairscore = pairdistance_array_type[na1][na2];
    if ((score = pairdistance_array_type[na1][na2_alt]) > pairscore) {
      pairscore = score;
    }
    scoreR += pairscore;
  }

  scoreL = 0;
  for (rL = 1, rR = rlength-1; rL < rlength; rL++, rR--) {
    na1 = rsequenceuc[rL-1];
    na2 = gsequenceL[rL-1];
    na2_alt = gsequenceL_alt[rL-1];
    pairscore = pairdistance_array_type[na1][na2];
    if ((score = pairdistance_array_type[na1][na2_alt]) > pairscore) {
      pairscore = score;
    }
    scoreL += pairscore;

    /* Read dinucleotides */
    left1 = gsequenceL[rL];
    left1_alt = gsequenceL_alt[rL];
    left2 = gsequenceL[rL+1];
    left2_alt = gsequenceL_alt[rL+1];

    if ((left1 == 'G' || left1_alt == 'G') && (left2 == 'T' || left2_alt == 'T')) {
      leftdi = LEFT_GT;
    } else if ((left1 == 'G' || left1_alt == 'G') && (left2 == 'C' || left2_alt == 'C')) {
      leftdi = LEFT_GC;
    } else if ((left1 == 'A' || left1_alt == 'A') && (left2 == 'T' || left2_alt == 'T')) {
      leftdi = LEFT_AT;
#ifndef PMAP
    } else if ((left1 == 'C' || left1_alt == 'C') && (left2 == 'T' || left2_alt == 'T')) {
      leftdi = LEFT_CT;
#endif
    } else {
      leftdi = 0x00;
    }

    right2 = rev_gsequenceR[-rR-1];
    right2_alt = rev_gsequenceR_alt[-rR-1];
    right1 = rev_gsequenceR[-rR];
    right1_alt = rev_gsequenceR_alt[-rR];
    if ((right2 == 'A' || right2_alt == 'A') && (right1 == 'G' || right1_alt == 'G')) {
      rightdi = RIGHT_AG;
    } else if ((right2 == 'A' || right2_alt == 'A') && (right1 == 'C' || right1_alt == 'C')) {
      rightdi = RIGHT_AC;
#ifndef PMAP
    } else if ((right2 == 'G' || right2_alt == 'G') && (right1 == 'C' || right1_alt == 'C')) {
      rightdi = RIGHT_GC;
    } else if ((right2 == 'A' || right2_alt == 'A') && (right1 == 'T' || right1_alt == 'T')) {
      rightdi = RIGHT_AT;
#endif
    } else {
      rightdi = 0x00;
    }

    scoreI = intron_score(&introntype,leftdi,rightdi,cdna_direction,canonical_reward,/*finalp*/false);
    if ((introntype != NONINTRON || left_known[rL] > 0 || right_known[rR] > 0) &&
	(score = scoreL + left_known[rL] + scoreI + right_known[rR] + scoreR) >= bestscore) {  /* Use >= for jump late */
      debug(printf("At %d left (%c%c) to %d right (%c%c), score is (%d)+(%d)+(%d)+(%d)+(%d) = %d (bestscore)\n",
		   rL,left1,left2,rR,right2,right1,scoreL,left_known[rL],scoreI,right_known[rR],scoreR,
		   scoreL+left_known[rL]+scoreI+right_known[rR]+scoreR));
      bestscore = score;
      bestscoreI = scoreI;
      bestrL = /* *bestcL = */ rL;
      bestrR = /* *bestcR = */ rR;
      *best_introntype = introntype;
    } else {
      debug(printf("At %d left (%c%c) to %d right (%c%c), score is (%d)+(%d)+(%d)+(%d)+(%d) = %d\n",
		   rL,left1,left2,rR,right2,right1,scoreL,left_known[rL],scoreI,right_known[rR],scoreR,
		   scoreL+left_known[rL]+scoreI+right_known[rR]+scoreR));
    }

    /* Subtract pairscore from cumulative scoreR */
    na1 = rev_rsequenceuc[1-rR];
    na2 = rev_gsequenceR[1-rR];
    na2_alt = rev_gsequenceR_alt[1-rR];
    pairscore = pairdistance_array_type[na1][na2];
    if ((score = pairdistance_array_type[na1][na2_alt]) > pairscore) {
      pairscore = score;
    }
    scoreR -= pairscore;
  }

  if (halfp == true) {
    *finalscore = (int) (bestscore - bestscoreI/2);
  } else {
    *finalscore = (int) bestscore;
  }

  if (*finalscore <= 0) {
    result = false;
  } else if (novelsplicingp == true) {
    result = true;
  } else if (splicing_iit == NULL) {
    result = true;
  } else if (donor_typeint >= 0 && acceptor_typeint >= 0) {
    /* If novelsplicingp is false and using splicing at splice site level, require sites to be known */
    if (left_known[bestrL] == 0 || right_known[bestrR] == 0) {
      debug(printf("Novel splicing not allowed, so bridge_intron_gap returning false\n"));
      result = false;
    } else {
      result = true;
    }
  } else {
    /* If novelsplicingp is false and using splicing at splice site level, result was already constrained */
    result = true;
  }

  if (result == true) {
    get_splicesite_probs(&(*left_prob),&(*right_prob),bestrL,bestrR,
			 left_known,right_known,leftoffset,rightoffset,
			 chroffset,chrhigh,cdna_direction,watsonp);
    debug(printf("Probabilities are %f and %f\n",*left_prob,*right_prob));
    if (*left_prob < 0.90 || *right_prob < 0.90) {
      result = false;
    }
  }

  if (result == true) {
    *nmatches = *nmismatches = 0;

    /* Push from left to right, so we don't need to do List_reverse() later */
    for (r = 1; r <= bestrL; r++) {
      querycoord = genomecoord = r-1;

      c1 = rsequence[querycoord];
      c1_uc = rsequenceuc[querycoord];
      c2 = gsequenceL[genomecoord];
      c2_alt = gsequenceL_alt[genomecoord];

      if (c2 == '*') {
	/* Don't push pairs past end of chromosome */
	debug(printf("Don't push pairs past end of chromosome: genomeoffset %u, genomecoord %u, chroffset %u, chrhigh %u, watsonp %d\n",
		     leftoffset,genomecoord,chroffset,chrhigh,watsonp));

      } else if (c1_uc == c2 || c1_uc == c2_alt) {
	debug(printf("Pushing simple %d,%d [%d,%d] (%c,%c) - match\n",
		     r,/*c*/r,roffset+querycoord,leftoffset+genomecoord,c1_uc,c2));
	*nmatches += 1;
	pairs = Pairpool_push(pairs,pairpool,roffset+querycoord,leftoffset+genomecoord,
			      c1,DYNPROG_MATCH_COMP,c2,c2_alt,dynprogindex);
	
      } else if (consistent_array[(int) c1_uc][(int) c2] == true || consistent_array[(int) c1_uc][(int) c2_alt] == true) {
	debug(printf("Pushing simple %d,%d [%d,%d] (%c,%c) - ambiguous\n",
		     r,/*c*/r,roffset+querycoord,leftoffset+genomecoord,c1_uc,c2));
	*nmatches += 1;
	pairs = Pairpool_push(pairs,pairpool,roffset+querycoord,leftoffset+genomecoord,
			      c1,AMBIGUOUS_COMP,c2,c2_alt,dynprogindex);
	
      } else {
	debug(printf("Pushing simple %d,%d [%d,%d] (%c,%c) - mismatch\n",
		     r,/*c*/r,roffset+querycoord,leftoffset+genomecoord,c1_uc,c2));
	*nmismatches += 1;
	pairs = Pairpool_push(pairs,pairpool,roffset+querycoord,leftoffset+genomecoord,
			      c1,MISMATCH_COMP,c2,c2_alt,dynprogindex);
      }
    }

    debug(printf("Pushing a gap\n"));
    *new_leftgenomepos = leftoffset+(bestrL-1);
    *new_rightgenomepos = *exonhead = rightoffset-(bestrR-1);
#ifndef NOGAPHOLDER
    pairs = Pairpool_push_gapholder(pairs,pairpool,/*queryjump*/0,/*genomejump*/(*new_rightgenomepos)-(*new_leftgenomepos)-1,
				    /*leftpair*/NULL,/*rightpair*/NULL,/*knownp*/false);
    gappair = (Pair_T) pairs->first;
    gappair->introntype = introntype;
    gappair->donor_prob = *left_prob;
    gappair->acceptor_prob = *right_prob;
#endif

    for (r = bestrR; r > 0; r--) {
      querycoord = genomecoord = 1-r;

      c1 = rev_rsequence[querycoord];
      c1_uc = rev_rsequenceuc[querycoord];
      c2 = rev_gsequenceR[genomecoord];
      c2_alt = rev_gsequenceR_alt[genomecoord];
      
      if (c2 == '*') {
	/* Don't push pairs past end of chromosome */
	debug(printf("Don't push pairs past end of chromosome: genomeoffset %u, genomecoord %u, chroffset %u, chrhigh %u, watsonp %d\n",
		     rightoffset,genomecoord,chroffset,chrhigh,watsonp));

      } else if (c1_uc == c2 || c1_uc == c2_alt) {
	debug(printf("Pushing simple %d,%d [%d,%d] (%c,%c) - match\n",
		     r,/*c*/r,rev_roffset+querycoord,rightoffset+genomecoord,c1_uc,c2));
	*nmatches += 1;
	pairs = Pairpool_push(pairs,pairpool,rev_roffset+querycoord,rightoffset+genomecoord,
			      c1,DYNPROG_MATCH_COMP,c2,c2_alt,dynprogindex);
	
      } else if (consistent_array[(int) c1_uc][(int) c2] == true || consistent_array[(int) c1_uc][(int) c2_alt] == true) {
	debug(printf("Pushing simple %d,%d [%d,%d] (%c,%c) - ambiguous\n",
		     r,/*c*/r,rev_roffset+querycoord,rightoffset+genomecoord,c1_uc,c2));
	*nmatches += 1;
	pairs = Pairpool_push(pairs,pairpool,rev_roffset+querycoord,rightoffset+genomecoord,
			      c1,AMBIGUOUS_COMP,c2,c2_alt,dynprogindex);
	
      } else {
	debug(printf("Pushing simple %d,%d [%d,%d] (%c,%c) - mismatch\n",
		     r,/*c*/r,rev_roffset+querycoord,rightoffset+genomecoord,c1_uc,c2));
	*nmismatches += 1;
	pairs = Pairpool_push(pairs,pairpool,rev_roffset+querycoord,rightoffset+genomecoord,
			      c1,MISMATCH_COMP,c2,c2_alt,dynprogindex);
      }
    }

  }

  FREEA(right_known);
  FREEA(left_known);

  return pairs;
}




/* A genome gap is usually an intron.  Sequence 2L and 2R represent
   the two genomic ends of the intron. */
List_T
Dynprog_genome_gap (int *dynprogindex, int *finalscore, int *new_leftgenomepos, int *new_rightgenomepos,
		    double *left_prob, double *right_prob,
		    int *nmatches, int *nmismatches, int *nopens, int *nindels, int *exonhead, int *introntype,
		    T dynprogL, T dynprogR,
		    char *rsequence, char *rsequenceuc, int rlength, int glengthL, int glengthR, 
		    int roffset, int goffsetL, int rev_goffsetR, 
		    Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh,
		    int cdna_direction, bool watsonp, bool jump_late_p, Pairpool_T pairpool, int extraband_paired,
		    double defect_rate, int maxpeelback, bool halfp, bool finalp,
		    bool splicingp) {
  List_T pairs = NULL;
  Pair_T gappair;
  char *gsequenceL, *gsequenceL_alt, *rev_gsequenceR, *rev_gsequenceR_alt;
  char *rev_rsequence, *rev_rsequenceuc;
  Mismatchtype_T mismatchtype;
  int canonical_reward;
  int rev_roffset, bestrL = -1, bestrR, bestcL, bestcR;
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
  /* int queryjump, genomejump; */

  debug(printf("\n"));
  debug(printf("%c:  ",*dynprogindex > 0 ? (*dynprogindex-1)%26+'a' : (-(*dynprogindex)-1)%26+'A'));
  debug(printf("Aligning genome gap with cdna_direction %d, defect_rate %f\n",cdna_direction,defect_rate));
#ifdef EXTRACT_GENOMICSEG
  debug(printf("At genomic offset %d-%d, %.*s\n",goffsetL,goffsetL+glengthL-1,glengthL,gsequenceL));
  debug(printf("At genomic offset %d-%d, %.*s\n",rev_goffsetR-glengthR+1,rev_goffsetR,glengthR,&(rev_gsequenceR[-glengthR+1])));
#endif
  debug(printf("\n"));

  /* ?check if offsets are too close.  But this eliminates a segment
     of the cDNA.  Should check in stage 3, and do single gap instead. */
  /*
  if (goffsetL+glengthL-1 >= rev_goffsetR-glengthR+1) {
    debug(printf("Bounds don't make sense\n"));
    *finalscore = NEG_INFINITY_16;
    return NULL;
  }
  */

  *nmatches = *nmismatches = *nopens = *nindels = 0;
  *left_prob = *right_prob = 0.0;
  *introntype = NONINTRON;
  if (rlength <= 1) {
    *finalscore = NEG_INFINITY_32;
    return (List_T) NULL;
  }

  if (defect_rate < DEFECT_HIGHQ) {
    mismatchtype = HIGHQ;
    if (rlength > maxpeelback * 4) {
      debug(printf("rlength %d is greater than maxpeelback %d * 4, so using single gap penalties\n",
		   rlength,maxpeelback));
      open = SINGLE_OPEN_HIGHQ;
      extend = SINGLE_EXTEND_HIGHQ;
    } else {
      open = PAIRED_OPEN_HIGHQ;
      extend = PAIRED_EXTEND_HIGHQ;
    }
    if (splicingp == false) {
      canonical_reward = 0;
    } else if (finalp == true) {
      canonical_reward = FINAL_CANONICAL_INTRON_HIGHQ;
    } else {
      canonical_reward = CANONICAL_INTRON_HIGHQ;
    }
  } else if (defect_rate < DEFECT_MEDQ) {
    mismatchtype = MEDQ;
    if (rlength > maxpeelback * 4) {
      debug(printf("rlength %d is greater than maxpeelback %d * 4, so using single gap penalties\n",
		   rlength,maxpeelback));
      open = SINGLE_OPEN_MEDQ;
      extend = SINGLE_EXTEND_MEDQ;
    } else {
      open = PAIRED_OPEN_MEDQ;
      extend = PAIRED_EXTEND_MEDQ;
    }
    if (splicingp == false) {
      canonical_reward = 0;
    } else if (finalp == true) {
      canonical_reward = FINAL_CANONICAL_INTRON_MEDQ;
    } else {
      canonical_reward = CANONICAL_INTRON_MEDQ;
    }
  } else {
    mismatchtype = LOWQ;
    if (rlength > maxpeelback * 4) {
      debug(printf("rlength %d is greater than maxpeelback %d * 4, so using single gap penalties\n",
		   rlength,maxpeelback));
      open = SINGLE_OPEN_LOWQ;
      extend = SINGLE_EXTEND_LOWQ;
    } else {
      open = PAIRED_OPEN_LOWQ;
      extend = PAIRED_EXTEND_LOWQ;
    }
    if (splicingp == false) {
      canonical_reward = 0;
    } else if (finalp == true) {
      canonical_reward = FINAL_CANONICAL_INTRON_LOWQ;
    } else {
      canonical_reward = CANONICAL_INTRON_LOWQ;
    }
  }

  if (rlength > dynprogL->max_rlength || glengthL > dynprogL->max_glength) {
    debug(printf("rlength %d or glengthL %d is too long.  Returning NULL\n",rlength,glengthL));
    *new_leftgenomepos = goffsetL-1;
    *new_rightgenomepos = rev_goffsetR+1;
    *exonhead = rev_roffset = roffset+rlength-1;
#ifndef NOGAPHOLDER
    /*
    queryjump = rev_roffset - roffset + 1;
    genomejump = rev_goffsetR - goffsetL + 1;
    pairs = Pairpool_push_gapholder(NULL,pairpool,queryjump,genomejump,false);
    */
#endif
    *dynprogindex += (*dynprogindex > 0 ? +1 : -1);
    *finalscore = NEG_INFINITY_32;
    *introntype = NONINTRON;
    return (List_T) NULL;
  }

  if (rlength > dynprogR->max_rlength || glengthR > dynprogR->max_glength) {
    debug(printf("rlength %d or glengthR %d is too long.  Returning NULL\n",rlength,glengthR));
    *new_leftgenomepos = goffsetL-1;
    *new_rightgenomepos = rev_goffsetR+1;
    *exonhead = rev_roffset = roffset+rlength-1;
#ifndef NOGAPHOLDER
    /*
    queryjump = rev_roffset - roffset + 1;
    genomejump = rev_goffsetR - goffsetL + 1;
    pairs = Pairpool_push_gapholder(NULL,pairpool,queryjump,genomejump,false);
    */
#endif
    *dynprogindex += (*dynprogindex > 0 ? +1 : -1);
    *finalscore = NEG_INFINITY_32;
    *introntype = NONINTRON;
    return (List_T) NULL;
  }

  rev_rsequence = &(rsequence[rlength-1]);
  rev_rsequenceuc = &(rsequenceuc[rlength-1]);
  debug(printf("At query offset %d-%d, %.*s\n",roffset,roffset+rlength-1,rlength,rsequence));

  rev_roffset = roffset+rlength-1;

  gsequenceL = (char *) MALLOCA((glengthL+1) * sizeof(char));
  gsequenceL_alt = (char *) MALLOCA((glengthL+1) * sizeof(char));
  rev_gsequenceR = (char *) MALLOCA((glengthR+1) * sizeof(char));
  rev_gsequenceR_alt = (char *) MALLOCA((glengthR+1) * sizeof(char));

  if (watsonp) {
    Genome_get_segment_blocks_right(gsequenceL,gsequenceL_alt,/*left*/chroffset+goffsetL,
                                    glengthL,chrhigh,/*revcomp*/false);
    Genome_get_segment_blocks_left(rev_gsequenceR,rev_gsequenceR_alt,/*right*/chroffset+rev_goffsetR+1,
                                   glengthR,chroffset,/*revcomp*/false);
  } else {
    Genome_get_segment_blocks_left(gsequenceL,gsequenceL_alt,/*right*/chrhigh-goffsetL+1,
                                   glengthL,chroffset,/*revcomp*/true);
    Genome_get_segment_blocks_right(rev_gsequenceR,rev_gsequenceR_alt,/*left*/chrhigh-rev_goffsetR,
                                    glengthR,chrhigh,/*revcomp*/true);
  }
  if (gsequenceL[0] == '\0' || rev_gsequenceR[0] == '\0') {
    FREEA(rev_gsequenceR_alt);
    FREEA(rev_gsequenceR);
    FREEA(gsequenceL_alt);
    FREEA(gsequenceL);
    *finalscore = NEG_INFINITY_32;
    return (List_T) NULL;
  }


  debug(printf("At genomic offset %d-%d, %.*s\n",goffsetL,goffsetL+glengthL-1,glengthL,gsequenceL));
  debug(printf("At genomic offset %d-%d, %.*s\n",rev_goffsetR-glengthR+1,rev_goffsetR,glengthR,rev_gsequenceR));

  /* In low-identity alignments, the simple procedure can
     lead to multiple mismatches, which will invalidate the intron
     because of its neighborhood */
  if (finalp == false && defect_rate < DEFECT_MEDQ &&
      (pairs = genome_gap_simple(&(*finalscore),&(*introntype),&(*new_leftgenomepos),&(*new_rightgenomepos),
				 &(*left_prob),&(*right_prob),&(*exonhead),&(*nmatches),&(*nmismatches),
				 rsequence,rsequenceuc,rev_rsequence,rev_rsequenceuc,
				 rlength,gsequenceL,gsequenceL_alt,
				 &(rev_gsequenceR[glengthR-1]),&(rev_gsequenceR_alt[glengthR-1]),
				 roffset,rev_roffset,goffsetL,rev_goffsetR,
				 chrnum,chroffset,chrhigh,pairpool,mismatchtype,canonical_reward,
				 cdna_direction,watsonp,*dynprogindex,halfp)) != NULL) {
    debug(printf("simple procedure worked\n"));

    FREEA(rev_gsequenceR_alt);
    FREEA(rev_gsequenceR);
    FREEA(gsequenceL_alt);
    FREEA(gsequenceL);

    *dynprogindex += (*dynprogindex > 0 ? +1 : -1);
    return pairs;		/* Already reversed */
  }


#if defined(HAVE_SSE4_1) || defined(HAVE_SSE2)
  /* Use || because we want the minimum length (which determines the diagonal length) to achieve a score less than 128 */
  if (rlength <= SIMD_MAXLENGTH_EPI8 || (glengthL <= SIMD_MAXLENGTH_EPI8 && glengthR <= SIMD_MAXLENGTH_EPI8)) {
    use8p = true;
  } else {
    use8p = false;
  }
#endif

#if defined(HAVE_SSE4_1) || defined(HAVE_SSE2)
  if (use8p == true) {
    Dynprog_compute_bands(&lbandL,&ubandL,rlength,glengthL,extraband_paired,/*widebandp*/true);
    debug3(printf("Computing matrix8L_upper\n"));
    matrix8L_upper = Dynprog_simd_8_upper(&directions8L_upper_nogap,&directions8L_upper_Egap,dynprogL,
					  rsequence,gsequenceL,gsequenceL_alt,rlength,glengthL,
#ifdef DEBUG14
					  goffsetL,chroffset,chrhigh,watsonp,
#endif
					  mismatchtype,open,extend,ubandL,jump_late_p,/*revp*/false);

    debug3(printf("Computing matrix8L_lower\n"));
    matrix8L_lower = Dynprog_simd_8_lower(&directions8L_lower_nogap,&directions8L_lower_Egap,dynprogL,
					  rsequence,gsequenceL,gsequenceL_alt,rlength,glengthL,
#ifdef DEBUG14
					  goffsetL,chroffset,chrhigh,watsonp,
#endif
					  mismatchtype,open,extend,lbandL,jump_late_p,/*revp*/false);
    

    Dynprog_compute_bands(&lbandR,&ubandR,rlength,glengthR,extraband_paired,/*widebandp*/true);
    debug3(printf("Computing matrix8R_upper\n"));
    matrix8R_upper = Dynprog_simd_8_upper(&directions8R_upper_nogap,&directions8R_upper_Egap,dynprogR,
					  rev_rsequence,&(rev_gsequenceR[glengthR-1]),&(rev_gsequenceR_alt[glengthR-1]),
					  rlength,glengthR,
#ifdef DEBUG14
					  rev_goffsetR,chroffset,chrhigh,watsonp,
#endif
					  mismatchtype,open,extend,ubandR,/*for revp true*/!jump_late_p,/*revp*/true);

    debug3(printf("Computing matrix8R_lower\n"));
    matrix8R_lower = Dynprog_simd_8_lower(&directions8R_lower_nogap,&directions8R_lower_Egap,dynprogR,
					  rev_rsequence,&(rev_gsequenceR[glengthR-1]),&(rev_gsequenceR_alt[glengthR-1]),
					  rlength,glengthR,
#ifdef DEBUG14
					  rev_goffsetR,chroffset,chrhigh,watsonp,
#endif
					  mismatchtype,open,extend,lbandR,/*for revp true*/!jump_late_p,/*revp*/true);

    if ((*finalscore = bridge_intron_gap_8_ud(&bestrL,&bestrR,&bestcL,&bestcR,
					      &(*introntype),&(*left_prob),&(*right_prob),
					      matrix8L_upper,matrix8L_lower,matrix8R_upper,matrix8R_lower,
					      directions8L_upper_nogap,directions8L_lower_nogap,
					      directions8R_upper_nogap,directions8R_lower_nogap,
					      gsequenceL,gsequenceL_alt,&(rev_gsequenceR[glengthR-1]),&(rev_gsequenceR_alt[glengthR-1]),
					      goffsetL,rev_goffsetR,rlength,glengthL,glengthR,
					      cdna_direction,watsonp,lbandL,ubandL,lbandR,ubandR,defect_rate,
					      canonical_reward,goffsetL,rev_goffsetR,
					      chrnum,chroffset,chrhigh,halfp,finalp,jump_late_p)) < 0) {
      FREEA(rev_gsequenceR_alt);
      FREEA(rev_gsequenceR);
      FREEA(gsequenceL_alt);
      FREEA(gsequenceL);

      return (List_T) NULL;

    } else {
      *new_leftgenomepos = goffsetL+(bestcL-1);
      *new_rightgenomepos = rev_goffsetR-(bestcR-1);
      debug(printf("New leftgenomepos = %d, New rightgenomepos = %d\n",*new_leftgenomepos,*new_rightgenomepos));

      *exonhead = rev_roffset-(bestrR-1);

      if (bestcR >= bestrR) {
	pairs = Dynprog_traceback_8_upper(NULL,&(*nmatches),&(*nmismatches),&(*nopens),&(*nindels),
					  directions8R_upper_nogap,directions8R_upper_Egap,
					  bestrR,bestcR,rev_rsequence,rev_rsequenceuc,
					  &(rev_gsequenceR[glengthR-1]),&(rev_gsequenceR_alt[glengthR-1]),
					  rev_roffset,rev_goffsetR,pairpool,/*revp*/true,
					  chroffset,chrhigh,cdna_direction,watsonp,*dynprogindex);
      } else {
	pairs = Dynprog_traceback_8_lower(NULL,&(*nmatches),&(*nmismatches),&(*nopens),&(*nindels),
					  directions8R_lower_nogap,directions8R_lower_Egap,
					  bestrR,bestcR,rev_rsequence,rev_rsequenceuc,
					  &(rev_gsequenceR[glengthR-1]),&(rev_gsequenceR_alt[glengthR-1]),
					  rev_roffset,rev_goffsetR,pairpool,/*revp*/true,
					  chroffset,chrhigh,cdna_direction,watsonp,*dynprogindex);
      }
      pairs = List_reverse(pairs);

      /* queryjump = (rev_roffset-bestrR) - (roffset+bestrL) + 1; */
      /* genomejump = (rev_goffsetR-bestcR) - (goffsetL+bestcL) + 1; */
      /* No need to revise queryjump or genomejump, because the above
	 coordinates are internal to the gap. */

      debug(printf("Pushing a gap with genomejump %d, introntype %s, prob %f and %f\n",
		   (*new_rightgenomepos)-(*new_leftgenomepos)-1,
		   Intron_type_string(*introntype),*left_prob,*right_prob));

#ifndef NOGAPHOLDER
      pairs = Pairpool_push_gapholder(pairs,pairpool,/*queryjump*/(rev_roffset-bestrR) - (roffset+bestrL) + 1,
				      /*genomejump*/(*new_rightgenomepos)-(*new_leftgenomepos)-1,
				      /*leftpair*/NULL,/*rightpair*/NULL,/*knownp*/false);
      gappair = (Pair_T) pairs->first;
      gappair->introntype = *introntype;
      gappair->donor_prob = *left_prob;
      gappair->acceptor_prob = *right_prob;
#endif

      if (bestcL >= bestrL) {
	pairs = Dynprog_traceback_8_upper(pairs,&(*nmatches),&(*nmismatches),&(*nopens),&(*nindels),
					  directions8L_upper_nogap,directions8L_upper_Egap,
					  bestrL,bestcL,rsequence,rsequenceuc,gsequenceL,gsequenceL_alt,
					  roffset,goffsetL,pairpool,/*revp*/false,
					  chroffset,chrhigh,cdna_direction,watsonp,*dynprogindex);
      } else {
	pairs = Dynprog_traceback_8_lower(pairs,&(*nmatches),&(*nmismatches),&(*nopens),&(*nindels),
					  directions8L_lower_nogap,directions8L_lower_Egap,
					  bestrL,bestcL,rsequence,rsequenceuc,gsequenceL,gsequenceL_alt,
					  roffset,goffsetL,pairpool,/*revp*/false,
					  chroffset,chrhigh,cdna_direction,watsonp,*dynprogindex);
      }

      if (List_length(pairs) == 1) {
	/* Only a gap inserted */
	pairs = (List_T) NULL;
      }

      FREEA(rev_gsequenceR_alt);
      FREEA(rev_gsequenceR);
      FREEA(gsequenceL_alt);
      FREEA(gsequenceL);

      debug(printf("End of dynprog genome gap\n"));

      *dynprogindex += (*dynprogindex > 0 ? +1 : -1);
      debug3(Pair_dump_list(pairs,true));
      debug3(printf("maxnegscore = %d\n",Pair_maxnegscore(pairs)));
      if (Pair_maxnegscore(pairs) < -10) {
	return (List_T) NULL;
      } else {
	return List_reverse(pairs);
      }
    }

  } else {
    /* Use 16-mers */
    Dynprog_compute_bands(&lbandL,&ubandL,rlength,glengthL,extraband_paired,/*widebandp*/true);
    debug3(printf("Computing matrix16L_upper\n"));
    matrix16L_upper = Dynprog_simd_16_upper(&directions16L_upper_nogap,&directions16L_upper_Egap,dynprogL,
					    rsequence,gsequenceL,gsequenceL_alt,rlength,glengthL,
#ifdef DEBUG14
					    goffsetL,chroffset,chrhigh,watsonp,
#endif
					    mismatchtype,open,extend,ubandL,jump_late_p,/*revp*/false);

    debug3(printf("Computing matrix16L_lower\n"));
    matrix16L_lower = Dynprog_simd_16_lower(&directions16L_lower_nogap,&directions16L_lower_Egap,dynprogL,
					    rsequence,gsequenceL,gsequenceL_alt,rlength,glengthL,
#ifdef DEBUG14
					    goffsetL,chroffset,chrhigh,watsonp,
#endif
					    mismatchtype,open,extend,lbandL,jump_late_p,/*revp*/false);

    Dynprog_compute_bands(&lbandR,&ubandR,rlength,glengthR,extraband_paired,/*widebandp*/true);
    debug3(printf("Computing matrix16R_upper\n"));
    matrix16R_upper = Dynprog_simd_16_upper(&directions16R_upper_nogap,&directions16R_upper_Egap,dynprogR,
					    rev_rsequence,&(rev_gsequenceR[glengthR-1]),&(rev_gsequenceR_alt[glengthR-1]),
					    rlength,glengthR,
#ifdef DEBUG14
					    rev_goffsetR,chroffset,chrhigh,watsonp,
#endif
					    mismatchtype,open,extend,ubandR,/*for revp true*/!jump_late_p,/*revp*/true);

    debug3(printf("Computing matrix16R_lower\n"));
    matrix16R_lower = Dynprog_simd_16_lower(&directions16R_lower_nogap,&directions16R_lower_Egap,dynprogR,
					    rev_rsequence,&(rev_gsequenceR[glengthR-1]),&(rev_gsequenceR_alt[glengthR-1]),
					    rlength,glengthR,
#ifdef DEBUG14
					    rev_goffsetR,chroffset,chrhigh,watsonp,
#endif
					    mismatchtype,open,extend,lbandR,/*for revp true*/!jump_late_p,/*revp*/true);
    
    if ((*finalscore = bridge_intron_gap_16_ud(&bestrL,&bestrR,&bestcL,&bestcR,
					       &(*introntype),&(*left_prob),&(*right_prob),
					       matrix16L_upper,matrix16L_lower,matrix16R_upper,matrix16R_lower,
					       directions16L_upper_nogap,directions16L_lower_nogap,
					       directions16R_upper_nogap,directions16R_lower_nogap,
					       gsequenceL,gsequenceL_alt,&(rev_gsequenceR[glengthR-1]),&(rev_gsequenceR_alt[glengthR-1]),
					       goffsetL,rev_goffsetR,rlength,glengthL,glengthR,
					       cdna_direction,watsonp,lbandL,ubandL,lbandR,ubandR,defect_rate,
					       canonical_reward,goffsetL,rev_goffsetR,
					       chrnum,chroffset,chrhigh,halfp,finalp,jump_late_p)) < 0) {

      FREEA(rev_gsequenceR_alt);
      FREEA(rev_gsequenceR);
      FREEA(gsequenceL_alt);
      FREEA(gsequenceL);

      return (List_T) NULL;

    } else {
      *new_leftgenomepos = goffsetL+(bestcL-1);
      *new_rightgenomepos = rev_goffsetR-(bestcR-1);
      debug(printf("New leftgenomepos = %d, New rightgenomepos = %d\n",*new_leftgenomepos,*new_rightgenomepos));

      *exonhead = rev_roffset-(bestrR-1);

      if (bestcR >= bestrR) {
	pairs = Dynprog_traceback_16_upper(NULL,&(*nmatches),&(*nmismatches),&(*nopens),&(*nindels),
					   directions16R_upper_nogap,directions16R_upper_Egap,
					   bestrR,bestcR,rev_rsequence,rev_rsequenceuc,
					   &(rev_gsequenceR[glengthR-1]),&(rev_gsequenceR_alt[glengthR-1]),
					   rev_roffset,rev_goffsetR,pairpool,/*revp*/true,
					   chroffset,chrhigh,cdna_direction,watsonp,*dynprogindex);
      } else {
	pairs = Dynprog_traceback_16_lower(NULL,&(*nmatches),&(*nmismatches),&(*nopens),&(*nindels),
					   directions16R_lower_nogap,directions16R_lower_Egap,
					   bestrR,bestcR,rev_rsequence,rev_rsequenceuc,
					   &(rev_gsequenceR[glengthR-1]),&(rev_gsequenceR_alt[glengthR-1]),
					   rev_roffset,rev_goffsetR,pairpool,/*revp*/true,
					   chroffset,chrhigh,cdna_direction,watsonp,*dynprogindex);
      }
      pairs = List_reverse(pairs);

      /* queryjump = (rev_roffset-bestrR) - (roffset+bestrL) + 1; */
      /* genomejump = (rev_goffsetR-bestcR) - (goffsetL+bestcL) + 1; */
      /* No need to revise queryjump or genomejump, because the above
	 coordinates are internal to the gap. */

      debug(printf("Pushing a gap with genomejump %d, introntype %s, prob %f and %f\n",
		   (*new_rightgenomepos)-(*new_leftgenomepos)-1,
		   Intron_type_string(*introntype),*left_prob,*right_prob));
#ifndef NOGAPHOLDER
      pairs = Pairpool_push_gapholder(pairs,pairpool,/*queryjump*/(rev_roffset-bestrR) - (roffset+bestrL) + 1,
				      /*genomejump*/(*new_rightgenomepos)-(*new_leftgenomepos)-1,
				      /*leftpair*/NULL,/*rightpair*/NULL,/*knownp*/false);
      gappair = (Pair_T) pairs->first;
      gappair->introntype = *introntype;
      gappair->donor_prob = *left_prob;
      gappair->acceptor_prob = *right_prob;
#endif

      if (bestcL >= bestrL) {
	pairs = Dynprog_traceback_16_upper(pairs,&(*nmatches),&(*nmismatches),&(*nopens),&(*nindels),
					   directions16L_upper_nogap,directions16L_upper_Egap,
					   bestrL,bestcL,rsequence,rsequenceuc,
					   gsequenceL,gsequenceL_alt,roffset,goffsetL,pairpool,/*revp*/false,
					   chroffset,chrhigh,cdna_direction,watsonp,*dynprogindex);
      } else {
	pairs = Dynprog_traceback_16_lower(pairs,&(*nmatches),&(*nmismatches),&(*nopens),&(*nindels),
					   directions16L_lower_nogap,directions16L_lower_Egap,
					   bestrL,bestcL,rsequence,rsequenceuc,
					   gsequenceL,gsequenceL_alt,roffset,goffsetL,pairpool,/*revp*/false,
					   chroffset,chrhigh,cdna_direction,watsonp,*dynprogindex);
      }

      if (List_length(pairs) == 1) {
	/* Only a gap inserted */
	pairs = (List_T) NULL;
      }

      FREEA(rev_gsequenceR_alt);
      FREEA(rev_gsequenceR);
      FREEA(gsequenceL_alt);
      FREEA(gsequenceL);

      debug(printf("End of dynprog genome gap\n"));

      *dynprogindex += (*dynprogindex > 0 ? +1 : -1);
      debug3(Pair_dump_list(pairs,true));
      debug3(printf("maxnegscore = %d\n",Pair_maxnegscore(pairs)));
      if (Pair_maxnegscore(pairs) < -10) {
	return (List_T) NULL;
      } else {
	return List_reverse(pairs);
      }
    }

  }

#else
  /* Non-SIMD methods */
  Dynprog_compute_bands(&lbandL,&ubandL,rlength,glengthL,extraband_paired,/*widebandp*/true);
  debug3(printf("Computing matrixL\n"));
  matrixL = Dynprog_standard(&directionsL_nogap,&directionsL_Egap,&directionsL_Fgap,dynprogL,
			     rsequence,gsequenceL,gsequenceL_alt,rlength,glengthL,
			     goffsetL,chroffset,chrhigh,watsonp,
			     mismatchtype,open,extend,lbandL,ubandL,
			     jump_late_p,/*revp*/false,/*saturation*/NEG_INFINITY_INT);
  
  Dynprog_compute_bands(&lbandR,&ubandR,rlength,glengthR,extraband_paired,/*widebandp*/true);
  debug3(printf("Computing matrixR\n"));
  matrixR = Dynprog_standard(&directionsR_nogap,&directionsR_Egap,&directionsR_Fgap,dynprogR,
			     rev_rsequence,&(rev_gsequenceR[glengthR-1]),&(rev_gsequenceR_alt[glengthR-1]),
			     rlength,glengthR,rev_goffsetR,chroffset,chrhigh,watsonp,
			     mismatchtype,open,extend,lbandL,ubandR,
			     /*for revp true*/!jump_late_p,/*revp*/true,/*saturation*/NEG_INFINITY_INT);
  
  if ((*finalscore = bridge_intron_gap(&bestrL,&bestrR,&bestcL,&bestcR,
				       &(*introntype),&(*left_prob),&(*right_prob),
				       matrixL,matrixR,directionsL_nogap,directionsR_nogap,
				       gsequenceL,gsequenceL_alt,&(rev_gsequenceR[glengthR-1]),&(rev_gsequenceR_alt[glengthR-1]),
				       goffsetL,rev_goffsetR,rlength,glengthL,glengthR,
				       cdna_direction,watsonp,extraband_paired,defect_rate,
				       canonical_reward,goffsetL,rev_goffsetR,
				       chrnum,chroffset,chrhigh,halfp,finalp,jump_late_p)) < 0) {
    
    FREEA(gsequenceL_alt);
    FREEA(rev_gsequenceR_alt);
    FREEA(gsequenceL);
    FREEA(rev_gsequenceR);

    return (List_T) NULL;
    
  } else {
    *new_leftgenomepos = goffsetL+(bestcL-1);
    *new_rightgenomepos = rev_goffsetR-(bestcR-1);
    debug(printf("New leftgenomepos = %d, New rightgenomepos = %d\n",*new_leftgenomepos,*new_rightgenomepos));
    
    *exonhead = rev_roffset-(bestrR-1);

    pairs = Dynprog_traceback_std(NULL,&(*nmatches),&(*nmismatches),&(*nopens),&(*nindels),
				  directionsR_nogap,directionsR_Egap,directionsR_Fgap,bestrR,bestcR,
				  rev_rsequence,rev_rsequenceuc,
				  &(rev_gsequenceR[glengthR-1]),&(rev_gsequenceR_alt[glengthR-1]),
				  rev_roffset,rev_goffsetR,pairpool,/*revp*/true,
				  chroffset,chrhigh,cdna_direction,watsonp,*dynprogindex);
    pairs = List_reverse(pairs);

    /* queryjump = (rev_roffset-bestrR) - (roffset+bestrL) + 1; */
    /* genomejump = (rev_goffsetR-bestcR) - (goffsetL+bestcL) + 1; */
    /* No need to revise queryjump or genomejump, because the above
       coordinates are internal to the gap. */
    
    debug(printf("Pushing a gap with genomejump %d, introntype %s, prob %f and %f\n",
		 (*new_rightgenomepos)-(*new_leftgenomepos)-1,
		 Intron_type_string(*introntype),*left_prob,*right_prob));
#ifndef NOGAPHOLDER
    pairs = Pairpool_push_gapholder(pairs,pairpool,/*queryjump*/(rev_roffset-bestrR) - (roffset+bestrL) + 1,
				    /*genomejump*/(*new_rightgenomepos)-(*new_leftgenomepos)-1,
				    /*leftpair*/NULL,/*rightpair*/NULL,/*knownp*/false);
    gappair = (Pair_T) pairs->first;
    gappair->introntype = *introntype;
    gappair->donor_prob = *left_prob;
    gappair->acceptor_prob = *right_prob;
#endif

    pairs = Dynprog_traceback_std(pairs,&(*nmatches),&(*nmismatches),&(*nopens),&(*nindels),
				  directionsL_nogap,directionsL_Egap,directionsL_Fgap,bestrL,bestcL,
				  rsequence,rsequenceuc,
				  gsequenceL,gsequenceL_alt,roffset,goffsetL,pairpool,/*revp*/false,
				  chroffset,chrhigh,cdna_direction,watsonp,*dynprogindex);

    if (List_length(pairs) == 1) {
      /* Only a gap inserted */
      pairs = (List_T) NULL;
    }

    FREEA(gsequenceL_alt);
    FREEA(rev_gsequenceR_alt);
    FREEA(gsequenceL);
    FREEA(rev_gsequenceR);

    debug(printf("End of dynprog genome gap\n"));
    
    *dynprogindex += (*dynprogindex > 0 ? +1 : -1);
    debug3(Pair_dump_list(pairs,true));
    debug3(printf("maxnegscore = %d\n",Pair_maxnegscore(pairs)));
    if (Pair_maxnegscore(pairs) < -10) {
      return (List_T) NULL;
    } else {
      return List_reverse(pairs);
    }
  }
#endif

}


