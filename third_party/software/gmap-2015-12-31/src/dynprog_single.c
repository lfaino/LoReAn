static char rcsid[] = "$Id: dynprog_single.c 181922 2016-01-08 00:43:31Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "dynprog_single.h"
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
#include "boyer-moore.h"
#include "complement.h"
#include "intron.h"
#include "maxent.h"
#include "maxent_hr.h"
#include "fastlog.h"
#include "dynprog_simd.h"


/* Tests whether get_genomic_nt == genomicseg in compute_scores procedures */
/* #define EXTRACT_GENOMICSEG 1 */


/* Prints parameters and results */
#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif

/* Microexon search */
#ifdef DEBUG1
#define debug1(x) x
#else
#define debug1(x)
#endif

/* Getting genomic nt */
#ifdef DEBUG8
#define debug8(x) x
#else
#define debug8(x)
#endif

/* print_vector */
#ifdef DEBUG15
#define debug15(x) x
#else
#define debug15(x)
#endif

/* Homopolymer (e.g., PacBio) */
#ifdef DEBUG16
#define debug16(x) x
#else
#define debug16(x)
#endif

/* Homopolymer details */
#ifdef DEBUG16A
#define debug16a(x) x
#else
#define debug16a(x)
#endif


#define MICROEXON_PVALUE_HIGHQ 0.01
#define MICROEXON_PVALUE_MEDQ 0.001
#define MICROEXON_PVALUE_LOWQ 0.0001
#define ENDSEQUENCE_PVALUE 0.001 /* Have stricter threshold for making end exons */

#define MIN_MICROEXON_LENGTH 3
#ifdef PMAP
#define MAX_MICROEXON_LENGTH 17	/* Should be oligomer length - 1 plus peelback */
#else
#define MAX_MICROEXON_LENGTH 12	/* Should be oligomer length - 1 plus peelback */
#endif
#define MICROINTRON_LENGTH 9



#define T Dynprog_T

bool homopolymerp;

void
Dynprog_single_setup (bool homopolymerp_in) {

  homopolymerp = homopolymerp_in;

  return;
}


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


/************************************************************************
 *  Homopolymer
 ************************************************************************/

static List_T
augment_pairs (List_T pairs, int *rsequence_nreps, int r_uniqlength, int roffset,
	       int *gsequence_nreps, int g_uniqlength, int goffset,
	       Pairpool_T pairpool, int dynprogindex) {
  List_T augmented = NULL, p;
  Pair_T pair;
  int r, c, r_nreps, c_nreps, i;
  int r_cum = 0, c_cum = 0;


#ifdef DEBUG16
  printf("r_nreps: ");
  for (i = 0; i < r_uniqlength; i++) {
    printf("%d",rsequence_nreps[i]);
  }
  printf("\n");

  printf("g_nreps: ");
  for (i = 0; i < g_uniqlength; i++) {
    printf("%d",gsequence_nreps[i]);
  }
  printf("\n");
#endif

  for (p = pairs; p != NULL; p = p->rest) {
    pair = (Pair_T) p->first;
    r = pair->querypos - roffset;
    c = pair->genomepos - goffset;

    if (pair->comp == INDEL_COMP) {
      augmented = Pairpool_push_existing(augmented,pairpool,pair);
      debug16a(Pair_dump_one(pair,true));
      pair->querypos += r_cum;
      pair->genomepos += c_cum;

      if (pair->cdna == ' ') {
	c_nreps = gsequence_nreps[c];
	debug16a(printf(" genomepos %u, c_cum %d, nreps %d\n",c,c_cum,c_nreps));
	if (c_nreps == 0) {
	  /* Do nothing */
	} else {
	  for (i = 1; i <= c_nreps; i++) {
	    augmented = Pairpool_push(augmented,pairpool,r+r_cum+roffset,c+c_cum+i+goffset,
				      /*cdna*/' ',INDEL_COMP,pair->genome,pair->genomealt,
				      dynprogindex);
	  }
	  c_cum += c_nreps;
	}

      } else if (pair->genome == ' ') {
	r_nreps = rsequence_nreps[r];
	debug16a(printf(" querypos %d, r_cum %d, nreps %d\n",r,r_cum,r_nreps));
	if (r_nreps == 0) {
	  /* Do nothing */
	} else {
	  for (i = 1; i <= r_nreps; i++) {
	    augmented = Pairpool_push(augmented,pairpool,r+r_cum+i+roffset,c+c_cum+goffset,
				      pair->cdna,INDEL_COMP,/*genome*/' ',/*genomealt*/' ',
				      dynprogindex);
	  }
	  r_cum += r_nreps;
	}

      } else {
	fprintf(stderr,"Indel pair is missing both cdna and genome nts\n");
	abort();
      }

    } else {
      r_nreps = rsequence_nreps[r];
      c_nreps = gsequence_nreps[c];
      debug16a(printf(" querypos %d, r_cum %d, nreps %d, genomepos %u, c_cum %d, nreps %d\n",
		      r,r_cum,r_nreps,c,c_cum,c_nreps));
      augmented = Pairpool_push_existing(augmented,pairpool,pair);
      pair->querypos += r_cum;
      pair->genomepos += c_cum;
      if (r_nreps == 0 && c_nreps == 0) {
	/* Do nothing */
      } else if (r_nreps == c_nreps) {
	for (i = 1; i <= r_nreps; i++) {
	  augmented = Pairpool_push_copy(augmented,pairpool,pair);
	  ((Pair_T) augmented->first)->querypos += i;
	  ((Pair_T) augmented->first)->genomepos += i;
	}
	r_cum += r_nreps;
	c_cum += c_nreps;

      } else if (r_nreps < c_nreps) {
	for (i = 1; i <= r_nreps; i++) {
	  augmented = Pairpool_push_copy(augmented,pairpool,pair);
	  ((Pair_T) augmented->first)->querypos += i;
	  ((Pair_T) augmented->first)->genomepos += i;
	}
	r_cum += r_nreps;

	for ( ; i <= c_nreps; i++) {
	  /* Add 1 to r to advance to next coordinate */
	  augmented = Pairpool_push(augmented,pairpool,r+r_cum+roffset + 1,c+c_cum+i+goffset,
				    /*cdna*/' ',INDEL_COMP,pair->genome,pair->genomealt,
				    dynprogindex);
	}
	c_cum += c_nreps;

      } else {
	for (i = 1; i <= c_nreps; i++) {
	  augmented = Pairpool_push_copy(augmented,pairpool,pair);
	  ((Pair_T) augmented->first)->querypos += i;
	  ((Pair_T) augmented->first)->genomepos += i;
	}
	c_cum += c_nreps;

	for ( ; i <= r_nreps; i++) {
	  /* Add 1 to c to advance to next coordinate */
	  augmented = Pairpool_push(augmented,pairpool,r+r_cum+i+roffset,c+c_cum+goffset + 1,
				    pair->cdna,INDEL_COMP,/*genome*/' ',/*genomealt*/' ',
				    dynprogindex);
	}
	r_cum += r_nreps;

      }
    }
  }

  debug16(Pair_dump_list(augmented,true));

  return List_reverse(augmented);
}

static char *
uniq_string (int **nreps, int *uniqlength, char *string, int length) {
  char *uniq, *p, nt, lastnt;
  int i, k, *a;

  *uniqlength = 1;
  lastnt = string[0];
  for (i = 1; i < length; i++) {
    if ((nt = string[i]) != lastnt) {
      (*uniqlength)++;
      lastnt = nt;
    }
  }

  p = uniq = (char *) MALLOC(((*uniqlength) + 1) * sizeof(char));
  a = *nreps = (int *) MALLOC((*uniqlength) * sizeof(int));
  k = 0;

  lastnt = string[0];
  for (i = 1; i < length; i++) {
    if ((nt = string[i]) != lastnt) {
      *p++ = lastnt;
      *a++ = k;
      lastnt = nt;
      k = 0;
    } else {
      k++;
    }
  }
  *p = lastnt;
  *a = k;

#ifdef DEBUG16
  printf("string: %.*s\n",length,string);
  printf("uniq:   %.*s\n",*uniqlength,uniq);
  printf("nreps   ");
  for (i = 0; i < *uniqlength; i++) {
    printf("%d",(*nreps)[i]);
  }
  printf("\n");
#endif

  return uniq;
}



static List_T
single_gap_simple (int *finalscore, int *nmatches, int *nmismatches,
		   char *rsequence, char *rsequenceuc, int rlength, char *gsequence, char *gsequence_alt,
		   int roffset, int goffset, Pairpool_T pairpool,
		   int mismatchtype, int dynprogindex) {
  int score;
  List_T pairs = NULL;
  int r;
  int querycoord, genomecoord;
  int c1, c1_uc, c2, c2_alt;
  Pairdistance_T **pairdistance_array_type;

  debug(printf("Starting single_gap_simple\n"));
  pairdistance_array_type = pairdistance_array[mismatchtype];

  *finalscore = 0;
  *nmatches = *nmismatches = 0;

  /* Push from left to right, so we don't need to do List_reverse() later */
  for (r = 1; r <= rlength; r++) {
    querycoord = genomecoord = r-1;

    c1 = rsequence[querycoord];
    c1_uc = rsequenceuc[querycoord];
    c2 = gsequence[genomecoord];
    c2_alt = gsequence_alt[genomecoord];

    if (c2 == '*') {
      /* Don't push pairs past end of chromosome */
      debug(printf("Don't push pairs past end of chromosome: genomeoffset %u, genomecoord %u\n",goffset,genomecoord));

    } else if (c1_uc == c2 || c1_uc == c2_alt) {
      debug(printf("Pushing simple %d,%d [%d,%d] (%c,%c) - match\n",
		   r,/*c*/r,roffset+querycoord,goffset+genomecoord,c1_uc,c2));
      score = pairdistance_array_type[c1_uc][c2];
      if (pairdistance_array_type[c1_uc][c2_alt] > score) {
	score = pairdistance_array_type[c1_uc][c2_alt];
      }
      *finalscore += score;
      *nmatches += 1;
      pairs = Pairpool_push(pairs,pairpool,roffset+querycoord,goffset+genomecoord,
			    c1,DYNPROG_MATCH_COMP,c2,c2_alt,dynprogindex);
	
    } else if (consistent_array[(int) c1_uc][(int) c2] == true || consistent_array[(int) c1_uc][(int) c2_alt] == true) {
      debug(printf("Pushing simple %d,%d [%d,%d] (%c,%c) - ambiguous\n",
		   r,/*c*/r,roffset+querycoord,goffset+genomecoord,c1_uc,c2));
      score = pairdistance_array_type[c1_uc][c2];
      if (pairdistance_array_type[c1_uc][c2_alt] > score) {
	score = pairdistance_array_type[c1_uc][c2_alt];
      }
      *finalscore += score;
      *nmatches += 1;
      pairs = Pairpool_push(pairs,pairpool,roffset+querycoord,goffset+genomecoord,
			    c1,AMBIGUOUS_COMP,c2,c2_alt,dynprogindex);
	
    } else {
      debug(printf("Pushing simple %d,%d [%d,%d] (%c,%c) - mismatch\n",
		   r,/*c*/r,roffset+querycoord,goffset+genomecoord,c1_uc,c2));
      score = pairdistance_array_type[c1_uc][c2];
      if (pairdistance_array_type[c1_uc][c2_alt] > score) {
	score = pairdistance_array_type[c1_uc][c2_alt];
      }
      *finalscore += score;
      *nmismatches += 1;
      pairs = Pairpool_push(pairs,pairpool,roffset+querycoord,goffset+genomecoord,
			    c1,MISMATCH_COMP,c2,c2_alt,dynprogindex);
    }
  }

  if (*nmismatches > 1) {
    return (List_T) NULL;
  } else {
    return pairs;
  }
}


List_T
Dynprog_single_gap (int *dynprogindex, int *finalscore, int *nmatches, int *nmismatches, int *nopens, int *nindels,
		    T dynprog, char *rsequence, char *rsequenceuc,
		    int rlength, int glength, int roffset, int goffset,
		    Univcoord_T chroffset, Univcoord_T chrhigh,
		    bool watsonp, bool jump_late_p, Pairpool_T pairpool,
		    int extraband_single, double defect_rate, bool widebandp) {
  List_T pairs = NULL;
  char *gsequence, *gsequence_alt;

  char *gsequence_orig, *rsequence_orig;
  int *gsequence_nreps, *rsequence_nreps;
  int glength_orig, rlength_orig;

  Mismatchtype_T mismatchtype;
  int lband, uband;
  int open, extend;
#if defined(HAVE_SSE4_1) || defined(HAVE_SSE2)
  Score8_T **matrix8;
  Direction8_T **directions8_nogap, **directions8_Egap, **directions8_Fgap;

  Score16_T **matrix16;
  Direction16_T **directions16_nogap, **directions16_Egap, **directions16_Fgap;
#else
  Score32_T **matrix;
  Direction32_T **directions_nogap, **directions_Egap, **directions_Fgap;
#endif
  /* bool onesidegapp; */

  if (defect_rate < DEFECT_HIGHQ) {
    mismatchtype = HIGHQ;
    open = SINGLE_OPEN_HIGHQ;
    extend = SINGLE_EXTEND_HIGHQ;
    /* onesidegapp = false; */
  } else if (defect_rate < DEFECT_MEDQ) {
    mismatchtype = MEDQ;
    open = SINGLE_OPEN_MEDQ;
    extend = SINGLE_EXTEND_MEDQ;
    /* onesidegapp = true; */
  } else {
    mismatchtype = LOWQ;
    open = SINGLE_OPEN_LOWQ;
    extend = SINGLE_EXTEND_LOWQ;
    /* onesidegapp = true; */
  }

#if 0
  if (close_indels_mode == +1) {
    /* Allow close indels */
    onesidegapp = false;
  } else if (close_indels_mode == -1) {
    /* Disallow close indels */
    onesidegapp = true;
  } else {
    /* Allow close indels for high quality alignments, as determined above */
  }
#endif    

  /* Rlength: maxlookback+MAXPEELBACK.  Glength +EXTRAMATERIAL */
  debug(printf("%c:  ",*dynprogindex > 0 ? (*dynprogindex-1)%26+'a' : (-(*dynprogindex)-1)%26+'A'));
  debug(printf("Aligning single gap middle with wideband = %d and extraband %d\n",widebandp,extraband_single));
#ifdef EXTRACT_GENOMICSEG
  debug(printf("At genomic offset %d-%d, %.*s\n",goffset,goffset+glength-1,glength,gsequence));
#endif
  debug(printf("\n"));

#if 0
  /* Can happen in bad genomic regions.  May want to give up, though, if rlength and glength differ greatly. */
  assert(glength < 1000);
#endif
  assert(glength > 0);

  if (rlength > dynprog->max_rlength || glength > dynprog->max_glength) {
    debug(printf("rlength %d or glength %d is too long.  Returning NULL\n",rlength,glength));
    *finalscore = NEG_INFINITY_32;
    *nmatches = *nmismatches = *nopens = *nindels = 0;
#if 0
    /* Don't push a gapholder for single gap, because gapholder already exists */
    pairs = Pairpool_push_gapholder(NULL,pairpool,rlength,glength,
				    /*leftpair*/NULL,/*rightpair*/NULL,/*knownp*/false);
#endif
    *dynprogindex += (*dynprogindex > 0 ? +1 : -1);
    return (List_T) NULL;
  }

  debug(printf("At query offset %d-%d, %.*s\n",roffset,roffset+rlength-1,rlength,rsequence));
  
    /* If extraband_single is too large, then gaps may be inserted on
       both sides, like this

       CACCC   AGAGCAGGCACTGCCT
       |||||--- ||| ||---||||| 
       CACCCCAGGGAGGAG   CTGCCC

    */


  if (homopolymerp == true) {
    rsequence_orig = rsequence;
    rlength_orig = rlength;
    rsequence = uniq_string(&rsequence_nreps,&rlength,rsequence_orig,rlength_orig);

    gsequence_orig = (char *) MALLOCA((glength+1) * sizeof(char));
    gsequence_alt = (char *) MALLOCA((glength+1) * sizeof(char));

    if (watsonp) {
      Genome_get_segment_blocks_right(gsequence_orig,gsequence_alt,/*left*/chroffset+goffset,
                                      glength,chrhigh,/*revcomp*/false);
    } else {
      Genome_get_segment_blocks_left(gsequence_orig,gsequence_alt,/*right*/chrhigh-goffset+1,
                                     glength,chroffset,/*revcomp*/true);
    }
    if (gsequence_orig[0] == '\0') {
      *finalscore = NEG_INFINITY_32;
      *nmatches = *nmismatches = *nopens = *nindels = 0;
      FREEA(gsequence_alt);
      FREEA(gsequence_orig);
      return (List_T) NULL;
    } else {
      glength_orig = glength;
      rsequence = uniq_string(&gsequence_nreps,&glength,gsequence_orig,glength_orig);
    }

  } else {
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
      *finalscore = NEG_INFINITY_32;
      *nmatches = *nmismatches = *nopens = *nindels = 0;
      FREEA(gsequence_alt);
      FREEA(gsequence);
      return (List_T) NULL;
    } else if (glength == rlength &&
	       (pairs = single_gap_simple(&(*finalscore),&(*nmatches),&(*nmismatches),
					  rsequence,rsequenceuc,rlength,gsequence,gsequence_alt,roffset,goffset,
					  pairpool,mismatchtype,*dynprogindex)) != NULL) {
      FREEA(gsequence_alt);
      FREEA(gsequence);

      *nopens = *nindels = 0;
      *dynprogindex += (*dynprogindex > 0 ? +1 : -1);
      return pairs;
    }
  }


  Dynprog_compute_bands(&lband,&uband,rlength,glength,extraband_single,widebandp);
#if defined(HAVE_SSE4_1) || defined(HAVE_SSE2)
  /* Use || because we want the minimum length (which determines the diagonal length) to achieve a score less than 128 */
  /* Use && because we don't want to overflow in either direction */
  if (rlength <= SIMD_MAXLENGTH_EPI8 && glength <= SIMD_MAXLENGTH_EPI8) {
    matrix8 = Dynprog_simd_8(&directions8_nogap,&directions8_Egap,&directions8_Fgap,dynprog,
			     rsequence,gsequence,gsequence_alt,rlength,glength,
#ifdef DEBUG14
			     goffset,chroffset,chrhigh,watsonp,
#endif
			     mismatchtype,open,extend,
			     lband,uband,jump_late_p,/*revp*/false);
    *finalscore = (int) matrix8[glength][rlength];

    *nmatches = *nmismatches = *nopens = *nindels = 0;
    pairs = Dynprog_traceback_8(NULL,&(*nmatches),&(*nmismatches),&(*nopens),&(*nindels),
				directions8_nogap,directions8_Egap,directions8_Fgap,rlength,glength,
				rsequence,rsequenceuc,
				gsequence,gsequence_alt,roffset,goffset,pairpool,/*revp*/false,
				chroffset,chrhigh,/*cdna_direction*/0,watsonp,*dynprogindex);

  } else {
    matrix16 = Dynprog_simd_16(&directions16_nogap,&directions16_Egap,&directions16_Fgap,dynprog,
			       rsequence,gsequence,gsequence_alt,rlength,glength,
#ifdef DEBUG14
			       goffset,chroffset,chrhigh,watsonp,
#endif
			       mismatchtype,open,extend,
			       lband,uband,jump_late_p,/*revp*/false);
    *finalscore = (int) matrix16[glength][rlength];

    *nmatches = *nmismatches = *nopens = *nindels = 0;
    pairs = Dynprog_traceback_16(NULL,&(*nmatches),&(*nmismatches),&(*nopens),&(*nindels),
				 directions16_nogap,directions16_Egap,directions16_Fgap,rlength,glength,
				 rsequence,rsequenceuc,
				 gsequence,gsequence_alt,roffset,goffset,pairpool,/*revp*/false,
				 chroffset,chrhigh,/*cdna_direction*/0,watsonp,*dynprogindex);
  }

#else

  matrix = Dynprog_standard(&directions_nogap,&directions_Egap,&directions_Fgap,dynprog,
			    rsequence,gsequence,gsequence_alt,rlength,glength,
			    goffset,chroffset,chrhigh,watsonp,mismatchtype,open,extend,
			    lband,uband,jump_late_p,/*revp*/false,/*saturation*/NEG_INFINITY_INT);
  *finalscore = (int) matrix[glength][rlength];

  *nmatches = *nmismatches = *nopens = *nindels = 0;
  pairs = Dynprog_traceback_std(NULL,&(*nmatches),&(*nmismatches),&(*nopens),&(*nindels),
				directions_nogap,directions_Egap,directions_Fgap,rlength,glength,
				rsequence,rsequenceuc,
				gsequence,gsequence_alt,roffset,goffset,pairpool,/*revp*/false,
				chroffset,chrhigh,/*cdna_direction*/0,watsonp,*dynprogindex);
#endif

  if (homopolymerp == true) {
    pairs = augment_pairs(pairs,rsequence_nreps,rlength,roffset,
			  gsequence_nreps,glength,goffset,pairpool,*dynprogindex);
    FREE(gsequence_nreps);
    FREEA(gsequence_orig);
    FREEA(gsequence);

    FREE(rsequence_nreps);
    /* Do not free rsequence_orig */
    FREE(rsequence);
    
  } else {
    FREEA(gsequence_alt);
    FREEA(gsequence);
  }

  /*
  Directions_free(directions);
  Matrix_free(matrix);
  */

  debug(printf("End of dynprog single gap\n"));

  *dynprogindex += (*dynprogindex > 0 ? +1 : -1);
  return List_reverse(pairs);
}



static const Except_T microexon_error = {"Microexon error"};

static List_T
make_microexon_pairs_double (int roffsetL, int roffsetM, int roffsetR,
			     int goffsetL, int goffsetM, int goffsetR,
			     int lengthL, int lengthM, int lengthR,
			     char *rsequence, char *rsequenceuc,
			     Univcoord_T chroffset, Univcoord_T chrhigh, bool watsonp,
			     Pairpool_T pairpool, char gapchar, int dynprogindex) {
  List_T pairs = NULL;
  Pair_T gappair;
  char c1, c1_uc, c2, c2_alt;
  int i;

  /* Left segment */
  for (i = 0; i < lengthL; i++) {
    c1 = rsequence[roffsetL+i];
    c1_uc = rsequenceuc[roffsetL+i];

    c2 = get_genomic_nt(&c2_alt,goffsetL+i,chroffset,chrhigh,watsonp);
#ifdef EXTRACT_GENOMICSEG
    assert(c2 == genomicseg[goffsetL+i]);
#endif

    if (c1_uc == c2 || c1_uc == c2_alt) {
      pairs = Pairpool_push(pairs,pairpool,roffsetL+i,goffsetL+i,c1,DYNPROG_MATCH_COMP,c2,c2_alt,
			    dynprogindex);
    } else if (consistent_array[(int) c1_uc][(int) c2] == true || consistent_array[(int) c1_uc][(int) c2_alt] == true) {
      pairs = Pairpool_push(pairs,pairpool,roffsetL+i,goffsetL+i,c1,AMBIGUOUS_COMP,c2,c2_alt,
			    dynprogindex);
    } else {
      pairs = Pairpool_push(pairs,pairpool,roffsetL+i,goffsetL+i,c1,MISMATCH_COMP,c2,c2_alt,
			    dynprogindex);
    }
  }

  /* First gap */
  /* Don't have to adjust querypos/genomepos, since no cdna/genome skips allowed */
  pairs = Pairpool_push_gapholder(pairs,pairpool,/*queryjump*/0,/*genomejump*/goffsetM-(goffsetL+lengthL),
				  /*leftpair*/NULL,/*rightpair*/NULL,/*knownp*/false);
  
  /* Assign pair->comp, because might occur after assign_gap_types */
  gappair = (Pair_T) List_head(pairs);
  gappair->comp = gapchar;


  /* Microexon */
  for (i = 0; i < lengthM; i++) {
    c1 = rsequence[roffsetM+i];
    c1_uc = rsequenceuc[roffsetM+i];

    c2 = get_genomic_nt(&c2_alt,goffsetM+i,chroffset,chrhigh,watsonp);
#ifdef EXTRACT_GENOMICSEG
    assert(c2 == genomicseg[goffsetM+i]);
#endif

    if (c1_uc == c2 || c1_uc == c2_alt) {
      pairs = Pairpool_push(pairs,pairpool,roffsetM+i,goffsetM+i,c1,DYNPROG_MATCH_COMP,c2,c2_alt,
			    dynprogindex);
    } else if (consistent_array[(int) c1_uc][(int) c2] == true || consistent_array[(int) c1_uc][(int) c2_alt] == true) {
      pairs = Pairpool_push(pairs,pairpool,roffsetM+i,goffsetM+i,c1,AMBIGUOUS_COMP,c2,c2_alt,
			    dynprogindex);
    } else {
      pairs = Pairpool_push(pairs,pairpool,roffsetM+i,goffsetM+i,c1,MISMATCH_COMP,c2,c2_alt,
			    dynprogindex);
    }
  }

  /* Second gap */
  /* Don't have to adjust querypos/genomepos, since no cdna/genome skips allowed */
  if (lengthR == 0) {
    /* If lengthR is zero, then we will have a gap after a gap */
    Except_raise(&microexon_error,__FILE__,__LINE__);
  } else {
    pairs = Pairpool_push_gapholder(pairs,pairpool,/*queryjump*/0,/*genomejump*/goffsetR-(goffsetM+lengthM),
				    /*leftpair*/NULL,/*rightpair*/NULL,/*knownp*/false);
  }

  /* Assign pair->comp, because might occur after assign_gap_types */
  gappair = (Pair_T) List_head(pairs);
  gappair->comp = gapchar;

  
  /* Right segment */
  for (i = 0; i < lengthR; i++) {
    c1 = rsequence[roffsetR+i];
    c1_uc = rsequenceuc[roffsetR+i];

    c2 = get_genomic_nt(&c2_alt,goffsetR+i,chroffset,chrhigh,watsonp);
#ifdef EXTRACT_GENOMICSEG
    assert(c2 == genomicseg[goffsetR+i]);
#endif

    if (c1_uc == c2 || c1_uc == c2_alt) {
      pairs = Pairpool_push(pairs,pairpool,roffsetR+i,goffsetR+i,c1,DYNPROG_MATCH_COMP,c2,c2_alt,
			    dynprogindex);
    } else if (consistent_array[(int) c1_uc][(int) c2] == true || consistent_array[(int) c1_uc][(int) c2_alt] == true) {
      pairs = Pairpool_push(pairs,pairpool,roffsetR+i,goffsetR+i,c1,AMBIGUOUS_COMP,c2,c2_alt,
			    dynprogindex);
    } else {
      pairs = Pairpool_push(pairs,pairpool,roffsetR+i,goffsetR+i,c1,MISMATCH_COMP,c2,c2_alt,
			    dynprogindex);
    }
  }

  return pairs;
}


#if 0
static List_T
make_microexon_pairs_single (int roffsetL, int roffsetR,
			     int goffsetL, int goffsetR,
			     int lengthL, int lengthR, char *rsequence, char *rsequenceuc,
			     Univcoord_T chroffset, Univcoord_T chrhigh, bool watsonp,
			     Pairpool_T pairpool, char gapchar, int dynprogindex) {
  List_T pairs = NULL;
  Pair_T gappair;
  char c1, c1_uc, c2, c2_alt;
  int i;

  /* Microexon */
  for (i = 0; i < lengthL; i++) {
    c1 = rsequence[roffsetL+i];
    c1_uc = rsequenceuc[roffsetL+i];

    c2 = get_genomic_nt(&c2_alt,goffsetL+i,chroffset,chrhigh,watsonp);
#ifdef EXTRACT_GENOMICSEG
    assert(c2 == genomicseg[goffsetL+i]);
#endif

    if (c1_uc == c2 || c1_uc == c2_alt) {
      pairs = Pairpool_push(pairs,pairpool,roffsetL+i,goffsetL+i,c1,DYNPROG_MATCH_COMP,c2,c2_alt,
			    dynprogindex);
    } else if (consistent_array[(int) c1_uc][(int) c2] == true || consistent_array[(int) c1_uc][(int) c2_alt] == true) {
      pairs = Pairpool_push(pairs,pairpool,roffsetL+i,goffsetL+i,c1,AMBIGUOUS_COMP,c2,c2_alt,
			    dynprogindex);
    } else {
      pairs = Pairpool_push(pairs,pairpool,roffsetL+i,goffsetL+i,c1,MISMATCH_COMP,c2,c2_alt,
			    dynprogindex);
    }
  }

  /* Gap */
  /* Don't have to adjust querypos/genomepos, since no cdna/genome skips allowed */
  pairs = Pairpool_push_gapholder(pairs,pairpool,/*queryjump*/0,/*genomejump*/goffsetR-(goffsetL+lengthL),
				  /*leftpair*/NULL,/*rightpair*/NULL,/*knownp*/false);

  /* Assign pair->comp, because might occur after assign_gap_types */
  gappair = (Pair_T) List_head(pairs);
  gappair->comp = gapchar;
  

  /* Right segment */
  for (i = 0; i < lengthR; i++) {
    c1 = rsequence[roffsetR+i];
    c1_uc = rsequenceuc[roffsetR+i];

    c2 = get_genomic_nt(&c2_alt,goffsetR+i,chroffset,chrhigh,watsonp);
#ifdef EXTRACT_GENOMICSEG
    assert(c2 == genomicseg[goffsetR+i]);
#endif

    if (c1_uc == c2 || c1_uc == c2_alt) {
      pairs = Pairpool_push(pairs,pairpool,roffsetR+i,goffsetR+i,c1,DYNPROG_MATCH_COMP,c2,c2_alt,
			    dynprogindex);
    } else if (consistent_array[(int) c1_uc][(int) c2] == true || consistent_array[(int) c1_uc][(int) c2_alt] == true) {
      pairs = Pairpool_push(pairs,pairpool,roffsetR+i,goffsetR+i,c1,AMBIGUOUS_COMP,c2,c2_alt,
			    dynprogindex);
    } else {
      pairs = Pairpool_push(pairs,pairpool,roffsetR+i,goffsetR+i,c1,MISMATCH_COMP,c2,c2_alt,
			    dynprogindex);
    }
  }

  return pairs;
}
#endif


List_T
Dynprog_microexon_int (double *bestprob2, double *bestprob3, int *dynprogindex, int *microintrontype,
		       char *rsequence, char *rsequenceuc,
		       int rlength, int glengthL, int glengthR,
		       int roffset, int goffsetL, int rev_goffsetR, int cdna_direction,
		       char *queryseq, char *queryuc,
		       Univcoord_T chroffset, Univcoord_T chrhigh, bool watsonp,
		       Pairpool_T pairpool, double defect_rate) {
  List_T pairs = NULL;
  Intlist_T hits = NULL, p;
#ifdef EXTRACT_GENOMICSEG
  Intlist_T hits_old;
#endif
  int bestcL = -1, bestcR = -1, best_middlelength;
  int middlelength, cL, cR, mincR, maxcR, leftbound, rightbound, textleft, textright,
    best_candidate, candidate, i;
  int span, nmismatches;
  char left1, left2, right2, right1, left1_alt, left2_alt, right2_alt, right1_alt;
  char c, c_alt;
  char c1_alt, c2_alt, c3_alt, c4_alt;
  char intron1, intron2, intron3, intron4, gapchar;
  float pvalue, bestprob = 0.0, prob2, prob3;
  Univcoord_T splicesitepos;


  *bestprob2 = *bestprob3 = 0.0;

  if (defect_rate < DEFECT_HIGHQ) {
    pvalue = MICROEXON_PVALUE_HIGHQ;
  } else if (defect_rate < DEFECT_MEDQ) {
    pvalue = MICROEXON_PVALUE_MEDQ;
  } else {
    pvalue = MICROEXON_PVALUE_LOWQ;
  }

#ifdef PMAP
  intron1 = 'G';
  intron2 = 'T';
  intron3 = 'A';
  intron4 = 'G';
  gapchar = FWD_CANONICAL_INTRON_COMP;
  *microintrontype = GTAG_FWD;
#else
  if (cdna_direction > 0) {
    intron1 = 'G';
    intron2 = 'T';
    intron3 = 'A';
    intron4 = 'G';
    gapchar = FWD_CANONICAL_INTRON_COMP;
    *microintrontype = GTAG_FWD;
  } else if (cdna_direction < 0) {
    intron1 = 'C';
    intron2 = 'T';
    intron3 = 'A';
    intron4 = 'C';
    gapchar = REV_CANONICAL_INTRON_COMP;
    *microintrontype = GTAG_REV;
  } else {
    /* Can occur when called by Stage3_merge_local_splice */
    /* fprintf(stderr,"cdna_direction is 0 in Dynprog_microexon_int\n"); */
    *microintrontype = NONINTRON;
    return NULL;
  }
#endif

#ifdef EXTRACT_GENOMICSEG
  debug1(printf("Begin microexon search for %.*s and %.*s\n",
	       glengthL,gsequenceL,glengthR,&(rev_gsequenceR[-glengthR+1])));
#else
  debug1(printf("Begin microexon search\n"));
#endif

  debug1(printf("  Query sequence is %.*s\n",rlength,rsequence));

  span = rev_goffsetR-goffsetL;
  debug1(printf("  Genomic span is of length %d\n",span));

#if 0
  if (span <= 0) {
    fprintf(stderr,"Bug in Dynprog_microexon_int.  span %d <= 0.  Please report to twu@gene.com\n",span);
    abort();
  } else {
    min_microexon_length = ceilf(-fasterlog(1.0-powf(1.0-pvalue,1.0/(float) span)) / /*log(4)*/1.386294);
  }
  min_microexon_length -= 8;	/* Two donor-acceptor pairs */
  debug1(printf("  Min microexon length is %d\n",min_microexon_length));
  if (min_microexon_length > MAX_MICROEXON_LENGTH) {
    *microintrontype = NONINTRON;
    return NULL;
  } else if (min_microexon_length < MIN_MICROEXON_LENGTH) {
    min_microexon_length = MIN_MICROEXON_LENGTH;
  }
#elif 0
  min_microexon_length = 6;
#endif

  debug1(printf("\nFinding starting boundary on left\n"));
  leftbound = 0;
  nmismatches = 0;
  while (leftbound < rlength - 1 && nmismatches <= 1) {
    debug1(printf("  leftbound = %d, nmismatches = %d.",leftbound,nmismatches));
    c = get_genomic_nt(&c_alt,goffsetL+leftbound,chroffset,chrhigh,watsonp);
#ifdef EXTRACT_GENOMICSEG
    assert(c == gsequence_ucL[leftbound]);
#endif
    debug1(printf("  Comparing %c with %c\n",rsequence[leftbound],c));
#ifdef PMAP
    if (matchtable[rsequence[leftbound]-'A'][c-'A'] == false) {
      nmismatches++;
    }
#else
    if (rsequenceuc[leftbound] != c) {
      nmismatches++;
    }
#endif
    leftbound++;
  }
  leftbound--;			/* This is where the leftmost mismatch occurred */

  debug1(printf("\nFinding starting boundary on right\n"));
  rightbound = 0;
  i = rlength-1;
  nmismatches = 0;
  while (i >= 0 && nmismatches <= 1) {
    debug1(printf("  rightbound = %d, nmismatches = %d.",rightbound,nmismatches));
    c = get_genomic_nt(&c_alt,rev_goffsetR-rightbound,chroffset,chrhigh,watsonp);
#ifdef EXTRACT_GENOMICSEG
    assert(c == rev_gsequence_ucR[-rightbound]);
#endif
    debug1(printf("  Comparing %c with %c\n",rsequence[i],c));
#ifdef PMAP
    if (matchtable[rsequence[i]-'A'][c-'A'] == false) {
      nmismatches++;
    }
#else
    if (rsequenceuc[i] != c) {
      nmismatches++;
    }
#endif
    rightbound++;
    i--;
  }
  rightbound--;			/* This is where the rightmost mismatch occurred */

  debug1(printf("  Left must start before %d from left end of query.  Right must start after %d from right end of query\n",
	       leftbound,rightbound));

  /* We require that cL >= 1 and cR >= 1 so that lengthL and lengthR are >= 1 */
  for (cL = 1; cL <= leftbound; cL++) {
    left1 = get_genomic_nt(&left1_alt,goffsetL+cL,chroffset,chrhigh,watsonp);
    left2 = get_genomic_nt(&left2_alt,goffsetL+cL+1,chroffset,chrhigh,watsonp);
#ifdef EXTRACT_GENOMICSEG
    assert(left1 == gsequence_ucL[cL]);
    assert(left2 == gsequence_ucL[cL+1]);
#endif

    debug1(printf("  %d: %c%c\n",cL,left1,left2));
    if (left1 == intron1 && left2 == intron2) {
      mincR = rlength - MAX_MICROEXON_LENGTH - cL;
      debug1(printf("  mincR %d = rlength %d - MAX_MICROEXON_LENGTH %d - cL %d\n",
		    mincR,rlength,MAX_MICROEXON_LENGTH,cL));
      if (mincR < 1) {
	mincR = 1;
      }
      maxcR = rlength - MIN_MICROEXON_LENGTH - cL;
      debug1(printf("  maxcR %d = rlength %d - MIN_MICROEXON_LENGTH %d - cL %d\n",
		    maxcR,rlength,MIN_MICROEXON_LENGTH,cL));
      if (maxcR > rightbound) {
	maxcR = rightbound;
      }
      debug1(printf("  Found left GT at %d.  Scanning from %d - cL - (1-7), or %d to %d\n",
		   cL,rlength,mincR,maxcR));
      for (cR = mincR; cR <= maxcR; cR++) {
	right2 = get_genomic_nt(&right2_alt,rev_goffsetR-cR-1,chroffset,chrhigh,watsonp);
	right1 = get_genomic_nt(&right1_alt,rev_goffsetR-cR,chroffset,chrhigh,watsonp);
#ifdef EXTRACT_GENOMICSEG
	assert(right2 == rev_gsequence_ucR[-cR-1]);
	assert(right1 == rev_gsequence_ucR[-cR]);
#endif
	debug1(printf("   Checking %d: %c%c\n",cR,right2,right1));
	if (right2 == intron3 && right1 == intron4) {
	  middlelength = rlength - cL - cR;
	  debug1(printf("  Found pair at %d to %d, length %d.  Middle sequence is %.*s\n",
		       cL,cR,middlelength,middlelength,&(rsequence[cL])));
	  
	  textleft = goffsetL + cL + MICROINTRON_LENGTH;
	  textright = rev_goffsetR - cR - MICROINTRON_LENGTH;

	  if (textright >= textleft + middlelength) {
	    hits = BoyerMoore_nt(&(rsequence[cL]),middlelength,textleft,textright-textleft,
				 chroffset,chrhigh,watsonp);
#ifdef EXTRACT_GENOMICSEG
	    hits_old = BoyerMoore(&(rsequenceuc[cL]),middlelength,&(genomicuc[textleft]),textright-textleft);
	    assert(Intlist_equal(hits,hits_old));
	    Intlist_free(&hits_old);
#endif
	    for (p = hits; p != NULL; p = Intlist_next(p)) {
	      candidate = textleft + Intlist_head(p);
#ifdef EXTRACT_GENOMICSEG
	      assert(get_genomic_nt(candidate-2,chroffset,chrhigh,watsonp) == genomicuc[candidate - 2]);
	      assert(get_genomic_nt(candidate-1,chroffset,chrhigh,watsonp) == genomicuc[candidate - 1]);
	      assert(get_genomic_nt(candidate+middlelength,chroffset,chrhigh,watsonp) == genomicuc[candidate + middlelength]);
	      assert(get_genomic_nt(candidate+middlelength+1,chroffset,chrhigh,watsonp) == genomicuc[candidate + middlelength+1]);
#endif
	      if (/*genomicuc[candidate - 2]*/ get_genomic_nt(&c3_alt,candidate-2,chroffset,chrhigh,watsonp) == intron3 &&
		  /*genomicuc[candidate - 1]*/ get_genomic_nt(&c4_alt,candidate-1,chroffset,chrhigh,watsonp)  == intron4 &&
		  /*genomicuc[candidate + middlelength]*/ get_genomic_nt(&c1_alt,candidate+middlelength,chroffset,chrhigh,watsonp) == intron1 &&
		  /*genomicuc[candidate + middlelength + 1]*/ get_genomic_nt(&c2_alt,candidate+middlelength+1,chroffset,chrhigh,watsonp) == intron2) {
		debug1(printf("  Successful microexon at %d >>> %d..%d >>> %d\n",goffsetL+cL,candidate,candidate+middlelength,rev_goffsetR-cR));

		/* Not handling known splice sites yet */
		if (watsonp == true) {
		  if (cdna_direction > 0) {
		    splicesitepos = chroffset + (candidate-1) + 1;
		    prob2 = Maxent_hr_acceptor_prob(splicesitepos,chroffset);
		    splicesitepos = chroffset + candidate+middlelength;
		    prob3 = Maxent_hr_donor_prob(splicesitepos,chroffset);
		  } else {
		    splicesitepos = chroffset + (candidate-1) + 1;
		    prob2 = Maxent_hr_antidonor_prob(splicesitepos,chroffset);
		    splicesitepos = chroffset + candidate+middlelength;
		    prob3 = Maxent_hr_antiacceptor_prob(splicesitepos,chroffset);
		  }
		} else {
		  if (cdna_direction > 0) {
		    splicesitepos = chrhigh - (candidate-1);
		    prob2 = Maxent_hr_antiacceptor_prob(splicesitepos,chroffset);
		    splicesitepos = chrhigh - (candidate+middlelength) + 1;
		    prob3 = Maxent_hr_antidonor_prob(splicesitepos,chroffset);
		  } else {
		    splicesitepos = chrhigh - (candidate-1);
		    prob2 = Maxent_hr_donor_prob(splicesitepos,chroffset);
		    splicesitepos = chrhigh - (candidate+middlelength) + 1;
		    prob3 = Maxent_hr_acceptor_prob(splicesitepos,chroffset);
		  }
		}
	      
		debug1(printf("microexon probabilities: prob2 = %f, prob3 = %f\n",prob2,prob3));
		if (prob2 + prob3 > bestprob) {
		  bestcL = cL;
		  bestcR = cR;
		  best_candidate = candidate;
		  best_middlelength = middlelength;
		  *bestprob2 = prob2;
		  *bestprob3 = prob3;
		  bestprob = prob2 + prob3;
		}
	      }
	    }
	    Intlist_free(&hits);
	  }

	}
      }
    }
  }

  if (bestcL < 0 || bestcR < 0) {
    debug1(printf("End of dynprog microexon int\n"));

    *microintrontype = NONINTRON;
    return NULL;

  } else {
    debug1(printf("Making microexon pairs with candidate %u\n",best_candidate));
    pairs = make_microexon_pairs_double(roffset,/*roffsetM*/roffset+bestcL,/*roffsetR*/roffset+bestcL+best_middlelength,
					goffsetL,/*candidate*/best_candidate,/*goffsetR*/rev_goffsetR-bestcR+1,
					/*lengthL*/bestcL,/*lengthM*/best_middlelength,/*lengthR*/bestcR,
					queryseq,queryuc,
					chroffset,chrhigh,watsonp,pairpool,gapchar,*dynprogindex);
    *dynprogindex += (*dynprogindex > 0 ? +1 : -1);
    return pairs;
  }
}


#if 0
/* Based on probability of seeing a pattern of length n in L is
   1-(1-p1)^L, where p1 is 4^n.  We determine L so chance probability
   is less than ENDSEQUENCE_PVALUE */
static int
search_length (int endlength, int maxlength, bool end_microexons_p) {
  double p1;
  int effective_maxlength, extrant, result;

  if (end_microexons_p == true) {
    extrant = 4;		/* Count the four nucleotides at the intron bounds */
    effective_maxlength = maxlength;
  } else {
    extrant = 0;		/* Don't count the four nucleotides */
    effective_maxlength = 5000;
    if (maxlength < effective_maxlength) {
      effective_maxlength = maxlength;
    }
  }

  if (endlength + extrant > 12) {
    debug(printf("  Search length for endlength of %d is maxlength %d\n",endlength,effective_maxlength));
    return effective_maxlength;
  } else {
    p1 = 1.0/pow(4.0,(double) (endlength + extrant));
    result = (int) (fasterlog(1.0-ENDSEQUENCE_PVALUE)/fasterlog(1-p1));
    debug(printf("  Search length for endlength of %d plus extra nt of %d is %d\n",endlength,extrant,result));
    if (result > effective_maxlength) {
      return effective_maxlength;
    } else {
      return result;
    }
  }
}
#endif


#if 0
/* Not currently used */
List_T
Dynprog_microexon_5 (int *dynprogindex, int *microintrontype, int *microexonlength,
		     char *rev_rsequence, char *rev_rsequenceuc, char *rev_gsequence, char *rev_gsequence_uc,
		     int rlength, int glength, int rev_roffset, int rev_goffset, int cdna_direction,
		     char *queryseq, char *queryuc, char *genomicseg, char *genomicuc,
		     Pairpool_T pairpool, bool end_microexons_p) {
  List_T pairs = NULL;
  Intlist_T hits = NULL, p;
  int endlength, maxc, c, textleft, textright, candidate, nmismatches = 0;
  char right2, right1;
  char intron1, intron2, intron3, intron4, gapchar;

#ifdef PMAP
  intron1 = 'G';
  intron2 = 'T';
  intron3 = 'A';
  intron4 = 'G';
  gapchar = FWD_CANONICAL_INTRON_COMP;
  *microintrontype = GTAG_FWD;
#else
  if (cdna_direction > 0) {
    intron1 = 'G';
    intron2 = 'T';
    intron3 = 'A';
    intron4 = 'G';
    gapchar = FWD_CANONICAL_INTRON_COMP;
    *microintrontype = GTAG_FWD;
  } else if (cdna_direction < 0) {
    intron1 = 'C';
    intron2 = 'T';
    intron3 = 'A';
    intron4 = 'C';
    gapchar = REV_CANONICAL_INTRON_COMP;
    *microintrontype = GTAG_REV;
  } else {
    *microintrontype = NONINTRON;
    return (List_T) NULL;
    abort();
  }
#endif

#ifdef EXTRACT_GENOMICSEG
  debug(printf("Begin microexon search at 5' for %.*s\n",
	       glength,&(rev_gsequence[-glength+1])));
#else
  debug(printf("Begin microexon search at 5'\n"));
#endif

  debug(printf("  Query sequence is %.*s\n",rlength,&(rev_rsequence[-rlength+1])));

  *microexonlength = 0;
  if (glength < rlength) {
    maxc = glength - MIN_MICROEXON_LENGTH;
  } else {
    maxc = rlength - MIN_MICROEXON_LENGTH;
  }
  for (c = 0; c < maxc; c++) {
    right2 = rev_gsequence_uc[-c-1];
    right1 = rev_gsequence_uc[-c];
    debug(printf("   Checking %c%c\n",right2,right1));
#ifdef PMAP
    if (c > 0 && matchtable[rev_rsequence[-c+1]-'A'][rev_gsequence_uc[-c+1]-'A'] == false) {
      nmismatches++;
    }
#else
    if (c > 0 && rev_rsequenceuc[-c+1] != rev_gsequence_uc[-c+1]) {
      nmismatches++;
    }
#endif
    if (nmismatches > 1) {
      debug(printf("   Aborting at %c != %c\n",rev_rsequence[-c+1],rev_gsequence[-c+1]));
      *microintrontype = NONINTRON;
      return NULL;
    }
    if (right2 == intron3 && right1 == intron4) {
      endlength = rlength - c;
      debug(printf("  Found acceptor at %d, length %d.  End sequence is %.*s\n",
		       c,endlength,endlength,&(rev_rsequence[-endlength+1])));

      textright = rev_goffset - c - MICROINTRON_LENGTH;
      textleft = textright - search_length(endlength,textright,end_microexons_p) + MICROINTRON_LENGTH;

      if (textright >= textleft + endlength) {
	hits = BoyerMoore_nt(&(rev_rsequence[-c-endlength+1]),endlength,textleft,textright-textleft,
			     chroffset,chrhigh,watsonp);
	for (p = hits; p != NULL; p = Intlist_next(p)) {
	  candidate = textleft + Intlist_head(p);
	  if (genomicseg[candidate + endlength] == intron1 &&
	      genomicseg[candidate + endlength + 1] == intron2) {
	    debug(printf("  Successful microexon at %d\n",candidate));

	    Intlist_free(&hits);
	    *microexonlength = endlength;
	    pairs = make_microexon_pairs_single(rev_roffset-c-endlength+1,rev_roffset-c+1,
						candidate,rev_goffset-c+1,endlength,c,
						queryseq,queryuc,chroffset,watsonp,pairpool,gapchar,*dynprogindex);
	    *dynprogindex += (*dynprogindex > 0 ? +1 : -1);
	    return pairs;
	  }
	}
	Intlist_free(&hits);
      }

    }
  }

  debug(printf("End of dynprog microexon 5\n"));

  *microintrontype = NONINTRON;
  return NULL;
}
#endif


#if 0
/* Not currently used */
List_T
Dynprog_microexon_3 (int *dynprogindex, int *microintrontype, int *microexonlength, 
		     char *rsequence, char *rsequenceuc, char *gsequence, char *gsequence_uc,
		     int rlength, int glength, int roffset, int goffset, int cdna_direction,
		     char *queryseq, char *queryuc, char *genomicseg, char *genomicuc,
		     int genomiclength, Pairpool_T pairpool, bool end_microexons_p) {
  List_T pairs = NULL;
  Intlist_T hits = NULL, p;
  int endlength, maxc, c, textleft, textright, candidate, nmismatches = 0;
  char left1, left2;
  char intron1, intron2, intron3, intron4, gapchar;

#ifdef PMAP
  intron1 = 'G';
  intron2 = 'T';
  intron3 = 'A';
  intron4 = 'G';
  gapchar = FWD_CANONICAL_INTRON_COMP;
  *microintrontype = GTAG_FWD;
#else
  if (cdna_direction > 0) {
    intron1 = 'G';
    intron2 = 'T';
    intron3 = 'A';
    intron4 = 'G';
    gapchar = FWD_CANONICAL_INTRON_COMP;
    *microintrontype = GTAG_FWD;
  } else if (cdna_direction < 0) {
    intron1 = 'C';
    intron2 = 'T';
    intron3 = 'A';
    intron4 = 'C';
    gapchar = REV_CANONICAL_INTRON_COMP;
    *microintrontype = GTAG_REV;
  } else {
    *microintrontype = NONINTRON;
    return (List_T) NULL;
    abort();
  }
#endif

#ifdef EXTRACT_GENOMICSEG
  debug(printf("Begin microexon search at 3' for %.*s\n",glength,gsequence));
#else
  debug(printf("Begin microexon search at 3'\n"));
#endif

  debug(printf("  Query sequence is %.*s\n",rlength,rsequence));

  *microexonlength = 0;
  if (glength < rlength) {
    maxc = glength - MIN_MICROEXON_LENGTH;
  } else {
    maxc = rlength - MIN_MICROEXON_LENGTH;
  }
  for (c = 0; c < maxc; c++) {
    left1 = gsequence_uc[c];
    left2 = gsequence_uc[c+1];
    debug(printf("   Checking %c%c\n",left1,left2));
#ifdef PMAP
    if (c > 0 && matchtable[rsequence[c-1]-'A'][gsequence_uc[c-1]-'A'] == false) {
      nmismatches++;
    }
#else
    if (c > 0 && rsequenceuc[c-1] != gsequence_uc[c-1]) {
      nmismatches++;
    }
#endif
    if (nmismatches > 1) {
      debug(printf("   Aborting at %c != %c\n",rsequence[c-1],gsequence[c-1]));
      *microintrontype = NONINTRON;
      return NULL;
    }
    if (left1 == intron1 && left2 == intron2) {
      endlength = rlength - c;
      debug(printf("  Found donor at %d, length %d.  End sequence is %.*s\n",
		   c,endlength,endlength,&(rsequence[c])));

      textleft = goffset + c;
      textright = textleft + search_length(endlength,genomiclength-textleft,end_microexons_p);
      
      if (textright >= textleft + endlength) {
	hits = BoyerMoore_nt(&(rsequence[c]),endlength,textleft,textright-textleft,
			     chroffset,chrhigh,watsonp);
	for (p = hits; p != NULL; p = Intlist_next(p)) {
	  candidate = textleft + Intlist_head(p);
	  if (genomicseg[candidate - 2] == intron3 &&
	      genomicseg[candidate - 1] == intron4) {
	    debug(printf("  Successful microexon at %d\n",candidate));

	    Intlist_free(&hits);
	    *microexonlength = endlength;
	    pairs = make_microexon_pairs_single(roffset,roffset+c,
						goffset,candidate,c,endlength,
						queryseq,queryuc,genomicseg,genomicuc,
						pairpool,gapchar,*dynprogindex);
	    *dynprogindex += (*dynprogindex > 0 ? +1 : -1);
	    return pairs;
	  }
	}
	Intlist_free(&hits);
      }

    }
  }

  debug(printf("End of dynprog microexon 3\n"));

  *microintrontype = NONINTRON;
  return NULL;
}
#endif

