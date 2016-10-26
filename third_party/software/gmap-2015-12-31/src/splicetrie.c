static char rcsid[] = "$Id: splicetrie.c 145990 2014-08-25 21:47:32Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "splicetrie.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>		/* For qsort */

#include "assert.h"
#include "mem.h"
#include "iitdef.h"
#include "interval.h"
#include "genome128_hr.h"
#include "splicetrie_build.h"	/* For single_leaf_p and multiple_leaf_p macros */
#include "dynprog_end.h"


#if 0
#define MAX_BEST_NMISMATCHES 1
#endif

/* Finding short-overlap splicing.  Also may want to turn on DEBUG4H in stage1hr.c.  No longer working; try debug2 */
#ifdef DEBUG1
#define debug1(x) x
#else
#define debug1(x)
#endif


/* Finding short exons and short overlap.  Also may want to turn on DEBUG4K in stage1hr.c. */
#ifdef DEBUG2
#define debug2(x) x
#else
#define debug2(x)
#endif

/* Dump coords.  May want to turn on DEBUG7 in dynprog.c. */
#ifdef DEBUG3
#define debug3(x) x
#else
#define debug3(x)
#endif

/* Splicecomp */
#ifdef DEBUG4
#define debug4(x) x
#else
#define debug4(x)
#endif


/* Solve ends */
#ifdef DEBUG7
#define debug7(x) x
#else
#define debug7(x)
#endif


#define clear_start(diff,startdiscard) (diff & (~0U << (startdiscard)))
#define clear_end(diff,enddiscard) (diff & ~(~0U << (enddiscard)))


#ifdef GSNAP
static Genomecomp_T *splicecomp;
#endif
static Univcoord_T *splicesites;
static Genomecomp_T *splicefrags_ref;
static Genomecomp_T *splicefrags_alt;

static Trieoffset_T *trieoffsets_obs;
static Triecontent_T *triecontents_obs;
static Trieoffset_T *trieoffsets_max;
static Triecontent_T *triecontents_max;

static bool snpp;
static bool amb_closest_p;
static int min_shortend;
static bool amb_clip_p;


void
Splicetrie_setup (
#ifdef GSNAP
		  Genomecomp_T *splicecomp_in,
#endif
		  Univcoord_T *splicesites_in, Genomecomp_T *splicefrags_ref_in, Genomecomp_T *splicefrags_alt_in,
		  Trieoffset_T *trieoffsets_obs_in, Triecontent_T *triecontents_obs_in,
		  Trieoffset_T *trieoffsets_max_in, Triecontent_T *triecontents_max_in,
		  bool snpp_in, bool amb_closest_p_in, bool amb_clip_p_in, int min_shortend_in) {

#ifdef GSNAP
  splicecomp = splicecomp_in;
#endif
  splicesites = splicesites_in;
  splicefrags_ref = splicefrags_ref_in;
  splicefrags_alt = splicefrags_alt_in;

  trieoffsets_obs = trieoffsets_obs_in;
  triecontents_obs = triecontents_obs_in;
  trieoffsets_max = trieoffsets_max_in;
  triecontents_max = triecontents_max_in;

  snpp = snpp_in;
  amb_closest_p = amb_closest_p_in;
  amb_clip_p = amb_clip_p_in;
  min_shortend = min_shortend_in;

  return;
}


#ifdef GSNAP
/* Modified from count_mismatches_limit */
bool
Splicetrie_splicesite_p (Univcoord_T left, int pos5, int pos3) {
  int startdiscard, enddiscard;
  Univcoord_T startblocki, endblocki;
  Genomecomp_T *endblock, *ptr;
  Genomecomp_T splicesitep;

  startblocki = (left+pos5)/32U;
  endblocki = (left+pos3)/32U;
  startdiscard = (left+pos5) % 32;
  enddiscard = (left+pos3) % 32;
  
  if (endblocki == startblocki) {
    debug4(printf("** Single block **\n"));
    splicesitep = splicecomp[startblocki];
    splicesitep = clear_start(splicesitep,startdiscard);
    splicesitep = clear_end(splicesitep,enddiscard);

    return (splicesitep ? true : false);

  } else if (endblocki == startblocki + 1) {
    /* Only two blocks to check */

    if (32 - startdiscard >= enddiscard) {
      /* Two blocks to check and more bits counted in startblock */
      debug4(printf("* Two blocks, start block first **\n"));

      /* 1/2: Startblock */
      splicesitep = splicecomp[startblocki];
      splicesitep = clear_start(splicesitep,startdiscard);
      if (splicesitep) {
	return true;
      }
      
      /* 2/2: Endblock */
      splicesitep = splicecomp[endblocki];
      splicesitep = clear_end(splicesitep,enddiscard);

      return (splicesitep ? true : false);

    } else {
      /* Two blocks to check and more bits counted in endblock */
      debug4(printf("** Two blocks, end block first **\n"));

      /* 1/2: Endblock */
      splicesitep = splicecomp[endblocki];
      splicesitep = clear_end(splicesitep,enddiscard);
      if (splicesitep) {
	return true;
      }

      /* 2/2: Startblock */
      splicesitep = splicecomp[startblocki];
      splicesitep = clear_start(splicesitep,startdiscard);
      return (splicesitep ? true : false);
    }

  } else {
    /* More than 2 blocks to check */
    debug4(printf("** More than two blocks **\n"));

    ptr = &(splicecomp[startblocki+1]);
    endblock = &(splicecomp[endblocki]);

    while (ptr < endblock) {
      if (*ptr) {
	return true;
      }
      ptr++;
    }

    if (enddiscard >= 32 - startdiscard) {
      /* More bits in end block */
      debug4(printf("** Final block, end block first **\n"));

      /* n/n: Go first to end block */
      splicesitep = *ptr;
      splicesitep = clear_end(splicesitep,enddiscard);
      if (splicesitep) {
	return true;
      }
      
      /* 1/n: Go second to start block */
      splicesitep = splicecomp[startblocki];
      splicesitep = clear_start(splicesitep,startdiscard);
      debug4(printf("adding start mask %08X\n",clear_start_mask(startdiscard)));
      
      return (splicesitep ? true : false);

    } else {
      /* 1/n: Go first to start block */
      debug4(printf("** Final block, start block first **\n"));

      splicesitep = splicecomp[startblocki];
      splicesitep = clear_start(splicesitep,startdiscard);
      if (splicesitep) {
	return true;
      }

      /* n/n: Go second to end block */
      splicesitep = splicecomp[endblocki];
      splicesitep = clear_end(splicesitep,enddiscard);

      return (splicesitep ? true : false);
    }
  }
}
#endif



/************************************************************************
 *   Using splicetries
 ************************************************************************/

#ifdef USE_2BYTE_RELOFFSETS
static void
get_offsets (int *offseta, int *offsetc, int *offsetg, int *offsett,
	     Trieoffset_T offsets1, Trieoffset_T offsets2) {

  *offsetc = (int) (offsets1 & 0xffff);
  *offseta = (int) ((offsets1 >>= 16) & 0xffff);

  *offsett = (int) (offsets2 & 0xffff);
  *offsetg = (int) ((offsets2 >>= 16) & 0xffff);

  return;
}
#endif


/* Modified from Splicetrie_dump in splicetrie_build.c */
static int
splicetrie_size (Triecontent_T *triestart) {
  int size;
  Triecontent_T leaf;
  int nleaves;
  int offseta, offsetc, offsetg, offsett;

  if (single_leaf_p(leaf = triestart[0])) {
    return 1;

  } else if (multiple_leaf_p(leaf)) {
    nleaves = (int) (-leaf);
    return nleaves;

  } else {
#ifdef USE_2BYTE_RELOFFSETS
    get_offsets(&offseta,&offsetc,&offsetg,&offsett,triestart[1],triestart[2]);
#else
    offseta = (int) triestart[1];
    offsetc = (int) triestart[2];
    offsetg = (int) triestart[3];
    offsett = (int) triestart[4];
#endif

    size = 0;
    if (offseta > 0) {
      size += splicetrie_size(&(triestart[-offseta]));
    }
    if (offsetc > 0) {
      size += splicetrie_size(&(triestart[-offsetc]));
    }
    if (offsetg > 0) {
      size += splicetrie_size(&(triestart[-offsetg]));
    }
    if (offsett > 0) {
      size += splicetrie_size(&(triestart[-offsett]));
    }

    return size;
  }
}



static List_T
solve_end5_aux (Univcoord_T **coordsptr, Univcoord_T *coords,
		List_T best_pairs, Triecontent_T *triestart,
		Univcoord_T knownsplice_limit_low, Univcoord_T knownsplice_limit_high,
		       
		int *finalscore, int *nmatches, int *nmismatches,
		int *nopens, int *nindels, bool *knownsplicep, int *ambig_end_length,
		int *threshold_miss_score, int obsmax_penalty, int perfect_score,

		Univcoord_T anchor_splicesite, char *splicejunction, char *splicejunction_alt,
		int splicelength, int contlength, Splicetype_T far_splicetype,
		Univcoord_T chroffset, Univcoord_T chrhigh,
		int *dynprogindex, Dynprog_T dynprog, 
		char *revsequence1, char *revsequenceuc1,
		int length1, int length2, int revoffset1, int revoffset2,
		int cdna_direction, bool watsonp, bool jump_late_p, Pairpool_T pairpool,
		int extraband_end, double defect_rate) {
  Triecontent_T leaf;
  int nleaves, i;
  Univcoord_T splicecoord;
  Chrpos_T shortest_intron_length = -1U, intron_length;
  int offseta, offsetc, offsetg, offsett;
  int spliceoffset2_anchor, spliceoffset2_far;

  int score, miss_score, nmatches0, nmismatches0, nopens0, nindels0;
  int length_distal;
  List_T pairs;
    
  if (single_leaf_p(leaf = triestart[0])) {
    splicecoord = splicesites[leaf];
    debug7(printf("Checking leaf %d at %u: ",(int) leaf,splicecoord));
    if (splicecoord >= knownsplice_limit_low && splicecoord <= knownsplice_limit_high &&
	Dynprog_make_splicejunction_5(splicejunction,splicejunction_alt,splicecoord,splicelength,contlength,far_splicetype,watsonp) == true) {
      debug7(printf("intron length %d, ",splicecoord - anchor_splicesite));
      debug7(printf("length1 = %d, length2 = %d, chroffset = %u, splicecoord = %u\n",
		    length1,length2,chroffset,splicecoord));
      if (watsonp) {
	spliceoffset2_anchor = revoffset2;
	spliceoffset2_far = spliceoffset2_anchor - anchor_splicesite + splicecoord;
      } else {
	spliceoffset2_anchor = revoffset2;
	spliceoffset2_far = spliceoffset2_anchor + anchor_splicesite - splicecoord;
      }
      if ((pairs = Dynprog_end5_splicejunction(&(*dynprogindex),&score,&miss_score,&nmatches0,&nmismatches0,
					       &nopens0,&nindels0,dynprog,revsequence1,revsequenceuc1,
					       /*revsequence2*/&(splicejunction[length2-1]),/*revsequenceuc2*/&(splicejunction[length2-1]),
					       /*revsequencealt2*/&(splicejunction_alt[length2-1]),
					       length1,length2,revoffset1,spliceoffset2_anchor,spliceoffset2_far,
					       chroffset,chrhigh,cdna_direction,watsonp,jump_late_p,pairpool,
					       extraband_end,defect_rate,contlength)) != NULL) {

	/* miss_score = perfect_score - score; */
	assert(miss_score <= 0);
	debug7(printf("score %d - perfect score %d = miss %d expected vs %d returned.  ",
		      score,perfect_score,score-perfect_score,miss_score));
	debug7(printf("  comparing against threshold_miss %d + obsmax_penalty %d\n",*threshold_miss_score,obsmax_penalty));
	if (score > 0 && miss_score > *threshold_miss_score + obsmax_penalty) {
	  debug7(printf("miss %d > threshold %d + %d",miss_score,*threshold_miss_score,obsmax_penalty));
#if 0
	  /* Just use results from Dynprog_end5_splicejunction */
	  pairs = Dynprog_add_known_splice_5(&length_distal,pairs,anchor_splicesite,splicecoord,chroffset,
					     watsonp,pairpool);
#else
	  length_distal = length1 - contlength;
#endif
	  best_pairs = pairs;
	  *finalscore = score;
	  *nmatches = nmatches0;
	  *nmismatches = nmismatches0;
	  *nopens = nopens0;
	  *nindels = nindels0;
	  *knownsplicep = true;
	  *ambig_end_length = length_distal;
	  *threshold_miss_score = miss_score - obsmax_penalty;
	  shortest_intron_length = (splicecoord > anchor_splicesite) ? splicecoord - anchor_splicesite : anchor_splicesite - splicecoord;
	  *coordsptr = coords;
	  *(*coordsptr)++ = splicecoord;

	} else if (miss_score == *threshold_miss_score + obsmax_penalty
#if 0
		   && Genomicposlist_find(*coords,splicecoord) == false
#endif
		   ) {
	  if (amb_closest_p == false) {
	    debug7(printf("miss %d == threshold %d + %d, so ambiguous",miss_score,*threshold_miss_score,obsmax_penalty));
	    /* best_pairs = (List_T) NULL; */
	    *(*coordsptr)++ = splicecoord;
	  } else {
	    intron_length = (splicecoord > anchor_splicesite) ? splicecoord - anchor_splicesite : anchor_splicesite - splicecoord;
	    if (intron_length > shortest_intron_length) {
	      debug7(printf("miss %d == threshold %d + %d, but intron_length %d > shortest %d, so ignore",
			    miss_score,*threshold_miss_score,obsmax_penalty,intron_length,shortest_intron_length));
	    } else {
	      debug7(printf("miss %d == threshold %d + %d, but intron_length %d < shortest %d, so new best",
			    miss_score,*threshold_miss_score,obsmax_penalty,intron_length,shortest_intron_length));
#if 0
	      /* Just use results from Dynprog_end5_splicejunction */
	      pairs = Dynprog_add_known_splice_5(&length_distal,pairs,anchor_splicesite,splicecoord,chroffset,
						 watsonp,pairpool);
#else
	      length_distal = length1 - contlength;
#endif
	      best_pairs = pairs;
	      *finalscore = score;
	      *nmatches = nmatches0;
	      *nmismatches = nmismatches0;
	      *nopens = nopens0;
	      *nindels = nindels0;
	      *knownsplicep = true;
	      *ambig_end_length = length_distal;
	      *threshold_miss_score = miss_score - obsmax_penalty;
	      shortest_intron_length = intron_length;
	      *coordsptr = coords;
	      *(*coordsptr)++ = splicecoord;
	    }
	  }
	}
      }
    }
    debug7(printf("\n"));

  } else if (multiple_leaf_p(leaf)) {
    nleaves = (int) (-leaf);
    for (i = 1; i <= nleaves; i++) {
      leaf = triestart[i];
      splicecoord = splicesites[leaf];
      debug7(printf("Checking leaf %d at %u: ",(int) leaf,splicecoord));
      if (splicecoord >= knownsplice_limit_low && splicecoord <= knownsplice_limit_high &&
	  Dynprog_make_splicejunction_5(splicejunction,splicejunction_alt,splicecoord,splicelength,contlength,far_splicetype,watsonp) == true) {
	debug7(printf("intron length %d, ",splicecoord - anchor_splicesite));
	debug7(printf("length1 = %d, length2 = %d, chroffset = %u, splicecoord = %u\n",
		      length1,length2,chroffset,splicecoord));
	if (watsonp) {
	  spliceoffset2_anchor = revoffset2;
	  spliceoffset2_far = spliceoffset2_anchor - anchor_splicesite + splicecoord;
	} else {
	  spliceoffset2_anchor = revoffset2;
	  spliceoffset2_far = spliceoffset2_anchor + anchor_splicesite - splicecoord;
	}
	if ((pairs = Dynprog_end5_splicejunction(&(*dynprogindex),&score,&miss_score,&nmatches0,&nmismatches0,
						 &nopens0,&nindels0,dynprog,revsequence1,revsequenceuc1,
						 /*revsequence2*/&(splicejunction[length2-1]),/*revsequenceuc2*/&(splicejunction[length2-1]),
						 /*revsequencealt2*/&(splicejunction_alt[length2-1]),
						 length1,length2,revoffset1,spliceoffset2_anchor,spliceoffset2_far,
						 chroffset,chrhigh,cdna_direction,watsonp,jump_late_p,pairpool,
						 extraband_end,defect_rate,contlength)) != NULL) {

	  /* miss_score = perfect_score - score; */
	  assert(miss_score <= 0);
	  debug7(printf("score %d - perfect score %d = miss %d expected vs %d returned.  ",
			score,perfect_score,score-perfect_score,miss_score));
	  debug7(printf("  comparing against threshold_miss %d + obsmax_penalty %d\n",*threshold_miss_score,obsmax_penalty));
	  if (score > 0 && miss_score > *threshold_miss_score + obsmax_penalty) {
	    debug7(printf("miss %d > threshold %d + %d",miss_score,*threshold_miss_score,obsmax_penalty));
#if 0
	    /* Just use results from Dynprog_end5_splicejunction */
	    pairs = Dynprog_add_known_splice_5(&length_distal,pairs,anchor_splicesite,splicecoord,chroffset,watsonp,pairpool);
#else
	    length_distal = length1 - contlength;
#endif
	    best_pairs = pairs;
	    *finalscore = score;
	    *nmatches = nmatches0;
	    *nmismatches = nmismatches0;
	    *nopens = nopens0;
	    *nindels = nindels0;
	    *knownsplicep = true;
	    *ambig_end_length = length_distal;
	    *threshold_miss_score = miss_score - obsmax_penalty;
	    shortest_intron_length = (splicecoord > anchor_splicesite) ? splicecoord - anchor_splicesite : anchor_splicesite - splicecoord;
	    *coordsptr = coords;
	    *(*coordsptr)++ = splicecoord;

	  } else if (miss_score == *threshold_miss_score + obsmax_penalty
#if 0
		     /* Filter for duplicates later */
		     && Genomicposlist_find(*coords,splicecoord) == false
#endif
		     ) {
	    if (amb_closest_p == false) {
	      debug7(printf("miss %d == threshold %d + %d, so ambiguous",miss_score,*threshold_miss_score,obsmax_penalty));
	      /* best_pairs = (List_T) NULL; */
	      *(*coordsptr)++ = splicecoord;
	    } else {
	      intron_length = (splicecoord > anchor_splicesite) ? splicecoord - anchor_splicesite : anchor_splicesite - splicecoord;
	      if (intron_length > shortest_intron_length) {
		debug7(printf("miss %d == threshold %d + %d, but intron_length %d > shortest %d, so ignore",
			      miss_score,*threshold_miss_score,obsmax_penalty,intron_length,shortest_intron_length));
	      } else {
		debug7(printf("miss %d == threshold %d + %d, but intron_length %d < shortest %d, so new best",
			      miss_score,*threshold_miss_score,obsmax_penalty,intron_length,shortest_intron_length));
#if 0
		/* Just use results from Dynprog_end5_splicejunction */
		pairs = Dynprog_add_known_splice_5(&length_distal,pairs,anchor_splicesite,splicecoord,chroffset,watsonp,pairpool);
#else
		length_distal = length1 - contlength;
#endif
		best_pairs = pairs;
		*finalscore = score;
		*nmatches = nmatches0;
		*nmismatches = nmismatches0;
		*nopens = nopens0;
		*nindels = nindels0;
		*knownsplicep = true;
		*ambig_end_length = length_distal;
		*threshold_miss_score = miss_score - obsmax_penalty;
		shortest_intron_length = intron_length;
		*coordsptr = coords;
		*(*coordsptr)++ = splicecoord;
	      }
	    }
	  }
	}
      }
      debug7(printf("\n"));

    }

  } else {
    offseta = (int) triestart[1];
    offsetc = (int) triestart[2];
    offsetg = (int) triestart[3];
    offsett = (int) triestart[4];

    if (offseta > 0) {
      best_pairs = solve_end5_aux(&(*coordsptr),coords,best_pairs,&(triestart[-offseta]),
				  knownsplice_limit_low,knownsplice_limit_high,
				  &(*finalscore),&(*nmatches),&(*nmismatches),
				  &(*nopens),&(*nindels),&(*knownsplicep),&(*ambig_end_length),
				  &(*threshold_miss_score),obsmax_penalty,perfect_score,
				  anchor_splicesite,splicejunction,splicejunction_alt,
				  splicelength,contlength,far_splicetype,
				  chroffset,chrhigh,&(*dynprogindex),dynprog,
				  revsequence1,revsequenceuc1,length1,length2,revoffset1,revoffset2,
				  cdna_direction,watsonp,jump_late_p,pairpool,extraband_end,defect_rate);
    }
    if (offsetc > 0) {
      best_pairs = solve_end5_aux(&(*coordsptr),coords,best_pairs,&(triestart[-offsetc]),
				  knownsplice_limit_low,knownsplice_limit_high,
				  &(*finalscore),&(*nmatches),&(*nmismatches),
				  &(*nopens),&(*nindels),&(*knownsplicep),&(*ambig_end_length),
				  &(*threshold_miss_score),obsmax_penalty,perfect_score,
				  anchor_splicesite,splicejunction,splicejunction_alt,
				  splicelength,contlength,far_splicetype,
				  chroffset,chrhigh,&(*dynprogindex),dynprog,
				  revsequence1,revsequenceuc1,length1,length2,revoffset1,revoffset2,
				  cdna_direction,watsonp,jump_late_p,pairpool,extraband_end,defect_rate);
    }
    if (offsetg > 0) {
      best_pairs = solve_end5_aux(&(*coordsptr),coords,best_pairs,&(triestart[-offsetg]),
				  knownsplice_limit_low,knownsplice_limit_high,
				  &(*finalscore),&(*nmatches),&(*nmismatches),
				  &(*nopens),&(*nindels),&(*knownsplicep),&(*ambig_end_length),
				  &(*threshold_miss_score),obsmax_penalty,perfect_score,
				  anchor_splicesite,splicejunction,splicejunction_alt,
				  splicelength,contlength,far_splicetype,
				  chroffset,chrhigh,&(*dynprogindex),dynprog,
				  revsequence1,revsequenceuc1,length1,length2,revoffset1,revoffset2,
				  cdna_direction,watsonp,jump_late_p,pairpool,extraband_end,defect_rate);
    }
    if (offsett > 0) {
      best_pairs = solve_end5_aux(&(*coordsptr),coords,best_pairs,&(triestart[-offsett]),
				  knownsplice_limit_low,knownsplice_limit_high,
				  &(*finalscore),&(*nmatches),&(*nmismatches),
				  &(*nopens),&(*nindels),&(*knownsplicep),&(*ambig_end_length),
				  &(*threshold_miss_score),obsmax_penalty,perfect_score,
				  anchor_splicesite,splicejunction,splicejunction_alt,
				  splicelength,contlength,far_splicetype,
				  chroffset,chrhigh,&(*dynprogindex),dynprog,
				  revsequence1,revsequenceuc1,length1,length2,revoffset1,revoffset2,
				  cdna_direction,watsonp,jump_late_p,pairpool,extraband_end,defect_rate);
    }
  }

  return best_pairs;
}


static List_T
solve_end3_aux (Univcoord_T **coordsptr, Univcoord_T *coords,
		List_T best_pairs, Triecontent_T *triestart,
		Univcoord_T knownsplice_limit_low, Univcoord_T knownsplice_limit_high,

		int *finalscore, int *nmatches, int *nmismatches,
		int *nopens, int *nindels, bool *knownsplicep, int *ambig_end_length,
		int *threshold_miss_score, int obsmax_penalty, int perfect_score,

		Univcoord_T anchor_splicesite, char *splicejunction, char *splicejunction_alt,
		int splicelength, int contlength, Splicetype_T far_splicetype,
		Univcoord_T chroffset, Univcoord_T chrhigh,
		int *dynprogindex, Dynprog_T dynprog, 
		char *sequence1, char *sequenceuc1,
		int length1, int length2, int offset1, int offset2,
		int cdna_direction, bool watsonp, bool jump_late_p, Pairpool_T pairpool,
		int extraband_end, double defect_rate) {
  Triecontent_T leaf;
  int nleaves, i;
  Univcoord_T splicecoord;
  Chrpos_T shortest_intron_length = -1U, intron_length;
  int offseta, offsetc, offsetg, offsett;
  int spliceoffset2_anchor, spliceoffset2_far;

  int score, miss_score, nmatches0, nmismatches0, nopens0, nindels0;
  int length_distal;
  List_T pairs;
    
  if (single_leaf_p(leaf = triestart[0])) {
    splicecoord = splicesites[leaf];
    debug7(printf("Checking leaf %d at %u: ",(int) leaf,splicecoord));
    if (splicecoord >= knownsplice_limit_low && splicecoord <= knownsplice_limit_high &&
	Dynprog_make_splicejunction_3(splicejunction,splicejunction_alt,splicecoord,splicelength,contlength,far_splicetype,watsonp) == true) {
      debug7(printf("intron length %d, ",splicecoord - anchor_splicesite));
      debug7(printf("length1 = %d, length2 = %d, chroffset = %u, splicecoord = %u\n",
		    length1,length2,chroffset,splicecoord));
      if (watsonp) {
	spliceoffset2_anchor = offset2;
	spliceoffset2_far = spliceoffset2_anchor - anchor_splicesite + splicecoord;
      } else {
	spliceoffset2_anchor = offset2;
	spliceoffset2_far = spliceoffset2_anchor + anchor_splicesite - splicecoord;
      }
      if ((pairs = Dynprog_end3_splicejunction(&(*dynprogindex),&score,&miss_score,&nmatches0,&nmismatches0,
					       &nopens0,&nindels0,dynprog,sequence1,sequenceuc1,
					       /*sequence2*/splicejunction,/*sequenceuc2*/splicejunction,
					       /*sequencealt2*/splicejunction_alt,
					       length1,length2,offset1,spliceoffset2_anchor,spliceoffset2_far,
					       chroffset,chrhigh,cdna_direction,watsonp,jump_late_p,pairpool,
					       extraband_end,defect_rate,contlength)) != NULL) {

	/* miss_score = perfect_score - score; */
	assert(miss_score <= 0);
	debug7(printf("score %d - perfect score %d = miss %d expected vs %d returned.  ",
		      score,perfect_score,score-perfect_score,miss_score));
	debug7(printf("  comparing against threshold_miss %d + obsmax_penalty %d\n",*threshold_miss_score,obsmax_penalty));
	if (score > 0 && miss_score > *threshold_miss_score + obsmax_penalty) {
	  debug7(printf("miss %d > threshold %d + %d",miss_score,*threshold_miss_score,obsmax_penalty));
#if 0
	  /* Just results of Dynprog_end3_splicejunction */
	  pairs = Dynprog_add_known_splice_3(&length_distal,pairs,anchor_splicesite,splicecoord,chroffset,watsonp,pairpool);
#else
	  length_distal = length1 - contlength;
#endif
	  best_pairs = pairs;
	  *finalscore = score;
	  *nmatches = nmatches0;
	  *nmismatches = nmismatches0;
	  *nopens = nopens0;
	  *nindels = nindels0;
	  *knownsplicep = true;
	  *ambig_end_length = length_distal;
	  *threshold_miss_score = miss_score - obsmax_penalty;
	  shortest_intron_length = (splicecoord > anchor_splicesite) ? splicecoord - anchor_splicesite : anchor_splicesite - splicecoord;
	  *coordsptr = coords;
	  *(*coordsptr)++ = splicecoord;

	} else if (miss_score == *threshold_miss_score + obsmax_penalty
#if 0
		   && Genomicposlist_find(*coords,splicecoord) == false
#endif
		   ) {
	  if (amb_closest_p == false) {
	    debug7(printf("miss %d == threshold %d + %d, so ambiguous",miss_score,*threshold_miss_score,obsmax_penalty));
	    /* best_pairs = (List_T) NULL; */
	    *(*coordsptr)++ = splicecoord;
	  } else {
	    intron_length = (splicecoord > anchor_splicesite) ? splicecoord - anchor_splicesite : anchor_splicesite - splicecoord;
	    if (intron_length > shortest_intron_length) {
	      debug7(printf("miss %d == threshold %d + %d, but intron_length %d > shortest %d, so ignore",
			    miss_score,*threshold_miss_score,obsmax_penalty,intron_length,shortest_intron_length));
	    } else {
	      debug7(printf("miss %d == threshold %d + %d, but intron_length %d < shortest %d, so new best",
			    miss_score,*threshold_miss_score,obsmax_penalty,intron_length,shortest_intron_length));
#if 0
	      /* Just use results of Dynprog_end3_splicejunction */
	      pairs = Dynprog_add_known_splice_3(&length_distal,pairs,anchor_splicesite,splicecoord,chroffset,watsonp,pairpool);
#else
	      length_distal = length1 - contlength;
#endif
	      best_pairs = pairs;
	      *finalscore = score;
	      *nmatches = nmatches0;
	      *nmismatches = nmismatches0;
	      *nopens = nopens0;
	      *nindels = nindels0;
	      *knownsplicep = true;
	      *ambig_end_length = length_distal;
	      *threshold_miss_score = miss_score - obsmax_penalty;
	      shortest_intron_length = intron_length;
	      *coordsptr = coords;
	      *(*coordsptr)++ = splicecoord;
	    }
	  }
	}
      }
    }
    debug7(printf("\n"));

  } else if (multiple_leaf_p(leaf)) {
    nleaves = (int) (-leaf);
    for (i = 1; i <= nleaves; i++) {
      leaf = triestart[i];
      splicecoord = splicesites[leaf];
      debug7(printf("Checking leaf %d at %u: ",(int) leaf,splicecoord));
      if (splicecoord >= knownsplice_limit_low && splicecoord <= knownsplice_limit_high &&
	  Dynprog_make_splicejunction_3(splicejunction,splicejunction_alt,splicecoord,splicelength,contlength,far_splicetype,watsonp) == true) {
	debug7(printf("intron length %d, ",splicecoord - anchor_splicesite));
	debug7(printf("length1 = %d, length2 = %d, chroffset = %u, splicecoord = %u\n",
		      length1,length2,chroffset,splicecoord));
	if (watsonp) {
	  spliceoffset2_anchor = offset2;
	  spliceoffset2_far = spliceoffset2_anchor - anchor_splicesite + splicecoord;
	} else {
	  spliceoffset2_anchor = offset2;
	  spliceoffset2_far = spliceoffset2_anchor + anchor_splicesite - splicecoord;
	}
	if ((pairs = Dynprog_end3_splicejunction(&(*dynprogindex),&score,&miss_score,&nmatches0,&nmismatches0,
						 &nopens0,&nindels0,dynprog,sequence1,sequenceuc1,
						 /*sequence2*/splicejunction,/*sequenceuc2*/splicejunction,
						 /*sequencealt2*/splicejunction_alt,
						 length1,length2,offset1,spliceoffset2_anchor,spliceoffset2_far,
						 chroffset,chrhigh,cdna_direction,watsonp,jump_late_p,pairpool,
						 extraband_end,defect_rate,contlength)) != NULL) {

	  /* miss_score = perfect_score - score; */
	  assert(miss_score <= 0);
	  debug7(printf("score %d - perfect score %d = miss %d expected vs %d returned.  ",
			score,perfect_score,score-perfect_score,miss_score));
	  debug7(printf("  comparing against threshold_miss %d + obsmax_penalty %d\n",*threshold_miss_score,obsmax_penalty));
	  if (score > 0 && miss_score > *threshold_miss_score + obsmax_penalty) {
	    debug7(printf("miss %d > threshold %d + %d",miss_score,*threshold_miss_score,obsmax_penalty));
#if 0
	    /* Just use results of Dynprog_end3_splicejunction */
	    pairs = Dynprog_add_known_splice_3(&length_distal,pairs,anchor_splicesite,splicecoord,chroffset,watsonp,pairpool);
#else
	    length_distal = length1 - contlength;
#endif
	    best_pairs = pairs;
	    *finalscore = score;
	    *nmatches = nmatches0;
	    *nmismatches = nmismatches0;
	    *nopens = nopens0;
	    *nindels = nindels0;
	    *knownsplicep = true;
	    *ambig_end_length = length_distal;
	    *threshold_miss_score = miss_score - obsmax_penalty;
	    shortest_intron_length = (splicecoord > anchor_splicesite) ? splicecoord - anchor_splicesite : anchor_splicesite - splicecoord;
	    *coordsptr = coords;
	    *(*coordsptr)++ = splicecoord;
	  } else if (miss_score == *threshold_miss_score + obsmax_penalty
#if 0
		     && Genomicposlist_find(*coords,splicecoord) == false
#endif
		     ) {
	    if (amb_closest_p == false) {
	      debug7(printf("miss %d == threshold %d + %d, so ambiguous",miss_score,*threshold_miss_score,obsmax_penalty));
	      /* best_pairs = (List_T) NULL; */
	      *(*coordsptr)++ = splicecoord;
	    } else {
	      intron_length = (splicecoord > anchor_splicesite) ? splicecoord - anchor_splicesite : anchor_splicesite - splicecoord;
	      if (intron_length > shortest_intron_length) {
		debug7(printf("miss %d == threshold %d + %d, but intron_length %d > shortest %d, so ignore",
			      miss_score,*threshold_miss_score,obsmax_penalty,intron_length,shortest_intron_length));
	      } else {
		debug7(printf("miss %d == threshold %d + %d, but intron_length %d < shortest %d, so new best",
			      miss_score,*threshold_miss_score,obsmax_penalty,intron_length,shortest_intron_length));
#if 0
		/* Just use results of Dynprog_end3_splicejunction */
		pairs = Dynprog_add_known_splice_3(&length_distal,pairs,anchor_splicesite,splicecoord,chroffset,watsonp,pairpool);
#else
		length_distal = length1 - contlength;
#endif
		best_pairs = pairs;
		*finalscore = score;
		*nmatches = nmatches0;
		*nmismatches = nmismatches0;
		*nopens = nopens0;
		*nindels = nindels0;
		*knownsplicep = true;
		*ambig_end_length = length_distal;
		*threshold_miss_score = miss_score - obsmax_penalty;
		shortest_intron_length = intron_length;
		*coordsptr = coords;
		*(*coordsptr)++ = splicecoord;
	      }
	    }
	  }
	}
      }
      debug7(printf("\n"));

    }

  } else {
    offseta = (int) triestart[1];
    offsetc = (int) triestart[2];
    offsetg = (int) triestart[3];
    offsett = (int) triestart[4];

    if (offseta > 0) {
      best_pairs = solve_end3_aux(&(*coordsptr),coords,best_pairs,&(triestart[-offseta]),
				  knownsplice_limit_low,knownsplice_limit_high,
				  &(*finalscore),&(*nmatches),&(*nmismatches),
				  &(*nopens),&(*nindels),&(*knownsplicep),&(*ambig_end_length),
				  &(*threshold_miss_score),obsmax_penalty,perfect_score,
				  anchor_splicesite,splicejunction,splicejunction_alt,
				  splicelength,contlength,far_splicetype,
				  chroffset,chrhigh,&(*dynprogindex),dynprog,
				  sequence1,sequenceuc1,length1,length2,offset1,offset2,
				  cdna_direction,watsonp,jump_late_p,pairpool,extraband_end,defect_rate);
    }
    if (offsetc > 0) {
      best_pairs = solve_end3_aux(&(*coordsptr),coords,best_pairs,&(triestart[-offsetc]),
				  knownsplice_limit_low,knownsplice_limit_high,
				  &(*finalscore),&(*nmatches),&(*nmismatches),
				  &(*nopens),&(*nindels),&(*knownsplicep),&(*ambig_end_length),
				  &(*threshold_miss_score),obsmax_penalty,perfect_score,
				  anchor_splicesite,splicejunction,splicejunction_alt,
				  splicelength,contlength,far_splicetype,
				  chroffset,chrhigh,&(*dynprogindex),dynprog,
				  sequence1,sequenceuc1,length1,length2,offset1,offset2,
				  cdna_direction,watsonp,jump_late_p,pairpool,extraband_end,defect_rate);
    }
    if (offsetg > 0) {
      best_pairs = solve_end3_aux(&(*coordsptr),coords,best_pairs,&(triestart[-offsetg]),
				  knownsplice_limit_low,knownsplice_limit_high,
				  &(*finalscore),&(*nmatches),&(*nmismatches),
				  &(*nopens),&(*nindels),&(*knownsplicep),&(*ambig_end_length),
				  &(*threshold_miss_score),obsmax_penalty,perfect_score,
				  anchor_splicesite,splicejunction,splicejunction_alt,
				  splicelength,contlength,far_splicetype,
				  chroffset,chrhigh,&(*dynprogindex),dynprog,
				  sequence1,sequenceuc1,length1,length2,offset1,offset2,
				  cdna_direction,watsonp,jump_late_p,pairpool,extraband_end,defect_rate);
    }
    if (offsett > 0) {
      best_pairs = solve_end3_aux(&(*coordsptr),coords,best_pairs,&(triestart[-offsett]),
				  knownsplice_limit_low,knownsplice_limit_high,
				  &(*finalscore),&(*nmatches),&(*nmismatches),
				  &(*nopens),&(*nindels),&(*knownsplicep),&(*ambig_end_length),
				  &(*threshold_miss_score),obsmax_penalty,perfect_score,
				  anchor_splicesite,splicejunction,splicejunction_alt,
				  splicelength,contlength,far_splicetype,
				  chroffset,chrhigh,&(*dynprogindex),dynprog,
				  sequence1,sequenceuc1,length1,length2,offset1,offset2,
				  cdna_direction,watsonp,jump_late_p,pairpool,extraband_end,defect_rate);
    }
  }

  return best_pairs;
}

List_T
Splicetrie_solve_end5 (List_T best_pairs, Triecontent_T *triecontents, Trieoffset_T *trieoffsets, int j,
		       Univcoord_T knownsplice_limit_low, Univcoord_T knownsplice_limit_high,

		       int *finalscore, int *nmatches, int *nmismatches,
		       int *nopens, int *nindels, bool *knownsplicep, int *ambig_end_length,
		       int *threshold_miss_score, int obsmax_penalty, int perfect_score,

		       Univcoord_T anchor_splicesite, char *splicejunction, char *splicejunction_alt,
		       int splicelength, int contlength, Splicetype_T far_splicetype,
		       Univcoord_T chroffset, Univcoord_T chrhigh,
		       int *dynprogindex, Dynprog_T dynprog, 
		       char *revsequence1, char *revsequenceuc1,
		       int length1, int length2, int revoffset1, int revoffset2,
		       int cdna_direction, bool watsonp, bool jump_late_p, Pairpool_T pairpool,
		       int extraband_end, double defect_rate) {
  Univcoord_T *coordsptr, *coords, splicecoord0;
  unsigned int *triestart;
  int size, i;

  debug7(printf("Entering Splicetrie_solve_end5 with limits %u..%u, anchor splicesite %u (%u)\n",
		knownsplice_limit_low,knownsplice_limit_high,anchor_splicesite,anchor_splicesite - chroffset));
  if (trieoffsets[j] == NULL_POINTER) {
    return best_pairs;
  } else {
    triestart = &(triecontents[trieoffsets[j]]);
  }

  if ((size = splicetrie_size(triestart)) == 0) {
    return best_pairs;
  } else {
    coordsptr = coords = (Univcoord_T *) MALLOCA(size * sizeof(Univcoord_T));

    best_pairs = solve_end5_aux(&coordsptr,coords,best_pairs,triestart,
				knownsplice_limit_low,knownsplice_limit_high,
				&(*finalscore),&(*nmatches),&(*nmismatches),
				&(*nopens),&(*nindels),&(*knownsplicep),&(*ambig_end_length),
				&(*threshold_miss_score),obsmax_penalty,perfect_score,
				anchor_splicesite,splicejunction,splicejunction_alt,
				splicelength,contlength,far_splicetype,
				chroffset,chrhigh,&(*dynprogindex),dynprog,
				revsequence1,revsequenceuc1,length1,length2,revoffset1,revoffset2,
				cdna_direction,watsonp,jump_late_p,pairpool,extraband_end,defect_rate);
    debug7(printf("\n"));

    /* Check for unique splicecoord */
    size = (int) (coordsptr - coords);
    if (size > 1) {
      splicecoord0 = coords[0];
      i = 1;
      while (i < size && coords[i] == splicecoord0) {
	i++;
      }
      if (i < size) {
	/* Signal non-uniqueness or ambiguity */
	best_pairs = (List_T) NULL;
      }
    }
    
    FREEA(coords);
    return best_pairs;
  }
}


List_T
Splicetrie_solve_end3 (List_T best_pairs, Triecontent_T *triecontents, Trieoffset_T *trieoffsets, int j,
		       Univcoord_T knownsplice_limit_low, Univcoord_T knownsplice_limit_high,

		       int *finalscore, int *nmatches, int *nmismatches,
		       int *nopens, int *nindels, bool *knownsplicep, int *ambig_end_length,
		       int *threshold_miss_score, int obsmax_penalty, int perfect_score,

		       Univcoord_T anchor_splicesite, char *splicejunction, char *splicejunction_alt,
		       int splicelength, int contlength, Splicetype_T far_splicetype,
		       Univcoord_T chroffset, Univcoord_T chrhigh,
		       int *dynprogindex, Dynprog_T dynprog, 
		       char *sequence1, char *sequenceuc1,
		       int length1, int length2, int offset1, int offset2,
		       int cdna_direction, bool watsonp, bool jump_late_p, Pairpool_T pairpool,
		       int extraband_end, double defect_rate) {
  Univcoord_T *coordsptr, *coords, splicecoord0;
  unsigned int *triestart;
  int size, i;

  debug7(printf("Entering Splicetrie_solve_end3 with limits %u..%u, anchor splicesite %u (%u)\n",
		knownsplice_limit_low,knownsplice_limit_high,anchor_splicesite,anchor_splicesite - chroffset));
  if (trieoffsets[j] == NULL_POINTER) {
    return best_pairs;
  } else {
    triestart = &(triecontents[trieoffsets[j]]);
  }

  if ((size = splicetrie_size(triestart)) == 0) {
    return best_pairs;
  } else {
    coordsptr = coords = (Univcoord_T *) MALLOCA(size * sizeof(Univcoord_T));

    best_pairs = solve_end3_aux(&coordsptr,coords,best_pairs,triestart,
				knownsplice_limit_low,knownsplice_limit_high,
				&(*finalscore),&(*nmatches),&(*nmismatches),
				&(*nopens),&(*nindels),&(*knownsplicep),&(*ambig_end_length),
				&(*threshold_miss_score),obsmax_penalty,perfect_score,
				anchor_splicesite,splicejunction,splicejunction_alt,
				splicelength,contlength,far_splicetype,
				chroffset,chrhigh,&(*dynprogindex),dynprog,
				sequence1,sequenceuc1,length1,length2,offset1,offset2,
				cdna_direction,watsonp,jump_late_p,pairpool,extraband_end,defect_rate);
    debug7(printf("\n"));

    /* Check for unique splicecoord */
    size = (int) (coordsptr - coords);
    if (size > 1) {
      splicecoord0 = coords[0];
      i = 1;
      while (i < size && coords[i] == splicecoord0) {
	i++;
      }
      if (i < size) {
	/* Signal non-uniqueness or ambiguity */
	best_pairs = (List_T) NULL;
      }
    }

    FREEA(coords);
    return best_pairs;
  }
}



#ifdef GSNAP
#ifdef END_KNOWNSPLICING_SHORTCUT

static Genomicposlist_T
dump_left_aux (int *best_nmismatches, Genomicposlist_T coords, Triecontent_T *triecontents,
	       char *queryptr, Compress_T query_compress,
	       int pos5, int pos3, bool plusp, int nmismatches, int charpos,
	       Univcoord_T knownsplice_limit_low, Univcoord_T knownsplice_limit_high) {
  Triecontent_T leaf;
  Univcoord_T segment_left;
  int nleaves, i;
  Univcoord_T position;
  int offseta, offsetc, offsetg, offsett;
  char c;

  if (nmismatches > *best_nmismatches) {
    return coords;

  } else if (single_leaf_p(leaf = triecontents[0])) {
    position = splicesites[leaf];
    if (position < knownsplice_limit_low || position > knownsplice_limit_high) {
      debug3(printf("Found leaf %u at %u, but outside knownsplice limits %u..%u\n",
		    leaf,position,knownsplice_limit_low,knownsplice_limit_high));
      return coords;
    } else if (charpos - 1 >= pos5) {
      if (pos3 - pos5 <= 16) {
	/* Recomputes entire segment to determine mismatches (necessary because of splicefrags) */
	nmismatches = Genome_count_mismatches_fragment_left(query_compress,pos5,pos3,
							    splicefrags_ref[leaf],splicefrags_alt[leaf]);
	debug3(printf("Found leaf %u at %u, but still have characters to check against Genome: %.*s => %d mismatches by fragment\n",
		      leaf,position,charpos + 1,&(queryptr[0]),nmismatches));
      } else {
	/* Can happen */
	segment_left = position - pos3;
	nmismatches =
	  Genome_count_mismatches_substring(query_compress,segment_left,pos5,pos3,plusp);
	debug3(printf("Found leaf %u at %u, but still have characters to check against Genome: %.*s => %d mismatches by substring\n",
		      leaf,position,pos5 - pos3,&(queryptr[pos5]),nmismatches));
      }

      if (nmismatches < *best_nmismatches) {
	debug3(printf("  nmismatches %d < best_nmismatches %d => setting splicesites_i to be %d (new best_nmismatches %d)\n",
		      nmismatches,*best_nmismatches,(int) leaf,nmismatches));
	*best_nmismatches = nmismatches;
	Genomicposlist_free(&coords);
	return Genomicposlist_push(NULL,position);
      } else if (nmismatches == *best_nmismatches) {
	debug3(printf("  nmismatches %d == best_nmismatches %d => pushing leaf %d (same best_nmismatches %d)\n",
		      nmismatches,*best_nmismatches,(int) leaf,nmismatches));
	return Genomicposlist_push(coords,position);
      } else {
	return coords;
      }

    } else {
      debug3(printf("Found leaf %u at %u, and completely finished query, so unique\n",leaf,position));
      if (nmismatches < *best_nmismatches) {
	debug3(printf("  nmismatches %d < best_nmismatches %d => setting splicesites_i to be %d (new best_nmismatches %d)\n",
		      nmismatches,*best_nmismatches,(int) leaf,nmismatches));
	*best_nmismatches = nmismatches;
	Genomicposlist_free(&coords);
	return Genomicposlist_push(NULL,position);
      } else if (nmismatches == *best_nmismatches) {
	debug3(printf("  nmismatches %d == best_nmismatches %d => pushing leaf %d (same best_nmismatches %d)\n",
		      nmismatches,*best_nmismatches,(int) leaf,nmismatches));
	return Genomicposlist_push(coords,position);
      } else {
	return coords;
      }
    }

  } else if (multiple_leaf_p(leaf)) {
    nleaves = (int) (-leaf);
    for (i = 1; i <= nleaves; i++) {
      leaf = triecontents[i];
      position = splicesites[leaf];
      if (position < knownsplice_limit_low || position > knownsplice_limit_high) {
	/* Skip */
      } else if (charpos - 1 >= pos5) {
	if (pos3 - pos5 <= 16) {
	  /* Recomputes entire segment to determine mismatches (necessary because of splicefrags) */
	  nmismatches = Genome_count_mismatches_fragment_left(query_compress,pos5,pos3,
							      splicefrags_ref[leaf],splicefrags_alt[leaf]);
	} else {
	  /* Can happen */
	  segment_left = position - pos3;
	  nmismatches =
	    Genome_count_mismatches_substring(query_compress,segment_left,pos5,pos3,plusp);
	}

	debug3(printf("Found leaf %u at %u => %d mismatches\n",leaf,position,nmismatches));
	if (nmismatches < *best_nmismatches) {
	  debug3(printf("  nmismatches %d < best_nmismatches %d => setting splicesites_i to be %d (new best_nmismatches %d)\n",
			nmismatches,*best_nmismatches,(int) leaf,nmismatches));
	  *best_nmismatches = nmismatches;
	  Genomicposlist_free(&coords);
	  coords = Genomicposlist_push(NULL,position);
	} else if (nmismatches == *best_nmismatches) {
	  debug3(printf("  nmismatches %d == best_nmismatches %d => pushing leaf %d (same best_nmismatches %d)\n",
			nmismatches,*best_nmismatches,(int) leaf,nmismatches));
	  coords = Genomicposlist_push(coords,position);
	}
      
      } else {
	debug3(printf("Found leaf %u at %u, and completely finished query, so unique\n",leaf,position));
	if (nmismatches < *best_nmismatches) {
	  debug3(printf("  nmismatches %d < best_nmismatches %d => setting splicesites_i to be %d (new best_nmismatches %d)\n",
			nmismatches,*best_nmismatches,(int) leaf,nmismatches));
	  *best_nmismatches = nmismatches;
	  Genomicposlist_free(&coords);
	  coords = Genomicposlist_push(NULL,position);
	} else if (nmismatches == *best_nmismatches) {
	  debug3(printf("  nmismatches %d == best_nmismatches %d => pushing leaf %d (same best_nmismatches %d)\n",
			nmismatches,*best_nmismatches,(int) leaf,nmismatches));
	  coords = Genomicposlist_push(coords,position);
	}
      }
    }

    return coords;

  } else {
#ifdef USE_2BYTE_RELOFFSETS
    get_offsets(&offseta,&offsetc,&offsetg,&offsett,triecontents[1],triecontents[2]);
#else
    offseta = (int) triecontents[1];
    offsetc = (int) triecontents[2];
    offsetg = (int) triecontents[3];
    offsett = (int) triecontents[4];
#endif

    if (charpos - 1 >= pos5) {
      c = queryptr[charpos - 1];
      debug3(printf("Character at pos %d is %c\n",charpos - 1,c));
    } else {
      c = '\0';
    }

    if (offseta > 0) {
      coords = dump_left_aux(&(*best_nmismatches),coords,&(triecontents[-offseta]),
			     queryptr,query_compress,
			     pos5,pos3,plusp,nmismatches+(c != '\0' && c != 'A'),charpos-1,
			     knownsplice_limit_low,knownsplice_limit_high);
    }
    if (offsetc > 0) {
      coords = dump_left_aux(&(*best_nmismatches),coords,&(triecontents[-offsetc]),
			     queryptr,query_compress,
			     pos5,pos3,plusp,nmismatches+(c != '\0' && c != 'C'),charpos-1,
			     knownsplice_limit_low,knownsplice_limit_high);
    }
    if (offsetg > 0) {
      coords = dump_left_aux(&(*best_nmismatches),coords,&(triecontents[-offsetg]),
			     queryptr,query_compress,
			     pos5,pos3,plusp,nmismatches+(c != '\0' && c != 'G'),charpos-1,
			     knownsplice_limit_low,knownsplice_limit_high);
    }
    if (offsett > 0) {
      coords = dump_left_aux(&(*best_nmismatches),coords,&(triecontents[-offsett]),
			     queryptr,query_compress,
			     pos5,pos3,plusp,nmismatches+(c != '\0' && c != 'T'),charpos-1,
			     knownsplice_limit_low,knownsplice_limit_high);
    }

    return coords;
  }
}


Genomicposlist_T
Splicetrie_dump_coords_left (int *best_nmismatches, Triecontent_T *triestart, int pos5, int pos3,
			     Compress_T query_compress, char *queryptr, bool plusp,
			     Univcoord_T knownsplice_limit_low, Univcoord_T knownsplice_limit_high) {
  Genomicposlist_T coords;
  Univcoord_T *array, lastvalue;
  int n, k;

  debug3(printf("Splicetrie_dump_coords_left called at pos5 %d, pos3 %d with query %s\n",
		pos5,pos3,queryptr));
  debug3(printf("  knownsplice limits are %u..%u\n",knownsplice_limit_low,knownsplice_limit_high));

  if (pos5 >= pos3) {
    debug3(printf("Splicetrie_dump_coords_left returning NULL, because pos5 >= pos3\n"));
    return (Genomicposlist_T) NULL;
  } else {
    *best_nmismatches = pos3 - pos5;
  }

  coords = dump_left_aux(&(*best_nmismatches),/*coords*/NULL,
			 triestart,queryptr,query_compress,
			 pos5,pos3,plusp,/*nmismatches*/0,/*charpos*/pos3,
			 knownsplice_limit_low,knownsplice_limit_high);

  assert(*best_nmismatches >= 0);

  if ((n = Genomicposlist_length(coords)) > 1 && snpp == true) {
    array = (Univcoord_T *) MALLOCA(n * sizeof(Univcoord_T));
    Genomicposlist_fill_array_and_free(array,&coords);
    qsort(array,n,sizeof(Univcoord_T),Univcoord_compare);

    coords = (Genomicposlist_T) NULL;
    coords = Genomicposlist_push(coords,array[0]);
    lastvalue = array[0];
    for (k = 1; k < n; k++) {
      if (array[k] == lastvalue) {
	/* Skip */
      } else {
	coords = Genomicposlist_push(coords,array[k]);
	lastvalue = array[k];
      }
    }
    FREEA(array);
  }

  debug3(printf("Splicetrie_dump_left returning %s\n",Genomicposlist_to_string(coords)));
  return coords;
}




static Genomicposlist_T
dump_right_aux (int *best_nmismatches, Genomicposlist_T coords,
		Triecontent_T *triecontents, char *queryptr, Compress_T query_compress,
		int pos5, int pos3, bool plusp, int nmismatches, int charpos,
		Univcoord_T knownsplice_limit_low, Univcoord_T knownsplice_limit_high) {
  Triecontent_T leaf;
  Univcoord_T segment_left;
  int nleaves, i;
  Univcoord_T position;
  int offseta, offsetc, offsetg, offsett;
  char c;

  if (nmismatches > *best_nmismatches) {
    return coords;

  } else if (single_leaf_p(leaf = triecontents[0])) {
    position = splicesites[leaf];
    if (position < knownsplice_limit_low || position > knownsplice_limit_high) {
      debug3(printf("Found leaf %u at %u, but outside knownsplice limits %u..%u\n",
		    leaf,position,knownsplice_limit_low,knownsplice_limit_high));
      return coords;
    } else if (charpos + 1 < pos3) {
      if (pos3 - pos5 <= 16) {
	/* Recomputes entire segment to determine mismatches (necessary because of splicefrags) */
	nmismatches = Genome_count_mismatches_fragment_right(query_compress,pos5,pos3,
							     splicefrags_ref[leaf],splicefrags_alt[leaf]);
	debug3(printf("Found leaf %u at %u, but still have characters to check against Genome: %.*s => %d mismatches by fragment\n",
		      leaf,position,charpos + 1,&(queryptr[0]),nmismatches));
      } else {
	/* Can happen */
	segment_left = position - pos5;
	nmismatches =
	  Genome_count_mismatches_substring(query_compress,segment_left,pos5,pos3,plusp);
	debug3(printf("Found leaf %u at %u, but still have characters to check against Genome: %.*s => %d mismatches by substring\n",
		      leaf,position,pos5 - pos3,&(queryptr[pos5]),nmismatches));
      }

      if (nmismatches < *best_nmismatches) {
	debug3(printf("  nmismatches %d < best_nmismatches %d => setting splicesites_i to be %d (new best_nmismatches %d)\n",
		      nmismatches,*best_nmismatches,(int) leaf,nmismatches));
	*best_nmismatches = nmismatches;
	Genomicposlist_free(&coords);
	return Genomicposlist_push(NULL,position);
      } else if (nmismatches == *best_nmismatches) {
	debug3(printf("  nmismatches %d == best_nmismatches %d => pushing leaf %d (same best_nmismatches %d)\n",
		      nmismatches,*best_nmismatches,(int) leaf,nmismatches));
	return Genomicposlist_push(coords,position);
      } else {
	return coords;
      }

    } else {
      debug3(printf("Found leaf %u at %u, and completely finished query, so unique\n",leaf,position));
      if (nmismatches < *best_nmismatches) {
	debug3(printf("  nmismatches %d < best_nmismatches %d => setting splicesites_i to be %d (new best_nmismatches %d)\n",
		      nmismatches,*best_nmismatches,(int) leaf,nmismatches));
	*best_nmismatches = nmismatches;
	Genomicposlist_free(&coords);
	return Genomicposlist_push(NULL,position);
      } else if (nmismatches == *best_nmismatches) {
	debug3(printf("  nmismatches %d == best_nmismatches %d => pushing leaf %d (same best_nmismatches %d)\n",
		      nmismatches,*best_nmismatches,(int) leaf,nmismatches));
	return Genomicposlist_push(coords,position);
      } else {
	return coords;
      }
    }

  } else if (multiple_leaf_p(leaf)) {
    nleaves = (int) (-leaf);
    for (i = 1; i <= nleaves; i++) {
      leaf = triecontents[i];
      position = splicesites[leaf];
      if (position < knownsplice_limit_low || position > knownsplice_limit_high) {
	/* Skip */
      } else if (charpos + 1 < pos3) {
	if (pos3 - pos5 <= 16) {
	  /* Recomputes entire segment to determine mismatches (necessary because of splicefrags) */
	  nmismatches = Genome_count_mismatches_fragment_right(query_compress,pos5,pos3,
							       splicefrags_ref[leaf],splicefrags_alt[leaf]);
	} else {
	  /* Can happen */
	  segment_left = position - pos5;
	  nmismatches =
	    Genome_count_mismatches_substring(query_compress,segment_left,pos5,pos3,plusp);
	}

	debug3(printf("Found leaf %u at %u => %d mismatches\n",leaf,position,nmismatches));
	if (nmismatches < *best_nmismatches) {
	  debug3(printf("  nmismatches %d < best_nmismatches %d => setting splicesites_i to be %d (new best_nmismatches %d)\n",
			nmismatches,*best_nmismatches,(int) leaf,nmismatches));
	  *best_nmismatches = nmismatches;
	  Genomicposlist_free(&coords);
	  coords = Genomicposlist_push(NULL,position);
	} else if (nmismatches == *best_nmismatches) {
	  debug3(printf("  nmismatches %d == best_nmismatches %d => pushing leaf %d (same best_nmismatches %d)\n",
			nmismatches,*best_nmismatches,(int) leaf,nmismatches));
	  coords = Genomicposlist_push(coords,position);
	}
      
      } else {
	debug3(printf("Found leaf %u at %u, and completely finished query, so unique\n",leaf,position));
	if (nmismatches < *best_nmismatches) {
	  debug3(printf("  nmismatches %d < best_nmismatches %d => setting splicesites_i to be %d (new best_nmismatches %d)\n",
			nmismatches,*best_nmismatches,(int) leaf,nmismatches));
	  *best_nmismatches = nmismatches;
	  Genomicposlist_free(&coords);
	  coords = Genomicposlist_push(NULL,position);
	} else if (nmismatches == *best_nmismatches) {
	  debug3(printf("  nmismatches %d == best_nmismatches %d => pushing leaf %d (same best_nmismatches %d)\n",
			nmismatches,*best_nmismatches,(int) leaf,nmismatches));
	  coords = Genomicposlist_push(coords,position);
	}
      }
    }

    return coords;

  } else {
#ifdef USE_2BYTE_RELOFFSETS
    get_offsets(&offseta,&offsetc,&offsetg,&offsett,triecontents[1],triecontents[2]);
#else
    offseta = (int) triecontents[1];
    offsetc = (int) triecontents[2];
    offsetg = (int) triecontents[3];
    offsett = (int) triecontents[4];
#endif

    if (charpos + 1 < pos3) {
      c = queryptr[charpos + 1];
      debug3(printf("Character at pos %d is %c\n",charpos + 1,c));
    } else {
      c = '\0';
    }

    if (offseta > 0) {
      coords = dump_right_aux(&(*best_nmismatches),coords,&(triecontents[-offseta]),
			      queryptr,query_compress,
			      pos5,pos3,plusp,nmismatches+(c != '\0' && c != 'A'),charpos+1,
			      knownsplice_limit_low,knownsplice_limit_high);
    }
    if (offsetc > 0) {
      coords = dump_right_aux(&(*best_nmismatches),coords,&(triecontents[-offsetc]),
			      queryptr,query_compress,
			      pos5,pos3,plusp,nmismatches+(c != '\0' && c != 'C'),charpos+1,
			      knownsplice_limit_low,knownsplice_limit_high);
    }
    if (offsetg > 0) {
      coords = dump_right_aux(&(*best_nmismatches),coords,&(triecontents[-offsetg]),
			      queryptr,query_compress,
			      pos5,pos3,plusp,nmismatches+(c != '\0' && c != 'G'),charpos+1,
			      knownsplice_limit_low,knownsplice_limit_high);
    }
    if (offsett > 0) {
      coords = dump_right_aux(&(*best_nmismatches),coords,&(triecontents[-offsett]),
			      queryptr,query_compress,
			      pos5,pos3,plusp,nmismatches+(c != '\0' && c != 'T'),charpos+1,
			      knownsplice_limit_low,knownsplice_limit_high);
    }

    return coords;
  }
}


Genomicposlist_T
Splicetrie_dump_coords_right (int *best_nmismatches, Triecontent_T *triestart, int pos5, int pos3,
			      Compress_T query_compress, char *queryptr, bool plusp,
			      Univcoord_T knownsplice_limit_low, Univcoord_T knownsplice_limit_high) {
  Genomicposlist_T coords;
  unsigned int *array, lastvalue;
  int n, k;

  debug3(printf("Splicetrie_dump_coords_right called at pos5 %d, pos3 %d with query %s\n",
		pos5,pos3,queryptr));
  debug3(printf("  knownsplice limits are %u..%u\n",knownsplice_limit_low,knownsplice_limit_high));

  if (pos5 >= pos3) {
    debug3(printf("Splicetrie_dump_coords_right returning NULL, because pos5 >= pos3\n"));
    return (Genomicposlist_T) NULL;
  } else {
    *best_nmismatches = pos3 - pos5;
  }

  coords = dump_right_aux(&(*best_nmismatches),/*coords*/NULL,
			  triestart,queryptr,query_compress,
			  pos5,pos3,plusp,/*nmismatches*/0,/*charpos*/pos5-1,
			  knownsplice_limit_low,knownsplice_limit_high);

  assert(*best_nmismatches >= 0);

  if ((n = Genomicposlist_length(coords)) > 1 && snpp == true) {
    array = (Univcoord_T *) MALLOCA(n * sizeof(Univcoord_T));
    Genomicposlist_fill_array_and_free(array,&coords);
    qsort(array,n,sizeof(Univcoord_T),Univcoord_compare);

    coords = (Genomicposlist_T) NULL;
    coords = Genomicposlist_push(coords,array[0]);
    lastvalue = array[0];
    for (k = 1; k < n; k++) {
      if (array[k] == lastvalue) {
	/* Skip */
      } else {
	coords = Genomicposlist_push(coords,array[k]);
	lastvalue = array[k];
      }
    }
    FREEA(array);
  }

  debug3(printf("Splicetrie_dump_right returning %s\n",Genomicposlist_to_string(coords)));
  return coords;
}

#endif
#endif



/************************************************************************
 *   General search procedures
 ************************************************************************/

#ifdef GSNAP

static Intlist_T
search_left (int *best_nmismatches, Intlist_T *nmismatches_list, Intlist_T splicesites_i,
	     Triecontent_T *triecontents, char *queryptr, Compress_T query_compress,
	     int pos5, int pos3, bool plusp, int genestrand, bool first_read_p, bool collect_all_p,
	     int max_mismatches_allowed, int nmismatches, int charpos, Univcoord_T limit_low) {
  Triecontent_T leaf;
  Univcoord_T segment_left, position;
  int nleaves, i;
  int offseta, offsetc, offsetg, offsett;
  char c;
  
#if 0
  debug2(printf("Entered search_left with nmismatches %d, charpos %d, best_nmismatches %d\n",
		nmismatches,charpos,*best_nmismatches));
#endif

  if (nmismatches > *best_nmismatches && !collect_all_p) {
    return splicesites_i;

  } else if (nmismatches > max_mismatches_allowed) {
    return splicesites_i;

  } else if (single_leaf_p(leaf = triecontents[0])) {
    position = splicesites[leaf];
    if (position <= limit_low) {
      debug2(printf("Found leaf %u at %u, but below limit_low %u\n",leaf,position,limit_low));
      return splicesites_i;

    } else if (charpos - 1 >= pos5) {
      if (pos3 - pos5 <= 16) {
	/* Recomputes entire segment to determine mismatches (necessary because of splicefrags) */
	nmismatches = Genome_count_mismatches_fragment_left(query_compress,pos5,pos3,
							    splicefrags_ref[leaf],splicefrags_alt[leaf]);
      } else {
	/* Can happen in search for short middle exon */
	segment_left = splicesites[leaf] - pos3;
	nmismatches =
	  Genome_count_mismatches_substring(query_compress,segment_left,pos5,pos3,plusp,genestrand,first_read_p);
      }

      debug2(printf("Found leaf %u at %u, but still have characters to check against Genome: %.*s => %d mismatches\n",
		    leaf,position,charpos + 1,&(queryptr[0]),nmismatches));
      if (nmismatches < *best_nmismatches) {
	debug2(printf("  nmismatches %d < best_nmismatches %d => setting splicesites_i to be %d (new best_nmismatches %d)\n",
		      nmismatches,*best_nmismatches,(int) leaf,nmismatches));
	*best_nmismatches = nmismatches;
	Intlist_free(&(*nmismatches_list));
	*nmismatches_list = Intlist_push(NULL,nmismatches);
	Intlist_free(&splicesites_i);
	return Intlist_push(NULL,(int) leaf);
      } else if (nmismatches == *best_nmismatches) {
	debug2(printf("  nmismatches %d == best_nmismatches %d => pushing leaf %d (same best_nmismatches %d)\n",
		      nmismatches,*best_nmismatches,(int) leaf,nmismatches));
	*nmismatches_list = Intlist_push(*nmismatches_list,nmismatches);
	return Intlist_push(splicesites_i,(int) leaf);
      } else if (nmismatches <= max_mismatches_allowed && collect_all_p) {
	debug2(printf("  nmismatches %d > best_nmismatches %d, but collecting all => pushing leaf %d\n",
		      nmismatches,*best_nmismatches,(int) leaf));
	*nmismatches_list = Intlist_push(*nmismatches_list,nmismatches);
	return Intlist_push(splicesites_i,(int) leaf);
      } else {
	return splicesites_i;
      }
      
    } else {
      debug2(printf("Found leaf %u at %u, and completely finished query, so unique\n",leaf,position));
      assert(nmismatches <= max_mismatches_allowed);
      if (nmismatches < *best_nmismatches) {
	debug2(printf("  nmismatches %d < best_nmismatches %d => setting splicesites_i to be %d (new best_nmismatches %d)\n",
		      nmismatches,*best_nmismatches,(int) leaf,nmismatches));
	*best_nmismatches = nmismatches;
	Intlist_free(&(*nmismatches_list));
	*nmismatches_list = Intlist_push(NULL,nmismatches);
	Intlist_free(&splicesites_i);
	return Intlist_push(NULL,(int) leaf);
      } else if (nmismatches == *best_nmismatches) {
	debug2(printf("  nmismatches %d == best_nmismatches %d => pushing leaf %d (same best_nmismatches %d)\n",
		      nmismatches,*best_nmismatches,(int) leaf,nmismatches));
	*nmismatches_list = Intlist_push(*nmismatches_list,nmismatches);
	return Intlist_push(splicesites_i,(int) leaf);
      } else if (/*nmismatches <= max_mismatches_allowed &&*/ collect_all_p) {
	debug2(printf("  nmismatches %d > best_nmismatches %d, but collecting all => pushing leaf %d\n",
		      nmismatches,*best_nmismatches,(int) leaf));
	*nmismatches_list = Intlist_push(*nmismatches_list,nmismatches);
	return Intlist_push(splicesites_i,(int) leaf);
      } else {
	return splicesites_i;
      }
    }

  } else if (multiple_leaf_p(leaf)) {
    nleaves = (int) (-leaf);

    /* Need to compute nmismatches multiple times if we have SNPs */
    for (i = 1; i <= nleaves; i++) {
      leaf = triecontents[i];
      position = splicesites[leaf];
      if (position <= limit_low) {
	/* Skip */

      } else if (charpos - 1 >= pos5) {
	if (pos3 - pos5 <= 16) {
	  /* Recomputes entire segment to determine mismatches (necessary because of splicefrags) */
	  nmismatches = Genome_count_mismatches_fragment_left(query_compress,pos5,pos3,
							      splicefrags_ref[leaf],splicefrags_alt[leaf]);
	} else {
	  /* Can happen in search for short middle exon */
	  segment_left = splicesites[leaf] - pos3;
	  nmismatches =
	    Genome_count_mismatches_substring(query_compress,segment_left,pos5,pos3,plusp,genestrand,first_read_p);
	}

	debug2(printf("Found leaf %u at %u => %d mismatches\n",leaf,position,nmismatches));
	if (nmismatches < *best_nmismatches) {
	  debug2(printf("  nmismatches %d < best_nmismatches %d => setting splicesites_i to be %d (new best_nmismatches %d)\n",
			nmismatches,*best_nmismatches,(int) leaf,nmismatches));
	  *best_nmismatches = nmismatches;
	  Intlist_free(&(*nmismatches_list));
	  *nmismatches_list = Intlist_push(NULL,nmismatches);
	  Intlist_free(&splicesites_i);
	  splicesites_i = Intlist_push(NULL,(int) leaf);
	} else if (nmismatches == *best_nmismatches) {
	  debug2(printf("  nmismatches %d == best_nmismatches %d => pushing leaf %d (same best_nmismatches %d)\n",
			nmismatches,*best_nmismatches,(int) leaf,nmismatches));
	  *nmismatches_list = Intlist_push(*nmismatches_list,nmismatches);
	  splicesites_i = Intlist_push(splicesites_i,(int) leaf);
	} else if (nmismatches <= max_mismatches_allowed && collect_all_p) {
	  debug2(printf("  nmismatches %d > best_nmismatches %d, but collecting all => pushing leaf %d\n",
			nmismatches,*best_nmismatches,(int) leaf));
	  *nmismatches_list = Intlist_push(*nmismatches_list,nmismatches);
	  splicesites_i = Intlist_push(splicesites_i,(int) leaf);
	}
      
      } else {
	debug2(printf("Found leaf %u at %u, and completely finished query, so unique\n",leaf,position));
	assert(nmismatches <= max_mismatches_allowed);
	if (nmismatches < *best_nmismatches) {
	  debug2(printf("  nmismatches %d < best_nmismatches %d => setting splicesites_i to be %d (new best_nmismatches %d)\n",
			nmismatches,*best_nmismatches,(int) leaf,nmismatches));
	  *best_nmismatches = nmismatches;
	  Intlist_free(&(*nmismatches_list));
	  *nmismatches_list = Intlist_push(NULL,nmismatches);
	  Intlist_free(&splicesites_i);
	  splicesites_i = Intlist_push(NULL,(int) leaf);
	} else if (nmismatches == *best_nmismatches) {
	  debug2(printf("  nmismatches %d == best_nmismatches %d => pushing leaf %d (same best_nmismatches %d)\n",
			nmismatches,*best_nmismatches,(int) leaf,nmismatches));
	  *nmismatches_list = Intlist_push(*nmismatches_list,nmismatches);
	  splicesites_i = Intlist_push(splicesites_i,(int) leaf);
	} else if (/*nmismatches <= max_mismatches_allowed &&*/ collect_all_p) {
	  debug2(printf("  nmismatches %d > best_nmismatches %d, but collecting all => pushing leaf %d\n",
			nmismatches,*best_nmismatches,(int) leaf));
	  *nmismatches_list = Intlist_push(*nmismatches_list,nmismatches);
	  splicesites_i = Intlist_push(splicesites_i,(int) leaf);
	}
      }
    }

    return splicesites_i;

  } else {
    /* Non-leaf, and characters left, so recurse */
#ifdef USE_2BYTE_RELOFFSETS
    get_offsets(&offseta,&offsetc,&offsetg,&offsett,triecontents[1],triecontents[2]);
#else
    offseta = (int) triecontents[1];
    offsetc = (int) triecontents[2];
    offsetg = (int) triecontents[3];
    offsett = (int) triecontents[4];
#endif

    if (charpos - 1 >= pos5) {
      c = queryptr[charpos - 1];
      debug2(printf("Character at pos %d is %c\n",charpos - 1,c));
    } else {
      c = '\0';
    }

    if (offseta > 0) {
      splicesites_i = search_left(&(*best_nmismatches),&(*nmismatches_list),splicesites_i,
				  &(triecontents[-offseta]),queryptr,query_compress,
				  pos5,pos3,plusp,genestrand,first_read_p,collect_all_p,max_mismatches_allowed,
				  nmismatches+(c != '\0' && c != 'A'),charpos-1,limit_low);
    }
      
    if (offsetc > 0) {
      splicesites_i = search_left(&(*best_nmismatches),&(*nmismatches_list),splicesites_i,
				  &(triecontents[-offsetc]),queryptr,query_compress,
				  pos5,pos3,plusp,genestrand,first_read_p,collect_all_p,max_mismatches_allowed,
				  nmismatches+(c != '\0' && c != 'C'),charpos-1,limit_low);
    }

    if (offsetg > 0) {
      splicesites_i = search_left(&(*best_nmismatches),&(*nmismatches_list),splicesites_i,
				  &(triecontents[-offsetg]),queryptr,query_compress,
				  pos5,pos3,plusp,genestrand,first_read_p,collect_all_p,max_mismatches_allowed,
				  nmismatches+(c != '\0' && c != 'G'),charpos-1,limit_low);
    }

    if (offsett > 0) {
      splicesites_i = search_left(&(*best_nmismatches),&(*nmismatches_list),splicesites_i,
				  &(triecontents[-offsett]),queryptr,query_compress,
				  pos5,pos3,plusp,genestrand,first_read_p,collect_all_p,max_mismatches_allowed,
				  nmismatches+(c != '\0' && c != 'T'),charpos-1,limit_low);
    }

    return splicesites_i;
  }
}


static Intlist_T
search_right (int *best_nmismatches, Intlist_T *nmismatches_list, Intlist_T splicesites_i,
	      Triecontent_T *triecontents, char *queryptr, Compress_T query_compress,
	      int pos5, int pos3, bool plusp, int genestrand, bool first_read_p, bool collect_all_p,
	      int max_mismatches_allowed, int nmismatches, int charpos, Univcoord_T limit_high) {
  Triecontent_T leaf;
  Univcoord_T segment_left, position;
  int nleaves, i;
  int offseta, offsetc, offsetg, offsett;
  char c;
  
#if 0
  debug2(printf("Entered search_right with nmismatches %d, charpos %d, best_nmismatches %d\n",
		nmismatches,charpos,*best_nmismatches));
#endif

  if (nmismatches > *best_nmismatches && !collect_all_p) {
    return splicesites_i;

  } else if (nmismatches > max_mismatches_allowed) {
    return splicesites_i;

  } else if (single_leaf_p(leaf = triecontents[0])) {
    position = splicesites[leaf];
    if (position > limit_high) {
      debug2(printf("Found leaf %u at %u, but above limit_high %u\n",leaf,position,limit_high));
      return splicesites_i;

    } else if (charpos + 1 < pos3) {
      if (pos3 - pos5 <= 16) {
	/* Recomputes entire segment to determine mismatches (necessary because of splicefrags) */
	nmismatches = Genome_count_mismatches_fragment_right(query_compress,pos5,pos3,
							     splicefrags_ref[leaf],splicefrags_alt[leaf]);
      } else {
	/* Can happen in search for short middle exon */
	segment_left = splicesites[leaf] - pos5;
	nmismatches =
	  Genome_count_mismatches_substring(query_compress,segment_left,pos5,pos3,plusp,genestrand,first_read_p);
      }

      debug2(printf("Found leaf %u at %u => %d mismatches\n",leaf,position,nmismatches));
      if (nmismatches < *best_nmismatches) {
	debug2(printf("  nmismatches %d < best_nmismatches %d => setting splicesites_i to be %d (new best_nmismatches %d)\n",
		      nmismatches,*best_nmismatches,(int) leaf,nmismatches));
	*best_nmismatches = nmismatches;
	Intlist_free(&(*nmismatches_list));
	*nmismatches_list = Intlist_push(NULL,nmismatches);
	Intlist_free(&splicesites_i);
	return Intlist_push(NULL,(int) leaf);
      } else if (nmismatches == *best_nmismatches) {
	debug2(printf("  nmismatches %d == best_nmismatches %d => pushing leaf %d (same best_nmismatches %d)\n",
		      nmismatches,*best_nmismatches,(int) leaf,nmismatches));
	*nmismatches_list = Intlist_push(*nmismatches_list,nmismatches);
	return Intlist_push(splicesites_i,(int) leaf);
      } else if (nmismatches <= max_mismatches_allowed && collect_all_p) {
	debug2(printf("  nmismatches %d > best_nmismatches %d, but collecting all => pushing leaf %d)\n",
		      nmismatches,*best_nmismatches,(int) leaf));
	*nmismatches_list = Intlist_push(*nmismatches_list,nmismatches);
	return Intlist_push(splicesites_i,(int) leaf);
      } else {
	return splicesites_i;
      }
      
    } else {
      debug2(printf("Found leaf %u at %u, and completely finished query, so unique\n",leaf,position));
      assert(nmismatches <= max_mismatches_allowed);
      if (nmismatches < *best_nmismatches) {
	debug2(printf("  nmismatches %d < best_nmismatches %d => setting splicesites_i to be %d (new best_nmismatches %d)\n",
		      nmismatches,*best_nmismatches,(int) leaf,nmismatches));
	*best_nmismatches = nmismatches;
	Intlist_free(&(*nmismatches_list));
	*nmismatches_list = Intlist_push(NULL,nmismatches);
	Intlist_free(&splicesites_i);
	return Intlist_push(NULL,(int) leaf);
      } else if (nmismatches == *best_nmismatches) {
	debug2(printf("  nmismatches %d == best_nmismatches %d => pushing leaf %d (same best_nmismatches %d)\n",
		      nmismatches,*best_nmismatches,(int) leaf,nmismatches));
	*nmismatches_list = Intlist_push(*nmismatches_list,nmismatches);
	return Intlist_push(splicesites_i,(int) leaf);
      } else if (/*nmismatches <= max_mismatches_allowed &&*/ collect_all_p) {
	debug2(printf("  nmismatches %d > best_nmismatches %d => pushing leaf %d\n",
		      nmismatches,*best_nmismatches,(int) leaf));
	*nmismatches_list = Intlist_push(*nmismatches_list,nmismatches);
	return Intlist_push(splicesites_i,(int) leaf);
      } else {
	return splicesites_i;
      }
    }

  } else if (multiple_leaf_p(leaf)) {
    nleaves = (int) (-leaf);

    /* Need to compute nmismatches multiple times if we have SNPs */
    for (i = 1; i <= nleaves; i++) {
      leaf = triecontents[i];
      position = splicesites[leaf];
      if (position > limit_high) {
	/* Skip */

      } else if (charpos + 1 < pos3) {
	if (pos3 - pos5 <= 16) {
	  /* Recomputes entire segment to determine mismatches (necessary because of splicefrags) */
	  nmismatches = Genome_count_mismatches_fragment_right(query_compress,pos5,pos3,
							       splicefrags_ref[leaf],splicefrags_alt[leaf]);
	} else {
	  /* Can happen in search for short middle exon */
	  segment_left = splicesites[leaf] - pos5;
	  nmismatches =
	    Genome_count_mismatches_substring(query_compress,segment_left,pos5,pos3,plusp,genestrand,first_read_p);
	}

	debug2(printf("Found leaf %u at %u => %d mismatches\n",leaf,position,nmismatches));
	if (nmismatches < *best_nmismatches) {
	  debug2(printf("  nmismatches %d < best_nmismatches %d => setting splicesites_i to be %d (new best_nmismatches %d)\n",
			nmismatches,*best_nmismatches,(int) leaf,nmismatches));
	  *best_nmismatches = nmismatches;
	  Intlist_free(&(*nmismatches_list));
	  *nmismatches_list = Intlist_push(NULL,nmismatches);
	  Intlist_free(&splicesites_i);
	  splicesites_i =  Intlist_push(NULL,(int) leaf);
	} else if (nmismatches == *best_nmismatches) {
	  debug2(printf("  nmismatches %d == best_nmismatches %d => pushing leaf %d (same best_nmismatches %d)\n",
			nmismatches,*best_nmismatches,(int) leaf,nmismatches));
	  *nmismatches_list = Intlist_push(*nmismatches_list,nmismatches);
	  splicesites_i = Intlist_push(splicesites_i,(int) leaf);
	} else if (nmismatches <= max_mismatches_allowed && collect_all_p) {
	  debug2(printf("  nmismatches %d > best_nmismatches %d, but collecting all => pushing leaf %d\n",
			nmismatches,*best_nmismatches,(int) leaf));
	  *nmismatches_list = Intlist_push(*nmismatches_list,nmismatches);
	  splicesites_i = Intlist_push(splicesites_i,(int) leaf);
	}

      } else {
	debug2(printf("Found leaf %u at %u, and completely finished query, so unique\n",leaf,position));
	assert(nmismatches <= max_mismatches_allowed);
	if (nmismatches < *best_nmismatches) {
	  debug2(printf("  nmismatches %d < best_nmismatches %d => setting splicesites_i to be %d (new best_nmismatches %d)\n",
			nmismatches,*best_nmismatches,(int) leaf,nmismatches));
	  *best_nmismatches = nmismatches;
	  Intlist_free(&(*nmismatches_list));
	  *nmismatches_list = Intlist_push(NULL,nmismatches);
	  Intlist_free(&splicesites_i);
	  splicesites_i = Intlist_push(NULL,(int) leaf);
	} else if (nmismatches == *best_nmismatches) {
	  debug2(printf("  nmismatches %d == best_nmismatches %d => pushing leaf %d (same best_nmismatches %d)\n",
			nmismatches,*best_nmismatches,(int) leaf,nmismatches));
	  *nmismatches_list = Intlist_push(*nmismatches_list,nmismatches);
	  splicesites_i = Intlist_push(splicesites_i,(int) leaf);
	} else if (/*nmismatches <= max_mismatches_allowed &&*/ collect_all_p) {
	  debug2(printf("  nmismatches %d > best_nmismatches %d, but collecting all => pushing leaf %d\n",
			nmismatches,*best_nmismatches,(int) leaf));
	  *nmismatches_list = Intlist_push(*nmismatches_list,nmismatches);
	  splicesites_i = Intlist_push(splicesites_i,(int) leaf);
	}
      }
    }

    return splicesites_i;

  } else {
    /* Non-leaf, and characters left, so recurse */
#ifdef USE_2BYTE_RELOFFSETS
    get_offsets(&offseta,&offsetc,&offsetg,&offsett,triecontents[1],triecontents[2]);
#else
    offseta = (int) triecontents[1];
    offsetc = (int) triecontents[2];
    offsetg = (int) triecontents[3];
    offsett = (int) triecontents[4];
#endif

    if (charpos + 1 < pos3) {
      c = queryptr[charpos + 1];
      debug2(printf("Character at pos %d is %c\n",charpos + 1,c));
    } else {
      c = '\0';
    }

    if (offseta > 0) {
      splicesites_i = search_right(&(*best_nmismatches),&(*nmismatches_list),splicesites_i,
				   &(triecontents[-offseta]),queryptr,query_compress,
				   pos5,pos3,plusp,genestrand,first_read_p,collect_all_p,max_mismatches_allowed,
				   nmismatches+(c != '\0' && c != 'A'),charpos+1,limit_high);
    }
      
    if (offsetc > 0) {
      splicesites_i = search_right(&(*best_nmismatches),&(*nmismatches_list),splicesites_i,
				   &(triecontents[-offsetc]),queryptr,query_compress,
				   pos5,pos3,plusp,genestrand,first_read_p,collect_all_p,max_mismatches_allowed,
				   nmismatches+(c != '\0' && c != 'C'),charpos+1,limit_high);
    }

    if (offsetg > 0) {
      splicesites_i = search_right(&(*best_nmismatches),&(*nmismatches_list),splicesites_i,
				   &(triecontents[-offsetg]),queryptr,query_compress,
				   pos5,pos3,plusp,genestrand,first_read_p,collect_all_p,max_mismatches_allowed,
				   nmismatches+(c != '\0' && c != 'G'),charpos+1,limit_high);
    }

    if (offsett > 0) {
      splicesites_i = search_right(&(*best_nmismatches),&(*nmismatches_list),splicesites_i,
				   &(triecontents[-offsett]),queryptr,query_compress,
				   pos5,pos3,plusp,genestrand,first_read_p,collect_all_p,max_mismatches_allowed,
				   nmismatches+(c != '\0' && c != 'T'),charpos+1,limit_high);
    }

    return splicesites_i;
  }
}

#define OBSMAX_PENALTY 1

/* Allowing 1 mismatches in the internal (exon) region of 6 */
#define SPLICE_REGION 6
#define SPLICE_REGION_MISMATCHES 1 /* was 2 */

Intlist_T
Splicetrie_find_left (int *best_nmismatches, Intlist_T *nmismatches_list, int i,
		      Univcoord_T origleft, int pos5, int pos3, Univcoord_T chroffset,
		      Compress_T query_compress, char *queryptr, int querylength,
		      int max_mismatches_allowed, bool plusp, int genestrand, bool first_read_p,
		      bool collect_all_p) {
  Intlist_T splicesites_i = NULL, p, q;
  int closesti;
  int spliceregion;
#if 0
  unsigned int *triebranch_obs = NULL, *triebranch_max = NULL;
#endif
  Triecontent_T *triestart_obs, *triestart_max;
  int best_nmismatches_obs, best_nmismatches_max, nmismatches_int, nmismatches;
  int obsmax_penalty;
#ifdef LOOSE_ALLOWANCE
  int nmismatches_ext;
#endif

  int *array_i, *array_mm, lastvalue;
  int n, k;

  debug2(printf("Splicetrie_find_left called with #%d at origleft %u, pos5 %d, pos3 %d, collect_all_p %d\n",
		i,origleft,pos5,pos3,collect_all_p));

  if (pos5 == pos3) {
    debug2(printf("Splicetrie_find_left returning NULL, because pos5 == pos3\n"));
    /* *penalty = 0; */
    return (Intlist_T) NULL;
  } else if (amb_clip_p == false && pos3 - pos5 < min_shortend) {
    return (Intlist_T) NULL;
  }

  debug2(printf("max_mismatches_allowed %d",max_mismatches_allowed));
  if (max_mismatches_allowed > (pos3 - pos5)/3) {
    max_mismatches_allowed = (pos3 - pos5)/3;
    debug2(printf(" => limited by %d mismatches for length of %d\n",
		  max_mismatches_allowed,pos3-pos5));
  }
  debug2(printf("\n"));


  /* Check internal to splice site */
  if ((spliceregion = pos3 + SPLICE_REGION) > querylength) {
    spliceregion = querylength;
  }

  nmismatches_int =
    Genome_count_mismatches_substring(query_compress,origleft,pos3,spliceregion,plusp,genestrand,first_read_p);
  debug2(printf("  internal at %u => %d nmismatches/%d\n",
		origleft,nmismatches_int,SPLICE_REGION));
  if (nmismatches_int > SPLICE_REGION_MISMATCHES) {
    debug2(printf("  too many mismatches internal to splice site => returning NULL\n"));
    *nmismatches_list = (Intlist_T) NULL;
    return (Intlist_T) NULL;
  }


#ifdef LOOSE_ALLOWANCE
  nmismatches_ext = 
    Genome_count_mismatches_substring(query_compress,origleft,pos5,pos3,plusp,genestrand);
  debug2(printf("  extension at %u => %d nmismatches\n",origleft,nmismatches_ext));

  if (max_mismatches_allowed > nmismatches_ext) {
    max_mismatches_allowed = nmismatches_ext;
  }
#else
  /* Want strict allowance to avoid reliance on bad gene models */
  max_mismatches_allowed = 1;
#endif


#if 0
  if (triecontents_max == NULL) {
    Splicetrie_build_one(&triebranch_obs,&triestart_obs,
			 &triebranch_max,&triestart_max,
			 nsplicepartners_skip,nsplicepartners_obs,nsplicepartners_max,
			 i,splicetypes,splicestrings);
  }
#endif


  *nmismatches_list = NULL;
  *best_nmismatches = max_mismatches_allowed;

  best_nmismatches_obs = max_mismatches_allowed;
  if (trieoffsets_obs == NULL) {
    obsmax_penalty = 0;

  } else if (trieoffsets_obs[i] == NULL_POINTER) {
    obsmax_penalty = 0;

  } else {
    triestart_obs = &(triecontents_obs[trieoffsets_obs[i]]);
#if 0
    debug2(printf("Trie for observed splice lengths:\n"));
    debug2(Splicetrie_dump(triestart_obs,splicesites,splicetypes[i],splicefrags_ref));
#endif

    splicesites_i = search_left(&best_nmismatches_obs,&(*nmismatches_list),/*splicesites_i*/NULL,
				triestart_obs,queryptr,query_compress,
				pos5,pos3,plusp,genestrand,first_read_p,collect_all_p,max_mismatches_allowed,
				/*nmismatches*/0,/*charpos*/pos3,/*limit_low*/chroffset+pos3-1);
    /* *penalty = 0; */
    *best_nmismatches = best_nmismatches_obs;
    obsmax_penalty = OBSMAX_PENALTY;
  }

  best_nmismatches_max = best_nmismatches_obs - obsmax_penalty;
  if (trieoffsets_max == NULL) {
    /* Skip */
    
  } else if (trieoffsets_max[i] == NULL_POINTER) {
    /* Skip */

  } else {
    triestart_max = &(triecontents_max[trieoffsets_max[i]]);
#if 0
    debug2(printf("Trie for maximum allowed splice lengths:\n"));
    debug2(Splicetrie_dump(triestart_max,splicesites,splicetypes[i],splicefrags_ref));
#endif

    splicesites_i = search_left(&best_nmismatches_max,&(*nmismatches_list),splicesites_i,
				triestart_max,queryptr,query_compress,
				pos5,pos3,plusp,genestrand,first_read_p,collect_all_p,
				/*max_mismatches_allowed*/best_nmismatches_obs - obsmax_penalty,
				/*nmismatches*/0,/*charpos*/pos3,/*limit_low*/chroffset+pos3-1);
    debug2(printf("best_nmismatches_obs = %d, best_nmismatches_max = %d\n",
		  best_nmismatches_obs,best_nmismatches_max));
    if (best_nmismatches_max + obsmax_penalty < best_nmismatches_obs) {
      /* *penalty = OBSMAX_PENALTY; */
      *best_nmismatches = best_nmismatches_max;
    }
  }

#if 0
  if (triecontents_max == NULL) {
    FREE(triebranch_max);
    FREE(triebranch_obs);
  }
#endif

  debug2(printf("Final number of splices: %d\n",Intlist_length(splicesites_i)));
  if ((n = Intlist_length(splicesites_i)) > 1) {
    if (snpp == true) {
      array_mm = (int *) MALLOCA(n * sizeof(int)); /* elts */
      array_i = (int *) MALLOCA(n * sizeof(int)); /* keyvalues */
      Intlist_array_dual_ascending_by_key(array_mm,array_i,n,*nmismatches_list,splicesites_i);
      Intlist_free(&splicesites_i);
      Intlist_free(&(*nmismatches_list));

      splicesites_i = (Intlist_T) NULL;
      *nmismatches_list = (Intlist_T) NULL;

      splicesites_i = Intlist_push(splicesites_i,array_i[0]);
      *nmismatches_list = Intlist_push(*nmismatches_list,array_mm[0]);
      lastvalue = array_i[0];
      for (k = 1; k < n; k++) {
	if (array_i[k] == lastvalue) {
	  /* Skip */
	} else {
	  splicesites_i = Intlist_push(splicesites_i,array_i[k]);
	  *nmismatches_list = Intlist_push(*nmismatches_list,array_mm[k]);
	  lastvalue = array_i[k];
	}
      }
      FREEA(array_mm);
      FREEA(array_i);
    }

    if (amb_closest_p == true) {
      debug2(printf("Before checking for closest, splicesites are %s\n",Intlist_to_string(splicesites_i)));
      closesti = Intlist_head(splicesites_i);
      nmismatches = Intlist_head(*nmismatches_list);
      for (p = Intlist_next(splicesites_i), q = Intlist_next(*nmismatches_list); p != NULL; p = Intlist_next(p), q = Intlist_next(q)) {
	if (Intlist_head(p) > closesti) {
	  closesti = Intlist_head(p);
	  nmismatches = Intlist_head(q);
	}
      }
      Intlist_free(&splicesites_i);
      Intlist_free(&(*nmismatches_list));
      splicesites_i = Intlist_push(NULL,closesti);
      *nmismatches_list = Intlist_push(NULL,nmismatches);

    } else if (amb_clip_p == false && Intlist_length(splicesites_i) > 1) {
      Intlist_free(&splicesites_i);
      Intlist_free(&(*nmismatches_list));
      splicesites_i = (Intlist_T) NULL;
    }
  }

  debug2(printf("Splicetrie_find_left returning %s (best_nmismatches %d)\n",
		Intlist_to_string(splicesites_i),*best_nmismatches));

  return splicesites_i;
}


Intlist_T
Splicetrie_find_right (int *best_nmismatches, Intlist_T *nmismatches_list, int i,
		       Univcoord_T origleft, int pos5, int pos3, Univcoord_T chrhigh,
		       Compress_T query_compress, char *queryptr, int max_mismatches_allowed,
		       bool plusp, int genestrand, bool first_read_p, bool collect_all_p) {
  Intlist_T splicesites_i = NULL, p, q;
  int closesti;
  int spliceregion;
#if 0
  unsigned int *triebranch_obs = NULL, *triebranch_max = NULL;
#endif
  Triecontent_T *triestart_obs, *triestart_max;
  int best_nmismatches_obs, best_nmismatches_max, nmismatches_int, nmismatches;
  int obsmax_penalty;
#ifdef LOOSE_ALLOWANCE
  int nmismatches_ext;
#endif

  int *array_i, *array_mm, lastvalue;
  int n, k;

  debug2(printf("Splicetrie_find_right called with #%d at origleft %u, pos5 %d, pos3 %d, collect_all_p %d\n",
		i,origleft,pos5,pos3,collect_all_p));
  debug2(printf("chrhigh is %u\n",chrhigh));

  if (pos5 == pos3) {
    debug2(printf("Splicetrie_find_right returning NULL, because pos5 == pos3\n"));
    /* *penalty = 0; */
    return (Intlist_T) NULL;
  } else if (amb_clip_p == false && pos3 - pos5 < min_shortend) {
    return (Intlist_T) NULL;
  }

  debug2(printf("max_mismatches_allowed %d",max_mismatches_allowed));
  if (max_mismatches_allowed > (pos3 - pos5)/3) {
    max_mismatches_allowed = (pos3 - pos5)/3;
    debug2(printf(" => limited by %d mismatches for length of %d\n",
		  max_mismatches_allowed,pos3-pos5));
  }
  debug2(printf("\n"));


  /* Check internal to splice site */
  if ((spliceregion = pos5 - SPLICE_REGION) < 0) {
    spliceregion = 0;
  }

  nmismatches_int = 
    Genome_count_mismatches_substring(query_compress,origleft,spliceregion,pos5,plusp,genestrand,first_read_p);
  debug2(printf("  internal at %u => %d nmismatches/%d\n",
		origleft,nmismatches_int,SPLICE_REGION));
  if (nmismatches_int > SPLICE_REGION_MISMATCHES) {
    debug2(printf("  too many mismatches internal to splice site => returning NULL\n"));
    *nmismatches_list = (Intlist_T) NULL;
    return (Intlist_T) NULL;
  }


#ifdef LOOSE_ALLOWANCE
  nmismatches_ext = 
    Genome_count_mismatches_substring(query_compress,origleft,pos5,pos3,plusp,genestrand);
  debug2(printf("  extension at %u => %d nmismatches\n",origleft,nmismatches_ext));
    
  if (max_mismatches_allowed > nmismatches_ext) {
    max_mismatches_allowed = nmismatches_ext;
  }
#else
  /* Want strict allowance to avoid reliance on bad gene models */
  max_mismatches_allowed = 1;
#endif


#if 0
  if (triecontents_max == NULL) {
    Splicetrie_build_one(&triebranch_obs,&triestart_obs,
			 &triebranch_max,&triestart_max,
			 nsplicepartners_skip,nsplicepartners_obs,nsplicepartners_max,
			 i,splicetypes,splicestrings);
  }
#endif


  *nmismatches_list = NULL;
  *best_nmismatches = max_mismatches_allowed;

  best_nmismatches_obs = max_mismatches_allowed;
  if (trieoffsets_obs == NULL) {
    obsmax_penalty = 0;

  } else if (trieoffsets_obs[i] == NULL_POINTER) {
    obsmax_penalty = 0;

  } else {
    triestart_obs = &(triecontents_obs[trieoffsets_obs[i]]);
#if 0
    debug2(printf("Trie for observed splice lengths:\n"));
    debug2(Splicetrie_dump(triestart_obs,splicesites,splicetypes[i],splicefrags_ref));
#endif

    splicesites_i = search_right(&best_nmismatches_obs,&(*nmismatches_list),/*splicesites_i*/NULL,
				 triestart_obs,queryptr,query_compress,
				 pos5,pos3,plusp,genestrand,first_read_p,collect_all_p,max_mismatches_allowed,
				 /*nmismatches*/0,/*charpos*/pos5-1,/*limit_high*/chrhigh+pos5-pos3);
    /* *penalty = 0; */
    *best_nmismatches = best_nmismatches_obs;
    obsmax_penalty = OBSMAX_PENALTY;
  }

  best_nmismatches_max = best_nmismatches_obs - obsmax_penalty;
  if (trieoffsets_max == NULL) {
    /* Skip */

  } else if (trieoffsets_max[i] == NULL_POINTER) {
    /* Skip */

  } else {
    triestart_max = &(triecontents_max[trieoffsets_max[i]]);
#if 0
    debug2(printf("Trie for maximum allowed splice lengths:\n"));
    debug2(Splicetrie_dump(triestart_max,splicesites,splicetypes[i],splicefrags_ref));
#endif

    splicesites_i = search_right(&best_nmismatches_max,&(*nmismatches_list),splicesites_i,
				 triestart_max,queryptr,query_compress,
				 pos5,pos3,plusp,genestrand,first_read_p,collect_all_p,
				 /*max_mismatches_allowed*/best_nmismatches_obs - obsmax_penalty,
				 /*nmismatches*/0,/*charpos*/pos5-1,/*limit_high*/chrhigh+pos5-pos3);
    debug2(printf("best_nmismatches_obs = %d, best_nmismatches_max = %d\n",
		  best_nmismatches_obs,best_nmismatches_max));
    if (best_nmismatches_max + obsmax_penalty < best_nmismatches_obs) {
      /* *penalty = OBSMAX_PENALTY; */
      *best_nmismatches = best_nmismatches_max;
    }
  }

#if 0
  if (triecontents_max == NULL) {
    FREE(triebranch_max);
    FREE(triebranch_obs);
  }
#endif

  debug2(printf("Final number of splices: %d\n",Intlist_length(splicesites_i)));
  if ((n = Intlist_length(splicesites_i)) > 1) {
    if (snpp == true) {
      array_mm = (int *) MALLOCA(n * sizeof(int)); /* elts */
      array_i = (int *) MALLOCA(n * sizeof(int)); /* keyvalues */
      Intlist_array_dual_ascending_by_key(array_mm,array_i,n,*nmismatches_list,splicesites_i);
      Intlist_free(&splicesites_i);
      Intlist_free(&(*nmismatches_list));

      splicesites_i = (Intlist_T) NULL;
      *nmismatches_list = (Intlist_T) NULL;

      splicesites_i = Intlist_push(splicesites_i,array_i[0]);
      *nmismatches_list = Intlist_push(*nmismatches_list,array_mm[0]);
      lastvalue = array_i[0];
      for (k = 1; k < n; k++) {
	if (array_i[k] == lastvalue) {
	  /* Skip */
	} else {
	  splicesites_i = Intlist_push(splicesites_i,array_i[k]);
	  *nmismatches_list = Intlist_push(*nmismatches_list,array_mm[k]);
	  lastvalue = array_i[k];
	}
      }
      FREEA(array_mm);
      FREEA(array_i);
    }

    if (amb_closest_p == true) {
      debug2(printf("Before checking for closest, splicesites are %s\n",Intlist_to_string(splicesites_i)));
      closesti = Intlist_head(splicesites_i);
      nmismatches = Intlist_head(*nmismatches_list);
      for (p = Intlist_next(splicesites_i), q = Intlist_next(*nmismatches_list); p != NULL; p = Intlist_next(p), q = Intlist_next(q)) {
	if (Intlist_head(p) < closesti) {
	  closesti = Intlist_head(p);
	  nmismatches = Intlist_head(q);
	}
      }
      Intlist_free(&splicesites_i);
      Intlist_free(&(*nmismatches_list));
      splicesites_i = Intlist_push(NULL,closesti);
      *nmismatches_list = Intlist_push(NULL,nmismatches);

    } else if (amb_clip_p == false && Intlist_length(splicesites_i) > 1) {
      Intlist_free(&splicesites_i);
      Intlist_free(&(*nmismatches_list));
      splicesites_i = (Intlist_T) NULL;
    }

  }

    

  debug2(printf("Splicetrie_find_right returning %s (best_nmismatches %d)\n",
		Intlist_to_string(splicesites_i),*best_nmismatches));

  return splicesites_i;
}

#endif

