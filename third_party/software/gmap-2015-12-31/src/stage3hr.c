static char rcsid[] = "$Id: stage3hr.c 184529 2016-02-18 19:28:56Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "stage3hr.h"

#include <stdlib.h>		/* For qsort */
#include <string.h>
#include <strings.h>
#include <ctype.h>		/* For islower */
#include <math.h>		/* For exp() and log10() */
#include "assert.h"
#include "mem.h"
#include "chrnum.h"
#include "complement.h"
#include "interval.h"
#include "listdef.h"
#include "substring.h"
#include "junction.h"
#include "genome128_hr.h"
#include "mapq.h"
#include "pair.h"		/* For Pair_print_gsnap and Pair_compute_mapq */
#include "comp.h"		/* For Stage3end_run_gmap */
#include "maxent_hr.h"
#include "fastlog.h"


/* Scores for amb_status_inside */
#define AMB_RESOLVED_BYLENGTH 0
#define AMB_RESOLVED_BYMATCHES 0
#define AMB_NOT_AMBIGUOUS 1
#define AMB_UNRESOLVED_MULTIPLE 2 /* Not so preferable, since none had the expected length */
#define AMB_UNRESOLVED_TOOCLOSE 3 /* Worse than multiple, since no options work */



/* Originally added to avoid CIGAR strings like 1S99H, but results in
   errors if first chromosome is circular.  Now checking in samprint.c
   whether the CIGAR string is bad.  But now needed again because we
   are allowing alignments that are out of bounds.  However, can cause
   left to be negative, so cannot use. */
/* #define SOFT_CLIPS_AVOID_CIRCULARIZATION 1 -- Turn off here, but on in pair.c */

#define TRANSLOC_SPECIAL 1	/* Eliminates translocations if non-translocations are found */

#define MAX_HITS 100000


#define CONCORDANT_TEXT "concordant"
#define PAIRED_TEXT "paired"
#define UNPAIRED_TEXT "unpaired"

#ifdef USE_TALLY_RATIO
#define TALLY_RATIO 2.0
#endif

#define SUBSUMPTION_SLOP 10	/* Should allow for short insert lengths */
/* #define TERMINAL_SECOND_CLASS 1 -- enabling this leads to poor results */

#define TERMINAL_COMPUTE_MINLENGTH 40
/* #define SCORE_INDELS 1 -- Needed to compare genomic positions with and without indels */

#define OUTERLENGTH_SLOP 100

/* #define USE_ALLOCA_FOR_HITS 1 -- can lead to stack overflow */


#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif

/* Stage3end_new */
#ifdef DEBUG0
#define debug0(x) x
#else
#define debug0(x)
#endif

/* printing */
#ifdef DEBUG1 
#define debug1(x) x
#else
#define debug1(x)
#endif


/* new_insertion */
#ifdef DEBUG2
#define debug2(x) x
#else
#define debug2(x)
#endif

/* new_deletion */
#ifdef DEBUG3
#define debug3(x) x
#else
#define debug3(x)
#endif

/* Stage3end_optimal_score */
#ifdef DEBUG4
#define debug4(x) x
#else
#define debug4(x)
#endif


/* Stage3_pair_up_concordant */
#ifdef DEBUG5
#define debug5(x) x
#else
#define debug5(x)
#endif

/* Stage3_pair_up_concordant, details */
#ifdef DEBUG5A
#define debug5a(x) x
#else
#define debug5a(x)
#endif


/* Stage3pair_optimal_score */
#ifdef DEBUG6
#define debug6(x) x
#else
#define debug6(x)
#endif


/* Stage3end_remove_duplicates and Stage3end_remove_overlaps */
#ifdef DEBUG7
#define debug7(x) x
#else
#define debug7(x)
#endif

/* Stage3pair_remove_overlaps */
#ifdef DEBUG8
#define debug8(x) x
#else
#define debug8(x)
#endif

/* Resolving ambiguous splice sites in Stage3pair_new */
#ifdef DEBUG9
#define debug9(x) x
#else
#define debug9(x)
#endif

/* insert length calculation */
#ifdef DEBUG10
#define debug10(x) x
#else
#define debug10(x)
#endif

/* genomicbound */
#ifdef DEBUG11
#define debug11(x) x
#else
#define debug11(x)
#endif

/* circular chromosomes */
#ifdef DEBUG12
#define debug12(x) x
#else
#define debug12(x)
#endif

/* substring_gmap */
#ifdef DEBUG13
#define debug13(x) x
#else
#define debug13(x)
#endif

/* Stage3_determine_pairtype */
#ifdef DEBUG14
#define debug14(x) x
#else
#define debug14(x)
#endif

/* Stage3pair_overlap */
#ifdef DEBUG15
#define debug15(x) x
#else
#define debug15(x)
#endif



#define MAPQ_MAXIMUM_SCORE 40


static bool want_random_p;
static bool invert_first_p;
static bool invert_second_p;
static Genome_T genome;

static Univ_IIT_T chromosome_iit;
static int nchromosomes;
static int circular_typeint;

static IIT_T genes_iit;
static int *genes_divint_crosstable;
static IIT_T tally_iit;
static int *tally_divint_crosstable;
static IIT_T runlength_iit;
static int *runlength_divint_crosstable;

static int pairmax;

#if 0
static int expected_pairlength;
static int pairlength_deviation;
#else
static int expected_pairlength_low;
static int expected_pairlength_high;
static int expected_pairlength_very_high;
#endif

static int amb_penalty = 2;
static int localsplicing_penalty;
static int indel_penalty_middle;
static int antistranded_penalty;
static bool favor_multiexon_p;
static int gmap_min_nconsecutive;

static int ambig_end_interval;	/* For penalizing large ambiguous ends
				   in GMAP alignments, since such ends
				   should have been found */
static bool novelsplicingp;
static bool merge_samechr_p;
static bool *circularp;

static char *failedinput_root;
static bool print_m8_p;


/* Probably not good to use in certain genomic regions, unless we also
   use known splicesites with distance information. */
/* But sometimes need to use to get correct mapping */
static bool favor_ambiguous_p;


void
Stage3hr_setup (bool invert_first_p_in, bool invert_second_p_in, Genome_T genome_in,
		Univ_IIT_T chromosome_iit_in, int nchromosomes_in, int circular_typeint_in,
		IIT_T genes_iit_in, int *genes_divint_crosstable_in,
		IIT_T tally_iit_in, int *tally_divint_crosstable_in,
		IIT_T runlength_iit_in, int *runlength_divint_crosstable_in,
		bool distances_observed_p, int pairmax_in,
		Chrpos_T expected_pairlength, Chrpos_T pairlength_deviation,
		int localsplicing_penalty_in, int indel_penalty_middle_in,
		int antistranded_penalty_in, bool favor_multiexon_p_in,
		int gmap_min_nconsecutive_in, int index1part,
		int index1interval, bool novelsplicingp_in, bool merge_samechr_p_in,
		bool *circularp_in, char *failedinput_root_in,
		bool print_m8_p_in, bool want_random_p_in) {
  invert_first_p = invert_first_p_in;
  invert_second_p = invert_second_p_in;
  genome = genome_in;

  chromosome_iit = chromosome_iit_in;
  nchromosomes = nchromosomes_in;
  circular_typeint = circular_typeint_in;
  genes_iit = genes_iit_in;
  genes_divint_crosstable = genes_divint_crosstable_in;
  tally_iit = tally_iit_in;
  tally_divint_crosstable = tally_divint_crosstable_in;
  runlength_iit = runlength_iit_in;
  runlength_divint_crosstable = runlength_divint_crosstable_in;
  localsplicing_penalty = localsplicing_penalty_in;
  indel_penalty_middle = indel_penalty_middle_in;
  antistranded_penalty = antistranded_penalty_in;
  favor_multiexon_p = favor_multiexon_p_in;
  gmap_min_nconsecutive = gmap_min_nconsecutive_in;

  pairmax = pairmax_in;
  if (pairlength_deviation > expected_pairlength) {
    expected_pairlength_low = 0;
  } else {
    expected_pairlength_low = expected_pairlength - pairlength_deviation;
  }
  expected_pairlength_high = expected_pairlength + pairlength_deviation;
  expected_pairlength_very_high = expected_pairlength + 10*pairlength_deviation;

  if (distances_observed_p == true) {
    favor_ambiguous_p = false;
  } else {
    favor_ambiguous_p = true;
  }

#if 0
  ambig_end_interval = index1part + (index1interval - 1);
#else
  ambig_end_interval = 8;	/* Since GMAP uses 8-mers */
#endif

  novelsplicingp = novelsplicingp_in;
  merge_samechr_p = merge_samechr_p_in;
  circularp = circularp_in;

  failedinput_root = failedinput_root_in;

  print_m8_p = print_m8_p_in;
  want_random_p = want_random_p_in;

  return;
}



static char *
print_sense (int sense) {
  if (sense == SENSE_NULL) {
    return "sense:null";
  } else if (sense == SENSE_ANTI) {
    return "sense:anti";
  } else if (sense == SENSE_FORWARD) {
    return "sense:fwd";
  } else {
    abort();
  }
}



/* Note: Substring_T has genomiclength, but not Stage3end_T */
#define T Stage3end_T
struct T {
  Hittype_T hittype;
  int genestrand;
  bool sarrayp;			/* true if alignment found by suffix array */
  GMAP_source_T gmap_source;
  bool improved_by_gmap_p;	/* true if GMAP alignment based on this hit is better */

  Chrnum_T chrnum; /* Needed for printing paired-end results.  A chrnum of 0 indicates a distant splice. */
  Chrnum_T effective_chrnum;	/* For determining concordance */
  Chrnum_T other_chrnum;	/* 0 for non-translocations, and other chrnum besides effective_chrnum for translocations */
  Univcoord_T chroffset;
  Univcoord_T chrhigh;
  Chrpos_T chrlength;

  int querylength;		/* Needed for overlap and pairlength calculations */
  int querylength_adj;		/* Adjusted for insertions */

  Univcoord_T genomicstart;
  Univcoord_T genomicend;
  bool plusp;

  Univcoord_T low;
  Univcoord_T high;
  Chrpos_T genomiclength;
  Chrpos_T guided_insertlength; /* Used only by Stage3end_eval_and_sort_guided */

  float mapq_loglik;
  int mapq_score;
  int absmq_score;		/* Absolute MAPQ, for XQ and X2 flags */

  int nsegments;
  int score;			/* Includes colordiffs and penalties */
  int ntscore;			/* Includes penalties */
  int nmatches;
  int nmatches_posttrim;
  int gmap_max_match_length;		/* Used only by GMAP */
  double gmap_min_splice_prob;		/* Used only by GMAP */

  /* trim_left and trim_right should really be named trim_start and trim_end */
  int trim_left; /* Used by Stage3end_optimal_score for comparing terminals and non-terminals */
  int trim_right;
  bool trim_left_splicep;
  bool trim_right_splicep;

#if 0
  int penalties;		/* Indel penalties */
#endif
  int score_eventrim;		/* Temporary storage used by Stage3end_optimal_score */

  Overlap_T gene_overlap;
  long int tally;

  int nmismatches_whole;
  int nmismatches_bothdiff;
  int nmismatches_refdiff;	/* Set only for display */

  int nindels;			/* for indels */

  Chrpos_T distance;	/* for splicing or shortexon (sum of two distances) */
  Chrpos_T shortexonA_distance; /* for shortexon */
  Chrpos_T shortexonD_distance;	  /* for shortexon */

  int gmap_nindelbreaks;
  int gmap_cdna_direction;
  int gmap_nintrons;
  int sensedir;			/* for splicing */

  int nsplices;

#if 0
  bool start_ambiguous_p;
  bool end_ambiguous_p;
  int amb_length_donor;	/* For shortexon only */
  int amb_length_acceptor;	/* For shortexon only */
  double amb_prob_donor;	/* For shortexon */
  double amb_prob_acceptor;	/* For shortexon */
#endif

  int gmap_start_amb_length;	/* Needed because GMAP doesn't have substrings */
  int gmap_end_amb_length;	/* Needed because GMAP doesn't have substrings */
  Endtype_T gmap_start_endtype;	/* For GMAP, which has no substrings */
  Endtype_T gmap_end_endtype;	/* For GMAP, which has no substrings */

#if 0
  Univcoord_T *start_ambcoords;	/* Pointer to either ambcoords_donor or ambcoords_acceptor */
  Univcoord_T *end_ambcoords;	/* Pointer to either ambcoords_donor or ambcoords_acceptor */
  int start_nambcoords;		/* Equal to either nambcoords_donor or nambcoords_acceptor */
  int end_nambcoords;           /* Equal to either nambcoords_donor or nambcoords_acceptor */
  Univcoord_T *ambcoords_donor;
  Univcoord_T *ambcoords_acceptor;
  int nambcoords_donor;
  int nambcoords_acceptor;

  int *start_amb_knowni;        /* Pointer to either amb_knowni_donor or amb_knowni_acceptor */
  int *end_amb_knowni;          /* Pointer to either amb_knowni_donor or amb_knowni_acceptor */
  int *amb_knowni_donor;
  int *amb_knowni_acceptor;

  int *start_amb_nmismatches;	/* Pointer to either amb_nmismatches_donor or amb_nmismatches_acceptor */
  int *end_amb_nmismatches;     /* Pointer to either amb_nmismatches_donor or amb_nmismatches_acceptor */
  double *start_amb_probs;	/* Pointer to either amb_probs_donor or amb_probs_acceptor */
  double *end_amb_probs;	/* Pointer to either amb_probs_donor or amb_probs_acceptor */

  int *amb_nmismatches_donor;
  int *amb_nmismatches_acceptor;
  double *amb_probs_donor;
  double *amb_probs_acceptor;
#endif


  /* For GMAP alignment */
  struct Pair_T *pairarray;
  int npairs;
  List_T cigar_tokens;
  bool gmap_intronp;

  List_T substrings_1toN;	/* query position 1 to N */
  List_T substrings_Nto1;	/* query position N to 1.  Keeps only pointers to the substrings. */
  List_T substrings_LtoH;	/* Chromosomal low-to-high.  Keeps only pointers to the substrings. */

  List_T junctions_LtoH;
  List_T junctions_1toN;
  List_T junctions_Nto1;

  bool paired_usedp;
  bool paired_seenp;   /* for paired-end.  set to true by Stage3_pair_up(). */
  bool concordantp;    /* for paired-end.  set to true by Stage3_pair_up(). */

  int alias;			/* -1 if all below chrlength, 0 if straddles or NA (e.g., transloc), and +1 if above */
                                /* -2 if extends below beginning of circular chromosome, +2 if extends beyond end of second copy */
  int circularpos;		/* if alias == 0, then amount of queryseq below chrlength */
};


struct Stage3pair_T {
  Pairtype_T pairtype;
  int genestrand;

  T hit5;
  T hit3;
  bool private5p;			/* A private copy separate from hits5 and hits3, and not a pointer */
  bool private3p;

  Univcoord_T low;
  Univcoord_T high;
  int insertlength;
  int insertlength_expected_sign;	/* 1 if in (expected_pairlength_low, expected_pairlength_high),
					   0 if in (expected_pairlength_low, expected_pairlength_very_high), and
					   -1 if < expected_pairlength_low or > expected_pairlength_very_high */

  Chrpos_T outerlength;

  float mapq_loglik;
  int mapq_score;
  int absmq_score;

  int score;
  int nmatches;
  int nmatches_posttrim;

  int score_eventrim;

  Overlap_T gene_overlap;
  long int tally;

#ifdef USE_ABSDIFFLENGTH
  Chrpos_T absdifflength;
#endif
#ifdef USE_BINGO
  bool absdifflength_bingo_p;
#endif
  int dir;			/* -1, 0, or +1 */
  bool sense_consistent_p;

  int nsplices;

  bool circularp;		/* If either hit5 or hit3 are circular */
  int amb_resolve_5; /* Resolution of ambiguous end for this particular pair */
  int amb_resolve_3; /* Resolution of ambiguous end for this particular pair */
  int amb_status_inside;
};



char *
Stage3end_deletion_string (T this) {
  abort();
}


Hittype_T
Stage3end_hittype (T this) {
  return this->hittype;
}

static char *
hittype_string (Hittype_T hittype) {
  switch (hittype) {
  case EXACT: return "exact";
  case SUB: return "sub";
  case INSERTION: return "insertion";
  case DELETION: return "deletion";
  case HALFSPLICE_DONOR: return "donor";
  case HALFSPLICE_ACCEPTOR: return "acceptor";
  case SPLICE: return "splice";
  case SAMECHR_SPLICE: return "samechr_splice";
  case TRANSLOC_SPLICE: return "transloc_splice";
  case ONE_THIRD_SHORTEXON: return "one-third-shortexon";
  case TWO_THIRDS_SHORTEXON: return "two-thirds-shortexon";
  case SHORTEXON: return "shortexon";
  case SUBSTRINGS: return "substrings";
  case GMAP: return "gmap";
  case TERMINAL: return "terminal";
  default: abort();
  }
}

char *
Stage3end_hittype_string (T this) {
  return hittype_string(this->hittype);
}

int
Stage3end_genestrand (T this) {
  return this->genestrand;
}

bool
Stage3end_sarrayp (T this) {
  if (this == NULL) {
    /* Can happen if we call upon a mate in a halfmapping */
    return false;
  } else {
    return this->sarrayp;
  }
}

bool
Stage3end_improved_by_gmap_p (T this) {
  return this->improved_by_gmap_p;
}

void
Stage3end_set_improved_by_gmap (T this) {
  this->improved_by_gmap_p = true;
  return;
}

bool
Stage3end_anomalous_splice_p (T this) {
  if (this->chrnum == 0) {
    return true;
  } else if (this->hittype == SAMECHR_SPLICE) {
    return true;
  } else {
    return false;
  }
}


Chrnum_T
Stage3end_chrnum (T this) {
  if (this == NULL) {
    /* Can happen if we call upon a mate in a halfmapping */
    return 0;
  } else {
    return this->chrnum;
  }
}

Chrnum_T
Stage3end_effective_chrnum (T this) {
  if (this == NULL) {
    /* Can happen if we call upon a mate in a halfmapping */
    return 0;
  } else {
    return this->effective_chrnum;
  }
}

Univcoord_T
Stage3end_chroffset (T this) {
  return this->chroffset;
}

Univcoord_T
Stage3end_chrhigh (T this) {
  return this->chrhigh;
}

Chrpos_T
Stage3end_chrlength (T this) {
  if (this == NULL) {
    /* Can happen if we call upon a mate in a halfmapping */
    return 0;
  } else {
    return this->chrlength;
  }
}

Univcoord_T
Stage3end_genomicstart (T this) {
  return this->genomicstart;
}

Univcoord_T
Stage3end_genomicend (T this) {
  return this->genomicend;
}

/* For Goby */
int
Stage3end_query_alignment_length (T this) {
  int length = 0;
  List_T p;
  Substring_T substring;
  Junction_T junction;

  for (p = this->substrings_LtoH; p != NULL; p = List_next(p)) {
    substring = (Substring_T) List_head(p);
    length += Substring_match_length(substring);
  }
  for (p = this->junctions_LtoH; p != NULL; p = List_next(p)) {
    junction = (Junction_T) List_head(p);
    if (Junction_type(junction) == INS_JUNCTION) {
      length += Junction_nindels(junction);
    }
  }

  return length;
}

Chrpos_T
Stage3end_genomic_alignment_length (T this) {
  Chrpos_T length = 0;
  List_T p;
  Substring_T substring;
  Junction_T junction;

  for (p = this->substrings_LtoH; p != NULL; p = List_next(p)) {
    substring = (Substring_T) List_head(p);
    length += Substring_genomic_alignment_length(substring);
  }
  for (p = this->junctions_LtoH; p != NULL; p = List_next(p)) {
    junction = (Junction_T) List_head(p);
    if (Junction_type(junction) == DEL_JUNCTION) {
      length += (Chrpos_T) Junction_nindels(junction);
    }
  }

  return length;
}


Chrpos_T
Stage3end_chrpos_low_trim (T this) {
  Substring_T substring_low;

  substring_low = (Substring_T) List_head(this->substrings_LtoH);
  if (this->plusp == true) {
    return Substring_alignstart_trim(substring_low) - Substring_chroffset(substring_low);
  } else {
    return Substring_alignend_trim(substring_low) - Substring_chroffset(substring_low);
  }
}

int
Stage3end_mapq_score (T this) {
  return this->mapq_score;
}

int
Stage3end_absmq_score (T this) {
  return this->absmq_score;
}

int
Stage3end_score (T this) {
  return this->score;
}

int
Stage3end_gmap_max_match_length (T this) {
  return this->gmap_max_match_length;
}

double
Stage3end_gmap_min_splice_prob (T this) {
  return this->gmap_min_splice_prob;
}


int
Stage3end_best_score (List_T hits) {
  List_T p;
  T hit;
  int best_score = 1000000;

  for (p = hits; p != NULL; p = p->rest) {
    hit = (T) p->first;
    if (hit->score < best_score) {
      best_score = hit->score;
    }
  }

  return best_score;
}

bool
Stage3end_equiv_score_unpaired_p (List_T hits, int best_score) {
  List_T p;
  T hit;

  for (p = hits; p != NULL; p = p->rest) {
    hit = (T) p->first;
    if (hit->paired_usedp == false && hit->score <= best_score) {
      return true;
    }
  }

  return false;
}

int
Stage3end_best_score_paired (List_T hits) {
  List_T p;
  T hit;
  int best_score = 1000000;

  for (p = hits; p != NULL; p = p->rest) {
    hit = (T) p->first;
    if (hit->paired_usedp == true) {
      if (hit->score < best_score) {
	best_score = hit->score;
      }
    }
  }

  return best_score;
}


int
Stage3end_nmatches_posttrim (T this) {
  return this->nmatches_posttrim;
}

int
Stage3end_nmismatches_whole (T this) {
  return this->nmismatches_whole;
}

int
Stage3end_nmismatches_bothdiff (T this) {
  return this->nmismatches_bothdiff;
}

int
Stage3end_nmismatches_refdiff (T this) {
  return this->nmismatches_refdiff;
}

/* Called only for terminals */
Endtype_T
Stage3end_start_endtype (T this) {
  Substring_T substring;

  substring = (Substring_T) List_head(this->substrings_LtoH);
  return Substring_start_endtype(substring);
}

/* Called only for terminals */
Endtype_T
Stage3end_end_endtype (T this) {
  Substring_T substring;

  substring = (Substring_T) List_head(this->substrings_LtoH);
  return Substring_end_endtype(substring);
}

Endtype_T
Stage3end_gmap_start_endtype (T this) {
  return this->gmap_start_endtype;
}

Endtype_T
Stage3end_gmap_end_endtype (T this) {
  return this->gmap_end_endtype;
}

int
Stage3end_nindels (T this) {
  return this->nindels;
}

int
Stage3end_querylength (T this) {
  return this->querylength;
}

bool
Stage3end_plusp (T this) {
  return this->plusp;
}

bool
Stage3end_paired_usedp (T this) {
  return this->paired_usedp;
}

/* Accounts for ambig ends */
int
Stage3end_trim_left (T this) {
  return this->trim_left;
}

/* Accounts for ambig ends */
int
Stage3end_trim_right (T this) {
  return this->trim_right;
}

static int
start_amb_length (T this) {
  if (this->hittype == GMAP) {
    return this->gmap_start_amb_length;
  } else {
    return Substring_match_length_amb((Substring_T) List_head(this->substrings_1toN));
  }
}

static int
end_amb_length (T this) {
  if (this->hittype == GMAP) {
    return this->gmap_end_amb_length;
  } else {
    return Substring_match_length_amb((Substring_T) List_head(this->substrings_Nto1));
  }
}

int
Stage3end_trim_left_raw (T this) {
  return this->trim_left + start_amb_length(this);
}

int
Stage3end_trim_right_raw (T this) {
  return this->trim_right + end_amb_length(this);
}

int
Stage3end_circularpos (T this) {
  return this->circularpos;
}


Junction_T
Stage3end_junctionD (T this) {
  if (this->sensedir == SENSE_ANTI) {
    return (Junction_T) List_head(this->junctions_Nto1);
  } else {
    return (Junction_T) List_head(this->junctions_1toN);
  }
}

Junction_T
Stage3end_junctionA (T this) {
  if (this->sensedir == SENSE_ANTI) {
    return (Junction_T) List_head(this->junctions_1toN);
  } else {
    return (Junction_T) List_head(this->junctions_Nto1);
  }
}

List_T
Stage3end_substrings_LtoH (T this) {
  return this->substrings_LtoH;
}

List_T
Stage3end_junctions_LtoH (T this) {
  return this->junctions_LtoH;
}


/* Called only by samprint currently */
Substring_T
Stage3end_substring1 (T this) {
  return (Substring_T) List_head(this->substrings_1toN);
}

/* Called only by samprint currently */
Substring_T
Stage3end_substring2 (T this) {
  return (Substring_T) List_head(this->substrings_Nto1);
}


Substring_T
Stage3end_substring_donor (T this) {
  if (this->sensedir == SENSE_ANTI) {
    return (Substring_T) List_head(this->substrings_Nto1);
  } else if (this->sensedir == SENSE_FORWARD) {
    return (Substring_T) List_head(this->substrings_1toN);
  } else {
    abort();
  }
}

Substring_T
Stage3end_substring_acceptor (T this) {
  if (this->sensedir == SENSE_ANTI) { 
    return (Substring_T) List_head(this->substrings_1toN);
  } else if (this->sensedir == SENSE_FORWARD) {
    return (Substring_T) List_head(this->substrings_Nto1);
  } else {
    abort();
  }
}

/* Now same as Stage3end_substring_donor */
Substring_T
Stage3end_substringD (T this) {
  if (this->sensedir == SENSE_ANTI) {
    return (Substring_T) List_head(this->substrings_Nto1);
  } else {
    return (Substring_T) List_head(this->substrings_1toN);
  }
}

/* Now same as Stage3end_substring_acceptor */
Substring_T
Stage3end_substringA (T this) {
  if (this->sensedir == SENSE_ANTI) {
    return (Substring_T) List_head(this->substrings_1toN);
  } else {
    return (Substring_T) List_head(this->substrings_Nto1);
  }
}


Substring_T
Stage3end_substringS (T this) {
  return (Substring_T) List_head(List_next(this->substrings_1toN));
}


/* Same logic as in print_substrings in samprint.c to get the first substring for CIGAR or MD string */
Substring_T
Stage3end_substring_low (T this, int hardclip_low) {
  List_T p;

  if (this == NULL) {
    return (Substring_T) NULL;

  } else if (this->plusp == true) {
    p = this->substrings_LtoH;
    if (Substring_ambiguous_p((Substring_T) List_head(p)) == true) {
      p = List_next(p);
    }
    while (p != NULL && Substring_queryend((Substring_T) List_head(p)) <= hardclip_low) {
      debug15(printf("Plus: Skippping substring %d..%d against hardclip_low %d\n",
		     Substring_querystart((Substring_T) List_head(p)),Substring_queryend((Substring_T) List_head(p)),
		     hardclip_low));
      p = List_next(p);
    }
    assert(p != NULL);
    if (p == NULL) {
      return (Substring_T) NULL;
    } else {
      debug15(printf("Plus: Returning substring %d..%d against hardclip_low %d\n",
		     Substring_querystart((Substring_T) List_head(p)),Substring_queryend((Substring_T) List_head(p)),
		     hardclip_low));
      return (Substring_T) List_head(p);
    }

  } else {
#ifdef DEBUG15
    for (p = this->substrings_LtoH; p != NULL; p = List_next(p)) {
      printf("LtoH: %d..%d\n",
	     Substring_querystart((Substring_T) List_head(p)),Substring_queryend((Substring_T) List_head(p)));
    }
#endif

    p = this->substrings_LtoH;
    if (Substring_ambiguous_p((Substring_T) List_head(p)) == true) {
      p = List_next(p);
    }

    while (p != NULL && Substring_querystart((Substring_T) List_head(p)) >= this->querylength - hardclip_low) {
      debug15(printf("Minus: Skipping substring %d..%d against %d = querylength %d - hardclip_low %d\n",
		     Substring_querystart((Substring_T) List_head(p)),Substring_queryend((Substring_T) List_head(p)),
		     this->querylength - hardclip_low,this->querylength,hardclip_low));
      p = List_next(p);
    }
    assert(p != NULL);
    if (p == NULL) {
      return (Substring_T) NULL;
    } else {
      debug15(printf("Minus: Returning substring %d..%d against %d = querylength %d - hardclip_low %d\n",
		     Substring_querystart((Substring_T) List_head(p)),Substring_queryend((Substring_T) List_head(p)),
		     this->querylength - hardclip_low,this->querylength,hardclip_low));
      return (Substring_T) List_head(p);
    }
  }
}


Substring_T
Stage3end_substring_containing (T this, int querypos) {
  Substring_T substring;
  List_T p;

  for (p = this->substrings_LtoH; p != NULL; p = List_next(p)) {
    substring = (Substring_T) List_head(p);
    if (Substring_contains_p(substring,querypos) == true) {
      return substring;
    }
  }
  return (Substring_T) NULL;
}



struct Pair_T *
Stage3end_pairarray (T this) {
  return this->pairarray;
}

int
Stage3end_npairs (T this) {
  return this->npairs;
}

List_T
Stage3end_cigar_tokens (T this) {
  return this->cigar_tokens;
}

bool
Stage3end_gmap_intronp (T this) {
  return this->gmap_intronp;
}


Chrpos_T
Stage3end_distance (T this) {
  return this->distance;
}

Chrpos_T
Stage3end_shortexonA_distance (T this) {
  return this->shortexonA_distance;
}

Chrpos_T
Stage3end_shortexonD_distance (T this) {
  return this->shortexonD_distance;
}

double
Stage3end_chimera_prob (T this) {
  List_T p;
  Junction_T junction;

  for (p = this->junctions_1toN; p != NULL; p = List_next(p)) {
    junction = (Junction_T) List_head(p);
    if (Junction_type(junction) == CHIMERA_JUNCTION) {
      return Junction_prob(junction);
    }
  }

  return 0.0;
}

double
Stage3end_shortexon_prob (T this) {
  double prob = 0.0;
  List_T p;
  Junction_T junction;

  for (p = this->junctions_LtoH; p != NULL; p = List_next(p)) {
    junction = (Junction_T) List_head(p);
    prob += Junction_prob(junction);
  }

  return prob;
}

static double
Stage3end_prob (T this) {
  double prob = 0.0;
  List_T p;
  Junction_T junction;

  for (p = this->junctions_LtoH; p != NULL; p = List_next(p)) {
    junction = (Junction_T) List_head(p);
    prob += Junction_prob(junction);
  }

  return prob;
}


/* Should eventually look for substrings adjacent to the chimeric junction */
Univcoord_T
Stage3end_chimera_segmenti_left (T this) {
  Univcoord_T x_segmenti, x_segmentj;
  Substring_T substring_donor, substring_acceptor;

  if (this->sensedir == SENSE_ANTI) {
    substring_donor = (Substring_T) List_head(this->substrings_Nto1);
    substring_acceptor = (Substring_T) List_head(this->substrings_1toN);
  } else {
    substring_donor = (Substring_T) List_head(this->substrings_1toN);
    substring_acceptor = (Substring_T) List_head(this->substrings_Nto1);
  }

  x_segmenti = Substring_left_genomicseg(substring_donor);
  x_segmentj = Substring_left_genomicseg(substring_acceptor);
  if (x_segmenti < x_segmentj) {
    return x_segmenti;
  } else {
    return x_segmentj;
  }
}  

/* Should eventually look for substrings adjacent to the chimeric junction */
Univcoord_T
Stage3end_chimera_segmentj_left (T this) {
  Univcoord_T x_segmenti, x_segmentj;
  Substring_T substring_donor, substring_acceptor;

  if (this->sensedir == SENSE_ANTI) {
    substring_donor = (Substring_T) List_head(this->substrings_Nto1);
    substring_acceptor = (Substring_T) List_head(this->substrings_1toN);
  } else {
    substring_donor = (Substring_T) List_head(this->substrings_1toN);
    substring_acceptor = (Substring_T) List_head(this->substrings_Nto1);
  }

  x_segmenti = Substring_left_genomicseg(substring_donor);
  x_segmentj = Substring_left_genomicseg(substring_acceptor);
  if (x_segmenti > x_segmentj) {
    return x_segmenti;
  } else {
    return x_segmentj;
  }
}  


int
Stage3end_chimera_segmenti_cmp (const void *a, const void *b) {
  T x = * (T *) a;
  T y = * (T *) b;
  Univcoord_T x_segmenti, x_segmentj, y_segmenti, y_segmentj, temp;
  Substring_T x_substring_donor, x_substring_acceptor,
    y_substring_donor, y_substring_acceptor;

  if (x->sensedir == SENSE_ANTI) {
    x_substring_donor = (Substring_T) List_head(x->substrings_Nto1);
    x_substring_acceptor = (Substring_T) List_head(x->substrings_1toN);
  } else {
    x_substring_donor = (Substring_T) List_head(x->substrings_1toN);
    x_substring_acceptor = (Substring_T) List_head(x->substrings_Nto1);
  }

  if (y->sensedir == SENSE_ANTI) {
    y_substring_donor = (Substring_T) List_head(y->substrings_Nto1);
    y_substring_acceptor = (Substring_T) List_head(y->substrings_1toN);
  } else {
    y_substring_donor = (Substring_T) List_head(y->substrings_1toN);
    y_substring_acceptor = (Substring_T) List_head(y->substrings_Nto1);
  }

  x_segmenti = Substring_left_genomicseg(x_substring_donor);
  x_segmentj = Substring_left_genomicseg(x_substring_acceptor);
  if (x_segmentj < x_segmenti) {
    temp = x_segmentj;
    x_segmentj = x_segmenti;
    x_segmenti = temp;
  }

  y_segmenti = Substring_left_genomicseg(y_substring_donor);
  y_segmentj = Substring_left_genomicseg(y_substring_acceptor);
  if (y_segmentj < y_segmenti) {
    temp = y_segmentj;
    y_segmentj = y_segmenti;
    y_segmenti = temp;
  }

  if (x_segmenti < y_segmenti) {
    return -1;
  } else if (y_segmenti < x_segmenti) {
    return +1;
  } else if (x_segmentj > y_segmentj) {
    return -1;
  } else if (y_segmentj > x_segmentj) {
    return +1;
  } else {
    return 0;
  }
}



int
Stage3end_chimera_segmentj_cmp (const void *a, const void *b) {
  T x = * (T *) a;
  T y = * (T *) b;
  Univcoord_T x_segmenti, x_segmentj, y_segmenti, y_segmentj, temp;
  Substring_T x_substring_donor, x_substring_acceptor,
    y_substring_donor, y_substring_acceptor;

  if (x->sensedir == SENSE_ANTI) {
    x_substring_donor = (Substring_T) List_head(x->substrings_Nto1);
    x_substring_acceptor = (Substring_T) List_head(x->substrings_1toN);
  } else {
    x_substring_donor = (Substring_T) List_head(x->substrings_1toN);
    x_substring_acceptor = (Substring_T) List_head(x->substrings_Nto1);
  }

  if (y->sensedir == SENSE_ANTI) {
    y_substring_donor = (Substring_T) List_head(y->substrings_Nto1);
    y_substring_acceptor = (Substring_T) List_head(y->substrings_1toN);
  } else {
    y_substring_donor = (Substring_T) List_head(y->substrings_1toN);
    y_substring_acceptor = (Substring_T) List_head(y->substrings_Nto1);
  }


  x_segmenti = Substring_left_genomicseg(x_substring_donor);
  x_segmentj = Substring_left_genomicseg(x_substring_acceptor);
  if (x_segmentj < x_segmenti) {
    temp = x_segmentj;
    x_segmentj = x_segmenti;
    x_segmenti = temp;
  }

  y_segmenti = Substring_left_genomicseg(y_substring_donor);
  y_segmentj = Substring_left_genomicseg(y_substring_acceptor);
  if (y_segmentj < y_segmenti) {
    temp = y_segmentj;
    y_segmentj = y_segmenti;
    y_segmenti = temp;
  }

  if (x_segmentj < y_segmentj) {
    return -1;
  } else if (y_segmentj < x_segmentj) {
    return +1;
  } else if (x_segmenti > y_segmenti) {
    return -1;
  } else if (y_segmenti > x_segmenti) {
    return +1;
  } else {
    return 0;
  }
}


int
Stage3end_shortexon_substringD_cmp (const void *a, const void *b) {
  T x = * (T *) a;
  T y = * (T *) b;
  Univcoord_T x_left, y_left;
  Substring_T x_substring_donor, y_substring_donor;

  if (x->sensedir == SENSE_ANTI) {
    x_substring_donor = (Substring_T) List_head(x->substrings_Nto1);
  } else {
    x_substring_donor = (Substring_T) List_head(x->substrings_1toN);
  }

  if (y->sensedir == SENSE_ANTI) {
    y_substring_donor = (Substring_T) List_head(y->substrings_Nto1);
  } else {
    y_substring_donor = (Substring_T) List_head(y->substrings_1toN);
  }


  x_left = Substring_left_genomicseg(x_substring_donor);
  y_left = Substring_left_genomicseg(y_substring_donor);
  if (x_left < y_left) {
    return -1;
  } else if (y_left < x_left) {
    return +1;
  } else {
    return 0;
  }
}

int
Stage3end_shortexon_substringA_cmp (const void *a, const void *b) {
  T x = * (T *) a;
  T y = * (T *) b;
  Univcoord_T x_left, y_left;
  Substring_T x_substring_acceptor, y_substring_acceptor;

  if (x->sensedir == SENSE_ANTI) {
    x_substring_acceptor = (Substring_T) List_head(x->substrings_1toN);
  } else {
    x_substring_acceptor = (Substring_T) List_head(x->substrings_Nto1);
  }

  if (y->sensedir == SENSE_ANTI) {
    y_substring_acceptor = (Substring_T) List_head(y->substrings_1toN);
  } else {
    y_substring_acceptor = (Substring_T) List_head(y->substrings_Nto1);
  }


  x_left = Substring_left_genomicseg(x_substring_acceptor);
  y_left = Substring_left_genomicseg(y_substring_acceptor);
  if (x_left < y_left) {
    return -1;
  } else if (y_left < x_left) {
    return +1;
  } else {
    return 0;
  }
}




int
Stage3end_sensedir (T this) {
  if (this == NULL) {
    /* Can happen if we call upon a mate in a halfmapping */
    return SENSE_NULL;
  } else {
    return this->sensedir;
  }
}


int
Stage3end_cdna_direction (T this) {
  if (this == NULL) {
    return SENSE_NULL;
#if 0
  } else if (this->pairarray != NULL) {
    return this->gmap_cdna_direction;
#endif
  } else if (this->sensedir == SENSE_FORWARD) {
    return +1;
  } else if (this->sensedir == SENSE_ANTI) {
    return -1;
  } else {
#if 0
    /* Leads to non-canonical XS:A:? output in SAM format */
    return 0;
#else
    return this->gmap_cdna_direction;
#endif
  }
}

int
Stage3end_nintrons (T this) {
  return this->gmap_nintrons;
}

bool
Stage3end_start_ambiguous_p (T this) {
  Substring_T substring;

  substring = (Substring_T) List_head(this->substrings_1toN);
  return Substring_ambiguous_p(substring);
}

bool
Stage3end_end_ambiguous_p (T this) {
  Substring_T substring;

  substring = (Substring_T) List_head(this->substrings_Nto1);
  return Substring_ambiguous_p(substring);
}

Univcoord_T *
Stage3end_start_ambcoords (T this) {
  Substring_T substring;

  substring = (Substring_T) List_head(this->substrings_1toN);
  if (Substring_ambiguous_p(substring) == false) {
    return (Univcoord_T *) NULL;
  } else {
    return Substring_ambcoords(substring);
  }
}

Univcoord_T *
Stage3end_end_ambcoords (T this) {
  Substring_T substring;

  substring = (Substring_T) List_head(this->substrings_Nto1);
  if (Substring_ambiguous_p(substring) == false) {
    return (Univcoord_T *) NULL;
  } else {
    return Substring_ambcoords(substring);
  }
}

int
Stage3end_start_nambcoords (T this) {
  Substring_T substring;

  substring = (Substring_T) List_head(this->substrings_1toN);
  if (Substring_ambiguous_p(substring) == false) {
    return 0;
  } else {
    return Substring_nambcoords(substring);
  }
}

int
Stage3end_end_nambcoords (T this) {
  Substring_T substring;

  substring = (Substring_T) List_head(this->substrings_Nto1);
  if (Substring_ambiguous_p(substring) == false) {
    return 0;
  } else {
    return Substring_nambcoords(substring);
  }
}


int
Stage3end_substrings_querystart (T this) {
  Substring_T substring;

  substring = (Substring_T) List_head(this->substrings_1toN);
  return Substring_querystart(substring);
}

int
Stage3end_substrings_queryend (T this) {
  Substring_T substring;

  substring = (Substring_T) List_head(this->substrings_Nto1);
  return Substring_queryend(substring);
}

int
Stage3end_gmap_querystart (T this) {
  return this->pairarray[0].querypos;
}

int
Stage3end_gmap_queryend (T this) {
  return this->pairarray[this->npairs - 1].querypos;
}

int
Stage3end_terminal_trim (T this) {
  Substring_T substring;

  if (this->hittype != TERMINAL) {
    return 0;
  } else {
    substring = (Substring_T) List_head(this->substrings_LtoH);
    return Substring_trim_left(substring) + Substring_trim_right(substring);
  }
}

int
Stage3end_trimlength (T this) {
  return this->trim_left + this->trim_right;
}


static Overlap_T
Stage3end_gene_overlap (T this) {
  Overlap_T overlap;
  bool foundp = false;
  Substring_T substring;
  List_T p;

  if (this->hittype == GMAP) {
    return Pair_gene_overlap(this->pairarray,this->npairs,genes_iit,
			     genes_divint_crosstable[this->chrnum],favor_multiexon_p);
  } else {
    for (p = this->substrings_LtoH; p != NULL; p = List_next(p)) {
      substring = (Substring_T) List_head(p);
      if ((overlap = Substring_gene_overlap(substring,favor_multiexon_p)) == KNOWN_GENE_MULTIEXON) {
	return KNOWN_GENE_MULTIEXON;
      } else if (overlap == KNOWN_GENE) {
	if (favor_multiexon_p == false) {
	  return KNOWN_GENE;
	} else {
	  foundp = true;
	}
      }
    }

    if (foundp == true) {
      return KNOWN_GENE;
    } else {
      return NO_KNOWN_GENE;
    }
  }
}


bool
Stage3end_contains_known_splicesite (T this) {
  List_T p;
  Substring_T substring;

  /* assert(this->hittype != GMAP); */

  /* indel + splice => requires gmap
     doublesplice + splice => requires gmap
     other cases should have already been covered by gsnap
  */

  if (this->hittype == GMAP) {
    /* Possible now because performing redo of GMAP for sense inconsistency */
    return false;
  } else {
    for (p = this->substrings_LtoH; p != NULL; p = List_next(p)) {
      substring = (Substring_T) List_head(p);
      if (Substring_contains_known_splicesite(substring) == true) {
	return true;
      }
    }
    return false;
  }
}


bool
Stage3end_indel_contains_known_splicesite (bool *leftp, bool *rightp, T this) {
  Substring_T substring1, substring2;

  /* indel + splice => requires gmap */
  substring1 = (Substring_T) List_head(this->substrings_1toN);
  substring2 = (Substring_T) List_head(this->substrings_Nto1);

  *leftp = Substring_contains_known_splicesite(substring1);
  *rightp = Substring_contains_known_splicesite(substring2);
  if (*leftp == true || *rightp == true) {
    return true;
  } else {
    return false;
  }
}


static long int
Stage3end_compute_tally (T this) {
  long int tally = 0L;
  List_T p;
  Substring_T substring;

  for (p = this->substrings_LtoH; p != NULL; p = List_next(p)) {
    substring = (Substring_T) List_head(p);
    tally += Substring_tally(substring,tally_iit,tally_divint_crosstable);
  }

  return tally;
}

static bool
Stage3end_runlength_p (T this) {
  List_T p;
  Substring_T substring;

  for (p = this->substrings_LtoH; p != NULL; p = List_next(p)) {
    substring = (Substring_T) List_head(p);
    if (Substring_runlength_p(substring,runlength_iit,runlength_divint_crosstable) == true) {
      return true;
    }
  }

  return false;
}


bool
Stage3end_genomicbound_from_start (Univcoord_T *genomicbound, T this, int overlap, Univcoord_T chroffset) {
  int substring_length;
  List_T p;
  Substring_T substring;

  debug11(printf("Stage3end_genomicbound_from_start with overlap %d\n",overlap));
  if (this->hittype == GMAP) {
    debug11(printf("  Computing on GMAP\n"));
    *genomicbound = chroffset + Pairarray_genomicbound_from_start(this->pairarray,this->npairs,overlap);
    return true;
  } else {
    debug11(printf("  Computing on substrings\n"));
    for (p = this->substrings_1toN; p != NULL; p = List_next(p)) {
      substring = (Substring_T) List_head(p);
      debug11(printf("  Substring as length %d\n",Substring_match_length_orig(substring)));
      if ((substring_length = Substring_match_length_orig(substring)) >= overlap) {
	if (this->plusp == true) {
	  *genomicbound = Substring_alignstart(substring) + overlap;
	} else {
	  *genomicbound = Substring_alignstart(substring) - overlap;
	}
	return true;
      } else {
	overlap -= substring_length;
      }
    }

    debug11(printf("Still have %d of overlap\n",overlap));
    return false;
  }
}

bool
Stage3end_genomicbound_from_end (Univcoord_T *genomicbound, T this, int overlap, Univcoord_T chroffset) {
  int substring_length;
  List_T p;
  Substring_T substring;

  debug11(printf("Stage3end_genomicbound_from_end with overlap %d\n",overlap));
  if (this->hittype == GMAP) {
    debug11(printf("  Computing on GMAP\n"));
    *genomicbound = chroffset + Pairarray_genomicbound_from_end(this->pairarray,this->npairs,overlap);
    return true;
  } else {
    debug11(printf("  Computing on substrings\n"));
    for (p = this->substrings_Nto1; p != NULL; p = List_next(p)) {
      substring = (Substring_T) List_head(p);
      debug11(printf("  Substring has length %d\n",Substring_match_length_orig(substring)));
      if ((substring_length = Substring_match_length_orig(substring)) >= overlap) {
	if (this->plusp == true) {
	  *genomicbound = Substring_alignend(substring) - overlap;
	} else {
	  *genomicbound = Substring_alignend(substring) + overlap;
	}
	return true;
      } else {
	overlap -= substring_length;
      }
    }

    debug11(printf("Still have %d of overlap\n",overlap));
    return false;
  }
}


void
Stage3end_free (T *old) {
  List_T p;
  Substring_T substring;
  Junction_T junction;

  debug0(printf("Freeing Stage3end %p of type %s\n",*old,hittype_string((*old)->hittype)));

#if 0
  FREE_OUT((*old)->ambcoords_donor);
  FREE_OUT((*old)->ambcoords_acceptor);
  FREE_OUT((*old)->amb_knowni_donor);
  FREE_OUT((*old)->amb_knowni_acceptor);
  FREE_OUT((*old)->amb_nmismatches_donor);
  FREE_OUT((*old)->amb_nmismatches_acceptor);
  FREE_OUT((*old)->amb_probs_donor);
  FREE_OUT((*old)->amb_probs_acceptor);
#endif

  if ((*old)->cigar_tokens != NULL) {
    Pair_tokens_free(&(*old)->cigar_tokens);
  }

  if ((*old)->pairarray != NULL) {
    FREE_OUT((*old)->pairarray);
  }

  for (p = (*old)->substrings_1toN; p != NULL; p = List_next(p)) {
    substring = (Substring_T) List_head(p);
    Substring_free(&substring);
  }
  List_free(&(*old)->substrings_1toN);
  List_free(&(*old)->substrings_Nto1);
  List_free(&(*old)->substrings_LtoH);

  for (p = (*old)->junctions_1toN; p != NULL; p = List_next(p)) {
    junction = (Junction_T) List_head(p);
    Junction_free(&junction);
  }
  List_free(&(*old)->junctions_1toN);
  List_free(&(*old)->junctions_Nto1);
  List_free(&(*old)->junctions_LtoH);

  FREE_OUT(*old);
  return;
}


/* Used for freeing gmap_history_values in stage1hr.c */
void
Stage3end_list_free (List_T *values) {
  List_T p;
  T hit;

  for (p = *values; p != NULL; p = p->rest) {
    if ((hit = (T) p->first) != NULL) {
      Stage3end_free(&hit);
    }
  }
  List_free(&(*values));
  return;
}



bool
Stage3pair_anomalous_splice_p (Stage3pair_T this) {
  if (this->hit5 != NULL && (this->hit5->chrnum == 0 || this->hit5->hittype == SAMECHR_SPLICE)) {
    return true;
  } else if (this->hit3 != NULL && (this->hit3->chrnum == 0 || this->hit3->hittype == SAMECHR_SPLICE)) {
    return true;
  } else {
    return false;
  }
}


int
Stage3pair_genestrand (Stage3pair_T this) {
  return this->genestrand;
}

Stage3end_T
Stage3pair_hit5 (Stage3pair_T this) {
  return this->hit5;
}

Stage3end_T
Stage3pair_hit3 (Stage3pair_T this) {
  return this->hit3;
}

int
Stage3pair_mapq_score (Stage3pair_T this) {
  return this->mapq_score;
}

int
Stage3pair_absmq_score (Stage3pair_T this) {
  return this->absmq_score;
}

Chrpos_T
Stage3pair_pairlength (Stage3pair_T this) {
  return this->insertlength;
}

int
Stage3pair_nmatches_posttrim (int *nmatches5, int *nmatches3, Stage3pair_T this) {
  *nmatches5 = this->hit5->nmatches_posttrim;
  *nmatches3 = this->hit3->nmatches_posttrim;
  return this->nmatches_posttrim;
}


bool
Stage3pair_concordantp (List_T hitpairs) {
  List_T p;
  Stage3pair_T hitpair;

  for (p = hitpairs; p != NULL; p = List_next(p)) {
    hitpair = (Stage3pair_T) List_head(p);
    if (hitpair->pairtype == CONCORDANT) {
      return true;
    }
  }
  return false;
}

List_T
Stage3pair_filter_nonconcordant (List_T hitpairs) {
  List_T filtered = NULL, p;
  Stage3pair_T hitpair;

  for (p = hitpairs; p != NULL; p = List_next(p)) {
    hitpair = (Stage3pair_T) List_head(p);
    if (hitpair->pairtype != CONCORDANT) {
      Stage3pair_free(&hitpair);
    } else {
      filtered = List_push(filtered,(void *) hitpair);
    }
  }
  List_free(&hitpairs);
  return filtered;
}


static Univcoord_T
gmap5_substring3_common_genomicpos (Stage3end_T hit5, Stage3end_T hit3, Substring_T substring) {
  Univcoord_T chroffset;
  Chrpos_T start, end;
  int i;

  chroffset = hit3->chroffset;
  if (hit5->plusp == true) {
    start = Substring_alignstart_trim(substring) - chroffset;
    end = Substring_alignend_trim(substring) - 1U - chroffset; /* inclusive */
    debug15(printf("plus goal: %u up to %u\n",start,end));
    i = 0;
    while (i < hit5->npairs) {
      debug15(printf("  pair %d: genomepos %u\n",i,hit5->pairarray[i].genomepos));
      if (hit5->pairarray[i].gapp == true) {
	/* Skip intron */
	i++;
      } else if (hit5->pairarray[i].cdna == ' ' || hit5->pairarray[i].genome == ' ') {
	/* Skip indel */
	i++;
      } else if (hit5->pairarray[i].genomepos < start) {
	i++;
      } else if (hit5->pairarray[i].genomepos > end) {
	i++;
      } else {
	debug15(printf("Returning common point at %llu\n",(unsigned long long) hit5->pairarray[i].genomepos));
	return hit5->pairarray[i].genomepos + chroffset;
      }
    }
    return 0;

  } else {
    start = Substring_alignstart_trim(substring) - 1U - chroffset; /* inclusive */
    end = Substring_alignend_trim(substring) - chroffset;
    debug15(printf("minus goal: %u up to %u\n",end,start));
    i = hit5->npairs - 1;
    while (i >= 0) {
      debug15(printf("  pair %d: genomepos %u\n",i,hit5->pairarray[i].genomepos));
      if (hit5->pairarray[i].gapp == true) {
	/* Skip intron */
	i--;
      } else if (hit5->pairarray[i].cdna == ' ' || hit5->pairarray[i].genome == ' ') {
	/* Skip indel */
	i--;
      } else if (hit5->pairarray[i].genomepos > start) {
	i--;
      } else if (hit5->pairarray[i].genomepos < end) {
	i--;
      } else {
	debug15(printf("Returning common point at %llu\n",(unsigned long long) hit5->pairarray[i].genomepos));
	return hit5->pairarray[i].genomepos + chroffset;
      }
    }
    return 0;
  }
}

static Univcoord_T
substring5_gmap3_common_genomicpos (Stage3end_T hit5, Stage3end_T hit3, Substring_T substring) {
  Univcoord_T chroffset;
  Chrpos_T start, end;
  int j;

  chroffset = hit5->chroffset;
  if (hit5->plusp == true) {
    start = Substring_alignstart_trim(substring) - chroffset;
    end = Substring_alignend_trim(substring) - 1U - chroffset; /* inclusive */
    debug15(printf("plus goal: %u up to %u\n",start,end));

    j = 0;
    while (j < hit3->npairs) {
      debug15(printf("  pair %d: genomepos %u\n",j,hit3->pairarray[j].genomepos));
      if (hit3->pairarray[j].gapp == true) {
	/* Skip intron */
	j++;
      } else if (hit3->pairarray[j].cdna == ' ' || hit3->pairarray[j].genome == ' ') {
	/* Skip indel */
	j++;
      } else if (hit3->pairarray[j].genomepos < start) {
	j++;
      } else if (hit3->pairarray[j].genomepos > end) {
	j++;
      } else {
	debug15(printf("Returning common point at %llu\n",(unsigned long long) hit3->pairarray[j].genomepos));
	return hit3->pairarray[j].genomepos + chroffset;
      }
    }
    return 0;

  } else {
    start = Substring_alignstart_trim(substring) - 1U - chroffset; /* inclusive */
    end = Substring_alignend_trim(substring) - chroffset;
    debug15(printf("minus goal: %u up to %u\n",end,start));
    j = hit3->npairs - 1;
    while (j >= 0) {
      debug15(printf("  pair %d: genomepos %u\n",j,hit3->pairarray[j].genomepos));
      if (hit3->pairarray[j].gapp == true) {
	/* Skip intron */
	j--;
      } else if (hit3->pairarray[j].cdna == ' ' || hit3->pairarray[j].genome == ' ') {
	/* Skip indel */
	j--;
      } else if (hit3->pairarray[j].genomepos > start) {
	j--;
      } else if (hit3->pairarray[j].genomepos < end) {
	j--;
      } else {
	debug15(printf("Returning common point at %llu\n",(unsigned long long) hit3->pairarray[j].genomepos));
	return hit3->pairarray[j].genomepos + chroffset;
      }
    }
    return 0;
  }
}


/* Returns true if ilengths are valid */
static bool
find_ilengths (int *ilength_low, int *ilength_high, Stage3end_T hit, Univcoord_T common_genomicpos, Univcoord_T chroffset) {
  int i;
  List_T p, q;
  Substring_T substring;
  Junction_T junction;


  debug15(printf("Finding ilengths for common_genomicpos %u\n",(Chrpos_T) (common_genomicpos - chroffset)));
  if (hit->hittype == GMAP) {
    debug15(printf("Type is GMAP\n"));
    debug15(Pair_dump_array(hit->pairarray,hit->npairs,true));
    i = 0;
    while (i < hit->npairs && hit->pairarray[i].genomepos != common_genomicpos - chroffset) {
      i++;
    }
    if (i >= hit->npairs) {
      return false;
    } else if (hit->plusp == true) {
      *ilength_low = hit->pairarray[i].querypos - hit->pairarray[0].querypos + 1;
      *ilength_high = hit->pairarray[hit->npairs - 1].querypos - hit->pairarray[i].querypos + 1;
    } else {
      *ilength_low = hit->pairarray[hit->npairs - 1].querypos - hit->pairarray[i].querypos + 1;
      *ilength_high = hit->pairarray[i].querypos - hit->pairarray[0].querypos + 1;
    }
    debug15(printf("GMAP: Have ilength_low %d and ilength_high %d\n",*ilength_low,*ilength_high));
    return true;

  } else if (hit->plusp == true) {
#ifdef DEBUG15
    printf("plus.  Checking common genomicpos %llu against\n",
	   common_genomicpos - hit->chroffset);
    for (p = hit->substrings_1toN; p != NULL; p = List_next(p)) {
      substring = (Substring_T) List_head(p);
      printf("substring %p: %u..%u\n",
	     substring,Substring_alignstart_trim(substring) - hit->chroffset,
	     Substring_alignend_trim(substring) - 1U - hit->chroffset);
    }
    printf("\n");
#endif
    /* Plus: Subtract 1 from alignend */
    *ilength_low = 0;
    for (p = hit->substrings_1toN, q = hit->junctions_1toN; p != NULL; p = List_next(p), q = List_next(q)) {
      substring = (Substring_T) List_head(p);
      debug15(printf("substring %p: %u..%u\n",substring,
		     Substring_alignstart_trim(substring) - hit->chroffset,
		     Substring_alignend_trim(substring) - 1U - hit->chroffset));
      if (Substring_overlap_point_trimmed_p(substring,common_genomicpos) == false) {
	*ilength_low += Substring_genomic_alignment_length(substring);
	if (q != NULL) {
	  junction = (Junction_T) List_head(q);
	  if (Junction_type(junction) == INS_JUNCTION) {
	    *ilength_low += Junction_nindels(junction);
	  }
	}

      } else {
	*ilength_low += (common_genomicpos - Substring_alignstart_trim(substring) + 1);
	*ilength_high = ((Substring_alignend_trim(substring) - 1) - common_genomicpos + 1);
	p = List_next(p);
	while (p != NULL) {
	  substring = (Substring_T) List_head(p);
	  *ilength_high += Substring_genomic_alignment_length(substring);
	  p = List_next(p);
	}
	while (q != NULL) {
	  junction = (Junction_T) List_head(q);
	  if (Junction_type(junction) == INS_JUNCTION) {
	    *ilength_high += Junction_nindels(junction);
	  }
	  q = List_next(q);
	}
	debug15(printf("Plus: Have ilength_low %d and ilength_high %d\n",*ilength_low,*ilength_high));
	return true;
      }
    }
  } else {
#ifdef DEBUG15
    printf("minus.  Checking common genomicpos %llu against\n",
      common_genomicpos - hit->chroffset);
    for (p = hit->substrings_1toN; p != NULL; p = List_next(p)) {
      substring = (Substring_T) List_head(p);
      printf("substring %p: %u..%u\n",
	     substring,Substring_alignstart_trim(substring) - hit->chroffset,
	     Substring_alignend_trim(substring) - 1U - hit->chroffset);
    }
    printf("\n");
#endif
    /* Minus: Subtract 1 from alignstart */
    *ilength_high = 0;
    for (p = hit->substrings_1toN, q = hit->junctions_1toN; p != NULL; p = List_next(p), q = List_next(q)) {
      substring = (Substring_T) List_head(p);
      debug15(printf("substring: %u..%u\n",
		     Substring_alignstart_trim(substring) - 1U - hit->chroffset,
		     Substring_alignend_trim(substring) - hit->chroffset));
      if (Substring_overlap_point_trimmed_p(substring,common_genomicpos) == false) {
	*ilength_high += Substring_genomic_alignment_length(substring);
	if (q != NULL) {
	  junction = (Junction_T) List_head(q);
	  if (Junction_type(junction) == INS_JUNCTION) {
	    *ilength_high += Junction_nindels(junction);
	  }
	}

      } else {
	*ilength_high += ((Substring_alignstart_trim(substring) - 1) - common_genomicpos + 1);
	*ilength_low = (common_genomicpos - (Substring_alignend_trim(substring) /*+ 1*/) + 1);
	p = List_next(p);
	while (p != NULL) {
	  substring = (Substring_T) List_head(p);
	  *ilength_low += Substring_genomic_alignment_length(substring);
	  p = List_next(p);
	}
	while (q != NULL) {
	  junction = (Junction_T) List_head(q);
	  if (Junction_type(junction) == INS_JUNCTION) {
	    *ilength_low += Junction_nindels(junction);
	  }
	  q = List_next(q);
	}
	debug15(printf("Minus: Have ilength_low %d and ilength_high %d\n",*ilength_low,*ilength_high));
	return true;
      }
    }
  }

  return false;
}



/* Needed to compute overlap properly.  Based on pair_insert_length below, plus code for handling GMAP. */
static Univcoord_T
pair_common_genomicpos (Stage3end_T hit5, Stage3end_T hit3) {
  Univcoord_T common_genomicpos;
  int i, j;
  Univcoord_T start5, end5, start3, end3;
  List_T p, q;
  Substring_T substring, substring5, substring3;

  if (hit5->hittype == GMAP && hit3->hittype == GMAP) {
    debug15(printf("Computing overlap using dual GMAP\n"));
    if (hit5->plusp == true) {
      i = j = 0;
      while (i < hit5->npairs && j < hit3->npairs) {
	if (hit5->pairarray[i].gapp == true) {
	  /* Skip intron */
	  i++;
	} else if (hit5->pairarray[i].cdna == ' ' || hit5->pairarray[i].genome == ' ') {
	  /* Skip indel */
	  i++;
	} else if (hit3->pairarray[j].gapp == true) {
	  /* Skip intron */
	  j++;
	} else if (hit3->pairarray[j].cdna == ' ' || hit3->pairarray[j].genome == ' ') {
	  /* Skip indel */
	  j++;
	} else if (hit5->pairarray[i].genomepos < hit3->pairarray[j].genomepos) {
	  i++;
	} else if (hit5->pairarray[i].genomepos > hit3->pairarray[j].genomepos) {
	  j++;
	} else {
	  debug15(printf("GMAP and GMAP show overlap at position %d, querypos %d and %d\n",
			 hit5->pairarray[i].genomepos,hit5->pairarray[i].querypos,hit3->pairarray[j].querypos));
	  return hit5->pairarray[i].genomepos + hit5->chroffset;
	}
      }
      debug15(printf("GMAP and GMAP show no overlap\n"));
      return 0U;

    } else {
      i = j = 0;
      while (i < hit5->npairs && j < hit3->npairs) {
	if (hit5->pairarray[i].gapp == true) {
	  /* Skip intron */
	  i++;
	} else if (hit5->pairarray[i].cdna == ' ' || hit5->pairarray[i].genome == ' ') {
	  /* Skip indel */
	  i++;
	} else if (hit3->pairarray[j].gapp == true) {
	  /* Skip intron */
	  j++;
	} else if (hit3->pairarray[j].cdna == ' ' || hit3->pairarray[j].genome == ' ') {
	  /* Skip indel */
	  j++;
	} else if (hit5->pairarray[i].genomepos > hit3->pairarray[j].genomepos) {
	  i++;
	} else if (hit5->pairarray[i].genomepos < hit3->pairarray[j].genomepos) {
	  j++;
	} else {
	  debug15(printf("GMAP and GMAP show overlap at position %d, querypos %d and %d\n",
			 hit5->pairarray[i].genomepos,hit5->pairarray[i].querypos,hit3->pairarray[j].querypos));
	  return hit5->pairarray[i].genomepos + hit5->chroffset;
	}
      }
      debug15(printf("GMAP and GMAP show no overlap\n"));
      return 0U;
    }

  } else if (hit5->hittype == GMAP) {
    debug15(printf("Computing common point using 5' GMAP\n"));
    for (p = hit3->substrings_LtoH; p != NULL; p = List_next(p)) {
      substring = (Substring_T) List_head(p);
      if ((common_genomicpos = gmap5_substring3_common_genomicpos(hit5,hit3,substring)) != 0) {
	return common_genomicpos;
      }
    }
    return 0U;

  } else if (hit3->hittype == GMAP) {
    debug15(printf("Computing common point using 3' GMAP\n"));
    for (p = hit5->substrings_LtoH; p != NULL; p = List_next(p)) {
      substring = (Substring_T) List_head(p);
      if ((common_genomicpos = substring5_gmap3_common_genomicpos(hit5,hit3,substring)) != 0) {
	return common_genomicpos;
      }
    }
    return 0U;

  } else if (hit5->plusp == true && hit3->plusp == true) {
    /* plus/plus */
    debug15(printf("Computing overlap using substrings plus/plus\n"));

    start5 = hit5->genomicstart + hit5->trim_left + start_amb_length(hit5);
    end5 = (hit5->genomicend - 1) - hit5->trim_right - end_amb_length(hit5);
    start3 = hit3->genomicstart + hit3->trim_left + start_amb_length(hit3);
    end3 = (hit3->genomicend - 1) - hit3->trim_right - end_amb_length(hit3);
    debug15(printf("hit5 endpoints are %u..%u.  hit3 endpoints are %u..%u\n",
		   start5-hit5->chroffset,end5-hit5->chroffset,start3-hit3->chroffset,end3-hit3->chroffset));

    if (end3 < start5) {
      /* Case 1 */
      return false;
    } else if (end5 < start3) {
      /* Case 6 */
      return false;
    } else if (start3 < start5) {
      if (end3 < end5) {
	/* Case 2: Tails overlap.  Go from start5 to end3 */
	debug15(printf("plus/plus case 2a: start5 %u\n",start5 - hit5->chroffset));
	for (p = hit3->substrings_1toN; p != NULL; p = List_next(p)) {
	  substring = (Substring_T) List_head(p);
	  if (Substring_overlap_point_trimmed_p(substring,start5)) {
	    return start5;
	  }
	}

	/* Case 2: Tails overlap.  Go from start5 to end3 */
	debug15(printf("plus/plus case 2b: end3 %u\n",end3 - hit3->chroffset));
	for (p = hit5->substrings_Nto1; p != NULL; p = List_next(p)) {
	  substring = (Substring_T) List_head(p);
	  if (Substring_overlap_point_trimmed_p(substring,end3)) {
	    return end3;
	  }
	}
	/* Fall through to general algorithm */

      } else {
	/* Case 3: hit3 subsumes hit5 */
	debug15(printf("plus/plus case 3\n"));
	for (p = hit3->substrings_Nto1; p != NULL; p = List_next(p)) {
	  substring = (Substring_T) List_head(p);
	  if (Substring_overlap_point_trimmed_p(substring,end5)) {
	    return end5;
	  }
	}
	/* Fall through to general algorithm */
      }

    } else {
      if (end3 < end5) {
	/* Case 4: hit5 subsumes hit3 */
	debug15(printf("plus/plus case 4\n"));
	for (p = hit5->substrings_1toN; p != NULL; p = List_next(p)) {
	  substring = (Substring_T) List_head(p);
	  if (Substring_overlap_point_trimmed_p(substring,start3)) {
	    return start3;
	  }
	}
	/* Fall through to general algorithm */

      } else {
	/* Case 5: Based on hit3_trimmed_length */
	debug15(printf("plus/plus case 5a\n"));
	for (p = hit5->substrings_1toN; p != NULL; p = List_next(p)) {
	  substring = (Substring_T) List_head(p);
	  if (Substring_overlap_point_trimmed_p(substring,start3)) {
	    return start3;
	  }
	}

	/* Case 5: Based on hit5_trimmed_length */
	debug15(printf("plus/plus case 5b\n"));
	for (p = hit3->substrings_Nto1; p != NULL; p = List_next(p)) {
	  substring = (Substring_T) List_head(p);
	  if (Substring_overlap_point_trimmed_p(substring,end5)) {
	    return end5;
	  }
	}
	/* Fall through to general algorithm */
      }
    }

    /* General algorithm */
    debug15(printf("plus/plus general\n"));
    for (p = hit3->substrings_1toN; p != NULL; p = List_next(p)) {
      substring3 = (Substring_T) List_head(p);
      for (q = hit5->substrings_1toN; q != NULL; q = List_next(q)) {
	substring5 = (Substring_T) List_head(q);
	if ((common_genomicpos = Substring_overlap_segment_trimmed(substring5,substring3)) != 0) {
	  return common_genomicpos;
	}
      }
    }

    return 0;

  } else if (hit5->plusp == true && hit3->plusp == false) {
    /* plus/minus */
    debug15(printf("Computing overlap using substrings plus/minus\n"));
    return 0;

#if 0
    start5 = hit5->genomicstart + hit5->trim_left + start_amb_length(hit5);
    end5 = hit5->genomicend - hit5->trim_right - end_amb_length(hit5);
    start3 = hit3->genomicstart - hit3->trim_left - start_amb_length(hit3);
    end3 = hit3->genomicend + hit3->trim_right + end_amb_length(hit3);

    if (start3 < start5) {
      /* Case 1 */
      return 0;
    } else if (end5 < end3) {
      /* Case 6 */
      return 0;
    } else if (end3 < start5) {
      if (start3 < end5) {
	/* Case 2: Tails overlap.  Go from start5 to start3 */
	debug15(printf("plus case 2a: start5 %u\n",start5 - hit5->chroffset));
	if (Substring_overlap_point_trimmed_p(hit3->substring0,start5)) {
	  return start5;
	} else if (Substring_overlap_point_trimmed_p(hit3->substring1,start5)) {
	  return start5;
	} else if (Substring_overlap_point_trimmed_p(hit3->substring2,start5)) {
	  return start5;
	}

	/* Case 2: Tails overlap.  Go from start5 to start3 */
	debug15(printf("plus case 2b: start3 %u\n",start3 - hit3->chroffset));
	if (Substring_overlap_point_trimmed_p(hit5->substring2,start3)) {
	  return start3;
	} else if (Substring_overlap_point_trimmed_p(hit5->substring1,start3)) {
	  return start3;
	} else if (Substring_overlap_point_trimmed_p(hit5->substring0,start3)) {
	  return start3;
	}
	/* Fall through to general algorithm */

      } else {
	/* Case 3: hit3 subsumes hit5 */
	debug15(printf("plus case 3\n"));
	if (Substring_overlap_point_trimmed_p(hit3->substring2,end5)) {
	  return end5;
	} else if (Substring_overlap_point_trimmed_p(hit3->substring1,end5)) {
	  return end5;
	} else if (Substring_overlap_point_trimmed_p(hit3->substring0,end5)) {
	  return end5;
	}
	/* Fall through to general algorithm */
      }

    } else {
      if (start3 < end5) {
	/* Case 4: hit5 subsumes hit3 */
	debug15(printf("plus case 4\n"));
	if (Substring_overlap_point_trimmed_p(hit5->substring0,end3)) {
	  return end3;
	} else if (Substring_overlap_point_trimmed_p(hit5->substring1,end3)) {
	  return end3;
	} else if (Substring_overlap_point_trimmed_p(hit5->substring2,end3)) {
	  return end3;
	}
	/* Fall through to general algorithm */

      } else {
	/* Case 5: Based on hit3_trimmed_length */
	debug15(printf("plus case 5a\n"));
	if (Substring_overlap_point_trimmed_p(hit5->substring0,end3)) {
	  return end3;
	} else if (Substring_overlap_point_trimmed_p(hit5->substring1,end3)) {
	  return end3;
	} else if (Substring_overlap_point_trimmed_p(hit5->substring2,end3)) {
	  return end3;
	}

	/* Case 5: Based on hit5_trimmed_length */
	debug15(printf("plus case 5b\n"));
	if (Substring_overlap_point_trimmed_p(hit3->substring2,end5)) {
	  return end5;
	} else if (Substring_overlap_point_trimmed_p(hit3->substring1,end5)) {
	  return end5;
	} else if (Substring_overlap_point_trimmed_p(hit3->substring0,end5)) {
	  return end5;
	}
	/* Fall through to general algorithm */
      }
    }

    /* General algorithm */
    debug15(printf("plus general: hit3->substring1\n"));
    if ((common_genomicpos = Substring_overlap_segment_trimmed(hit5->substring1,hit3->substring1)) != 0) {
      return common_genomicpos;
    } else if (hit5->substring2 != NULL &&
	       (common_genomicpos = Substring_overlap_segment_trimmed(hit5->substring2,hit3->substring1)) != 0) {
      return common_genomicpos;
    } else if (hit5->substring0 != NULL &&
	       (common_genomicpos = Substring_overlap_segment_trimmed(hit5->substring0,hit3->substring1)) != 0) {
      return common_genomicpos;
    }

    if (hit3->substring2 != NULL) {
      debug15(printf("plus general: hit3->substring2\n"));
      if ((common_genomicpos = Substring_overlap_segment_trimmed(hit5->substring1,hit3->substring2)) != 0) {
	return common_genomicpos;
      } else if (hit5->substring2 != NULL &&
		 (common_genomicpos = Substring_overlap_segment_trimmed(hit5->substring2,hit3->substring2)) != 0) {
	return common_genomicpos;
      } else if (hit5->substring0 != NULL &&
		 (common_genomicpos = Substring_overlap_segment_trimmed(hit5->substring0,hit3->substring2)) != 0) {
	return common_genomicpos;
      }
    }

    if (hit3->substring0 != NULL) {
      debug15(printf("plus general: hit3->substring0\n"));
      if ((common_genomicpos = Substring_overlap_segment_trimmed(hit5->substring1,hit3->substring0)) != 0) {
	return common_genomicpos;
      } else if (hit5->substring2 != NULL &&
		 (common_genomicpos = Substring_overlap_segment_trimmed(hit5->substring2,hit3->substring0)) != 0) {
	return common_genomicpos;
      } else if (hit5->substring0 != NULL &&
		 (common_genomicpos = Substring_overlap_segment_trimmed(hit5->substring0,hit3->substring0)) != 0) {
	return common_genomicpos;
      }
    }

    return 0U;
#endif

  } else if (hit5->plusp == false && hit3->plusp == true) {
    /* minus/plus */
    debug15(printf("Computing overlap using substrings minus/plus\n"));
    return 0;

#if 0
    start5 = hit5->genomicstart - hit5->trim_left - start_amb_length(hit5);
    end5 = hit5->genomicend + hit5->trim_right + end_amb_length(hit5);
    start3 = hit3->genomicstart + hit3->trim_left + start_amb_length(hit3);
    end3 = hit3->genomicend - hit3->trim_right - end_amb_length(hit3);

    if (end3 < end5) {
      /* Case 1 */
      return 0;
    } else if (start5 < start3) {
      /* Case 6 */
      return 0;
    } else if (start3 < end5) {
      if (end3 < start5) {
	/* Case 2: Tails overlap.  Go from end5 to end3 */
	debug15(printf("plus case 2a: end5 %u\n",end5 - hit5->chroffset));
	if (Substring_overlap_point_trimmed_p(hit3->substring0,end5)) {
	  return end5;
	} else if (Substring_overlap_point_trimmed_p(hit3->substring1,end5)) {
	  return end5;
	} else if (Substring_overlap_point_trimmed_p(hit3->substring2,end5)) {
	  return end5;
	}

	/* Case 2: Tails overlap.  Go from end5 to end3 */
	debug15(printf("plus case 2b: end3 %u\n",end3 - hit3->chroffset));
	if (Substring_overlap_point_trimmed_p(hit5->substring2,end3)) {
	  return end3;
	} else if (Substring_overlap_point_trimmed_p(hit5->substring1,end3)) {
	  return end3;
	} else if (Substring_overlap_point_trimmed_p(hit5->substring0,end3)) {
	  return end3;
	}
	/* Fall through to general algorithm */

      } else {
	/* Case 3: hit3 subsumes hit5 */
	debug15(printf("plus case 3\n"));
	if (Substring_overlap_point_trimmed_p(hit3->substring2,start5)) {
	  return start5;
	} else if (Substring_overlap_point_trimmed_p(hit3->substring1,start5)) {
	  return start5;
	} else if (Substring_overlap_point_trimmed_p(hit3->substring0,start5)) {
	  return start5;
	}
	/* Fall through to general algorithm */
      }

    } else {
      if (end3 < start5) {
	/* Case 4: hit5 subsumes hit3 */
	debug15(printf("plus case 4\n"));
	if (Substring_overlap_point_trimmed_p(hit5->substring0,start3)) {
	  return start3;
	} else if (Substring_overlap_point_trimmed_p(hit5->substring1,start3)) {
	  return start3;
	} else if (Substring_overlap_point_trimmed_p(hit5->substring2,start3)) {
	  return start3;
	}
	/* Fall through to general algorithm */

      } else {
	/* Case 5: Based on hit3_trimmed_length */
	debug15(printf("plus case 5a\n"));
	if (Substring_overlap_point_trimmed_p(hit5->substring0,start3)) {
	  return start3;
	} else if (Substring_overlap_point_trimmed_p(hit5->substring1,start3)) {
	  return start3;
	} else if (Substring_overlap_point_trimmed_p(hit5->substring2,start3)) {
	  return start3;
	}

	/* Case 5: Based on hit5_trimmed_length */
	debug15(printf("plus case 5b\n"));
	if (Substring_overlap_point_trimmed_p(hit3->substring2,start5)) {
	  return start5;
	} else if (Substring_overlap_point_trimmed_p(hit3->substring1,start5)) {
	  return start5;
	} else if (Substring_overlap_point_trimmed_p(hit3->substring0,start5)) {
	  return start5;
	}
	/* Fall through to general algorithm */
      }
    }

    /* General algorithm */
    debug15(printf("plus general: hit3->substring1\n"));
    if ((common_genomicpos = Substring_overlap_segment_trimmed(hit5->substring1,hit3->substring1)) != 0) {
      return common_genomicpos;
    } else if (hit5->substring2 != NULL &&
	       (common_genomicpos = Substring_overlap_segment_trimmed(hit5->substring2,hit3->substring1)) != 0) {
      return common_genomicpos;
    } else if (hit5->substring0 != NULL &&
	       (common_genomicpos = Substring_overlap_segment_trimmed(hit5->substring0,hit3->substring1)) != 0) {
      return common_genomicpos;
    }

    if (hit3->substring2 != NULL) {
      debug15(printf("plus general: hit3->substring2\n"));
      if ((common_genomicpos = Substring_overlap_segment_trimmed(hit5->substring1,hit3->substring2)) != 0) {
	return common_genomicpos;
      } else if (hit5->substring2 != NULL &&
		 (common_genomicpos = Substring_overlap_segment_trimmed(hit5->substring2,hit3->substring2)) != 0) {
	return common_genomicpos;
      } else if (hit5->substring0 != NULL &&
		 (common_genomicpos = Substring_overlap_segment_trimmed(hit5->substring0,hit3->substring2)) != 0) {
	return common_genomicpos;
      }
    }

    if (hit3->substring0 != NULL) {
      debug15(printf("plus general: hit3->substring0\n"));
      if ((common_genomicpos = Substring_overlap_segment_trimmed(hit5->substring1,hit3->substring0)) != 0) {
	return common_genomicpos;
      } else if (hit5->substring2 != NULL &&
		 (common_genomicpos = Substring_overlap_segment_trimmed(hit5->substring2,hit3->substring0)) != 0) {
	return common_genomicpos;
      } else if (hit5->substring0 != NULL &&
		 (common_genomicpos = Substring_overlap_segment_trimmed(hit5->substring0,hit3->substring0)) != 0) {
	return common_genomicpos;
      }
    }

    return 0;
#endif

  } else if (hit5->plusp == false && hit3->plusp == false) {
    /* minus/minus */
    debug15(printf("Computing overlap using substrings minus/minus\n"));

    start5 = (hit5->genomicstart - 1) - hit5->trim_left - start_amb_length(hit5);
    end5 = hit5->genomicend + hit5->trim_right + end_amb_length(hit5);
    start3 = (hit3->genomicstart - 1) - hit3->trim_left - start_amb_length(hit3);
    end3 = hit3->genomicend + hit3->trim_right + end_amb_length(hit3);
    debug15(printf("hit5 endpoints are %u..%u.  hit3 endpoints are %u..%u\n",
		   start5-hit5->chroffset,end5-hit5->chroffset,start3-hit3->chroffset,end3-hit3->chroffset));

    if (end3 > start5) {
      /* Case 1 */
      return 0;
    } else if (end5 > start3) {
      /* Case 6 */
      return 0;
    } else if (start3 > start5) {
      if (end3 > end5) {
	/* Case 2: Tails overlap.  Go from start5 to end3 */
	debug15(printf("minus/minus case 2a: start5 %llu (%u)\n",start5,start5 - hit5->chroffset));
	for (p = hit3->substrings_1toN; p != NULL; p = List_next(p)) {
	  substring = (Substring_T) List_head(p);
	  if (Substring_overlap_point_trimmed_p(substring,start5)) {
	    return start5;
	  }
	}

	/* Case 2: Tails overlap.  Go from start5 to end3 */
	debug15(printf("plus case 2b: end3 %u\n",end3 - hit3->chroffset));
	for (p = hit5->substrings_Nto1; p != NULL; p = List_next(p)) {
	  substring = (Substring_T) List_head(p);
	  if (Substring_overlap_point_trimmed_p(substring,end3)) {
	    return end3;
	  }
	}
	/* Fall through to general algorithm */

      } else {
	/* Case 3: hit3 subsumes hit5 */
	debug15(printf("minus/minus case 3: end5 %u\n",end5 - hit5->chroffset));
	for (p = hit3->substrings_1toN; p != NULL; p = List_next(p)) {
	  substring = (Substring_T) List_head(p);
	  if (Substring_overlap_point_trimmed_p(substring,end5)) {
	    return end5;
	  }
	}

	/* Fall through to general algorithm */
      }

    } else {
      if (end3 > end5) {
	/* Case 4: hit5 subsumes hit3 */
	debug15(printf("minus/minus case 4: start3 %u\n",(Chrpos_T) (start3 - hit3->chroffset)));
	for (p = hit5->substrings_1toN; p != NULL; p = List_next(p)) {
	  substring = (Substring_T) List_head(p);
	  if (Substring_overlap_point_trimmed_p(substring,start3)) {
	    return start3;
	  }
	}
	/* Fall through to general algorithm */

      } else {
	/* Case 5: Based on hit3_trimmed_length */
	debug15(printf("minus case 5a: start3 %u\n",start3 - hit3->chroffset));
	for (p = hit5->substrings_1toN; p != NULL; p = List_next(p)) {
	  substring = (Substring_T) List_head(p);
	  if (Substring_overlap_point_trimmed_p(substring,start3)) {
	    return start3;
	  }
	}

	/* Case 5: Based on hit5_trimmed_length */
	debug15(printf("minus case 5b: end5 %u\n",end5 - hit5->chroffset));
	for (p = hit3->substrings_Nto1; p != NULL; p = List_next(p)) {
	  substring = (Substring_T) List_head(p);
	  if (Substring_overlap_point_trimmed_p(substring,end5)) {
	    return end5;
	  }
	}
	/* Fall through to general algorithm */
      }
    }

    /* General algorithm */
    debug15(printf("minus/minus general\n"));
    for (p = hit3->substrings_1toN; p != NULL; p = List_next(p)) {
      substring3 = (Substring_T) List_head(p);
      for (q = hit5->substrings_1toN; q != NULL; q = List_next(q)) {
	substring5 = (Substring_T) List_head(q);
	if ((common_genomicpos = Substring_overlap_segment_trimmed(substring5,substring3)) != 0) {
	  return common_genomicpos;
	}
      }
    }

    return 0;

  } else {
    abort();
    return 0;
  }
}


static bool
test_hardclips (Univcoord_T *common_genomicpos, int hardclip_low, Stage3end_T hit_low,
		int hardclip_high, Stage3end_T hit_high, Univcoord_T chroffset) {
  Substring_T low_substring, high_substring;
  struct Pair_T *low_pairarray, *high_pairarray;
  int low_npairs, high_npairs;
  int low_querypos, high_querypos;
  int low_querylength, high_querylength;
  bool plusp;

  low_querylength = hit_low->querylength;
  high_querylength = hit_high->querylength;

  debug15(printf("Entering test_hardclips with hardclip_low %d, hardclip_high %d\n",
		 hardclip_low,hardclip_high));
  debug15(printf("querylength_low %d, querylength_high %d\n",low_querylength,high_querylength));

  plusp = Stage3end_plusp(hit_low);

  if (Stage3end_hittype(hit_low) == GMAP && Stage3end_hittype(hit_high) == GMAP) {
    low_pairarray = Stage3end_pairarray(hit_low);
    low_npairs = Stage3end_npairs(hit_low);
    high_pairarray = Stage3end_pairarray(hit_high);
    high_npairs = Stage3end_npairs(hit_high);

    if (plusp == true) {
      low_querypos = hardclip_low;
      high_querypos = high_querylength /*- 1*/ - hardclip_high;
    } else {
      low_querypos = low_querylength /*- 1*/ - hardclip_low;
      high_querypos = hardclip_high;
    }
    debug15(printf("Dual GMAP.  low_querypos %d, high_querypos %d\n",low_querypos,high_querypos));

    if (Pairarray_contains_p(low_pairarray,low_npairs,low_querypos) == false) {
      debug15(printf("Fails because low_querypos %d is not in low_pairarray\n",low_querypos));
      return false;
    } else if (Pairarray_contains_p(low_pairarray,low_npairs,low_querypos-1) == false) {
      debug15(printf("Fails because low_querypos %d - 1 is not in low_pairarray\n",low_querypos));
      return false;
    } else if (Pairarray_contains_p(low_pairarray,low_npairs,low_querypos+1) == false) {
      debug15(printf("Fails because low_querypos %d + 1 is not in low_pairarray\n",low_querypos));
      return false;
    } else if (Pairarray_contains_p(high_pairarray,high_npairs,high_querypos) == false) {
      debug15(printf("Fails because high_querypos %d is not in high_pairarray\n",low_querypos));
      return false;
    } else if (Pairarray_contains_p(high_pairarray,high_npairs,high_querypos-1) == false) {
      debug15(printf("Fails because high_querypos %d - 1 is not in high_pairarray\n",low_querypos));
      return false;
    } else if (Pairarray_contains_p(high_pairarray,high_npairs,high_querypos+1) == false) {
      debug15(printf("Fails because high_querypos %d + 1 is not in high_pairarray\n",low_querypos));
      return false;
    } else if (Pairarray_lookup(low_pairarray,low_npairs,low_querypos) != Pairarray_lookup(high_pairarray,high_npairs,high_querypos)) {
      debug15(printf("Fails because low genomicpos %u != high genomicpos %u\n",
		     Pairarray_lookup(low_pairarray,low_npairs,low_querypos),
		     Pairarray_lookup(high_pairarray,high_npairs,high_querypos)));
      return false;
    } else {
      *common_genomicpos = Pairarray_lookup(low_pairarray,low_npairs,low_querypos) + chroffset;
      debug15(printf("Succeeds with common point %u\n",*common_genomicpos - chroffset));
      return true;
    }

  } else if (Stage3end_hittype(hit_low) == GMAP) {
    low_pairarray = Stage3end_pairarray(hit_low);
    low_npairs = Stage3end_npairs(hit_low);

    if (plusp == true) {
      low_querypos = hardclip_low;
      high_querypos = high_querylength /*- 1*/ - hardclip_high;
    } else {
      low_querypos = low_querylength /*- 1*/ - hardclip_low;
      high_querypos = hardclip_high;
    }
    debug15(printf("Low GMAP.  low_querypos %d, high_querypos %d\n",low_querypos,high_querypos));

    if (Pairarray_contains_p(low_pairarray,low_npairs,low_querypos) == false) {
      debug15(printf("Fails because low_querypos %d is not in low_pairarray\n",low_querypos));
      return false;
    } else if (Pairarray_contains_p(low_pairarray,low_npairs,low_querypos-1) == false) {
      debug15(printf("Fails because low_querypos %d - 1 is not in low_pairarray\n",low_querypos));
      return false;
    } else if (Pairarray_contains_p(low_pairarray,low_npairs,low_querypos+1) == false) {
      debug15(printf("Fails because low_querypos %d + 1 is not in low_pairarray\n",low_querypos));
      return false;
    } else if ((high_substring = Stage3end_substring_containing(hit_high,high_querypos)) == NULL) {
      debug15(printf("Fails because high_querypos %d gives a NULL substring\n",high_querypos));
      return false;
    } else if (Stage3end_substring_containing(hit_high,high_querypos-1) != high_substring) {
      debug15(printf("Fails because high_querypos %d - 1 gives substring %p\n",
		     high_querypos,Stage3end_substring_containing(hit_high,high_querypos-1)));
      return false;
    } else if (Stage3end_substring_containing(hit_high,high_querypos+1) != high_substring) {
      debug15(printf("Fails because high_querypos %d + 1 gives substring %p\n",
		     high_querypos,Stage3end_substring_containing(hit_high,high_querypos+1)));
      return false;
    } else if (plusp == true) {
      if (Pairarray_lookup(low_pairarray,low_npairs,low_querypos) != Substring_genomicstart(high_substring) + high_querypos - chroffset) {
	debug15(printf("Fails because low chrpos %u != high chrpos %u\n",
		       Pairarray_lookup(low_pairarray,low_npairs,low_querypos),
		       Substring_genomicstart(high_substring) + high_querypos - chroffset));
	return false;
      }
    } else {
      if (Pairarray_lookup(low_pairarray,low_npairs,low_querypos) != (Substring_genomicstart(high_substring) - 1) - high_querypos - chroffset) {
	debug15(printf("Fails because low chrpos %u != high chrpos %u\n",
		       Pairarray_lookup(low_pairarray,low_npairs,low_querypos),
		       (Substring_genomicstart(high_substring) - 1) - high_querypos - chroffset));
	return false;
      }
    }

    *common_genomicpos = Pairarray_lookup(low_pairarray,low_npairs,low_querypos) + chroffset;
    debug15(printf("Succeeds with common point %u\n",*common_genomicpos - chroffset));
    return true;

  } else if (Stage3end_hittype(hit_high) == GMAP) {
    high_pairarray = Stage3end_pairarray(hit_high);
    high_npairs = Stage3end_npairs(hit_high);

    if (plusp == true) {
      low_querypos = hardclip_low;
      high_querypos = high_querylength /*- 1*/ - hardclip_high;
    } else {
      low_querypos = low_querylength /*- 1*/ - hardclip_low;
      high_querypos = hardclip_high;
    }
    debug15(printf("High GMAP.  low_querypos %d, high_querypos %d\n",low_querypos,high_querypos));

    if ((low_substring = Stage3end_substring_containing(hit_low,low_querypos)) == NULL) {
      debug15(printf("Fails because low_querypos %d gives a NULL substring\n",low_querypos));
      return false;
    } else if (Stage3end_substring_containing(hit_low,low_querypos-1) != low_substring) {
      debug15(printf("Fails because low_querypos %d - 1 gives substring %p\n",
		     low_querypos,Stage3end_substring_containing(hit_low,low_querypos-1)));
      return false;
    } else if (Stage3end_substring_containing(hit_low,low_querypos+1) != low_substring) {
      debug15(printf("Fails because low_querypos %d + 1 gives substring %p\n",
		     low_querypos,Stage3end_substring_containing(hit_low,low_querypos+1)));
      return false;
    } else if (Pairarray_contains_p(high_pairarray,high_npairs,high_querypos) == false) {
      debug15(printf("Fails because high_querypos %d is not in high_pairarray\n",low_querypos));
      return false;
    } else if (Pairarray_contains_p(high_pairarray,high_npairs,high_querypos-1) == false) {
      debug15(printf("Fails because high_querypos %d - 1 is not in high_pairarray\n",low_querypos));
      return false;
    } else if (Pairarray_contains_p(high_pairarray,high_npairs,high_querypos+1) == false) {
      debug15(printf("Fails because high_querypos %d + 1 is not in high_pairarray\n",low_querypos));
      return false;
    } else if (plusp == true) {
      if (Pairarray_lookup(high_pairarray,high_npairs,high_querypos) != Substring_genomicstart(low_substring) + low_querypos - chroffset) {
	debug15(printf("Fails because low chrpos %u != high chrpos %u\n",
		       Substring_genomicstart(low_substring) + low_querypos - chroffset,
		       Pairarray_lookup(high_pairarray,high_npairs,high_querypos)));
	return false;
      }
    } else {
      if (Pairarray_lookup(high_pairarray,high_npairs,high_querypos) != (Substring_genomicstart(low_substring) - 1) - low_querypos - chroffset) {
	debug15(printf("Fails because low chrpos %u != high chrpos %u\n",
		       (Substring_genomicstart(low_substring) - 1) - low_querypos - chroffset,
		       Pairarray_lookup(high_pairarray,high_npairs,high_querypos)));
	return false;
      }
    }

    *common_genomicpos = Pairarray_lookup(high_pairarray,high_npairs,high_querypos) + chroffset;
    debug15(printf("Succeeds with common point %u\n",*common_genomicpos - chroffset));
    return true;

  } else {
    if (plusp == true) {
      low_querypos = hardclip_low;
      high_querypos = high_querylength /*- 1*/ - hardclip_high;
      debug15(printf("Both substrings, plus.  low_querypos %d, high_querypos %d\n",low_querypos,high_querypos));

      if ((low_substring = Stage3end_substring_containing(hit_low,low_querypos)) == NULL) {
	debug15(printf("Fails because low_querypos %d gives a NULL substring\n",low_querypos));
	return false;
      } else if (Stage3end_substring_containing(hit_low,low_querypos-1) != low_substring) {
	debug15(printf("Fails because low_querypos %d - 1 gives substring %p\n",
		       low_querypos,Stage3end_substring_containing(hit_low,low_querypos-1)));
	return false;
      } else if (Stage3end_substring_containing(hit_low,low_querypos+1) != low_substring) {
	debug15(printf("Fails because low_querypos %d + 1 gives substring %p\n",
		       low_querypos,Stage3end_substring_containing(hit_low,low_querypos+1)));
	return false;
      } else if ((high_substring = Stage3end_substring_containing(hit_high,high_querypos)) == NULL) {
	debug15(printf("Fails because high_querypos %d gives a NULL substring\n",high_querypos));
	return false;
      } else if (Stage3end_substring_containing(hit_high,high_querypos-1) != high_substring) {
	debug15(printf("Fails because high_querypos %d - 1 gives substring %p\n",
		       high_querypos,Stage3end_substring_containing(hit_high,high_querypos-1)));
	return false;
      } else if (Stage3end_substring_containing(hit_high,high_querypos+1) != high_substring) {
	debug15(printf("Fails because high_querypos %d + 1 gives substring %p\n",
		       high_querypos,Stage3end_substring_containing(hit_high,high_querypos+1)));
	return false;
      } else if (Substring_genomicstart(low_substring) + low_querypos - chroffset != Substring_genomicstart(high_substring) + high_querypos - chroffset) {
	debug15(printf("Fails because low chrpos %u != high chrpos %u\n",
		       Substring_genomicstart(low_substring) + low_querypos - chroffset,
		       Substring_genomicstart(high_substring) + high_querypos - chroffset));
	return false;
      } else {
	*common_genomicpos = Substring_genomicstart(low_substring) + low_querypos; /* Want univcoord */
	debug15(printf("Succeeds with common point %u\n",*common_genomicpos - chroffset));
	return true;
      }

    } else {
      low_querypos = low_querylength /*- 1*/ - hardclip_low;
      high_querypos = hardclip_high;
      debug15(printf("Both substrings, minus.  low_querypos %d, high_querypos %d\n",low_querypos,high_querypos));

      if ((low_substring = Stage3end_substring_containing(hit_low,low_querypos)) == NULL) {
	debug15(printf("Fails because low_querypos %d gives a NULL substring\n",low_querypos));
	return false;
      } else if (Stage3end_substring_containing(hit_low,low_querypos-1) != low_substring) {
	debug15(printf("Fails because low_querypos %d - 1 gives substring %p\n",
		       low_querypos,Stage3end_substring_containing(hit_low,low_querypos-1)));
	return false;
      } else if (Stage3end_substring_containing(hit_low,low_querypos+1) != low_substring) {
	debug15(printf("Fails because low_querypos %d + 1 gives substring %p\n",
		       low_querypos,Stage3end_substring_containing(hit_low,low_querypos+1)));
	return false;
      } else if ((high_substring = Stage3end_substring_containing(hit_high,high_querypos)) == NULL) {
	debug15(printf("Fails because high_querypos %d gives a NULL substring\n",high_querypos));
	return false;
      } else if (Stage3end_substring_containing(hit_high,high_querypos-1) != high_substring) {
	debug15(printf("Fails because high_querypos %d - 1 gives substring %p\n",
		       high_querypos,Stage3end_substring_containing(hit_high,high_querypos-1)));
	return false;
      } else if (Stage3end_substring_containing(hit_high,high_querypos+1) != high_substring) {
	debug15(printf("Fails because high_querypos %d + 1 gives substring %p\n",
		       high_querypos,Stage3end_substring_containing(hit_high,high_querypos+1)));
	return false;
      }  else if ((Substring_genomicstart(low_substring) - 1) - low_querypos - chroffset != (Substring_genomicstart(high_substring) - 1) - high_querypos - chroffset) {
	debug15(printf("Fails because low chrpos %u != high chrpos %u\n",
		       (Substring_genomicstart(low_substring) - 1) - low_querypos - chroffset,
		       (Substring_genomicstart(high_substring) - 1) - high_querypos - chroffset));
	return false;
      } else {
	*common_genomicpos = (Substring_genomicstart(low_substring) - 1) - low_querypos; /* Want univcoord */
	debug15(printf("Succeeds with common point %u\n",*common_genomicpos - chroffset));
	return true;
      }
    }
  }
}



/* Replaces adjust_hardclips in samprint.c */
static Univcoord_T
adjust_hardclips_right (int *shift, int hardclip_low, Stage3end_T hit_low,
			int hardclip_high, Stage3end_T hit_high, Univcoord_T chroffset) {
  Substring_T low_substring, high_substring;
  struct Pair_T *low_pairarray, *high_pairarray;
  int low_npairs, high_npairs;
  int low_querypos, high_querypos;
  int low_querylength, high_querylength;
  Chrpos_T low_chrpos, high_chrpos;
  bool plusp;


  low_querylength = hit_low->querylength;
  high_querylength = hit_high->querylength;

  debug15(printf("Entering adjust_hardclips_right with hardclip_low %d, hardclip_high %d\n",
		 hardclip_low,hardclip_high));
  *shift = 1;			/* Making an initial move before each while loop */
  plusp = Stage3end_plusp(hit_low);

  if (Stage3end_hittype(hit_low) == GMAP && Stage3end_hittype(hit_high) == GMAP) {
    low_pairarray = Stage3end_pairarray(hit_low);
    low_npairs = Stage3end_npairs(hit_low);
    high_pairarray = Stage3end_pairarray(hit_high);
    high_npairs = Stage3end_npairs(hit_high);

    if (plusp == true) {
      low_querypos = hardclip_low;
      high_querypos = high_querylength /*- 1*/ - hardclip_high;
      debug15(printf("Dual GMAP, plus.  low_querypos %d, high_querypos %d\n",low_querypos,high_querypos));

      low_querypos++;
      high_querypos++;
      debug15(printf("right shift %d: Advancing to low_querypos %d and high_querypos %d\n",*shift,low_querypos,high_querypos));
      while ((low_querypos + 1) < low_querylength && (high_querypos + 1) < high_querylength &&
	     (Pairarray_contains_p(low_pairarray,low_npairs,low_querypos) == false ||
	      Pairarray_contains_p(low_pairarray,low_npairs,low_querypos-1) == false ||
	      Pairarray_contains_p(low_pairarray,low_npairs,low_querypos+1) ==  false ||
	      Pairarray_contains_p(high_pairarray,high_npairs,high_querypos) == false ||
	      Pairarray_contains_p(high_pairarray,high_npairs,high_querypos-1) == false ||
	      Pairarray_contains_p(high_pairarray,high_npairs,high_querypos+1) == false ||
	      Pairarray_lookup(low_pairarray,low_npairs,low_querypos) != Pairarray_lookup(high_pairarray,high_npairs,high_querypos))) {
	(*shift) += 1;
	if (Pairarray_contains_p(low_pairarray,low_npairs,low_querypos) == false) {
	  low_querypos++;
	} else if (Pairarray_contains_p(high_pairarray,high_npairs,high_querypos) == false) {
	  high_querypos++;
	} else {
	  low_chrpos = Pairarray_lookup(low_pairarray,low_npairs,low_querypos);
	  high_chrpos = Pairarray_lookup(high_pairarray,high_npairs,high_querypos);
	  if (low_chrpos < high_chrpos) {
	    debug15(printf("low_chrpos %u < high_chrpos %u, so advancing low_querypos\n",low_chrpos,high_chrpos));
	    low_querypos++;
	  } else if (high_chrpos < low_chrpos) {
	    debug15(printf("high_chrpos %u < low_chrpos %u, so advancing high_querypos\n",high_chrpos,low_chrpos));
	    high_querypos++;
	  } else {
	    low_querypos++;
	    high_querypos++;
	  }
	}
	debug15(printf("right shift %d: Advancing to low_querypos %d and high_querypos %d\n",*shift,low_querypos,high_querypos));
      }

      if ((low_querypos + 1) >= low_querylength || (high_querypos + 1) >= high_querylength) {
	*shift = 0;
	return 0;
      } else {
	assert(Pairarray_contains_p(low_pairarray,low_npairs,low_querypos) == true);
	assert(Pairarray_contains_p(low_pairarray,low_npairs,low_querypos-1) == true);
	assert(Pairarray_contains_p(low_pairarray,low_npairs,low_querypos+1) == true);
	assert(Pairarray_contains_p(high_pairarray,high_npairs,high_querypos) == true);
	assert(Pairarray_contains_p(high_pairarray,high_npairs,high_querypos-1) == true);
	assert(Pairarray_contains_p(high_pairarray,high_npairs,high_querypos+1) == true);
	return Pairarray_lookup(low_pairarray,low_npairs,low_querypos) + chroffset;
      }

    } else {
      low_querypos = low_querylength /*- 1*/ - hardclip_low;
      high_querypos = hardclip_high;
      debug15(printf("Dual GMAP, minus.  low_querypos %d, high_querypos %d\n",low_querypos,high_querypos));

      low_querypos--;
      high_querypos--;
      debug15(printf("right shift %d: Advancing to low_querypos %d and high_querypos %d\n",*shift,low_querypos,high_querypos));
      while ((low_querypos - 1) >= 0 && (high_querypos - 1) >= 0 &&
	     (Pairarray_contains_p(low_pairarray,low_npairs,low_querypos) == false ||
	      Pairarray_contains_p(low_pairarray,low_npairs,low_querypos-1) == false ||
	      Pairarray_contains_p(low_pairarray,low_npairs,low_querypos+1) ==  false ||
	      Pairarray_contains_p(high_pairarray,high_npairs,high_querypos) == false ||
	      Pairarray_contains_p(high_pairarray,high_npairs,high_querypos-1) == false ||
	      Pairarray_contains_p(high_pairarray,high_npairs,high_querypos+1) == false ||
	      Pairarray_lookup(low_pairarray,low_npairs,low_querypos) != Pairarray_lookup(high_pairarray,high_npairs,high_querypos))) {
	(*shift) += 1;
	if (Pairarray_contains_p(low_pairarray,low_npairs,low_querypos) == false) {
	  low_querypos--;
	} else if (Pairarray_contains_p(high_pairarray,high_npairs,high_querypos) == false) {
	  high_querypos--;
	} else {
	  low_chrpos = Pairarray_lookup(low_pairarray,low_npairs,low_querypos);
	  high_chrpos = Pairarray_lookup(high_pairarray,high_npairs,high_querypos);
	  if (low_chrpos < high_chrpos) {
	    debug15(printf("low_chrpos %u < high_chrpos %u, so decreasing low_querypos\n",low_chrpos,high_chrpos));
	    low_querypos--;
	  } else if (high_chrpos < low_chrpos) {
	    debug15(printf("high_chrpos %u < low_chrpos %u, so decreasing high_querypos\n",high_chrpos,low_chrpos));
	    high_querypos--;
	  } else {
	    low_querypos--;
	    high_querypos--;
	  }
	}
	debug15(printf("right shift %d: Advancing to low_querypos %d and high_querypos %d\n",*shift,low_querypos,high_querypos));
      }

      if ((low_querypos - 1) < 0 || (high_querypos - 1) < 0) {
	*shift = 0;
	return 0;
      } else {
	assert(Pairarray_contains_p(low_pairarray,low_npairs,low_querypos) == true);
	assert(Pairarray_contains_p(low_pairarray,low_npairs,low_querypos-1) == true);
	assert(Pairarray_contains_p(low_pairarray,low_npairs,low_querypos+1) == true);
	assert(Pairarray_contains_p(high_pairarray,high_npairs,high_querypos) == true);
	assert(Pairarray_contains_p(high_pairarray,high_npairs,high_querypos-1) == true);
	assert(Pairarray_contains_p(high_pairarray,high_npairs,high_querypos+1) == true);
	return Pairarray_lookup(low_pairarray,low_npairs,low_querypos) + chroffset;
      }
    }

  } else if (Stage3end_hittype(hit_low) == GMAP) {
    low_pairarray = Stage3end_pairarray(hit_low);
    low_npairs = Stage3end_npairs(hit_low);

    if (plusp == true) {
      low_querypos = hardclip_low;
      high_querypos = high_querylength /*- 1*/ - hardclip_high;
      debug15(printf("Low GMAP, plus.  low_querypos %d, high_querypos %d\n",low_querypos,high_querypos));

      low_querypos++;
      high_querypos++;
      debug15(printf("right shift %d: Advancing to low_querypos %d and high_querypos %d\n",*shift,low_querypos,high_querypos));
      while ((low_querypos + 1) < low_querylength && (high_querypos + 1) < high_querylength &&
	     (Pairarray_contains_p(low_pairarray,low_npairs,low_querypos) == false ||
	      Pairarray_contains_p(low_pairarray,low_npairs,low_querypos-1) == false ||
	      Pairarray_contains_p(low_pairarray,low_npairs,low_querypos+1) == false ||
	      (high_substring = Stage3end_substring_containing(hit_high,high_querypos)) == NULL ||
	      Stage3end_substring_containing(hit_high,high_querypos-1) != high_substring ||
	      Stage3end_substring_containing(hit_high,high_querypos+1) != high_substring ||
	      Pairarray_lookup(low_pairarray,low_npairs,low_querypos) != Substring_genomicstart(high_substring) + high_querypos - chroffset)) {
	(*shift) += 1;
	if (Pairarray_contains_p(low_pairarray,low_npairs,low_querypos) == false) {
	  low_querypos++;
	} else if ((high_substring = Stage3end_substring_containing(hit_high,high_querypos)) == NULL) {
	  high_querypos++;
	} else {
	  low_chrpos = Pairarray_lookup(low_pairarray,low_npairs,low_querypos);
	  high_chrpos = Substring_genomicstart(high_substring) + high_querypos - chroffset;
	  if (low_chrpos < high_chrpos) {
	    debug15(printf("low_chrpos %u < high_chrpos %u, so advancing low_querypos\n",low_chrpos,high_chrpos));
	    low_querypos++;
	  } else if (high_chrpos < low_chrpos) {
	    debug15(printf("high_chrpos %u < low_chrpos %u, so advancing high_querypos\n",high_chrpos,low_chrpos));
	    high_querypos++;
	  } else {
	    low_querypos++;
	    high_querypos++;
	  }
	}
	debug15(printf("right shift %d: Advancing to low_querypos %d and high_querypos %d\n",*shift,low_querypos,high_querypos));
      }

      if ((low_querypos + 1) >= low_querylength || (high_querypos + 1) >= high_querylength) {
	debug15(printf("Failing because low_querypos %d + 1 >= low_querylength %d\n",low_querypos,low_querylength));
	*shift = 0;
	return 0;
      } else if (Stage3end_substring_containing(hit_high,high_querypos) == NULL) {
	debug15(printf("Failing because no substring contains high_querypos %d\n",high_querypos));
	*shift = 0;
	return 0;
      } else {
	assert(Pairarray_contains_p(low_pairarray,low_npairs,low_querypos) == true);
	assert(Pairarray_contains_p(low_pairarray,low_npairs,low_querypos-1) == true);
	assert(Pairarray_contains_p(low_pairarray,low_npairs,low_querypos+1) == true);
	assert((high_substring = Stage3end_substring_containing(hit_high,high_querypos)) != NULL);
	assert(Stage3end_substring_containing(hit_high,high_querypos-1) == high_substring);
	assert(Stage3end_substring_containing(hit_high,high_querypos+1) == high_substring);
	return Pairarray_lookup(low_pairarray,low_npairs,low_querypos) + chroffset;
      }

    } else {
      low_querypos = low_querylength /*- 1*/ - hardclip_low;
      high_querypos = hardclip_high;
      debug15(printf("Low GMAP, minus.  low_querypos %d, high_querypos %d\n",low_querypos,high_querypos));

      low_querypos--;
      high_querypos--;
      debug15(printf("right shift %d: Advancing to low_querypos %d and high_querypos %d\n",*shift,low_querypos,high_querypos));
      while ((low_querypos - 1) >= 0 && (high_querypos - 1) >= 0 &&
	     (Pairarray_contains_p(low_pairarray,low_npairs,low_querypos) == false ||
	      Pairarray_contains_p(low_pairarray,low_npairs,low_querypos-1) == false ||
	      Pairarray_contains_p(low_pairarray,low_npairs,low_querypos+1) == false ||
	      (high_substring = Stage3end_substring_containing(hit_high,high_querypos)) == NULL ||
	      Stage3end_substring_containing(hit_high,high_querypos-1) != high_substring ||
	      Stage3end_substring_containing(hit_high,high_querypos+1) != high_substring ||
	      Pairarray_lookup(low_pairarray,low_npairs,low_querypos) != (Substring_genomicstart(high_substring) - 1) - high_querypos - chroffset)) {
	(*shift) += 1;
	if (Pairarray_contains_p(low_pairarray,low_npairs,low_querypos) == false) {
	  low_querypos--;
	} else if ((high_substring = Stage3end_substring_containing(hit_high,high_querypos)) == NULL) {
	  high_querypos--;
	} else {
	  low_chrpos = Pairarray_lookup(low_pairarray,low_npairs,low_querypos);
	  high_chrpos = (Substring_genomicstart(high_substring) - 1) - high_querypos - chroffset;
	  if (low_chrpos < high_chrpos) {
	    debug15(printf("low_chrpos %u < high_chrpos %u, so decreasing low_querypos\n",low_chrpos,high_chrpos));
	    low_querypos--;
	  } else if (high_chrpos < low_chrpos) {
	    debug15(printf("high_chrpos %u < low_chrpos %u, so decreasing high_querypos\n",high_chrpos,low_chrpos));
	    high_querypos--;
	  } else {
	    low_querypos--;
	    high_querypos--;
	  }
	}
	debug15(printf("right shift %d: Advancing to low_querypos %d and high_querypos %d\n",*shift,low_querypos,high_querypos));
      }

      if ((low_querypos - 1) < 0 || (high_querypos - 1) < 0) {
	debug15(printf("Failing because low_querypos %d - 1 < 0\n",low_querypos));
	*shift = 0;
	return 0;
      } else if (Stage3end_substring_containing(hit_high,high_querypos) == NULL) {
	debug15(printf("Failing because no substring contains high_querypos %d\n",high_querypos));
	*shift = 0;
	return 0;
      } else {
	assert(Pairarray_contains_p(low_pairarray,low_npairs,low_querypos) == true);
	assert(Pairarray_contains_p(low_pairarray,low_npairs,low_querypos-1) == true);
	assert(Pairarray_contains_p(low_pairarray,low_npairs,low_querypos+1) == true);
	assert((high_substring = Stage3end_substring_containing(hit_high,high_querypos)) != NULL);
	assert(Stage3end_substring_containing(hit_high,high_querypos-1) == high_substring);
	assert(Stage3end_substring_containing(hit_high,high_querypos+1) == high_substring);
	return Pairarray_lookup(low_pairarray,low_npairs,low_querypos) + chroffset;
      }
    }

  } else if (Stage3end_hittype(hit_high) == GMAP) {
    high_pairarray = Stage3end_pairarray(hit_high);
    high_npairs = Stage3end_npairs(hit_high);

    if (plusp == true) {
      low_querypos = hardclip_low;
      high_querypos = high_querylength /*- 1*/ - hardclip_high;
      debug15(printf("High GMAP.  low_querypos %d, high_querypos %d\n",low_querypos,high_querypos));

      low_querypos++;
      high_querypos++;
      debug15(printf("right shift %d: Advancing to low_querypos %d and high_querypos %d\n",*shift,low_querypos,high_querypos));
      while ((high_querypos + 1) < high_querylength && (low_querypos + 1) < low_querylength &&
	     (Pairarray_contains_p(high_pairarray,high_npairs,high_querypos) == false ||
	      Pairarray_contains_p(high_pairarray,high_npairs,high_querypos-1) == false ||
	      Pairarray_contains_p(high_pairarray,high_npairs,high_querypos+1) == false ||
	      (low_substring = Stage3end_substring_containing(hit_low,low_querypos)) == NULL ||
	      Stage3end_substring_containing(hit_low,low_querypos-1) != low_substring ||
	      Stage3end_substring_containing(hit_low,low_querypos+1) != low_substring ||
	      Pairarray_lookup(high_pairarray,high_npairs,high_querypos) != Substring_genomicstart(low_substring) + low_querypos - chroffset)) {
	(*shift) += 1;
	if ((low_substring = Stage3end_substring_containing(hit_low,low_querypos)) == NULL) {
	  low_querypos++;
	} else if (Pairarray_contains_p(high_pairarray,high_npairs,high_querypos) == false) {
	  high_querypos++;
	} else {
	  low_chrpos = Substring_genomicstart(low_substring) + low_querypos - chroffset;
	  high_chrpos = Pairarray_lookup(high_pairarray,high_npairs,high_querypos);
	  if (low_chrpos < high_chrpos) {
	    debug15(printf("low_chrpos %u < high_chrpos %u, so advancing low_querypos\n",low_chrpos,high_chrpos));
	    low_querypos++;
	  } else if (high_chrpos < low_chrpos) {
	    debug15(printf("high_chrpos %u < low_chrpos %u, so advancing high_querypos\n",high_chrpos,low_chrpos));
	    high_querypos++;
	  } else {
	    low_querypos++;
	    high_querypos++;
	  }
	}
	debug15(printf("right shift %d: Advancing to low_querypos %d and high_querypos %d\n",*shift,low_querypos,high_querypos));
      }

      if ((high_querypos + 1) >= high_querylength || (low_querypos + 1) >= low_querylength ||
	  Stage3end_substring_containing(hit_low,low_querypos) == NULL) {
	*shift = 0;
	return 0;
      } else {
	assert((low_substring = Stage3end_substring_containing(hit_low,low_querypos)) != NULL);
	assert(Stage3end_substring_containing(hit_low,low_querypos-1) == low_substring);
	assert(Stage3end_substring_containing(hit_low,low_querypos+1) == low_substring);
	assert(Pairarray_contains_p(high_pairarray,high_npairs,high_querypos) == true);
	assert(Pairarray_contains_p(high_pairarray,high_npairs,high_querypos-1) == true);
	assert(Pairarray_contains_p(high_pairarray,high_npairs,high_querypos+1) == true);
	return Pairarray_lookup(high_pairarray,high_npairs,high_querypos) + chroffset;
      }

    } else {
      low_querypos = low_querylength /*- 1*/ - hardclip_low;
      high_querypos = hardclip_high;
      debug15(printf("High GMAP, plus.  low_querypos %d, high_querypos %d\n",low_querypos,high_querypos));

      low_querypos--;
      high_querypos--;
      debug15(printf("right shift %d: Advancing to low_querypos %d and high_querypos %d\n",*shift,low_querypos,high_querypos));
      while ((high_querypos - 1) >= 0 && (low_querypos - 1) >= 0 &&
	     (Pairarray_contains_p(high_pairarray,high_npairs,high_querypos) == false ||
	      Pairarray_contains_p(high_pairarray,high_npairs,high_querypos-1) == false ||
	      Pairarray_contains_p(high_pairarray,high_npairs,high_querypos+1) == false ||
	      (low_substring = Stage3end_substring_containing(hit_low,low_querypos)) == NULL ||
	      Stage3end_substring_containing(hit_low,low_querypos-1) != low_substring ||
	      Stage3end_substring_containing(hit_low,low_querypos+1) != low_substring ||
	      Pairarray_lookup(high_pairarray,high_npairs,high_querypos) != (Substring_genomicstart(low_substring) - 1) - low_querypos - chroffset)) {
	(*shift) += 1;
	if ((low_substring = Stage3end_substring_containing(hit_low,low_querypos)) == NULL) {
	  low_querypos--;
	} else if (Pairarray_contains_p(high_pairarray,high_npairs,high_querypos) == false) {
	  high_querypos--;
	} else {
	  low_chrpos = (Substring_genomicstart(low_substring) - 1) - low_querypos - chroffset;
	  high_chrpos = Pairarray_lookup(high_pairarray,high_npairs,high_querypos);
	  if (low_chrpos < high_chrpos) {
	    debug15(printf("low_chrpos %u < high_chrpos %u, so decreasing low_querypos\n",low_chrpos,high_chrpos));
	    low_querypos--;
	  } else if (high_chrpos < low_chrpos) {
	    debug15(printf("high_chrpos %u < low_chrpos %u, so decreasing high_querypos\n",high_chrpos,low_chrpos));
	    high_querypos--;
	  } else {
	    low_querypos--;
	    high_querypos--;
	  }
	}
	debug15(printf("right shift %d: Advancing to low_querypos %d and high_querypos %d\n",*shift,low_querypos,high_querypos));
      }

      if ((high_querypos - 1) < 0 || (low_querypos - 1) < 0 || 
	  Stage3end_substring_containing(hit_low,low_querypos) == NULL) {
	*shift = 0;
	return 0;
      } else {
	assert((low_substring = Stage3end_substring_containing(hit_low,low_querypos)) != NULL);
	assert(Stage3end_substring_containing(hit_low,low_querypos-1) == low_substring);
	assert(Stage3end_substring_containing(hit_low,low_querypos+1) == low_substring);
	assert(Pairarray_contains_p(high_pairarray,high_npairs,high_querypos) == true);
	assert(Pairarray_contains_p(high_pairarray,high_npairs,high_querypos-1) == true);
	assert(Pairarray_contains_p(high_pairarray,high_npairs,high_querypos+1) == true);
	return Pairarray_lookup(high_pairarray,high_npairs,high_querypos) + chroffset;
      }
    }

  } else {
    if (plusp == true) {
      low_querypos = hardclip_low;
      high_querypos = high_querylength /*- 1*/ - hardclip_high;
      debug15(printf("Both substrings, plus.  low_querypos %d, high_querypos %d\n",low_querypos,high_querypos));

      low_querypos++;
      high_querypos++;
      debug15(printf("right shift %d: Advancing to low_querypos %d and high_querypos %d\n",*shift,low_querypos,high_querypos));
      while ((low_querypos + 1) < low_querylength && (high_querypos + 1) < high_querylength &&
	     ((low_substring = Stage3end_substring_containing(hit_low,low_querypos)) == NULL ||
	      Stage3end_substring_containing(hit_low,low_querypos-1) != low_substring ||
	      Stage3end_substring_containing(hit_low,low_querypos+1) != low_substring ||
	      (high_substring = Stage3end_substring_containing(hit_high,high_querypos)) == NULL ||
	      Stage3end_substring_containing(hit_high,high_querypos-1) != high_substring ||
	      Stage3end_substring_containing(hit_high,high_querypos+1) != high_substring ||
	      Substring_genomicstart(low_substring) + low_querypos - chroffset != Substring_genomicstart(high_substring) + high_querypos - chroffset)) {
	(*shift) += 1;
	if ((low_substring = Stage3end_substring_containing(hit_low,low_querypos)) == NULL) {
	  low_querypos++;
	} else if ((high_substring = Stage3end_substring_containing(hit_high,high_querypos)) == NULL) {
	  high_querypos++;
	} else {
	  low_chrpos = Substring_genomicstart(low_substring) + low_querypos - chroffset;
	  high_chrpos = Substring_genomicstart(high_substring) + high_querypos - chroffset;
	  if (low_chrpos < high_chrpos) {
	    debug15(printf("low_chrpos %u < high_chrpos %u, so advancing low_querypos\n",low_chrpos,high_chrpos));
	    low_querypos++;
	  } else if (high_chrpos < low_chrpos) {
	    debug15(printf("high_chrpos %u < low_chrpos %u, so advancing high_querypos\n",high_chrpos,low_chrpos));
	    high_querypos++;
	  } else {
	    low_querypos++;
	    high_querypos++;
	  }
	}
	debug15(printf("right shift %d: Advancing to low_querypos %d and high_querypos %d\n",*shift,low_querypos,high_querypos));
      }

      if ((low_querypos + 1) >= low_querylength ||
	  (high_querypos + 1) >= high_querylength ||
	  (low_substring = Stage3end_substring_containing(hit_low,low_querypos)) == NULL ||
	  Stage3end_substring_containing(hit_high,high_querypos) == NULL) {
	*shift = 0;
	return 0;
      } else {
	debug15(printf("Returning %u + %d\n",Substring_genomicstart(low_substring) - chroffset,
		       low_querypos));
	assert((low_substring = Stage3end_substring_containing(hit_low,low_querypos)) != NULL);
	assert((high_substring = Stage3end_substring_containing(hit_high,high_querypos)) != NULL);
	assert(Stage3end_substring_containing(hit_low,low_querypos-1) == low_substring);
	assert(Stage3end_substring_containing(hit_low,low_querypos+1) == low_substring);
	assert(Stage3end_substring_containing(hit_high,high_querypos-1) == high_substring);
	assert(Stage3end_substring_containing(hit_high,high_querypos+1) == high_substring);
	return Substring_genomicstart(low_substring) + low_querypos; /* Want univcoord */
      }

    } else {
      low_querypos = low_querylength /*- 1*/ - hardclip_low;
      high_querypos = hardclip_high;
      debug15(printf("Both substrings, minus.  low_querypos %d, high_querypos %d\n",low_querypos,high_querypos));

      low_querypos--;
      high_querypos--;
      debug15(printf("right shift %d: Advancing to low_querypos %d and high_querypos %d\n",*shift,low_querypos,high_querypos));
      while ((low_querypos - 1) >= 0 && (high_querypos - 1) >= 0 &&
	     ((low_substring = Stage3end_substring_containing(hit_low,low_querypos)) == NULL ||
	      Stage3end_substring_containing(hit_low,low_querypos-1) != low_substring ||
	      Stage3end_substring_containing(hit_low,low_querypos+1) != low_substring ||
	      (high_substring = Stage3end_substring_containing(hit_high,high_querypos)) == NULL ||
	      Stage3end_substring_containing(hit_high,high_querypos-1) != high_substring ||
	      Stage3end_substring_containing(hit_high,high_querypos+1) != high_substring ||
	      (Substring_genomicstart(low_substring) - 1) - low_querypos - chroffset != (Substring_genomicstart(high_substring) - 1) - high_querypos - chroffset)) {
	(*shift) += 1;
	if ((low_substring = Stage3end_substring_containing(hit_low,low_querypos)) == NULL) {
	  low_querypos--;
	} else if ((high_substring = Stage3end_substring_containing(hit_high,high_querypos)) == NULL) {
	  high_querypos--;
	} else {
	  low_chrpos = (Substring_genomicstart(low_substring) - 1) - low_querypos - chroffset;
	  high_chrpos = (Substring_genomicstart(high_substring) - 1) - high_querypos - chroffset;
	  if (low_chrpos < high_chrpos) {
	    debug15(printf("low_chrpos %u < high_chrpos %u, so decreasing low_querypos\n",low_chrpos,high_chrpos));
	    low_querypos--;
	  } else if (high_chrpos < low_chrpos) {
	    debug15(printf("high_chrpos %u < low_chrpos %u, so decreasing high_querypos\n",high_chrpos,low_chrpos));
	    high_querypos--;
	  } else {
	    low_querypos--;
	    high_querypos--;
	  }
	}
	debug15(printf("right shift %d: Advancing to low_querypos %d and high_querypos %d\n",*shift,low_querypos,high_querypos));
      }

      if ((low_querypos - 1) < 0 ||
	  (high_querypos - 1) < 0 ||
	  (low_substring = Stage3end_substring_containing(hit_low,low_querypos)) == NULL ||
	  Stage3end_substring_containing(hit_high,high_querypos) == NULL) {
	*shift = 0;
	return 0;
      } else {
	debug15(printf("Returning %u - %d\n",Substring_genomicstart(low_substring) - chroffset,
		       low_querypos));
	assert((low_substring = Stage3end_substring_containing(hit_low,low_querypos)) != NULL);
	assert((high_substring = Stage3end_substring_containing(hit_high,high_querypos)) != NULL);
	assert(Stage3end_substring_containing(hit_low,low_querypos-1) == low_substring);
	assert(Stage3end_substring_containing(hit_low,low_querypos+1) == low_substring);
	assert(Stage3end_substring_containing(hit_high,high_querypos-1) == high_substring);
	assert(Stage3end_substring_containing(hit_high,high_querypos+1) == high_substring);
	return (Substring_genomicstart(low_substring) - 1) - low_querypos; /* Want univcoord */
      }
    }
  }
}


/* Replaces adjust_hardclips in samprint.c */
static Univcoord_T
adjust_hardclips_left (int *shift, int hardclip_low, Stage3end_T hit_low,
		       int hardclip_high, Stage3end_T hit_high, Univcoord_T chroffset) {
  Substring_T low_substring, high_substring;
  struct Pair_T *low_pairarray, *high_pairarray;
  int low_npairs, high_npairs;
  int low_querypos, high_querypos;
  int low_querylength, high_querylength;
  Chrpos_T low_chrpos, high_chrpos;
  bool plusp;


  low_querylength = hit_low->querylength;
  high_querylength = hit_high->querylength;

  debug15(printf("Entering adjust_hardclips_left with hardclip_low %d, hardclip_high %d\n",
		 hardclip_low,hardclip_high));
  *shift = 1;			/* Making an initial move before each while loop */
  plusp = Stage3end_plusp(hit_low);

  if (Stage3end_hittype(hit_low) == GMAP && Stage3end_hittype(hit_high) == GMAP) {
    low_pairarray = Stage3end_pairarray(hit_low);
    low_npairs = Stage3end_npairs(hit_low);
    high_pairarray = Stage3end_pairarray(hit_high);
    high_npairs = Stage3end_npairs(hit_high);

    if (plusp == true) {
      low_querypos = hardclip_low;
      high_querypos = high_querylength /*- 1*/ - hardclip_high;
      debug15(printf("Dual GMAP, plus.  low_querypos %d, high_querypos %d\n",low_querypos,high_querypos));

      low_querypos--;
      high_querypos--;
      debug15(printf("left shift %d: Advancing to low_querypos %d and high_querypos %d\n",*shift,low_querypos,high_querypos));
      while ((low_querypos - 1) >= 0 && (high_querypos - 1) >= 0 &&
	     (Pairarray_contains_p(low_pairarray,low_npairs,low_querypos) == false ||
	      Pairarray_contains_p(low_pairarray,low_npairs,low_querypos-1) == false ||
	      Pairarray_contains_p(low_pairarray,low_npairs,low_querypos+1) ==  false ||
	      Pairarray_contains_p(high_pairarray,high_npairs,high_querypos) == false ||
	      Pairarray_contains_p(high_pairarray,high_npairs,high_querypos-1) == false ||
	      Pairarray_contains_p(high_pairarray,high_npairs,high_querypos+1) == false ||
	      Pairarray_lookup(low_pairarray,low_npairs,low_querypos) != Pairarray_lookup(high_pairarray,high_npairs,high_querypos))) {
	(*shift) += 1;
	if (Pairarray_contains_p(low_pairarray,low_npairs,low_querypos) == false) {
	  low_querypos--;
	} else if (Pairarray_contains_p(high_pairarray,high_npairs,high_querypos) == false) {
	  high_querypos--;
	} else {
	  low_chrpos = Pairarray_lookup(low_pairarray,low_npairs,low_querypos);
	  high_chrpos = Pairarray_lookup(high_pairarray,high_npairs,high_querypos);
	  if (low_chrpos > high_chrpos) {
	    debug15(printf("low_chrpos %u > high_chrpos %u, so decreasing low_querypos\n",low_chrpos,high_chrpos));
	    low_querypos--;
	  } else if (high_chrpos > low_chrpos) {
	    debug15(printf("high_chrpos %u > low_chrpos %u, so decreasing high_querypos\n",high_chrpos,low_chrpos));
	    high_querypos--;
	  } else {
	    low_querypos--;
	    high_querypos--;
	  }
	}
	debug15(printf("left shift %d: Advancing to low_querypos %d and high_querypos %d\n",*shift,low_querypos,high_querypos));
      }

      if ((low_querypos - 1) < 0 || (high_querypos - 1) < 0) {
	*shift = 0;
	return 0;
      } else {
	assert(Pairarray_contains_p(low_pairarray,low_npairs,low_querypos) == true);
	assert(Pairarray_contains_p(low_pairarray,low_npairs,low_querypos-1) == true);
	assert(Pairarray_contains_p(low_pairarray,low_npairs,low_querypos+1) == true);
	assert(Pairarray_contains_p(high_pairarray,high_npairs,high_querypos) == true);
	assert(Pairarray_contains_p(high_pairarray,high_npairs,high_querypos-1) == true);
	assert(Pairarray_contains_p(high_pairarray,high_npairs,high_querypos+1) == true);
	return Pairarray_lookup(low_pairarray,low_npairs,low_querypos) + chroffset;
      }

    } else {
      low_querypos = low_querylength /*- 1*/ - hardclip_low;
      high_querypos = hardclip_high;
      debug15(printf("Dual GMAP, minus.  low_querypos %d, high_querypos %d\n",low_querypos,high_querypos));

      low_querypos++;
      high_querypos++;
      debug15(printf("left shift %d: Advancing to low_querypos %d and high_querypos %d\n",*shift,low_querypos,high_querypos));
      while ((low_querypos + 1) < low_querylength && (high_querypos + 1) < high_querylength &&
	     (Pairarray_contains_p(low_pairarray,low_npairs,low_querypos) == false ||
	      Pairarray_contains_p(low_pairarray,low_npairs,low_querypos-1) == false ||
	      Pairarray_contains_p(low_pairarray,low_npairs,low_querypos+1) ==  false ||
	      Pairarray_contains_p(high_pairarray,high_npairs,high_querypos) == false ||
	      Pairarray_contains_p(high_pairarray,high_npairs,high_querypos-1) == false ||
	      Pairarray_contains_p(high_pairarray,high_npairs,high_querypos+1) == false ||
	      Pairarray_lookup(low_pairarray,low_npairs,low_querypos) != Pairarray_lookup(high_pairarray,high_npairs,high_querypos))) {
	(*shift) += 1;
	if (Pairarray_contains_p(low_pairarray,low_npairs,low_querypos) == false) {
	  low_querypos++;
	} else if (Pairarray_contains_p(high_pairarray,high_npairs,high_querypos) == false) {
	  high_querypos++;
	} else {
	  low_chrpos = Pairarray_lookup(low_pairarray,low_npairs,low_querypos);
	  high_chrpos = Pairarray_lookup(high_pairarray,high_npairs,high_querypos);
	  if (low_chrpos > high_chrpos) {
	    debug15(printf("low_chrpos %u > high_chrpos %u, so advancing low_querypos\n",low_chrpos,high_chrpos));
	    low_querypos++;
	  } else if (high_chrpos > low_chrpos) {
	    debug15(printf("high_chrpos %u > low_chrpos %u, so advancing high_querypos\n",high_chrpos,low_chrpos));
	    high_querypos++;
	  } else {
	    low_querypos++;
	    high_querypos++;
	  }
	}
	debug15(printf("left shift %d: Advancing to low_querypos %d and high_querypos %d\n",*shift,low_querypos,high_querypos));
      }

      if ((low_querypos + 1) >= low_querylength || (high_querypos + 1) >= high_querylength) {
	*shift = 0;
	return 0;
      } else {
	assert(Pairarray_contains_p(low_pairarray,low_npairs,low_querypos) == true);
	assert(Pairarray_contains_p(low_pairarray,low_npairs,low_querypos-1) == true);
	assert(Pairarray_contains_p(low_pairarray,low_npairs,low_querypos+1) == true);
	assert(Pairarray_contains_p(high_pairarray,high_npairs,high_querypos) == true);
	assert(Pairarray_contains_p(high_pairarray,high_npairs,high_querypos-1) == true);
	assert(Pairarray_contains_p(high_pairarray,high_npairs,high_querypos+1) == true);
	return Pairarray_lookup(low_pairarray,low_npairs,low_querypos) + chroffset;
      }
    }

  } else if (Stage3end_hittype(hit_low) == GMAP) {
    low_pairarray = Stage3end_pairarray(hit_low);
    low_npairs = Stage3end_npairs(hit_low);

    if (plusp == true) {
      low_querypos = hardclip_low;
      high_querypos = high_querylength /*- 1*/ - hardclip_high;
      debug15(printf("Low GMAP, plus.  low_querypos %d, high_querypos %d\n",low_querypos,high_querypos));

      low_querypos--;
      high_querypos--;
      debug15(printf("left shift %d: Advancing to low_querypos %d and high_querypos %d\n",*shift,low_querypos,high_querypos));
      while ((low_querypos - 1) >= 0 && (high_querypos - 1) >= 0 &&
	     (Pairarray_contains_p(low_pairarray,low_npairs,low_querypos) == false ||
	      Pairarray_contains_p(low_pairarray,low_npairs,low_querypos-1) == false ||
	      Pairarray_contains_p(low_pairarray,low_npairs,low_querypos+1) ==  false ||
	      (high_substring = Stage3end_substring_containing(hit_high,high_querypos)) == NULL ||
	      Stage3end_substring_containing(hit_high,high_querypos-1) != high_substring ||
	      Stage3end_substring_containing(hit_high,high_querypos+1) != high_substring ||
	      Pairarray_lookup(low_pairarray,low_npairs,low_querypos) != Substring_genomicstart(high_substring) + high_querypos - chroffset)) {
	(*shift) += 1;
	if (Pairarray_contains_p(low_pairarray,low_npairs,low_querypos) == false) {
	  low_querypos--;
	} else if ((high_substring = Stage3end_substring_containing(hit_high,high_querypos)) == NULL) {
	  high_querypos--;
	} else {
	  low_chrpos = Pairarray_lookup(low_pairarray,low_npairs,low_querypos);
	  high_chrpos = Substring_genomicstart(high_substring) + high_querypos - chroffset;
	  if (low_chrpos > high_chrpos) {
	    debug15(printf("low_chrpos %u > high_chrpos %u, so decreasing low_querypos\n",low_chrpos,high_chrpos));
	    low_querypos--;
	  } else if (high_chrpos > low_chrpos) {
	    debug15(printf("high_chrpos %u > low_chrpos %u, so decreasing high_querypos\n",high_chrpos,low_chrpos));
	    high_querypos--;
	  } else {
	    low_querypos--;
	    high_querypos--;
	  }
	}
	debug15(printf("left shift %d: Advancing to low_querypos %d and high_querypos %d\n",*shift,low_querypos,high_querypos));
      }

      if ((low_querypos - 1) < 0 || (high_querypos - 1) < 0 ||
	  Stage3end_substring_containing(hit_high,high_querypos) == NULL) {
	*shift = 0;
	return 0;
      } else {
	assert(Pairarray_contains_p(low_pairarray,low_npairs,low_querypos) == true);
	assert(Pairarray_contains_p(low_pairarray,low_npairs,low_querypos-1) == true);
	assert(Pairarray_contains_p(low_pairarray,low_npairs,low_querypos+1) == true);
	assert((high_substring = Stage3end_substring_containing(hit_high,high_querypos)) != NULL);
	assert(Stage3end_substring_containing(hit_high,high_querypos-1) == high_substring);
	assert(Stage3end_substring_containing(hit_high,high_querypos+1) == high_substring);
	return Pairarray_lookup(low_pairarray,low_npairs,low_querypos) + chroffset;
      }

    } else {
      low_querypos = low_querylength /*- 1*/ - hardclip_low;
      high_querypos = hardclip_high;
      debug15(printf("Low GMAP, minus.  low_querypos %d, high_querypos %d\n",low_querypos,high_querypos));

      low_querypos++;
      high_querypos++;
      debug15(printf("left shift %d: Advancing to low_querypos %d and high_querypos %d\n",*shift,low_querypos,high_querypos));
      while ((low_querypos + 1) < low_querylength && (high_querypos + 1) < high_querylength &&
	     (Pairarray_contains_p(low_pairarray,low_npairs,low_querypos) == false ||
	      Pairarray_contains_p(low_pairarray,low_npairs,low_querypos-1) == false ||
	      Pairarray_contains_p(low_pairarray,low_npairs,low_querypos+1) ==  false ||
	      (high_substring = Stage3end_substring_containing(hit_high,high_querypos)) == NULL ||
	      Stage3end_substring_containing(hit_high,high_querypos-1) != high_substring ||
	      Stage3end_substring_containing(hit_high,high_querypos+1) != high_substring ||
	      Pairarray_lookup(low_pairarray,low_npairs,low_querypos) != (Substring_genomicstart(high_substring) - 1) - high_querypos - chroffset)) {
	(*shift) += 1;
	if (Pairarray_contains_p(low_pairarray,low_npairs,low_querypos) == false) {
	  low_querypos++;
	} else if ((high_substring = Stage3end_substring_containing(hit_high,high_querypos)) == NULL) {
	  high_querypos++;
	} else {
	  low_chrpos = Pairarray_lookup(low_pairarray,low_npairs,low_querypos);
	  high_chrpos = (Substring_genomicstart(high_substring) - 1) - high_querypos - chroffset;
	  if (low_chrpos > high_chrpos) {
	    debug15(printf("low_chrpos %u > high_chrpos %u, so advancing low_querypos\n",low_chrpos,high_chrpos));
	    low_querypos++;
	  } else if (high_chrpos > low_chrpos) {
	    debug15(printf("high_chrpos %u > low_chrpos %u, so advancing high_querypos\n",high_chrpos,low_chrpos));
	    high_querypos++;
	  } else {
	    low_querypos++;
	    high_querypos++;
	  }
	}
	debug15(printf("left shift %d: Advancing to low_querypos %d and high_querypos %d\n",*shift,low_querypos,high_querypos));
      }

      if ((low_querypos + 1) >= low_querylength || (high_querypos + 1) >= high_querylength ||
	  Stage3end_substring_containing(hit_high,high_querypos) == NULL) {
	*shift = 0;
	return 0;
      } else {
	assert(Pairarray_contains_p(low_pairarray,low_npairs,low_querypos) == true);
	assert(Pairarray_contains_p(low_pairarray,low_npairs,low_querypos-1) == true);
	assert(Pairarray_contains_p(low_pairarray,low_npairs,low_querypos+1) == true);
	assert((high_substring = Stage3end_substring_containing(hit_high,high_querypos)) != NULL);
	assert(Stage3end_substring_containing(hit_high,high_querypos-1) == high_substring);
	assert(Stage3end_substring_containing(hit_high,high_querypos+1) == high_substring);
	return Pairarray_lookup(low_pairarray,low_npairs,low_querypos) + chroffset;
      }
    }

  } else if (Stage3end_hittype(hit_high) == GMAP) {
    high_pairarray = Stage3end_pairarray(hit_high);
    high_npairs = Stage3end_npairs(hit_high);

    if (plusp == true) {
      low_querypos = hardclip_low;
      high_querypos = high_querylength /*- 1*/ - hardclip_high;
      debug15(printf("High GMAP, plus.  low_querypos %d, high_querypos %d\n",low_querypos,high_querypos));

      low_querypos--;
      high_querypos--;
      debug15(printf("left shift %d: Advancing to low_querypos %d and high_querypos %d\n",*shift,low_querypos,high_querypos));
      while ((high_querypos - 1) >= 0 && (low_querypos - 1) >= 0 &&
	     (Pairarray_contains_p(high_pairarray,high_npairs,high_querypos) == false ||
	      Pairarray_contains_p(high_pairarray,high_npairs,high_querypos-1) == false ||
	      Pairarray_contains_p(high_pairarray,high_npairs,high_querypos+1) == false ||
	      (low_substring = Stage3end_substring_containing(hit_low,low_querypos)) == NULL ||
	      Stage3end_substring_containing(hit_low,low_querypos-1) != low_substring ||
	      Stage3end_substring_containing(hit_low,low_querypos+1) != low_substring ||
	      Pairarray_lookup(high_pairarray,high_npairs,high_querypos) != Substring_genomicstart(low_substring) + low_querypos - chroffset)) {
	(*shift) += 1;
	if ((low_substring = Stage3end_substring_containing(hit_low,low_querypos)) == NULL) {
	  low_querypos--;
	} else if (Pairarray_contains_p(high_pairarray,high_npairs,high_querypos) == false) {
	  high_querypos--;
	} else {
	  low_chrpos = Substring_genomicstart(low_substring) + low_querypos - chroffset;
	  high_chrpos = Pairarray_lookup(high_pairarray,high_npairs,high_querypos);
	  if (low_chrpos > high_chrpos) {
	    debug15(printf("low_chrpos %u > high_chrpos %u, so decreasing low_querypos\n",low_chrpos,high_chrpos));
	    low_querypos--;
	  } else if (high_chrpos > low_chrpos) {
	    debug15(printf("high_chrpos %u > low_chrpos %u, so decreasing high_querypos\n",high_chrpos,low_chrpos));
	    high_querypos--;
	  } else {
	    low_querypos--;
	    high_querypos--;
	  }
	}
	debug15(printf("left shift %d: Advancing to low_querypos %d and high_querypos %d\n",*shift,low_querypos,high_querypos));
      }

      if ((high_querypos - 1) < 0 || (low_querypos - 1) < 0 ||
	  Stage3end_substring_containing(hit_low,low_querypos) == NULL) {
	*shift = 0;
	return 0;
      } else {
	assert((low_substring = Stage3end_substring_containing(hit_low,low_querypos)) != NULL);
	assert(Stage3end_substring_containing(hit_low,low_querypos-1) == low_substring);
	assert(Stage3end_substring_containing(hit_low,low_querypos+1) == low_substring);
	assert(Pairarray_contains_p(high_pairarray,high_npairs,high_querypos) == true);
	assert(Pairarray_contains_p(high_pairarray,high_npairs,high_querypos-1) == true);
	assert(Pairarray_contains_p(high_pairarray,high_npairs,high_querypos+1) == true);
	return Pairarray_lookup(high_pairarray,high_npairs,high_querypos) + chroffset;
      }

    } else {
      low_querypos = low_querylength /*- 1*/ - hardclip_low;
      high_querypos = hardclip_high;
      debug15(printf("High GMAP, minus.  low_querypos %d, high_querypos %d\n",low_querypos,high_querypos));

      low_querypos++;
      high_querypos++;
      debug15(printf("left shift %d: Advancing to low_querypos %d and high_querypos %d\n",*shift,low_querypos,high_querypos));
      while ((high_querypos + 1) < high_querylength && (low_querypos + 1) < low_querylength &&
	     (Pairarray_contains_p(high_pairarray,high_npairs,high_querypos) == false ||
	      Pairarray_contains_p(high_pairarray,high_npairs,high_querypos-1) == false ||
	      Pairarray_contains_p(high_pairarray,high_npairs,high_querypos+1) == false ||
	      (low_substring = Stage3end_substring_containing(hit_low,low_querypos)) == NULL ||
	      Stage3end_substring_containing(hit_low,low_querypos-1) != low_substring ||
	      Stage3end_substring_containing(hit_low,low_querypos+1) != low_substring ||
	      Pairarray_lookup(high_pairarray,high_npairs,high_querypos) != (Substring_genomicstart(low_substring) - 1) - low_querypos - chroffset)) {
	(*shift) += 1;
	if ((low_substring = Stage3end_substring_containing(hit_low,low_querypos)) == NULL) {
	  low_querypos++;
	} else if (Pairarray_contains_p(high_pairarray,high_npairs,high_querypos) == false) {
	  high_querypos++;
	} else {
	  low_chrpos = (Substring_genomicstart(low_substring) - 1) - low_querypos - chroffset;
	  high_chrpos = Pairarray_lookup(high_pairarray,high_npairs,high_querypos);
	  if (low_chrpos > high_chrpos) {
	    debug15(printf("low_chrpos %u > high_chrpos %u, so advancing low_querypos\n",low_chrpos,high_chrpos));
	    low_querypos++;
	  } else if (high_chrpos > low_chrpos) {
	    debug15(printf("high_chrpos %u > low_chrpos %u, so advancing high_querypos\n",high_chrpos,low_chrpos));
	    high_querypos++;
	  } else {
	    low_querypos++;
	    high_querypos++;
	  }
	}
	debug15(printf("left shift %d: Advancing to low_querypos %d and high_querypos %d\n",*shift,low_querypos,high_querypos));
      }

      if ((high_querypos + 1) >= high_querylength || (low_querypos + 1) >= low_querylength ||
	  Stage3end_substring_containing(hit_low,low_querypos) == NULL) {
	*shift = 0;
	return 0;
      } else {
	assert((low_substring = Stage3end_substring_containing(hit_low,low_querypos)) != NULL);
	assert(Stage3end_substring_containing(hit_low,low_querypos-1) == low_substring);
	assert(Stage3end_substring_containing(hit_low,low_querypos+1) == low_substring);
	assert(Pairarray_contains_p(high_pairarray,high_npairs,high_querypos) == true);
	assert(Pairarray_contains_p(high_pairarray,high_npairs,high_querypos-1) == true);
	assert(Pairarray_contains_p(high_pairarray,high_npairs,high_querypos+1) == true);
	return Pairarray_lookup(high_pairarray,high_npairs,high_querypos) + chroffset;
      }
    }

  } else {
    if (plusp == true) {
      low_querypos = hardclip_low;
      high_querypos = high_querylength /*- 1*/ - hardclip_high;
      debug15(printf("Both substrings, plus.  low_querypos %d, high_querypos %d\n",low_querypos,high_querypos));

      low_querypos--;
      high_querypos--;
      debug15(printf("left shift %d: Advancing to low_querypos %d and high_querypos %d\n",*shift,low_querypos,high_querypos));
      while ((low_querypos - 1) >= 0 && (high_querypos - 1) >= 0 &&
	     ((low_substring = Stage3end_substring_containing(hit_low,low_querypos)) == NULL ||
	      Stage3end_substring_containing(hit_low,low_querypos-1) != low_substring ||
	      Stage3end_substring_containing(hit_low,low_querypos+1) != low_substring ||
	      (high_substring = Stage3end_substring_containing(hit_high,high_querypos)) == NULL ||
	      Stage3end_substring_containing(hit_high,high_querypos-1) != high_substring ||
	      Stage3end_substring_containing(hit_high,high_querypos+1) != high_substring ||
	      Substring_genomicstart(low_substring) + low_querypos - chroffset != Substring_genomicstart(high_substring) + high_querypos - chroffset)) {
	(*shift) += 1;
	if ((low_substring = Stage3end_substring_containing(hit_low,low_querypos)) == NULL) {
	  low_querypos--;
	} else if ((high_substring = Stage3end_substring_containing(hit_high,high_querypos)) == NULL) {
	  high_querypos--;
	} else {
	  low_chrpos = Substring_genomicstart(low_substring) + low_querypos - chroffset;
	  high_chrpos = Substring_genomicstart(high_substring) + high_querypos - chroffset;
	  if (low_chrpos > high_chrpos) {
	    debug15(printf("low_chrpos %u > high_chrpos %u, so decreasing low_querypos\n",low_chrpos,high_chrpos));
	    low_querypos--;
	  } else if (high_chrpos > low_chrpos) {
	    debug15(printf("high_chrpos %u > low_chrpos %u, so decreasing high_querypos\n",high_chrpos,low_chrpos));
	    high_querypos--;
	  } else {
	    low_querypos--;
	    high_querypos--;
	  }
	}
	debug15(printf("left shift %d: Advancing to low_querypos %d and high_querypos %d\n",*shift,low_querypos,high_querypos));
      }

      if ((low_querypos - 1) < 0 || (high_querypos - 1) < 0 ||
	  (low_substring = Stage3end_substring_containing(hit_low,low_querypos)) == NULL ||
	  Stage3end_substring_containing(hit_high,high_querypos) == NULL) {
	*shift = 0;
	return 0;
      } else {
	debug15(printf("Returning %u + %d\n",Substring_genomicstart(low_substring) - chroffset,
		       low_querypos));
	assert((low_substring = Stage3end_substring_containing(hit_low,low_querypos)) != NULL);
	assert((high_substring = Stage3end_substring_containing(hit_high,high_querypos)) != NULL);
	assert(Stage3end_substring_containing(hit_low,low_querypos-1) == low_substring);
	assert(Stage3end_substring_containing(hit_low,low_querypos+1) == low_substring);
	assert(Stage3end_substring_containing(hit_high,high_querypos-1) == high_substring);
	assert(Stage3end_substring_containing(hit_high,high_querypos+1) == high_substring);
	return Substring_genomicstart(low_substring) + low_querypos; /* Want univcoord */
      }

    } else {
      low_querypos = low_querylength /*- 1*/ - hardclip_low;
      high_querypos = hardclip_high;
      debug15(printf("Both substrings, minus.  low_querypos %d, high_querypos %d\n",low_querypos,high_querypos));

      low_querypos++;
      high_querypos++;
      debug15(printf("left shift %d: Advancing to low_querypos %d and high_querypos %d\n",*shift,low_querypos,high_querypos));
      while ((low_querypos + 1) < low_querylength && (high_querypos + 1) < high_querylength &&
	     ((low_substring = Stage3end_substring_containing(hit_low,low_querypos)) == NULL ||
	      Stage3end_substring_containing(hit_low,low_querypos-1) != low_substring ||
	      Stage3end_substring_containing(hit_low,low_querypos+1) != low_substring ||
	      (high_substring = Stage3end_substring_containing(hit_high,high_querypos)) == NULL ||
	      Stage3end_substring_containing(hit_high,high_querypos-1) != high_substring ||
	      Stage3end_substring_containing(hit_high,high_querypos+1) != high_substring ||
	      (Substring_genomicstart(low_substring) - 1) - low_querypos - chroffset != (Substring_genomicstart(high_substring) - 1) - high_querypos - chroffset)) {
	(*shift) += 1;
	if ((low_substring = Stage3end_substring_containing(hit_low,low_querypos)) == NULL) {
	  low_querypos++;
	} else if ((high_substring = Stage3end_substring_containing(hit_high,high_querypos)) == NULL) {
	  high_querypos++;
	} else {
	  low_chrpos = (Substring_genomicstart(low_substring) - 1) - low_querypos - chroffset;
	  high_chrpos = (Substring_genomicstart(high_substring) - 1) - high_querypos - chroffset;
	  if (low_chrpos > high_chrpos) {
	    debug15(printf("low_chrpos %u > high_chrpos %u, so advancing low_querypos\n",low_chrpos,high_chrpos));
	    low_querypos++;
	  } else if (high_chrpos > low_chrpos) {
	    debug15(printf("high_chrpos %u > low_chrpos %u, so advancing high_querypos\n",high_chrpos,low_chrpos));
	    high_querypos++;
	  } else {
	    low_querypos++;
	    high_querypos++;
	  }
	}
	debug15(printf("left shift %d: Advancing to low_querypos %d and high_querypos %d\n",*shift,low_querypos,high_querypos));
      }

      if ((low_querypos + 1) >= low_querylength || (high_querypos + 1) >= high_querylength ||
	  (low_substring = Stage3end_substring_containing(hit_low,low_querypos)) == NULL ||
	  Stage3end_substring_containing(hit_high,high_querypos) == NULL) {
	*shift = 0;
	return 0;
      } else {
	debug15(printf("Returning %u - %d\n",Substring_genomicstart(low_substring) - chroffset,
		       low_querypos));
	assert((low_substring = Stage3end_substring_containing(hit_low,low_querypos)) != NULL);
	assert((high_substring = Stage3end_substring_containing(hit_high,high_querypos)) != NULL);
	assert(Stage3end_substring_containing(hit_low,low_querypos-1) == low_substring);
	assert(Stage3end_substring_containing(hit_low,low_querypos+1) == low_substring);
	assert(Stage3end_substring_containing(hit_high,high_querypos-1) == high_substring);
	assert(Stage3end_substring_containing(hit_high,high_querypos+1) == high_substring);
	return (Substring_genomicstart(low_substring) - 1) - low_querypos; /* Want univcoord */
      }
    }
  }
}



/* Note: Do not alter this->insertlength, which is used for SAM
   output.  The insertlength computed here is used only for performing
   --clip-overlap or --merge-overlap */
int
Stage3pair_overlap (int *hardclip5_low, int *hardclip5_high, int *hardclip3_low, int *hardclip3_high, Stage3pair_T this) {
  Stage3end_T hit5, hit3;
  int overlap;
  int clipdir;
  int ilength53, ilength35, ilength5_low, ilength5_high, ilength3_low, ilength3_high;
  int common_shift, common_left, common_right;
  Univcoord_T common_genomicpos, common_genomicpos_right, common_genomicpos_left;
  int shift_right, shift_left;


  *hardclip5_low = *hardclip5_high = *hardclip3_low = *hardclip3_high = 0;

  hit5 = this->hit5;
  hit3 = this->hit3;

  debug15(printf("Entered Stage3pair_overlap with hittype %s and %s\n",
		 hittype_string(hit5->hittype),hittype_string(hit3->hittype)));
  if (hit5->hittype == SAMECHR_SPLICE || hit5->hittype == TRANSLOC_SPLICE) {
    return 0;
  } else if (hit3->hittype == SAMECHR_SPLICE || hit3->hittype == TRANSLOC_SPLICE) {
    return 0;
  } else if (hit5->plusp != hit3->plusp) {
    debug15(printf("The two ends are not on the same strand, so returning 0\n"));
    return 0;
  } else {
    debug15(printf("hit5 trim_left %d + amb_start %d, trim_right %d + amb_end %d, hit3 trim_left %d + amb_start %d, trim_right %d + amb_end %d\n",
		   hit5->trim_left,start_amb_length(hit5),hit5->trim_right,end_amb_length(hit5),
		   hit3->trim_left,start_amb_length(hit3),hit3->trim_right,end_amb_length(hit3)));
    if (hit5->plusp == true) {
      /* plus */
#if 0
      hit5_trimmed_length = hit5->querylength - hit5->trim_left - hit5->trim_right - start_amb_length(hit5) - end_amb_length(hit5);
      hit3_trimmed_length = hit3->querylength - hit3->trim_left - hit3->trim_right - start_amb_length(hit3) - end_amb_length(hit3);
      totallength = hit5_trimmed_length + hit3_trimmed_length;
      debug15(printf("totallength = %d, hit5 trimmed length = %d, hit3 trimmed length = %d\n",
		     totallength,hit5_trimmed_length,hit3_trimmed_length));
      debug15(printf("original insertlength: %d, trim+amb5: %d..%d, trim+amb3: %d..%d\n",
		     this->insertlength,hit5->trim_left + start_amb_length(hit5),
		     hit5->trim_right + end_amb_length(hit5),hit3->trim_left + start_amb_length(hit3),
		     hit3->trim_right + end_amb_length(hit3)));
#endif

      if ((common_genomicpos = pair_common_genomicpos(hit5,hit3)) == 0) {
	debug15(printf("Cannot determine a common point, so returning 0\n"));
	return 0;

      } else if (find_ilengths(&ilength5_low,&ilength5_high,hit5,common_genomicpos,hit5->chroffset) == false ||
		 find_ilengths(&ilength3_low,&ilength3_high,hit3,common_genomicpos,hit3->chroffset) == false) {
	debug15(printf("Cannot determine ilengths, so returning 0\n"));
	return 0;

      } else {
	debug15(printf("Inclusive: ilengths5: %d|%d.  ilengths3: %d|%d\n",ilength5_low,ilength5_high,ilength3_low,ilength3_high));
	debug15(printf("ilength53 is %d, ilength 35 is %d\n",ilength5_low + ilength3_high - 1,ilength3_low + ilength5_high - 1));

	common_left = (ilength5_low < ilength3_low) ? ilength5_low : ilength3_low;
	common_right = (ilength5_high < ilength3_high) ? ilength5_high : ilength3_high;
	if (common_right > common_left) {
	  common_shift = common_right/2 - (common_left - 1)/2;
	  debug15(printf("Common shift is %d = common_right %d/2 - (common_left %d - 1)/2\n",
			 common_shift,common_right,common_left));
	  ilength5_low -= 1;
	  ilength3_low -= 1;
	} else {
	  common_shift = (common_right - 1)/2 - common_left/2;
	  debug15(printf("Common shift is %d = (common_right %d - 1)/2 - common_left %d/2\n",
			 common_shift,common_right,common_left));
	  ilength5_high -= 1;
	  ilength3_high -= 1;
	}
	debug15(printf("Exclusive: ilengths5: %d|%d.  ilengths3: %d|%d\n",ilength5_low,ilength5_high,ilength3_low,ilength3_high));


	if ((ilength53 = ilength5_low + ilength3_high) >= (ilength35 = ilength3_low + ilength5_high)) {
	  /* Use >=, not >, so we favor clipping heads over clipping tails in case of a tie */
	  debug15(printf("plus, ilength53 is longer.  Clipping heads.\n"));
	  debug15(printf("Overlap is %d = common_left %d + common_right %d - 1\n",
			 common_left+common_right-1,common_left,common_right));
	  clipdir = +1;

	  /* Want to clip 5 high and 3 low */
	  *hardclip5_high = ilength5_high - common_shift;
	  *hardclip3_low = ilength3_low + common_shift;
	  debug15(printf("Overlap clip for ilength53 plus is hardclip5 %d..%d and hardclip3 %d..%d\n",
			 *hardclip5_low,*hardclip5_high,*hardclip3_low,*hardclip3_high));
	  *hardclip5_high += hit5->trim_right + end_amb_length(hit5);
	  *hardclip3_low += hit3->trim_left + start_amb_length(hit3);
	  debug15(printf("Ambig clip for ilength53 plus is hardclip5 %d..%d and hardclip3 %d..%d\n",
			 *hardclip5_low,*hardclip5_high,*hardclip3_low,*hardclip3_high));

	  if (common_shift != 0) {
	    if (test_hardclips(&common_genomicpos,*hardclip3_low,hit3,*hardclip5_high,hit5,hit3->chroffset) == true) {
	      /* No adjustment needed, but need to recompute ilengths for shifted common_genomicpos */
	    } else {
	      common_genomicpos_right = adjust_hardclips_right(&shift_right,*hardclip3_low,hit3,*hardclip5_high,hit5,hit3->chroffset);
	      common_genomicpos_left = adjust_hardclips_left(&shift_left,*hardclip3_low,hit3,*hardclip5_high,hit5,hit3->chroffset);
	      debug15(printf("shift_right %d, shift_left %d\n",shift_right,shift_left));
	      if (shift_right == 0 && shift_left == 0) {
		/* Try original position without a shift */
		*hardclip5_high = ilength5_high /*- common_shift*/;
		*hardclip3_low = ilength3_low /*+ common_shift*/;
		*hardclip5_high += hit5->trim_right + end_amb_length(hit5);
		*hardclip3_low += hit3->trim_left + start_amb_length(hit3);
		if (test_hardclips(&common_genomicpos,*hardclip3_low,hit3,*hardclip5_high,hit5,hit3->chroffset) == false) {
		  *hardclip5_low = *hardclip5_high = *hardclip3_low = *hardclip3_high = 0;
		  return 0;
		}
	      } else if (shift_left == 0) {
		common_genomicpos = common_genomicpos_right;
	      } else if (shift_right == 0) {
		common_genomicpos = common_genomicpos_left;
	      } else if (shift_right <= shift_left) {
		common_genomicpos = common_genomicpos_right;
	      } else {
		common_genomicpos = common_genomicpos_left;
	      }
	    }

	    debug15(printf("New common point is %u\n",common_genomicpos - hit3->chroffset));
	    /* Recompute hardclips */
	    if (find_ilengths(&ilength5_low,&ilength5_high,hit5,common_genomicpos,hit5->chroffset) == false ||
		find_ilengths(&ilength3_low,&ilength3_high,hit3,common_genomicpos,hit3->chroffset) == false) {
	      *hardclip5_low = *hardclip5_high = *hardclip3_low = *hardclip3_high = 0;
	      return 0;
	    } else if (ilength3_low > ilength5_high) {
	      debug15(printf("Uneven: ilengths5: %d|%d.  ilengths3: %d|%d\n",ilength5_low,ilength5_high,ilength3_low,ilength3_high));
	      ilength3_low -= 1;
	    } else {
	      debug15(printf("Uneven: ilengths5: %d|%d.  ilengths3: %d|%d\n",ilength5_low,ilength5_high,ilength3_low,ilength3_high));
	      ilength5_high -= 1;
	    }
	    debug15(printf("Even: ilengths5: %d|%d.  ilengths3: %d|%d\n",ilength5_low,ilength5_high,ilength3_low,ilength3_high));

	    *hardclip5_high = ilength5_high /*- common_shift*/;
	    *hardclip3_low = ilength3_low /*+ common_shift*/;
	    *hardclip5_high += hit5->trim_right + end_amb_length(hit5);
	    *hardclip3_low += hit3->trim_left + start_amb_length(hit3);

	    debug15(printf("Recomputed clip for ilength53 plus is hardclip5 %d..%d and hardclip3 %d..%d\n",
			   *hardclip5_low,*hardclip5_high,*hardclip3_low,*hardclip3_high));
	  }

	  if (hit5->hittype == GMAP || hit3->hittype == GMAP) {
	    /* Revise only for paired-ends involving GMAP and when successful.  Observed to be the correct action. */
	    this->insertlength = ilength53;
	  }

#if 0
	  if (*hardclip5_high < 0) {
	    *hardclip5_high = 0;
	  }
	  if (*hardclip3_low < 0) {
	    *hardclip3_low = 0;
	  }
	  debug15(printf("Positive clip for ilength53 plus is hardclip5 %d..%d and hardclip3 %d..%d\n",
			 *hardclip5_low,*hardclip5_high,*hardclip3_low,*hardclip3_high));
#endif

	} else {
	  debug15(printf("plus, ilength35 is longer.  Clipping tails.\n"));
	  debug15(printf("Overlap is %d = common_left %d + common_right %d - 1\n",
			 common_left+common_right-1,common_left,common_right));
	  clipdir = -1;

	  /* Want to clip 5 low and 3 high */
	  *hardclip5_low = ilength5_low + common_shift;
	  *hardclip3_high = ilength3_high - common_shift;
	  debug15(printf("Overlap clip for ilength35 plus is hardclip5 %d..%d and hardclip3 %d..%d\n",
			 *hardclip5_low,*hardclip5_high,*hardclip3_low,*hardclip3_high));
	  *hardclip5_low += hit5->trim_left + start_amb_length(hit5);
	  *hardclip3_high += hit3->trim_right + end_amb_length(hit3);
	  debug15(printf("Ambig clip for ilength35 plus is hardclip5 %d..%d and hardclip3 %d..%d\n",
			 *hardclip5_low,*hardclip5_high,*hardclip3_low,*hardclip3_high));

	  if (common_shift != 0) {
	    if (test_hardclips(&common_genomicpos,*hardclip5_low,hit5,*hardclip3_high,hit3,hit3->chroffset) == true) {
	      /* No adjustment needed, but need to recompute ilengths for shifted common_genomicpos */
	    } else {
	      common_genomicpos_right = adjust_hardclips_right(&shift_right,*hardclip5_low,hit5,*hardclip3_high,hit3,hit3->chroffset);
	      common_genomicpos_left = adjust_hardclips_left(&shift_left,*hardclip5_low,hit5,*hardclip3_high,hit3,hit3->chroffset);
	      debug15(printf("shift_right %d, shift_left %d\n",shift_right,shift_left));
	      if (shift_right == 0 && shift_left == 0) {
		/* Try original position without a shift */
		*hardclip5_low = ilength5_low /*+ common_shift*/;
		*hardclip3_high = ilength3_high /*- common_shift*/;
		*hardclip5_low += hit5->trim_left + start_amb_length(hit5);
		*hardclip3_high += hit3->trim_right + end_amb_length(hit3);
		if (test_hardclips(&common_genomicpos,*hardclip3_low,hit3,*hardclip5_high,hit5,hit3->chroffset) == false) {
		  *hardclip5_low = *hardclip5_high = *hardclip3_low = *hardclip3_high = 0;
		  return 0;
		}
	      } else if (shift_left == 0) {
		common_genomicpos = common_genomicpos_right;
	      } else if (shift_right == 0) {
		common_genomicpos = common_genomicpos_left;
	      } else if (shift_right <= shift_left) {
		common_genomicpos = common_genomicpos_right;
	      } else {
		common_genomicpos = common_genomicpos_left;
	      }
	    }

	    debug15(printf("New common point is %u\n",common_genomicpos - hit3->chroffset));
	    /* Recompute hardclips */
	    if (find_ilengths(&ilength5_low,&ilength5_high,hit5,common_genomicpos,hit5->chroffset) == false ||
		find_ilengths(&ilength3_low,&ilength3_high,hit3,common_genomicpos,hit3->chroffset) == false) {
	      *hardclip5_low = *hardclip5_high = *hardclip3_low = *hardclip3_high = 0;
	      return 0;
	    } else if (ilength5_low > ilength3_high) {
	      debug15(printf("Uneven: ilengths5: %d|%d.  ilengths3: %d|%d\n",ilength5_low,ilength5_high,ilength3_low,ilength3_high));
	      ilength5_low -= 1;
	    } else {
	      debug15(printf("Uneven: ilengths5: %d|%d.  ilengths3: %d|%d\n",ilength5_low,ilength5_high,ilength3_low,ilength3_high));
	      ilength3_high -= 1;
	    }
	    debug15(printf("Even: ilengths5: %d|%d.  ilengths3: %d|%d\n",ilength5_low,ilength5_high,ilength3_low,ilength3_high));

	    *hardclip5_low = ilength5_low /*+ common_shift*/;
	    *hardclip3_high = ilength3_high /*- common_shift*/;
	    *hardclip5_low += hit5->trim_left + start_amb_length(hit5);
	    *hardclip3_high += hit3->trim_right + end_amb_length(hit3);
	    debug15(printf("Recomputed clip for ilength35 plus is hardclip5 %d..%d and hardclip3 %d..%d\n",
			   *hardclip5_low,*hardclip5_high,*hardclip3_low,*hardclip3_high));
	  }

	  if (hit5->hittype == GMAP || hit3->hittype == GMAP) {
	    /* Revise only for paired-ends involving GMAP and when successful.  Observed to be the correct action. */
	    this->insertlength = ilength35;
	  }

#if 0
	  if (*hardclip5_low < 0) {
	    *hardclip5_low = 0;
	  }
	  if (*hardclip3_high < 0) {
	    *hardclip3_high = 0;
	  }
	  debug15(printf("Positive clip for ilength35 plus is hardclip5 %d..%d and hardclip3 %d..%d\n",
			 *hardclip5_low,*hardclip5_high,*hardclip3_low,*hardclip3_high));
#endif
	}

	debug15(printf("returning clipdir %d\n",clipdir));
	return clipdir;
      }

    } else {
      /* minus */
#if 0
      hit5_trimmed_length = hit5->querylength - hit5->trim_left - hit5->trim_right - start_amb_length(hit5) - end_amb_length(hit5);
      hit3_trimmed_length = hit3->querylength - hit3->trim_left - hit3->trim_right - start_amb_length(hit3) - end_amb_length(hit3);
      totallength = hit5_trimmed_length + hit3_trimmed_length;
      debug15(printf("totallength = %d, hit5 trimmed length = %d, hit3 trimmed length = %d\n",
		     totallength,hit5_trimmed_length,hit3_trimmed_length));
      debug15(printf("original insertlength: %d, trim+amb5: %d..%d, trim+amb3: %d..%d\n",
		     this->insertlength,hit5->trim_left + start_amb_length(hit5),
		     hit5->trim_right + hit5->end_amb_length,hit3->trim_left + start_amb_length(hit3),
		     hit3->trim_right + hit3->end_amb_length));
#endif

      if ((common_genomicpos = pair_common_genomicpos(hit5,hit3)) == 0) {
	debug15(printf("Cannot determine a common point, so returning 0\n"));
	return 0;

      } else if (find_ilengths(&ilength5_low,&ilength5_high,hit5,common_genomicpos,hit5->chroffset) == false ||
		 find_ilengths(&ilength3_low,&ilength3_high,hit3,common_genomicpos,hit3->chroffset) == false) {
	debug15(printf("Cannot determine ilengths, so returning 0\n"));
	return 0;

      } else {
	debug15(printf("Inclusive: ilengths5: %d|%d.  ilengths3: %d|%d\n",ilength5_low,ilength5_high,ilength3_low,ilength3_high));
	debug15(printf("ilength53lh is %d, ilength35lh is %d\n",ilength5_low + ilength3_high - 1,ilength3_low + ilength5_high - 1));

	common_left = (ilength5_low < ilength3_low) ? ilength5_low : ilength3_low;
	common_right = (ilength5_high < ilength3_high) ? ilength5_high : ilength3_high;
	if (common_right > common_left) {
	  common_shift = common_right/2 - (common_left - 1)/2;
	  debug15(printf("Common shift is %d = common_right %d/2 - (common_left %d - 1)/2\n",
			 common_shift,common_right,common_left));
	  ilength5_low -= 1;
	  ilength3_low -= 1;
	} else {
	  common_shift = (common_right - 1)/2 - common_left/2;
	  debug15(printf("Common shift is %d = (common_right %d - 1)/2 - common_left %d/2\n",
			 common_shift,common_right,common_left));
	  ilength5_high -= 1;
	  ilength3_high -= 1;
	}
	debug15(printf("Exclusive: ilengths5: %d|%d.  ilengths3: %d|%d\n",ilength5_low,ilength5_high,ilength3_low,ilength3_high));

	if ((ilength53 = ilength5_low + ilength3_high) > (ilength35 = ilength3_low + ilength5_high)) {
	  /* Use >, not >=, so we favor clipping heads over clipping tails in case of a tie */
	  debug15(printf("minus, ilength53 is longer.  Clipping tails.\n"));
	  debug15(overlap = common_left + common_right - 1);
	  debug15(printf("Overlap is %d = common_left %d + common_right %d - 1\n",
			 overlap,common_left,common_right));
	  clipdir = +1;


	  /* Want to clip 5 high and 3 low */
	  *hardclip5_high = ilength5_high - common_shift;
	  *hardclip3_low = ilength3_low + common_shift;
	  debug15(printf("Overlap clip for ilength53 minus is hardclip5 %d..%d and hardclip3 %d..%d\n",
			 *hardclip5_low,*hardclip5_high,*hardclip3_low,*hardclip3_high));
	  *hardclip5_high += hit5->trim_left + start_amb_length(hit5);
	  *hardclip3_low += hit3->trim_right + end_amb_length(hit3);
	  debug15(printf("Ambig clip for ilength53 minus is hardclip5 %d..%d and hardclip3 %d..%d\n",
			 *hardclip5_low,*hardclip5_high,*hardclip3_low,*hardclip3_high));

	  if (common_shift != 0) {
	    if (test_hardclips(&common_genomicpos,*hardclip3_low,hit3,*hardclip5_high,hit5,hit3->chroffset) == true) {
	      /* No adjustment needed, but need to recompute ilengths for shifted common_genomicpos */
	    } else {
	      common_genomicpos_right = adjust_hardclips_right(&shift_right,*hardclip3_low,hit3,*hardclip5_high,hit5,hit3->chroffset);
	      common_genomicpos_left = adjust_hardclips_left(&shift_left,*hardclip3_low,hit3,*hardclip5_high,hit5,hit3->chroffset);
	      debug15(printf("shift_right %d, shift_left %d\n",shift_right,shift_left));
	      if (shift_right == 0 && shift_left == 0) {
		/* Try original position without a shift */
		*hardclip5_high = ilength5_high /*- common_shift*/;
		*hardclip3_low = ilength3_low /*+ common_shift*/;
		*hardclip5_high += hit5->trim_left + start_amb_length(hit5);
		*hardclip3_low += hit3->trim_right + end_amb_length(hit3);
		if (test_hardclips(&common_genomicpos,*hardclip3_low,hit3,*hardclip5_high,hit5,hit3->chroffset) == false) {
		  *hardclip5_low = *hardclip5_high = *hardclip3_low = *hardclip3_high = 0;
		  return 0;
		}
	      } else if (shift_left == 0) {
		common_genomicpos = common_genomicpos_right;
	      } else if (shift_right == 0) {
		common_genomicpos = common_genomicpos_left;
	      } else if (shift_right <= shift_left) {
		common_genomicpos = common_genomicpos_right;
	      } else {
		common_genomicpos = common_genomicpos_left;
	      }
	    }

	    debug15(printf("New common point is %u\n",common_genomicpos - hit3->chroffset));
	    /* Recompute hardclips */
	    if (find_ilengths(&ilength5_low,&ilength5_high,hit5,common_genomicpos,hit5->chroffset) == false ||
		find_ilengths(&ilength3_low,&ilength3_high,hit3,common_genomicpos,hit3->chroffset) == false) {
	      *hardclip5_low = *hardclip5_high = *hardclip3_low = *hardclip3_high = 0;
	      return 0;
	    } else if (ilength3_low > ilength5_high) {
	      debug15(printf("Uneven: ilengths5: %d|%d.  ilengths3: %d|%d\n",ilength5_low,ilength5_high,ilength3_low,ilength3_high));
	      ilength3_low -= 1;
	    } else {
	      debug15(printf("Uneven: ilengths5: %d|%d.  ilengths3: %d|%d\n",ilength5_low,ilength5_high,ilength3_low,ilength3_high));
	      ilength5_high -= 1;
	    }
	    debug15(printf("Even: ilengths5: %d|%d.  ilengths3: %d|%d\n",ilength5_low,ilength5_high,ilength3_low,ilength3_high));

	    *hardclip5_high = ilength5_high /*- common_shift*/;
	    *hardclip3_low = ilength3_low /*+ common_shift*/;
	    *hardclip5_high += hit5->trim_left + start_amb_length(hit5);
	    *hardclip3_low += hit3->trim_right + end_amb_length(hit3);
	    debug15(printf("Recomputed clip for ilength53 minus is hardclip5 %d..%d and hardclip3 %d..%d\n",
			   *hardclip5_low,*hardclip5_high,*hardclip3_low,*hardclip3_high));
	  }

	  if (hit5->hittype == GMAP || hit3->hittype == GMAP) {
	    /* Revise only for paired-ends involving GMAP and when successful.  Observed to be the correct action. */
	    this->insertlength = ilength53;
	  }
#if 0
	  if (*hardclip5_high < 0) {
	    *hardclip5_high = 0;
	  }
	  if (*hardclip3_low < 0) {
	    *hardclip3_low = 0;
	  }
	  debug15(printf("Positive clip for ilength53 minus is hardclip5 %d..%d and hardclip3 %d..%d\n",
			 *hardclip5_low,*hardclip5_high,*hardclip3_low,*hardclip3_high));
#endif

	} else {
	  debug15(printf("minus, ilength35 is longer.  Clipping heads.\n"));
	  debug15(overlap = common_left + common_right - 1);
	  debug15(printf("Overlap is %d = common_left %d + common_right %d - 1\n",
			 overlap,common_left,common_right));
	  clipdir = -1;

	  /* Want to clip 5 low and 3 high */
	  *hardclip5_low = ilength5_low + common_shift;
	  *hardclip3_high = ilength3_high - common_shift;
	  debug15(printf("Overlap clip for ilength35 minus is hardclip5 %d..%d and hardclip3 %d..%d\n",
			 *hardclip5_low,*hardclip5_high,*hardclip3_low,*hardclip3_high));
	  *hardclip5_low += hit5->trim_right + end_amb_length(hit5);
	  *hardclip3_high += hit3->trim_left + start_amb_length(hit3);
	  debug15(printf("Ambig clip for ilength35 minus is hardclip5 %d..%d and hardclip3 %d..%d\n",
			 *hardclip5_low,*hardclip5_high,*hardclip3_low,*hardclip3_high));

	  if (common_shift != 0) {
	    if (test_hardclips(&common_genomicpos,*hardclip5_low,hit5,*hardclip3_high,hit3,hit3->chroffset) == true) {
	      /* No adjustment needed, but need to recompute ilengths for shifted common_genomicpos */
	    } else {
	      common_genomicpos_right = adjust_hardclips_right(&shift_right,*hardclip5_low,hit5,*hardclip3_high,hit3,hit3->chroffset);
	      common_genomicpos_left = adjust_hardclips_left(&shift_left,*hardclip5_low,hit5,*hardclip3_high,hit3,hit3->chroffset);
	      debug15(printf("shift_right %d, shift_left %d\n",shift_right,shift_left));
	      if (shift_right == 0 && shift_left == 0) {
		/* Try original position without a shift */
		*hardclip5_low = ilength5_low /*+ common_shift*/;
		*hardclip3_high = ilength3_high /*- common_shift*/;
		*hardclip5_low += hit5->trim_right + end_amb_length(hit5);
		*hardclip3_high += hit3->trim_left + start_amb_length(hit3);
		if (test_hardclips(&common_genomicpos,*hardclip3_low,hit3,*hardclip5_high,hit5,hit3->chroffset) == false) {
		  *hardclip5_low = *hardclip5_high = *hardclip3_low = *hardclip3_high = 0;
		  return 0;
		}
	      } else if (shift_left == 0) {
		common_genomicpos = common_genomicpos_right;
	      } else if (shift_right == 0) {
		common_genomicpos = common_genomicpos_left;
	      } else if (shift_right <= shift_left) {
		common_genomicpos = common_genomicpos_right;
	      } else {
		common_genomicpos = common_genomicpos_left;
	      }
	    }

	    debug15(printf("New common point is %u\n",common_genomicpos - hit3->chroffset));
	    /* Recompute hardclips */
	    if (find_ilengths(&ilength5_low,&ilength5_high,hit5,common_genomicpos,hit5->chroffset) == false ||
		find_ilengths(&ilength3_low,&ilength3_high,hit3,common_genomicpos,hit3->chroffset) == false) {
	      *hardclip5_low = *hardclip5_high = *hardclip3_low = *hardclip3_high = 0;
	      return 0;
	    } else if (ilength5_low > ilength3_high) {
	      debug15(printf("Uneven: ilengths5: %d|%d.  ilengths3: %d|%d\n",ilength5_low,ilength5_high,ilength3_low,ilength3_high));
	      ilength5_low -= 1;
	    } else {
	      debug15(printf("Uneven: ilengths5: %d|%d.  ilengths3: %d|%d\n",ilength5_low,ilength5_high,ilength3_low,ilength3_high));
	      ilength3_high -= 1;
	    }
	    debug15(printf("Even: ilengths5: %d|%d.  ilengths3: %d|%d\n",ilength5_low,ilength5_high,ilength3_low,ilength3_high));

	    *hardclip5_low = ilength5_low /*+ common_shift*/;
	    *hardclip3_high = ilength3_high /*- common_shift*/;
	    *hardclip5_low += hit5->trim_right + end_amb_length(hit5);
	    *hardclip3_high += hit3->trim_left + start_amb_length(hit3);
	    debug15(printf("Recomputed clip for ilength35 minus is hardclip5 %d..%d and hardclip3 %d..%d\n",
			   *hardclip5_low,*hardclip5_high,*hardclip3_low,*hardclip3_high));
	  }

	  if (hit5->hittype == GMAP || hit3->hittype == GMAP) {
	    /* Revise only for paired-ends involving GMAP and when successful.  Observed to be the correct action. */
	    this->insertlength = ilength35;
	  }

#if 0
	  if (*hardclip5_low < 0) {
	    *hardclip5_low = 0;
	  }
	  if (*hardclip3_high < 0) {
	    *hardclip3_high = 0;
	  }
	  debug15(printf("Positive clip for ilength35 minus is hardclip5 %d..%d and hardclip3 %d..%d\n",
			 *hardclip5_low,*hardclip5_high,*hardclip3_low,*hardclip3_high));
#endif
	}
      }

      debug15(printf("returning clipdir %d\n",clipdir));
      return clipdir;
    }
  }
}


void
Stage3pair_set_private5p (Stage3pair_T this) {
  this->private5p = true;
  return;
}

void
Stage3pair_clear_private5p (Stage3pair_T this) {
  this->private5p = false;
  return;
}

void
Stage3pair_set_private3p (Stage3pair_T this) {
  this->private3p = true;
  return;
}

void
Stage3pair_clear_private3p (Stage3pair_T this) {
  this->private3p = false;
  return;
}


void
Stage3pair_free (Stage3pair_T *old) {
  debug0(printf("Freeing pair %p with hits %p (privatep %d) and %p (privatep %d)\n",*old,(*old)->hit5,(*old)->private5p,(*old)->hit3,(*old)->private3p));
  if ((*old)->private3p == true) {
    assert((*old)->hit3 != NULL);
    debug0(printf("Freeing end3 at %p\n",(*old)->hit3));
    Stage3end_free(&(*old)->hit3);
  }
  if ((*old)->private5p == true) {
    assert((*old)->hit5 != NULL);
    debug0(printf("Freeing end5 at %p\n",(*old)->hit5));
    Stage3end_free(&(*old)->hit5);
  }
  FREE_OUT(*old);
  return;
}


static Overlap_T
Stage3pair_gene_overlap (Stage3pair_T this) {
  Overlap_T overlap;
  bool foundp = false;

  if (genes_iit == NULL) {
    return NO_KNOWN_GENE;
  } else {
    if ((overlap = Stage3end_gene_overlap(this->hit5)) == KNOWN_GENE_MULTIEXON) {
      return KNOWN_GENE_MULTIEXON;
    } else if (overlap == KNOWN_GENE) {
      if (favor_multiexon_p == false) {
	return KNOWN_GENE;
      } else {
	foundp = true;
      }
    }

    if ((overlap = Stage3end_gene_overlap(this->hit3)) == KNOWN_GENE_MULTIEXON) {
      return KNOWN_GENE_MULTIEXON;
    } else if (overlap == KNOWN_GENE) {
      if (favor_multiexon_p == false) {
	return KNOWN_GENE;
      } else {
	foundp = true;
      }
    }

    if (foundp == true) {
      return KNOWN_GENE;
    } else {
      return NO_KNOWN_GENE;
    }
  }
}

#if 0
static long int
Stage3pair_tally (Stage3pair_T this) {

  if (tally_iit == NULL) {
    return 0L;
  } else if (this->tally >= 0) {
    return this->tally;
  } else {
    this->tally = Stage3end_compute_tally(this->hit5) + Stage3end_compute_tally(this->hit3);
    return this->tally;
  }
}
#endif


static char complCode[128] = COMPLEMENT_LC;

static char *
make_complement_buffered (char *complement, char *sequence, unsigned int length) {
  int i, j;

  /* complement = (char *) CALLOC_OUT(length+1,sizeof(char)); */
  for (i = length-1, j = 0; i >= 0; i--, j++) {
    complement[j] = complCode[(int) sequence[i]];
  }
  complement[length] = '\0';
  return complement;
}


const Except_T Copy_Substring = { "Substring invalid during copy" };

T
Stage3end_copy (T old) {
  T new = (T) MALLOC_OUT(sizeof(*new));
  List_T p;
  Substring_T old_substring, new_substring;
  Junction_T old_junction, new_junction;

  debug0(printf("Copying Stage3end %p -> %p of type %s\n",
		old,new,hittype_string(old->hittype)));

  new->hittype = old->hittype;
  new->genestrand = old->genestrand;
  new->sarrayp = old->sarrayp;
  new->gmap_source = old->gmap_source;
  new->improved_by_gmap_p = old->improved_by_gmap_p;

  new->chrnum = old->chrnum;
  new->effective_chrnum = old->effective_chrnum;
  new->other_chrnum = old->other_chrnum;
  new->chroffset = old->chroffset;
  new->chrhigh = old->chrhigh;
  new->chrlength = old->chrlength;

  new->querylength = old->querylength;
  new->querylength_adj = old->querylength_adj;

  new->genomicstart = old->genomicstart;
  new->genomicend = old->genomicend;
  new->plusp = old->plusp;

  new->low = old->low;
  new->high = old->high;
  new->genomiclength = old->genomiclength;
  new->guided_insertlength = old->guided_insertlength;

  new->mapq_loglik = old->mapq_loglik;
  new->mapq_score = old->mapq_score;
  new->absmq_score = old->absmq_score;

  new->nsegments = old->nsegments;
  new->score = old->score;
  new->ntscore = old->ntscore;
  new->nmatches = old->nmatches;
  new->nmatches_posttrim = old->nmatches_posttrim;
  new->gmap_max_match_length = old->gmap_max_match_length;
  new->gmap_min_splice_prob = old->gmap_min_splice_prob;

  new->trim_left = old->trim_left;
  new->trim_right = old->trim_right;
  new->trim_left_splicep = old->trim_left_splicep;
  new->trim_right_splicep = old->trim_right_splicep;

  /* new->penalties = old->penalties; */
  new->score_eventrim = old->score_eventrim;

  new->gene_overlap = old->gene_overlap;
  new->tally = old->tally;

  new->nmismatches_whole = old->nmismatches_whole;
  new->nmismatches_bothdiff = old->nmismatches_bothdiff;
  new->nmismatches_refdiff = old->nmismatches_refdiff;

  new->nindels = old->nindels;

  new->distance = old->distance;
  new->shortexonA_distance = old->shortexonA_distance;
  new->shortexonD_distance = old->shortexonD_distance;

  new->gmap_nindelbreaks = old->gmap_nindelbreaks;
  new->gmap_cdna_direction = old->gmap_cdna_direction;
  new->gmap_nintrons = old->gmap_nintrons;
  new->sensedir = old->sensedir;

  new->gmap_start_amb_length = old->gmap_start_amb_length;
  new->gmap_end_amb_length = old->gmap_end_amb_length;
  new->gmap_start_endtype = old->gmap_start_endtype;
  new->gmap_end_endtype = old->gmap_end_endtype;

  new->nsplices = old->nsplices;

  new->substrings_1toN = (List_T) NULL;
  new->substrings_Nto1 = (List_T) NULL;
  new->substrings_LtoH = (List_T) NULL;

  new->junctions_1toN = (List_T) NULL;
  new->junctions_Nto1 = (List_T) NULL;
  new->junctions_LtoH = (List_T) NULL;

  if (old->hittype == GMAP) {
    new->pairarray = Pairpool_copy_array(old->pairarray,old->npairs);
    new->npairs = old->npairs;
    new->cigar_tokens = Pair_tokens_copy(old->cigar_tokens);
    new->gmap_intronp = old->gmap_intronp;

  } else {
    new->pairarray = (struct Pair_T *) NULL;
    new->npairs = 0;
    new->cigar_tokens = (List_T) NULL;
    new->gmap_intronp = false;

    for (p = old->substrings_1toN; p != NULL; p = List_next(p)) {
      old_substring = (Substring_T) List_head(p);
      new_substring = Substring_copy(old_substring);
      new->substrings_1toN = List_push(new->substrings_1toN,(void *) new_substring);
    }

    for (p = old->junctions_1toN; p != NULL; p = List_next(p)) {
      old_junction = (Junction_T) List_head(p);
      new_junction = Junction_copy(old_junction);
      new->junctions_1toN = List_push(new->junctions_1toN,(void *) new_junction);
    }

    new->substrings_Nto1 = List_copy(new->substrings_1toN); /* Before reversal of 1toN */
    new->junctions_Nto1 = List_copy(new->junctions_1toN);   /* Before reversal of 1toN */

    /* Reversals to handle builds of 1toN */
    new->substrings_1toN = List_reverse(new->substrings_1toN);
    new->junctions_1toN = List_reverse(new->junctions_1toN);

    if (old->chrnum == 0) {
      /* Translocation */
      if (old->sensedir == SENSE_FORWARD) {
	new->substrings_LtoH = List_copy(new->substrings_1toN);
	new->junctions_LtoH = List_copy(new->junctions_1toN);
      } else if (old->sensedir == SENSE_ANTI) {
	new->substrings_LtoH = List_copy(new->substrings_Nto1);
	new->junctions_LtoH = List_copy(new->junctions_Nto1);
      } else {
	abort();
      }

    } else {
      if (old->plusp == true) {
	new->substrings_LtoH = List_copy(new->substrings_1toN);
	new->junctions_LtoH = List_copy(new->junctions_1toN);
      } else {
	new->substrings_LtoH = List_copy(new->substrings_Nto1);
	new->junctions_LtoH = List_copy(new->junctions_Nto1);
      }
    }
    assert(Substring_querystart(List_head(new->substrings_1toN)) <= Substring_querystart(List_head(new->substrings_Nto1)));
  }

  new->paired_usedp = old->paired_usedp;
  new->paired_seenp = old->paired_seenp;
  new->concordantp = old->concordantp;

  new->alias = old->alias;
  new->circularpos = old->circularpos;

  return new;
}


static int
compute_circularpos (int *alias, T hit) {
  int circularpos;
  List_T p;
  Substring_T substring;


  debug12(printf("Computing circularpos on hit at %u..%u, plusp %d, with trim left %d and trim right %d\n",
		 hit->genomicstart - hit->chroffset,hit->genomicend - hit->chroffset,
		 hit->plusp,hit->trim_left,hit->trim_right));
  if (circularp[hit->chrnum] == false) {
    debug12(printf("Chromosome #%d is not circular\n",hit->chrnum));
    /* This also handles hit->chrnum == 0, where translocation cannot be circular */
    *alias = 0;
    return -1;

  } else if (hit->hittype == GMAP) {
    debug12(printf("Pair circularpos is %d\n",Pair_circularpos(&(*alias),hit->pairarray,hit->npairs,hit->chrlength,
							       hit->plusp,hit->querylength)));
    return Pair_circularpos(&(*alias),hit->pairarray,hit->npairs,hit->chrlength,
			    hit->plusp,hit->querylength);

  } else if (hit->plusp == true) {
    if (hit->low >= hit->chroffset + hit->chrlength) {
      /* hit->low + hit->trim_left >= hit->chroffset + hit->chrlength (SOFT_CLIPS_AVOID_CIRCULARIZATION) */
      /* hit->low >= hit->chroffset + hit->chrlength */

      /* All of read after trimming is in circular alias */
      debug12(printf("Soft clip of %d on left avoids circularization\n",hit->trim_left));
      if (hit->high >= hit->chrhigh) {
	*alias = +2;
      } else {
	*alias = +1;
      }
      return -1;

    } else if (hit->high < hit->chroffset + hit->chrlength) {
      /* hit->high - hit->trim_right <= hit->chroffset + hit->chrlength (SOFT_CLIPS_AVOID_CIRCULARIZATION) */
      /* hit->high <= hit->chroffset + hit->chrlength */

      /* All of read after trimming is in circular proper */
      debug12(printf("Soft clip of %d on right avoids circularization\n",hit->trim_right));
      if (hit->low < hit->chroffset) {
	*alias = -2;
      } else {
	*alias = -1;
      }
      return -1;

    } else {
      *alias = 0;
      for (p = hit->substrings_1toN; p != NULL; p = List_next(p)) {
	substring = (Substring_T) List_head(p);
	if ((circularpos = Substring_circularpos(substring)) > 0) {
	  debug12(printf("Returning circularpos %d from substring (plus)\n",circularpos));
	  return circularpos;
	}
      }
      return -1;
    }

  } else {
    /* Without SOFT_CLIPS_AVOID_CIRCULARIZATION, the branches are similar for plus/minus hits */
    if (hit->low >= hit->chroffset + hit->chrlength) {
      /* hit->low + hit->trim_right >= hit->chroffset + hit->chrlength (SOFT_CLIPS_AVOID_CIRCULARIZATION) */
      /* hit->low >= hit->chroffset + hit->chrlength */

      /* All of read after trimming is in circular alias */
      debug12(printf("Soft clip of %d on right avoids circularization\n",hit->trim_right));
      debug12(printf("All of read after trimming is in circular alias\n"));
      if (hit->high >= hit->chrhigh) {
	*alias = +2;
      } else {
	*alias = +1;
      }
      return -1;

    } else if (hit->high < hit->chroffset + hit->chrlength) {
      /* hit->high - hit->trim_left <= hit->chroffset + hit->chrlength (SOFT_CLIPS_AVOID_CIRCULARIZATION) */
      /* hit->high <= hit->chroffset + hit->chrlength */

      /* All of read after trimming is in circular proper */
      debug12(printf("Soft clip of %d on left avoids circularization\n",hit->trim_left));
      debug12(printf("All of read after trimming is in circular proper\n"));
      if (hit->low < hit->chroffset) {
	*alias = -2;
      } else {
	*alias = -1;
      }
      return -1;

    } else {
      *alias = 0;
      for (p = hit->substrings_Nto1; p != NULL; p = List_next(p)) {
	substring = (Substring_T) List_head(p);
	if ((circularpos = Substring_circularpos(substring)) > 0) {
	  debug12(printf("Returning circularpos %d from substring (minus)\n",circularpos));
	  return circularpos;
	}
      }
      return -1;
    }
  }
}


T
Stage3end_new_substrings (int *found_score, Intlist_T endpoints,
#ifdef LARGE_GENOMES
			  Uint8list_T lefts,
#else
			  Uintlist_T lefts,
#endif
			  Intlist_T nmismatches_list, List_T junctions, int querylength,
			  Compress_T query_compress,
			  Substring_T right_ambig, Substring_T left_ambig,
			  bool plusp, int genestrand, int sensedir, bool first_read_p,
			  Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh,
			  Chrpos_T chrlength, bool sarrayp) {
  T new;

  Univcoord_T genomicstart, genomicend, genomicstart_adj, genomicend_adj,
    alignstart, alignend, alignstart_trim, alignend_trim;
  int querylength_trimmed = 0;
  int querystart, queryend;
  Univcoord_T left;
  Intlist_T r, x;
#ifdef LARGE_GENOMES
  Uint8list_T q;
#else
  Uintlist_T q;
#endif
  Substring_T substring, substring1, substringN;
  Junction_T junction, junction_ambig = NULL;
  Junctiontype_T type;
  List_T substrings = NULL, p, j;
  bool trim_left_p = false, trim_right_p = false;
  int outofbounds_start = 0, outofbounds_end = 0;
  int adj = 0, adj0;			/* deletions - insertions */
  int nmismatches_whole = 0, nmismatches, indel_score = 0, nindels = 0;
  int nmismatches_bothdiff = 0;


  debug0(printf("%s read: Entered Stage3end_new_substrings at left %u, with plusp %d, sensedir %d, and endpoints %s\n",
		first_read_p ? "First" : "Second",Uintlist_head(lefts),plusp,sensedir,Intlist_to_string(endpoints)));
  debug0(printf("There are %d endpoints, %d lefts, %d nmismatches, and %d junctions\n",
		Intlist_length(endpoints),Uintlist_length(lefts),Intlist_length(nmismatches_list),List_length(junctions)));
  debug0(printf("Ambig left %p, right %p\n",left_ambig,right_ambig));
  debug0(printf("Endpoints: %s\n",Intlist_to_string(endpoints)));
  debug0(printf("Lefts: %s\n",Uintlist_to_string(lefts)));
  debug0(printf("Mismatches: %s\n",Intlist_to_string(nmismatches_list)));

  assert(Uintlist_length(lefts) == Intlist_length(endpoints) - 1);
  assert(Intlist_length(nmismatches_list) == Intlist_length(endpoints) - 1);
  assert(List_length(junctions) == Intlist_length(endpoints) - 2);


#ifdef DEBUG0
  for (p = junctions; p != NULL; p = List_next(p)) {
    Junction_print((Junction_T) List_head(p));
  }
  printf("\n");
#endif

  querystart = Intlist_head(endpoints);
  if (plusp == true) {
    j = junctions;		/* Put here before we handle left_ambig */
    if (left_ambig != NULL) {
      substrings = List_push(substrings,(void *) left_ambig);
      junctions = List_push(junctions,(void *) Junction_new_splice(/*distance*/0,sensedir,
								   Substring_amb_donor_prob(left_ambig),
								   Substring_amb_acceptor_prob(left_ambig)));
    } else {
      trim_left_p = true;
    }

    /* Add querypos to get alignstart/alignend */
    for (q = lefts, x = nmismatches_list, r = Intlist_next(endpoints); q != NULL;
#ifdef LARGE_GENOMES
	 q = Uint8list_next(q),
#else
	 q = Uintlist_next(q),
#endif
	   x = Intlist_next(x), r = Intlist_next(r), j = List_next(j)) {
      queryend = Intlist_head(r);
#ifdef LARGE_GENOMES
      left = Uint8list_head(q);
#else
      left = Uintlist_head(q);
#endif
      debug0(printf("Working on querystart %d..queryend %d at left %u\n",querystart,queryend,left));

      genomicstart = left;
      genomicend = left + querylength;
      genomicstart_adj = genomicstart + adj;
      genomicend_adj = genomicend + adj;

      alignstart = genomicstart + querystart;
      alignend = genomicstart + queryend;

      if (genomicstart < chroffset && genomicend > chrhigh) {
	/* Out of bounds on both sides */
	Junction_gc(&junctions);
	return (T) NULL;

      } else if (genomicstart < chroffset) {
	outofbounds_start = chroffset - genomicstart;
	outofbounds_end = genomicend - chroffset;
	debug0(printf("Out of bounds left (low) %d, out of bounds right (high) %d\n",outofbounds_start,outofbounds_end));
	Junction_gc(&junctions);
	return (T) NULL;
#if 0
	/* Could consider this for the lowest substring */
	if (outofbounds_start > outofbounds_end) {
	  /* Consider high part to be out of bounds and keep existing chromosome */
	  outofbounds_start = 0;
	  Junction_gc(&junctions);
	  return (T) NULL;
	} else {
	  /* Consider low part to be out of bounds and stay in this chromosome */
	  /* Univ_IIT_interval_bounds(&chroffset,&chrhigh,&chrlength,chromosome_iit,chrnum,circular_typeint); */
	  outofbounds_end = 0;
	}
#endif
      } else if (genomicend > chrhigh) {
	outofbounds_start = chrhigh - genomicstart;
	outofbounds_end = genomicend - chrhigh;
	debug0(printf("Out of bounds left (low) %d, out of bounds right (high) %d\n",outofbounds_start,outofbounds_end));
	Junction_gc(&junctions);
	return (T) NULL;
#if 0
	/* Could consider this for the highest substring */
	if (outofbounds_start > outofbounds_end) {
	  /* Consider high part to be out of bounds and keep existing chromosome */
	  outofbounds_start = 0;
	} else if (++chrnum > nchromosomes) {
	  debug0(printf("Returning NULL from Stage3end_new_substrings\n"));
	  Junction_gc(&junctions);
	  return (T) NULL;
	} else {
	  /* Consider low part to be out of bounds and move to next chromosome */
	  Univ_IIT_interval_bounds(&chroffset,&chrhigh,&chrlength,chromosome_iit,chrnum,circular_typeint);
	  outofbounds_end = 0;
	}
#endif
      }

      if ((nmismatches = Intlist_head(x)) < 0) {
	nmismatches = Genome_count_mismatches_substring(query_compress,left,/*pos5*/querystart,/*pos3*/queryend,
							/*plusp*/true,genestrand,first_read_p);
	debug0(printf("nmismatches %d from genome\n",nmismatches));
      }
      nmismatches_whole += nmismatches;
      debug0(printf("nmismatches %d from sarray\n",nmismatches));
#ifdef LARGE_GENOMES
      if (Uint8list_next(q) == NULL && right_ambig == NULL) {
	trim_right_p = true;
      }
#else
      if (Uintlist_next(q) == NULL && right_ambig == NULL) {
	trim_right_p = true;
      }
#endif
      if ((substring = Substring_new(/*nmismatches_whole*/nmismatches,chrnum,chroffset,chrhigh,chrlength,
				     query_compress,/*start_endtype*/END,/*end_endtype*/END,
				     querystart,queryend,querylength,alignstart,alignend,
				     /*genomiclength*/querylength,
				     /*exactp*/Intlist_head(x) == 0 ? true : false,plusp,genestrand,first_read_p,
				     trim_left_p,trim_right_p,outofbounds_start,outofbounds_end,
				     /*minlength*/0)) == NULL) {
	/* Don't know how to fix the junctions */
	debug0(printf("Don't know how to fix the junctions, so returning NULL from Stage3end_new_substrings\n"));
	for (p = substrings; p != NULL; p = List_next(p)) {
	  substring = (Substring_T) List_head(p);
	  if (substring == left_ambig) {
	    /* left_ambig freed by calling procedure.  Need to free junction created for left_ambig. */
	    junctions = List_pop(junctions,(void **) &junction);
	    Junction_free(&junction);
	  } else {
	    Substring_free(&substring);
	  }
	}
	List_free(&substrings);
	Junction_gc(&junctions);
	return (T) NULL;
      } else {
	substrings = List_push(substrings,(void *) substring);
	nmismatches_bothdiff += Substring_nmismatches_bothdiff(substring);
	querylength_trimmed += Substring_querylength(substring);
      }

      /* Prepare for next iteration */
      querystart = queryend;
      if (j != NULL) {
	junction = (Junction_T) List_head(j);
	if ((adj0 = Junction_adj(junction)) != 0) {
	  adj += adj0;
	  indel_score += indel_penalty_middle;
	  nindels += Junction_nindels(junction);
	  if (adj0 < 0) {
	    querystart -= adj0;	/* Insertion */
	  }
	}
      }
      trim_left_p = false;
    }

  } else {
    j = junctions;		/* Put here before we handle left_ambig */
    if (left_ambig != NULL) {
      substrings = List_push(substrings,(void *) left_ambig);
      junctions = List_push(junctions,(void *) Junction_new_splice(/*distance*/0,sensedir,
								   Substring_amb_donor_prob(left_ambig),
								   Substring_amb_acceptor_prob(left_ambig)));
    } else {
      trim_right_p = true;
    }

    /* Subtract querypos to get alignstart/alignend */
    for (q = lefts, x = nmismatches_list, r = Intlist_next(endpoints); q != NULL;
#ifdef LARGE_GENOMES
	 q = Uint8list_next(q),
#else
	 q = Uintlist_next(q),
#endif
	   x = Intlist_next(x), r = Intlist_next(r), j = List_next(j)) {
      queryend = Intlist_head(r);
#ifdef LARGE_GENOMES
      left = Uint8list_head(q);
#else
      left = Uintlist_head(q);
#endif
      debug0(printf("Working on querystart %d..queryend %d at left %u\n",querystart,queryend,left));

      genomicend = left;
      genomicstart = left + querylength;
      genomicend_adj = genomicend - adj;
      genomicstart_adj = genomicend - adj;

      alignstart = genomicstart - (querylength - queryend);
      alignend = genomicstart - (querylength - querystart);

      if (genomicend < chroffset && genomicstart > chrhigh) {
	/* Out of bounds on both sides */
	Junction_gc(&junctions);
	return (T) NULL;

      } else if (genomicend < chroffset) {
	outofbounds_end = chroffset - genomicend;
	outofbounds_start = genomicstart - chroffset;
	debug0(printf("Out of bounds left (high) %d, out of bounds right (low) %d\n",outofbounds_start,outofbounds_end));
	Junction_gc(&junctions);
	return (T) NULL;
#if 0
	/* Could consider this for the lowest substring */
	if (outofbounds_end > outofbounds_start) {
	  /* Consider high part to be out of bounds and keep existing chromosome */
	  outofbounds_end = 0;
	} else {
	  /* Consider low part to be out of bounds and stay in this chromosome */
	  /* Univ_IIT_interval_bounds(&chroffset,&chrhigh,&chrlength,chromosome_iit,chrnum,circular_typeint); */
	  outofbounds_start = 0;
	}
#endif

      } else if (genomicstart > chrhigh) {
	outofbounds_end = chrhigh - genomicend;
	outofbounds_start = genomicstart - chrhigh;
	debug0(printf("Out of bounds left (high) %d, out of bounds right (low) %d\n",outofbounds_start,outofbounds_end));
	Junction_gc(&junctions);
	return (T) NULL;
#if 0
	/* Could consider this for the highest substring */
	if (outofbounds_end > outofbounds_start) {
	  /* Consider high part to be out of bounds and keep existing chromosome */
	  outofbounds_end = 0;
	} else if (++chrnum > nchromosomes) {
	  debug0(printf("Returning NULL from Stage3end_new_substrings\n"));
	  Junction_gc(&junctions);
	  return (T) NULL;
	} else {
	  /* Consider low part to be out of bounds and move to next chromosome */
	  Univ_IIT_interval_bounds(&chroffset,&chrhigh,&chrlength,chromosome_iit,chrnum,circular_typeint);
	  outofbounds_start = 0;
	}
#endif
      }

      if ((nmismatches = Intlist_head(x)) < 0) {
	nmismatches = Genome_count_mismatches_substring(query_compress,left,/*pos5*/querystart,/*pos3*/queryend,
							/*plusp*/false,genestrand,first_read_p);
	debug0(printf("nmismatches %d from genome\n",nmismatches));
      }
      nmismatches_whole += nmismatches;
      debug0(printf("nmismatches %d from sarray\n",nmismatches));
#ifdef LARGE_GENOMES
      if (Uint8list_next(q) == NULL && right_ambig == NULL) {
	trim_left_p = true;
      }
#else
      if (Uintlist_next(q) == NULL && right_ambig == NULL) {
	trim_left_p = true;
      }
#endif
      if ((substring = Substring_new(/*nmismatches_whole*/nmismatches,chrnum,chroffset,chrhigh,chrlength,
				     query_compress,/*start_endtype*/END,/*end_endtype*/END,
				     /*querystart*/querylength - queryend,/*queryend*/querylength - querystart,querylength,
				     alignstart,alignend,/*genomiclength*/querylength,
				     /*exactp*/Intlist_head(x) == 0 ? true : false,plusp,genestrand,first_read_p,
				     trim_left_p,trim_right_p,outofbounds_start,outofbounds_end,
				     /*minlength*/0)) == NULL) {
	/* Don't know how to fix the junctions */
	debug0(printf("Don't know how to fix the junctions, so returning NULL from Stage3end_new_substrings\n"));
	for (p = substrings; p != NULL; p = List_next(p)) {
	  substring = (Substring_T) List_head(p);
	  if (substring == left_ambig) {
	    /* left_ambig freed by calling procedure.  Need to free junction created for left_ambig. */
	    junctions = List_pop(junctions,(void **) &junction);
	    Junction_free(&junction);
	  } else {
	    Substring_free(&substring);
	  }
	}
	List_free(&substrings);
	Junction_gc(&junctions);
	return (T) NULL;
      } else {
	substrings = List_push(substrings,(void *) substring);
	nmismatches_bothdiff += Substring_nmismatches_bothdiff(substring);
	querylength_trimmed += Substring_querylength(substring);
      }

      /* Prepare for next iteration */
      querystart = queryend;
      if (j != NULL) {
	junction = (Junction_T) List_head(j);
	if ((adj0 = Junction_adj(junction)) != 0) {
	  adj += adj0;
	  indel_score += indel_penalty_middle;
	  nindels += Junction_nindels(junction);
	  if (adj0 < 0) {
	    querystart -= adj0;	/* Insertion */
	  }
	}
      }
      trim_right_p = false;
    }
  }

  if (right_ambig != NULL) {
    substrings = List_push(substrings,(void *) right_ambig);
    junctions = List_reverse(junctions);
    junctions = List_push(junctions,(void *) Junction_new_splice(/*distance*/0,sensedir,
								 Substring_amb_donor_prob(right_ambig),
								 Substring_amb_acceptor_prob(right_ambig)));
    junctions = List_reverse(junctions);
  }

#ifdef DEBUG0
  printf("NEW JUNCTIONS\n");
  for (p = junctions; p != NULL; p = List_next(p)) {
    Junction_print(List_head(p));
  }
  printf("\n");
#endif

  new = (T) MALLOC(sizeof(*new));
  new->hittype = SUBSTRINGS;

  new->pairarray = (struct Pair_T *) NULL;
  new->cigar_tokens = (List_T) NULL;
  new->gmap_intronp = false;

  new->querylength = querylength;
  new->querylength_adj = querylength + adj;

  new->substrings_LtoH = substrings;
  new->substrings_1toN = List_copy(substrings); /* Takes over as primary holder of substrings */
  new->substrings_Nto1 = List_copy(substrings);

  new->junctions_LtoH = junctions;
  new->junctions_1toN = List_copy(junctions); /* Takes over as primary holder of substrings */
  new->junctions_Nto1 = List_copy(junctions);

  /* Note differences between substrings and junctions.  Substrings
     were pushed onto lists above, and junctions were created by the
     caller, so they are originally in opposite orders */

#if 0
  if (plusp == true) {
    new->substrings_LtoH = List_reverse(new->substrings_LtoH);
    new->substrings_1toN = List_reverse(new->substrings_1toN);
    new->junctions_Nto1 = List_reverse(new->junctions_Nto1);
  } else {
    new->junctions_LtoH = List_reverse(new->junctions_LtoH);
    new->substrings_Nto1 = List_reverse(new->substrings_Nto1);
    new->junctions_1toN = List_reverse(new->junctions_1toN);
  }
#else
  /* Correct for both plus and minus */
  new->substrings_LtoH = List_reverse(new->substrings_LtoH);
  if (plusp == true) {
    new->substrings_1toN = List_reverse(new->substrings_1toN);
    new->junctions_Nto1 = List_reverse(new->junctions_Nto1);
  } else {
    new->substrings_Nto1 = List_reverse(new->substrings_Nto1);
    new->junctions_1toN = List_reverse(new->junctions_1toN);
  }
#endif

#ifdef DEBUG0
  printf("NEW SUBSTRINGS\n");
  for (p = new->substrings_1toN; p != NULL; p = List_next(p)) {
    substring = List_head(p);
    printf("%d..%d\n",Substring_querystart(substring),Substring_queryend(substring));
  }
  printf("\n");
#endif
  

  substring1 = (Substring_T) List_head(new->substrings_1toN);
  substringN = (Substring_T) List_head(new->substrings_Nto1);

  genomicstart = Substring_genomicstart(substring1);
  genomicend = Substring_genomicend(substringN); /* DOESN'T WORK FOR AMBIGUOUS */
  new->genomicstart = genomicstart;
  new->genomicend = genomicend;

  if (genomicstart < genomicend) {
    new->low = genomicstart;
    new->high = genomicend;
  } else {
    new->low = genomicend;
    new->high = genomicstart;
  }
  new->genomiclength = new->high - new->low;
  new->guided_insertlength = 0U;

  new->genestrand = genestrand;
  new->sarrayp = sarrayp;
  new->gmap_source = GMAP_NOT_APPLICABLE;
  new->improved_by_gmap_p = false;

  new->chrnum = new->effective_chrnum = chrnum;
  new->other_chrnum = 0;
  new->chroffset = chroffset;
  new->chrhigh = chrhigh;
  new->chrlength = chrlength;
  new->plusp = plusp;
  new->sensedir = sensedir;

  new->nindels = nindels;
  new->nmismatches_whole = nmismatches_whole;
  new->nmismatches_bothdiff = nmismatches_bothdiff; /* Trimmed */
  /* new->nmismatches_refdiff = 0; */
  new->ntscore = nmismatches_whole + indel_score;
  new->score = nmismatches_whole + indel_score; /* Want untrimmed */
  new->nsegments = List_length(new->substrings_1toN);

  new->nmatches = querylength - nmismatches_whole;
  new->nmatches_posttrim = querylength_trimmed - nmismatches_whole;

  new->trim_left = Substring_trim_left(substring1);
  new->trim_right = Substring_trim_right(substringN);
  new->trim_left_splicep = Substring_trim_left_splicep(substring1);
  new->trim_right_splicep = Substring_trim_right_splicep(substringN);
  debug0(printf("substrings trim_left %d, trim_right %d\n",new->trim_left,new->trim_right));

  /* new->penalties = 0; */

  /* new->gene_overlap = NO_KNOWN_GENE; -- initialized later when resolving multimappers */
  new->tally = -1L;
  *found_score = new->score;

  new->nsplices = 0;
  for (p = junctions; p != NULL; p = List_next(p)) {
    junction = (Junction_T) List_head(p);
    if (Junction_type(junction) == SPLICE_JUNCTION) {
      new->nsplices += 1;
    }
  }

  new->distance = 0U;
  new->shortexonA_distance = new->shortexonD_distance = 0U;

  new->paired_usedp = false;
  new->paired_seenp = false;
  new->concordantp = false;

  new->circularpos = compute_circularpos(&new->alias,new);
  if (new->alias == +2 || new->alias == -2) {
    debug0(printf("Returning NULL from Stage3end_new_substrings because of alias of %d\n",new->alias));
    Stage3end_free(&new);
    /* Junction_gc(&junctions); -- Done by Stage3end_free */
    return (T) NULL;
  } else {
    debug0(printf("Returning %p from Stage3end_new_substrings with found_score %d\n",new,*found_score));
    return new;
  }
}


#define add_bounded(x,plusterm,highbound) ((x + (plusterm) >= highbound) ? (highbound - 1) : x + (plusterm))
#define subtract_bounded(x,minusterm,lowbound) ((x < lowbound + (minusterm)) ? lowbound : x - (minusterm))


/* Modified from run_gmap_plus in sarray-read.c */
T
Stage3end_substrings_run_gmap_plus (T this, char *queryuc_ptr, int querylength,
				    int genestrand, bool first_read_p,
				    int maxpeelback, Pairpool_T pairpool, Dynprog_T dynprogL, Dynprog_T dynprogM, Dynprog_T dynprogR,
				    Oligoindex_array_T oligoindices_minor, Diagpool_T diagpool, Cellpool_T cellpool) {
  T hit;
  List_T stage2pairs, all_stage2_starts, all_stage2_ends;
  List_T p, startp;

  int k, i;
  int querystart, queryend;
  Univcoord_T *ambcoords;

  int sensedir;

  struct Pair_T *pairarray;
  List_T pairs;
  Substring_T substring, first_ambig, last_ambig;
  int querypos, seglength;
  Chrpos_T genomepos;
  char c, g, g_alt, comp;

  int npairs, goodness, cdna_direction, matches, nmatches_posttrim,
    max_match_length, ambig_end_length_5, ambig_end_length_3,
    unknowns, mismatches, qopens, qindels, topens, tindels,
    ncanonical, nsemicanonical, nnoncanonical;
  double ambig_prob_5, ambig_prob_3, min_splice_prob;
  Splicetype_T ambig_splicetype_5, ambig_splicetype_3;
  Univcoord_T knownsplice_limit_low, knownsplice_limit_high;
  Univcoord_T start, end;
  int nsegments, nmismatches_whole, nindels, nintrons, nindelbreaks;
  char *gsequence_orig, *gsequence_alt;


  debug13(printf("Entered Stage3hr_substrings_run_gmap_plus\n"));

#ifdef HAVE_ALLOCA
    gsequence_orig = (char *) MALLOCA((querylength+1) * sizeof(char));
    gsequence_alt = (char *) MALLOCA((querylength+1) * sizeof(char));
#else
    gsequence_orig = (char *) MALLOC((querylength+1) * sizeof(char));
    gsequence_alt = (char *) MALLOC((querylength+1) * sizeof(char));
#endif

  startp = this->substrings_1toN;
  if (Substring_ambiguous_p((Substring_T) List_head(startp)) == true) {
    first_ambig = (Substring_T) List_head(startp);
    startp = List_next(startp);
  } else {
    first_ambig = (Substring_T) NULL;
  }

  p = this->substrings_Nto1;
  if (Substring_ambiguous_p((Substring_T) List_head(p)) == true) {
    last_ambig = (Substring_T) List_head(p);
  } else {
    last_ambig = (Substring_T) NULL;
  }


  /* D.  Make all_stage2_starts (paths) */
  all_stage2_starts = (List_T) NULL;
  if (first_ambig != NULL) {
    debug13(printf("Handling first ambig\n"));
    querystart = Substring_querystart(first_ambig);
    queryend = Substring_queryend(first_ambig);
    seglength = queryend - querystart;
    ambcoords = Substring_ambcoords(first_ambig);

    for (k = 0; k < Substring_nambcoords(first_ambig); k++) {
      debug13(printf("START, PLUS\n"));
      stage2pairs = (List_T) NULL;
      querypos = querystart;	/* Should be 0 */
      genomepos = (ambcoords[k] - seglength) - this->chroffset;
      Genome_get_segment_blocks_left(gsequence_orig,gsequence_alt,/*right*/ambcoords[k] /*- seglength*/,
				     seglength,this->chroffset,/*revcomp*/false);
      for (i = 0; i < seglength; i++) {
	c = queryuc_ptr[querypos];
	g = gsequence_orig[i];
	g_alt = gsequence_alt[i];
	if (g == c || g_alt == c) {
	  comp = MATCH_COMP;
	} else {
	  comp = MISMATCH_COMP;
	}
	debug13(printf("Pushing %c %c %c at %d,%u\n",c,comp,g,querypos,genomepos));
	stage2pairs = Pairpool_push(stage2pairs,pairpool,querypos,genomepos,
				    /*cdna*/c,comp,/*genome*/g,/*genomealt*/g_alt,
				    /*dynprogindex*/0);
	querypos++;
	genomepos++;
      }
      debug13(Pair_dump_list(stage2pairs,true));
      all_stage2_starts = List_push(all_stage2_starts,(void *) stage2pairs);
    }
  }


  /* E.  Make all_stage2_ends (pairs) */
  all_stage2_ends = (List_T) NULL;
  if (last_ambig != NULL) {
    debug13(printf("Handling last ambig\n"));
    querystart = Substring_querystart(last_ambig);
    queryend = Substring_queryend(last_ambig);
    seglength = queryend - querystart;
    ambcoords = Substring_ambcoords(last_ambig);

    for (k = 0; k < Substring_nambcoords(last_ambig); k++) {
      debug13(printf("END, PLUS\n"));
      stage2pairs = (List_T) NULL;
      querypos = querystart;
      genomepos = ambcoords[k] - this->chroffset;
      Genome_get_segment_blocks_right(gsequence_orig,gsequence_alt,/*left*/ambcoords[k],
                                      seglength,this->chrhigh,/*revcomp*/false);

      for (i = 0; i < seglength; i++) {
	c = queryuc_ptr[querypos];
	g = gsequence_orig[i];
	g_alt = gsequence_alt[i];
	if (g == c || g_alt == c) {
	  comp = MATCH_COMP;
	} else {
	  comp = MISMATCH_COMP;
	}
	debug13(printf("Pushing %c %c %c at %d,%u\n",c,comp,g,querypos,genomepos));
	stage2pairs = Pairpool_push(stage2pairs,pairpool,querypos,genomepos,
				    /*cdna*/c,comp,/*genome*/g,/*genomealt*/g_alt,
				    /*dynprogindex*/0);
	querypos++;
	genomepos++;
      }
      debug13(Pair_dump_list(stage2pairs,true));
      all_stage2_ends = List_push(all_stage2_ends,(void *) List_reverse(stage2pairs));
    }
  }

  /* F.  Make stage2pairs */
  stage2pairs = (List_T) NULL;
  for (p = startp; p != NULL; p = List_next(p)) {
    substring = (Substring_T) List_head(p);
    if (Substring_ambiguous_p(substring) == true) {
      /* Skip */
    } else {
      debug13(printf("Handling substring for %d..%d, %u..%u, chrlength %u\n",
		     Substring_querystart(substring),Substring_queryend(substring),
		     Substring_alignstart_trim_chr(substring),Substring_alignend_trim_chr(substring),
		     this->chrlength));
      querypos = Substring_querystart(substring);
      seglength = Substring_queryend(substring) - querypos;

      genomepos = Substring_alignstart_trim_chr(substring);
      Genome_get_segment_blocks_right(gsequence_orig,gsequence_alt,/*left*/Substring_alignstart_trim(substring),
                                      seglength,this->chrhigh,/*revcomp*/false);

      for (i = 0; i < seglength; i++) {
	c = queryuc_ptr[querypos];
	g = gsequence_orig[i];
	g_alt = gsequence_alt[i];
	if (g == c || g_alt == c) {
	  comp = MATCH_COMP;
	} else {
	  comp = MISMATCH_COMP;
	}
	debug13(printf("Pushing %c %c %c at %d,%u\n",c,comp,g,querypos,genomepos));
	stage2pairs = Pairpool_push(stage2pairs,pairpool,querypos,genomepos,
				    /*cdna*/c,comp,/*genome*/g,/*genomealt*/g_alt,
				    /*dynprogindex*/0);
	querypos++;
	genomepos++;
      }
      debug13(Pair_dump_list(stage2pairs,true));
      debug13(printf("\n"));
    }
  }

  if (stage2pairs == NULL) {
    hit = (T) NULL;
  } else {
    knownsplice_limit_high = ((Pair_T) stage2pairs->first)->genomepos + this->chroffset;
    stage2pairs = List_reverse(stage2pairs);
    knownsplice_limit_low = ((Pair_T) stage2pairs->first)->genomepos + this->chroffset;

    if ((pairarray = Stage3_compute(&pairs,&npairs,&goodness,&cdna_direction,&sensedir,
				    &matches,&nmatches_posttrim,&max_match_length,
				    &ambig_end_length_5,&ambig_end_length_3,
				    &ambig_splicetype_5,&ambig_splicetype_3,
				    &ambig_prob_5,&ambig_prob_3,
				    &unknowns,&mismatches,&qopens,&qindels,&topens,&tindels,
				    &ncanonical,&nsemicanonical,&nnoncanonical,&min_splice_prob,
				    stage2pairs,all_stage2_starts,all_stage2_ends,
#ifdef END_KNOWNSPLICING_SHORTCUT
				    cutoff_level,/*queryptr*/watsonp ? queryuc_ptr : queryrc,
				    watsonp ? query_compress_fwd : query_compress_rev,
#endif
				    /*queryseq_ptr*/queryuc_ptr,queryuc_ptr,querylength,/*skiplength*/0,
#ifdef EXTRACT_GENOMICSEG
				    /*query_subseq_offset*/0,
#else
				    /*query_subseq_offset*/0,
#endif
				    this->chrnum,this->chroffset,this->chrhigh,
				    knownsplice_limit_low,knownsplice_limit_high,/*plusp*/true,genestrand,
				    /*jump_late_p*/false,maxpeelback,pairpool,dynprogL,dynprogM,dynprogR,
				    /*sense_try*/0,/*sense_filter*/0,
				    oligoindices_minor,diagpool,cellpool)) == NULL) {
      hit = (T) NULL;

    } else {
      nsegments = Pair_gsnap_nsegments(&nmismatches_whole,&nindels,&nintrons,&nindelbreaks,
				       pairarray,npairs);
      start = subtract_bounded(this->chroffset + Pair_genomepos(&(pairarray[0])),
			       /*minusterm*/Pair_querypos(&(pairarray[0])),this->chroffset);
      end = add_bounded(this->chroffset + Pair_genomepos(&(pairarray[npairs-1])),
			/*plusterm*/querylength - 1 - Pair_querypos(&(pairarray[npairs-1])),this->chrhigh);

      if ((hit = Stage3end_new_gmap(nmismatches_whole,nmatches_posttrim,max_match_length,
				    ambig_end_length_5,ambig_end_length_3,
				    ambig_splicetype_5,ambig_splicetype_3,
				    ambig_prob_5,ambig_prob_3,min_splice_prob,
				    pairarray,npairs,nsegments,nintrons,nindelbreaks,
				    /*left*/start,/*genomiclength*/end - start + 1,
				    /*plusp*/true,genestrand,first_read_p,
				    /*accession*/NULL,querylength,this->chrnum,this->chroffset,this->chrhigh,this->chrlength,
				    cdna_direction,sensedir,/*gmap_source*/GMAP_VIA_SUBSTRINGS)) == NULL) {
	FREE_OUT(pairarray);
      }
    }
  }

  List_free(&all_stage2_ends);
  List_free(&all_stage2_starts);

#ifdef HAVE_ALLOCA
  FREEA(gsequence_alt);
  FREEA(gsequence_orig);
#else
  FREE(gsequence_alt);
  FREE(gsequence_orig);
#endif

  return hit;
}


/* Modified from run_gmap_minus in sarray-read.c */
T
Stage3end_substrings_run_gmap_minus (T this, char *queryuc_ptr, int querylength,
				     int genestrand, bool first_read_p,
				     int maxpeelback, Pairpool_T pairpool, Dynprog_T dynprogL, Dynprog_T dynprogM, Dynprog_T dynprogR,
				     Oligoindex_array_T oligoindices_minor, Diagpool_T diagpool, Cellpool_T cellpool) {
  T hit;
  List_T stage2pairs, all_stage2_starts, all_stage2_ends;
  List_T p, startp;

  int k, i;
  int querystart, queryend;
  Univcoord_T *ambcoords;

  int sensedir;

  struct Pair_T *pairarray;
  List_T pairs;
  Substring_T substring, first_ambig, last_ambig;
  int querypos, seglength;
  Chrpos_T genomepos;
  char c, g, g_alt, comp;

  int npairs, goodness, cdna_direction, matches, nmatches_posttrim,
    max_match_length, ambig_end_length_5, ambig_end_length_3,
    unknowns, mismatches, qopens, qindels, topens, tindels,
    ncanonical, nsemicanonical, nnoncanonical;
  double ambig_prob_5, ambig_prob_3, min_splice_prob;
  Splicetype_T ambig_splicetype_5, ambig_splicetype_3;
  Univcoord_T knownsplice_limit_low, knownsplice_limit_high;
  Univcoord_T start, end;
  int nsegments, nmismatches_whole, nindels, nintrons, nindelbreaks;

  char *gsequence_orig, *gsequence_alt;

  debug13(printf("Entered Stage3hr_substrings_run_gmap_minus\n"));

#ifdef HAVE_ALLOCA
    gsequence_orig = (char *) MALLOCA((querylength+1) * sizeof(char));
    gsequence_alt = (char *) MALLOCA((querylength+1) * sizeof(char));
#else
    gsequence_orig = (char *) MALLOC((querylength+1) * sizeof(char));
    gsequence_alt = (char *) MALLOC((querylength+1) * sizeof(char));
#endif

  startp = this->substrings_1toN;
  if (Substring_ambiguous_p((Substring_T) List_head(startp)) == true) {
    first_ambig = (Substring_T) List_head(startp);
    startp = List_next(startp);
  } else {
    first_ambig = (Substring_T) NULL;
  }

  p = this->substrings_Nto1;
  if (Substring_ambiguous_p((Substring_T) List_head(p)) == true) {
    last_ambig = (Substring_T) List_head(p);
  } else {
    last_ambig = (Substring_T) NULL;
  }


  /* D.  Make all_stage2_starts (paths) */
  all_stage2_starts = (List_T) NULL;
  if (first_ambig != NULL) {
    debug13(printf("Handling first ambig\n"));
    querystart = Substring_querystart(first_ambig);
    queryend = Substring_queryend(first_ambig);
    seglength = queryend - querystart;
    ambcoords = Substring_ambcoords(first_ambig);

    for (k = 0; k < Substring_nambcoords(first_ambig); k++) {
      debug13(printf("START, MINUS\n"));
      stage2pairs = (List_T) NULL;
      querypos = querystart;
      genomepos = (this->chrhigh + 1) - ambcoords[k] - seglength;
      Genome_get_segment_blocks_right(gsequence_orig,gsequence_alt,/*left*/ambcoords[k],
                                      seglength,this->chrhigh,/*revcomp*/true);

      for (i = 0; i < seglength; i++) {
	c = queryuc_ptr[querypos];
	g = gsequence_orig[i];
	g_alt = gsequence_alt[i];
	if (g == c || g_alt == c) {
	  comp = MATCH_COMP;
	} else {
	  comp = MISMATCH_COMP;
	}
	debug13(printf("Pushing %c %c %c at %d,%u\n",c,comp,g,querypos,genomepos));
	stage2pairs = Pairpool_push(stage2pairs,pairpool,querypos,genomepos,
				    /*cdna*/c,comp,/*genome*/g,/*genomealt*/g_alt,
				    /*dynprogindex*/0);
	querypos++;
	genomepos++;
      }
      debug13(Pair_dump_list(stage2pairs,true));
      all_stage2_starts = List_push(all_stage2_starts,(void *) stage2pairs);
    }
  }


  /* E.  Make all_stage2_ends (pairs) */
  all_stage2_ends = (List_T) NULL;
  if (last_ambig != NULL) {
    debug13(printf("Handling last ambig\n"));
    querystart = Substring_querystart(last_ambig);
    queryend = Substring_queryend(last_ambig);
    seglength = queryend - querystart;
    ambcoords = Substring_ambcoords(last_ambig);

    for (k = 0; k < Substring_nambcoords(last_ambig); k++) {
      debug13(printf("END, MINUS\n"));
      stage2pairs = (List_T) NULL;
      querypos = querystart;
      genomepos = (this->chrhigh + 1) - ambcoords[k];
      Genome_get_segment_blocks_left(gsequence_orig,gsequence_alt,/*right*/ambcoords[k] /*- seglength*/,
				     seglength,this->chroffset,/*revcomp*/true);

      for (i = 0; i < seglength; i++) {
	c = queryuc_ptr[querypos];
	g = gsequence_orig[i];
	g_alt = gsequence_alt[i];
	if (g == c || g_alt == c) {
	  comp = MATCH_COMP;
	} else {
	  comp = MISMATCH_COMP;
	}
	debug13(printf("Pushing %c %c %c at %d,%u\n",c,comp,g,querypos,genomepos));
	stage2pairs = Pairpool_push(stage2pairs,pairpool,querypos,genomepos,
				    /*cdna*/c,comp,/*genome*/g,/*genomealt*/g_alt,
				    /*dynprogindex*/0);
	querypos++;
	genomepos++;
      }
      debug13(Pair_dump_list(stage2pairs,true));
      all_stage2_ends = List_push(all_stage2_ends,(void *) List_reverse(stage2pairs));
    }
  }

  /* F.  Make stage2pairs */
  stage2pairs = (List_T) NULL;
  for (p = startp; p != NULL; p = List_next(p)) {
    substring = (Substring_T) List_head(p);
    if (Substring_ambiguous_p(substring) == true) {
      /* Skip */
    } else {
      debug13(printf("Handling substring for %d..%d, %u..%u, chrlength %u\n",
		     Substring_querystart(substring),Substring_queryend(substring),
		     Substring_alignstart_trim_chr(substring),Substring_alignend_trim_chr(substring),
		     this->chrlength));
      querypos = Substring_querystart(substring);
      seglength = Substring_queryend(substring) - querypos;

      /* Don't understand why it isn't this->chrhigh - 1, but it
	 looks like the minus coordinates for substrings are +1 higher
	 than they should be */
      genomepos = (this->chrhigh + 1) - Substring_alignstart_trim(substring);
      Genome_get_segment_blocks_right(gsequence_orig,gsequence_alt,/*left*/Substring_alignend_trim(substring),
                                      seglength,this->chrhigh,/*revcomp*/true);

      for (i = 0; i < seglength; i++) {
	c = queryuc_ptr[querypos];
	g = gsequence_orig[i];
	g_alt = gsequence_alt[i];
	if (g == c || g_alt == c) {
	  comp = MATCH_COMP;
	} else {
	  comp = MISMATCH_COMP;
	}
	debug13(printf("Pushing %c %c %c at %d,%u\n",c,comp,g,querypos,genomepos));
	stage2pairs = Pairpool_push(stage2pairs,pairpool,querypos,genomepos,
				    /*cdna*/c,comp,/*genome*/g,/*genomealt*/g_alt,
				    /*dynprogindex*/0);
	querypos++;
	genomepos++;
      }
      debug13(Pair_dump_list(stage2pairs,true));
      debug13(printf("\n"));
    }
  }

  if (stage2pairs == NULL) {
    hit = (T) NULL;
  } else {
    knownsplice_limit_low = ((Pair_T) stage2pairs->first)->genomepos + this->chroffset;
    stage2pairs = List_reverse(stage2pairs);
    knownsplice_limit_high = ((Pair_T) stage2pairs->first)->genomepos + this->chroffset;


    if ((pairarray = Stage3_compute(&pairs,&npairs,&goodness,&cdna_direction,&sensedir,
				    &matches,&nmatches_posttrim,&max_match_length,
				    &ambig_end_length_5,&ambig_end_length_3,
				    &ambig_splicetype_5,&ambig_splicetype_3,
				    &ambig_prob_5,&ambig_prob_3,
				    &unknowns,&mismatches,&qopens,&qindels,&topens,&tindels,
				    &ncanonical,&nsemicanonical,&nnoncanonical,&min_splice_prob,
				    stage2pairs,all_stage2_starts,all_stage2_ends,
#ifdef END_KNOWNSPLICING_SHORTCUT
				    cutoff_level,/*queryptr*/watsonp ? queryuc_ptr : queryrc,
				    watsonp ? query_compress_fwd : query_compress_rev,
#endif
				    /*queryseq_ptr*/queryuc_ptr,queryuc_ptr,querylength,/*skiplength*/0,
#ifdef EXTRACT_GENOMICSEG
				    /*query_subseq_offset*/0,
#else
				    /*query_subseq_offset*/0,
#endif
				    this->chrnum,this->chroffset,this->chrhigh,
				    knownsplice_limit_low,knownsplice_limit_high,/*plusp*/false,genestrand,
				    /*jump_late_p*/true,maxpeelback,pairpool,dynprogL,dynprogM,dynprogR,
				    /*sense_try*/0,/*sense_filter*/0,
				    oligoindices_minor,diagpool,cellpool)) == NULL) {
      hit = (T) NULL;

    } else {
      nsegments = Pair_gsnap_nsegments(&nmismatches_whole,&nindels,&nintrons,&nindelbreaks,
				       pairarray,npairs);
      start = add_bounded(this->chroffset + Pair_genomepos(&(pairarray[0])),
			  /*plusterm*/Pair_querypos(&(pairarray[0])),this->chrhigh);
      end = subtract_bounded(this->chroffset + Pair_genomepos(&(pairarray[npairs-1])),
			     /*minusterm*/querylength - 1 - Pair_querypos(&(pairarray[npairs-1])),this->chroffset);
      if ((hit = Stage3end_new_gmap(nmismatches_whole,nmatches_posttrim,max_match_length,
				    ambig_end_length_5,ambig_end_length_3,
				    ambig_splicetype_5,ambig_splicetype_3,
				    ambig_prob_5,ambig_prob_3,min_splice_prob,
				    pairarray,npairs,nsegments,nintrons,nindelbreaks,
				    /*left*/end,/*genomiclength*/start - end + 1,
				    /*plusp*/false,genestrand,first_read_p,
				    /*accession*/NULL,querylength,this->chrnum,this->chroffset,this->chrhigh,this->chrlength,
				    cdna_direction,sensedir,/*gmap_source*/GMAP_VIA_SUBSTRINGS)) == NULL) {
	FREE_OUT(pairarray);
      }
    }
  }

  List_free(&all_stage2_ends);
  List_free(&all_stage2_starts);

#ifdef HAVE_ALLOCA
  FREEA(gsequence_alt);
  FREEA(gsequence_orig);
#else
  FREE(gsequence_alt);
  FREE(gsequence_orig);
#endif

  return hit;
}




T
Stage3end_new_exact (int *found_score, Univcoord_T left, int genomiclength, Compress_T query_compress,
		     bool plusp, int genestrand, bool first_read_p,
		     Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh,
		     Chrpos_T chrlength, bool sarrayp) {
  T new;
  Substring_T substring;
  Univcoord_T genomicstart, genomicend;
  bool exactp = true;
  int outofbounds_start = 0, outofbounds_end = 0;

  if (plusp == true) {
    genomicstart = left;
    genomicend = left + genomiclength;
    if (genomicstart < chroffset && genomicend > chrhigh) {
      /* Out of bounds on both sides */
      return (T) NULL;
      
    } else if (genomicstart < chroffset) {
      outofbounds_start = chroffset - genomicstart;
      outofbounds_end = genomicend - chroffset;
      debug0(printf("Out of bounds left (low) %d, out of bounds right (high) %d\n",outofbounds_start,outofbounds_end));
      if (outofbounds_start > outofbounds_end) {
	/* Consider high part to be out of bounds and keep existing chromosome */
	outofbounds_start = 0;
      } else {
	/* Consider low part to be out of bounds and stay in this chromosome */
	/* Univ_IIT_interval_bounds(&chroffset,&chrhigh,&chrlength,chromosome_iit,chrnum,circular_typeint); */
	outofbounds_end = 0;
      }
      exactp = false;

    } else if (genomicend > chrhigh) {
      outofbounds_start = chrhigh - genomicstart;
      outofbounds_end = genomicend - chrhigh;
      debug0(printf("Out of bounds left (low) %d, out of bounds right (high) %d\n",outofbounds_start,outofbounds_end));
      if (outofbounds_start > outofbounds_end) {
	/* Consider high part to be out of bounds and keep existing chromosome */
	outofbounds_start = 0;
      } else if (++chrnum > nchromosomes) {
	return (T) NULL;
      } else {
	/* Consider low part to be out of bounds and move to next chromosome */
	Univ_IIT_interval_bounds(&chroffset,&chrhigh,&chrlength,chromosome_iit,chrnum,circular_typeint);
	outofbounds_end = 0;
      }
      exactp = false;
    }

  } else {
    genomicend = left;
    genomicstart = left + genomiclength;
    
    if (genomicend < chroffset && genomicstart > chrhigh) {
      /* Out of bounds on both ends */
      return (T) NULL;

    } else if (genomicend < chroffset) {
      outofbounds_end = chroffset - genomicend;
      outofbounds_start = genomicstart - chroffset;
      debug0(printf("Out of bounds left (high) %d, out of bounds right (low) %d\n",outofbounds_start,outofbounds_end));
      if (outofbounds_end > outofbounds_start) {
	/* Consider high part to be out of bounds and keep existing chromosome */
	outofbounds_end = 0;
      } else {
	/* Consider low part to be out of bounds and stay in this chromosome */
	/* Univ_IIT_interval_bounds(&chroffset,&chrhigh,&chrlength,chromosome_iit,chrnum,circular_typeint); */
	outofbounds_start = 0;
      }
      exactp = false;

    } else if (genomicstart > chrhigh) {
      outofbounds_end = chrhigh - genomicend;
      outofbounds_start = genomicstart - chrhigh;
      debug0(printf("Out of bounds left (high) %d, out of bounds right (low) %d\n",outofbounds_start,outofbounds_end));
      if (outofbounds_end > outofbounds_start) {
	/* Consider high part to be out of bounds and keep existing chromosome */
	outofbounds_end = 0;
      } else if (++chrnum > nchromosomes) {
	return (T) NULL;
      } else {
	/* Consider low part to be out of bounds and move to next chromosome */
	Univ_IIT_interval_bounds(&chroffset,&chrhigh,&chrlength,chromosome_iit,chrnum,circular_typeint);
	outofbounds_start = 0;
      }
      exactp = false;
    }
  }

  if ((substring = Substring_new(/*nmismatches*/0,chrnum,chroffset,chrhigh,chrlength,
				 query_compress,/*start_endtype*/END,/*end_endtype*/END,
				 /*querystart*/0,/*queryend*/genomiclength,/*querylength*/genomiclength,
				 /*alignstart*/genomicstart,/*alignend*/genomicend,
				 genomiclength,exactp,plusp,genestrand,first_read_p,/*trim_left_p*/false,/*trim_right_p*/false,
				 outofbounds_start,outofbounds_end,/*minlength*/0)) == NULL) {
    return (T) NULL;

  } else {
    new = (T) MALLOC_OUT(sizeof(*new));
    debug0(printf("Stage3end_new_exact %p: left %llu, chrnum %d, sarrayp %d\n",new,(unsigned long long) left,chrnum,sarrayp));

    new->substrings_LtoH = List_push(NULL,(void *) substring);
    new->substrings_1toN = List_push(NULL,(void *) substring);
    new->substrings_Nto1 = List_push(NULL,(void *) substring);

    new->junctions_LtoH = (List_T) NULL;
    new->junctions_1toN = (List_T) NULL;
    new->junctions_Nto1 = (List_T) NULL;

    new->pairarray = (struct Pair_T *) NULL;
    new->cigar_tokens = (List_T) NULL;
    new->gmap_intronp = false;

    new->querylength_adj = new->querylength = genomiclength;
    new->genomicstart = genomicstart;
    new->genomicend = genomicend;

    if (genomicstart < genomicend) {
      new->low = genomicstart;
      new->high = genomicend;
    } else {
      new->low = genomicend;
      new->high = genomicstart;
    }
    new->genomiclength = new->high - new->low;
    new->guided_insertlength = 0U;
    debug0(printf("Assigned %llu to low and %llu to high\n",(unsigned long long) new->low,(unsigned long long) new->high));


    if (exactp == true) {
      new->hittype = EXACT;
    } else {
      new->hittype = SUB;
    }
    new->genestrand = genestrand;
    new->sarrayp = sarrayp;
    new->gmap_source = GMAP_NOT_APPLICABLE;
    new->improved_by_gmap_p = false;

    new->chrnum = new->effective_chrnum = chrnum;
    new->other_chrnum = 0;
    new->chroffset = chroffset;
    new->chrhigh = chrhigh;
    new->chrlength = chrlength;
    new->plusp = plusp;
    new->sensedir = SENSE_NULL;

#if 0
    new->mapq_loglik = Substring_mapq_loglik(substring);
    new->mapq_score = 0;
    new->absmq_score = 0;
#endif

    new->nindels = 0;
    new->nmismatches_whole = 0;
    new->nmismatches_bothdiff = 0;
    /* new->nmismatches_refdiff = 0; */
    new->ntscore = 0;
    new->score = 0;
    new->nsegments = 1;
    new->nmatches = genomiclength;
    new->nmatches_posttrim = genomiclength;

    new->trim_left = 0;
    new->trim_right = 0;
    new->trim_left_splicep = false;
    new->trim_right_splicep = false;

    /* new->penalties = 0; */

    /* new->gene_overlap = NO_KNOWN_GENE; -- initialized later when resolving multimappers */
    new->tally = -1L;
    *found_score = 0;

    new->nsplices = 0;

    new->distance = 0U;
    new->shortexonA_distance = new->shortexonD_distance = 0U;

    new->paired_usedp = false;
    new->paired_seenp = false;
    new->concordantp = false;

    new->circularpos = compute_circularpos(&new->alias,new);
    if (new->alias == +2 || new->alias == -2) {
      Stage3end_free(&new);
      return (T) NULL;
    } else {
      return new;
    }
  }
}


T
Stage3end_new_substitution (int *found_score, int nmismatches_whole, Univcoord_T left,
			    int genomiclength, Compress_T query_compress,
			    bool plusp, int genestrand, bool first_read_p,
			    Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh,
			    Chrpos_T chrlength, bool sarrayp) {
  T new;
  Substring_T substring;
  Univcoord_T genomicstart, genomicend;
  int outofbounds_start = 0, outofbounds_end = 0;

  debug0(printf("Entered Stage3end_new_substitution at left %u and chrhigh %u, sarrayp %d\n",left,chrhigh,sarrayp));

  if (plusp == true) {
    genomicstart = left;
    genomicend = left + genomiclength;
    if (genomicstart < chroffset && genomicend > chrhigh) {
      /* Out of bounds on both sides */
      return (T) NULL;
      
    } else if (genomicstart < chroffset) {
      outofbounds_start = chroffset - genomicstart;
      outofbounds_end = genomicend - chroffset;
      debug0(printf("Out of bounds left (low) %d, out of bounds right (high) %d\n",outofbounds_start,outofbounds_end));
      if (outofbounds_start > outofbounds_end) {
	/* Consider high part to be out of bounds and keep existing chromosome */
	outofbounds_start = 0;
      } else {
	/* Consider low part to be out of bounds and stay in this chromosome */
	/* Univ_IIT_interval_bounds(&chroffset,&chrhigh,&chrlength,chromosome_iit,chrnum,circular_typeint); */
	outofbounds_end = 0;
      }

    } else if (genomicend > chrhigh) {
      outofbounds_start = chrhigh - genomicstart;
      outofbounds_end = genomicend - chrhigh;
      debug0(printf("Out of bounds left (low) %d, out of bounds right (high) %d\n",outofbounds_start,outofbounds_end));
      if (outofbounds_start > outofbounds_end) {
	/* Consider high part to be out of bounds and keep existing chromosome */
	outofbounds_start = 0;
      } else if (++chrnum > nchromosomes) {
	return (T) NULL;
      } else {
	/* Consider low part to be out of bounds and move to next chromosome */
	Univ_IIT_interval_bounds(&chroffset,&chrhigh,&chrlength,chromosome_iit,chrnum,circular_typeint);
	outofbounds_end = 0;
      }
    }

  } else {
    genomicend = left;
    genomicstart = left + genomiclength;

    if (genomicend < chroffset && genomicstart > chrhigh) {
      /* Out of bounds on both ends */
      return (T) NULL;

    } else if (genomicend < chroffset) {
      outofbounds_end = chroffset - genomicend;
      outofbounds_start = genomicstart - chroffset;
      debug0(printf("Out of bounds left (high) %d, out of bounds right (low) %d\n",outofbounds_start,outofbounds_end));
      if (outofbounds_end > outofbounds_start) {
	/* Consider high part to be out of bounds and keep existing chromosome */
	outofbounds_end = 0;
      } else {
	/* Consider low part to be out of bounds and stay in this chromosome */
	/* Univ_IIT_interval_bounds(&chroffset,&chrhigh,&chrlength,chromosome_iit,chrnum,circular_typeint); */
	outofbounds_start = 0;
      }

    } else if (genomicstart > chrhigh) {
      outofbounds_end = chrhigh - genomicend;
      outofbounds_start = genomicstart - chrhigh;
      debug0(printf("Out of bounds left (high) %d, out of bounds right (low) %d\n",outofbounds_start,outofbounds_end));
      if (outofbounds_end > outofbounds_start) {
	/* Consider high part to be out of bounds and keep existing chromosome */
	outofbounds_end = 0;
      } else if (++chrnum > nchromosomes) {
	return (T) NULL;
      } else {
	/* Consider low part to be out of bounds and move to next chromosome */
	Univ_IIT_interval_bounds(&chroffset,&chrhigh,&chrlength,chromosome_iit,chrnum,circular_typeint);
	outofbounds_start = 0;
      }
    }
  }

  if ((substring = Substring_new(nmismatches_whole,chrnum,chroffset,chrhigh,chrlength,
				 query_compress,/*start_endtype*/END,/*end_endtype*/END,
				 /*querystart*/0,/*queryend*/genomiclength,/*querylength*/genomiclength,
				 /*alignstart*/genomicstart,/*alignend*/genomicend,
				 genomiclength,/*exactp*/false,plusp,genestrand,first_read_p,/*trim_left_p*/true,/*trim_right_p*/true,
				 outofbounds_start,outofbounds_end,/*minlength*/genomiclength/2)) == NULL) {
    debug0(printf("Returning NULL\n"));
    return (T) NULL;

  } else {
    new = (T) MALLOC_OUT(sizeof(*new));
    debug0(printf("Stage3end_new_substitution %p: left %llu, chrnum %d, nmismatches %d, sarrayp %d\n",
		  new,(unsigned long long) left,chrnum,nmismatches_whole,sarrayp));

    new->substrings_LtoH = List_push(NULL,(void *) substring);
    new->substrings_1toN = List_push(NULL,(void *) substring);
    new->substrings_Nto1 = List_push(NULL,(void *) substring);

    new->junctions_LtoH = (List_T) NULL;
    new->junctions_1toN = (List_T) NULL;
    new->junctions_Nto1 = (List_T) NULL;

    new->pairarray = (struct Pair_T *) NULL;
    new->cigar_tokens = (List_T) NULL;
    new->gmap_intronp = false;

    new->querylength_adj = new->querylength = genomiclength;
    new->genomicstart = genomicstart;
    new->genomicend = genomicend;

    if (genomicstart < genomicend) {
      new->low = genomicstart;
      new->high = genomicend;
    } else {
      new->low = genomicend;
      new->high = genomicstart;
    }
    new->genomiclength = new->high - new->low;
    new->guided_insertlength = 0U;

    if (nmismatches_whole == 0) {
      /* Proper hittype needed so we can eliminate identical hits */
      new->hittype = EXACT;
    } else {
      new->hittype = SUB;
    }
    new->genestrand = genestrand;
    new->sarrayp = sarrayp;
    new->gmap_source = GMAP_NOT_APPLICABLE;
    new->improved_by_gmap_p = false;

    new->chrnum = new->effective_chrnum = chrnum;
    new->other_chrnum = 0;
    new->chroffset = chroffset;
    new->chrhigh = chrhigh;
    new->chrlength = chrlength;
    new->plusp = plusp;
    new->sensedir = SENSE_NULL;

#if 0
    new->mapq_loglik = Substring_mapq_loglik(substring);
    new->mapq_score = 0;
    new->absmq_score = 0;
#endif

    new->nindels = 0;
    new->nmismatches_whole = nmismatches_whole;
    new->ntscore = nmismatches_whole;
    new->score = nmismatches_whole;
    new->nsegments = 1;

    new->nmismatches_bothdiff = Substring_nmismatches_bothdiff(substring);
    /* new->nmismatches_refdiff = Substring_nmismatches_refdiff(new->substring1); */

#if 0
    /* This method was previously the only one correct for SNP-tolerant alignment */
    new->nmatches = Substring_match_length(new->substring1) - new->total_nmismatches;
#else
    /* This method is now correct for SNP-tolerant alignment */
    new->nmatches = Substring_nmatches(substring);
    new->nmatches_posttrim = Substring_nmatches_posttrim(substring);
#endif

    new->trim_left = Substring_trim_left(substring);
    new->trim_right = Substring_trim_right(substring);
    new->trim_left_splicep = Substring_trim_left_splicep(substring);
    new->trim_right_splicep = Substring_trim_right_splicep(substring);

    /* new->penalties = 0; */

    /* new->gene_overlap = NO_KNOWN_GENE; -- initialized later when resolving multimappers */
    new->tally = -1L;

    if (new->score < *found_score) {
      *found_score = new->score;
    }

    new->nsplices = 0;

    new->distance = 0U;
    new->shortexonA_distance = new->shortexonD_distance = 0U;

    new->paired_usedp = false;
    new->paired_seenp = false;
    new->concordantp = false;

    new->circularpos = compute_circularpos(&new->alias,new);
    if (new->alias == +2 || new->alias == -2) {
      Stage3end_free(&new);
      return (T) NULL;
    } else {
      debug(printf("Returning substitution %p\n",new));
      return new;
    }
  }
}



T
Stage3end_new_insertion (int *found_score, int nindels, int indel_pos, int nmismatches1_whole, int nmismatches2_whole,
			 Univcoord_T left, int genomiclength, Compress_T query_compress,
			 int querylength, bool plusp, int genestrand, bool first_read_p,
			 Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh, Chrpos_T chrlength,
			 int indel_penalty, bool sarrayp) {
  T new;
  Substring_T substring1, substring2;
  Junction_T junction;
  int querystart1, queryend1, querystart2, queryend2;
  Univcoord_T genomicstart1, genomicend1, genomicstart2, genomicend2;
  Univcoord_T genomicstart_adj_2, genomicend_adj_2;
  Univcoord_T alignstart1, alignend1, alignstart2, alignend2;
  int outofbounds_start = 0, outofbounds_end = 0;


  debug2(printf("Entered with left %llu, querylength %d, genomiclength %d, indel_pos %d\n",
		(unsigned long long) left,querylength,genomiclength,indel_pos));
#if 0
  debug2(printf("q: %s\n",query));
  debug2(printf("g: %s\n",genomicseg));
#endif

  assert(nindels > 0);

  querystart1 = 0;
  queryend1 = indel_pos;
  querystart2 = indel_pos + nindels;
  queryend2 = querylength;

  if (plusp == true) {
    alignstart1 = left /*+ querystart1 (0)*/;
    alignend1 = alignstart2 = left + /*queryend1*/indel_pos;
    alignend2 = (left - nindels) + /*queryend2*/querylength;

    genomicstart1 = alignstart1;
    genomicend1 = alignend1;
    genomicstart2 = alignstart2;
    genomicend2 = alignend2;

    if (genomicstart1 < chroffset && genomicend2 > chrhigh) {
      /* Out of bounds on both sides */
      return (T) NULL;

    } else if (genomicstart1 < chroffset) {
      outofbounds_start = chroffset - genomicstart1;
      outofbounds_end = genomicend2 - chroffset;
      debug0(printf("Out of bounds left (low) %d, out of bounds right (high) %d\n",outofbounds_start,outofbounds_end));
      if (outofbounds_start > outofbounds_end) {
	/* Consider high part to be out of bounds and keep existing chromosome */
	if (querylength - indel_pos - nindels < outofbounds_end) {
	  /* indel is in eliminated part, so abort */
	  return (T) NULL;
	}
	outofbounds_start = 0;
      } else {
	/* Consider low part to be out of bounds and stay in this chromosome */
	if (indel_pos < outofbounds_start) {
	  /* indel is in eliminated part, so abort */
	  return (T) NULL;
	}
	/* Univ_IIT_interval_bounds(&chroffset,&chrhigh,&chrlength,chromosome_iit,chrnum,circular_typeint); */
	outofbounds_end = 0;
      }

    } else if (genomicend2 > chrhigh) {
      outofbounds_start = chrhigh - genomicstart1;
      outofbounds_end = genomicend2 - chrhigh;
      debug0(printf("Out of bounds left (low) %d, out of bounds right (high) %d\n",outofbounds_start,outofbounds_end));
      if (outofbounds_start > outofbounds_end) {
	/* Consider high part to be out of bounds and keep existing chromosome */
	if (querylength - indel_pos - nindels < outofbounds_end) {
	  /* indel is in eliminated part, so abort */
	  return (T) NULL;
	}
	outofbounds_start = 0;
      } else if (++chrnum > nchromosomes) {
	return (T) NULL;
      } else {
	/* Consider low part to be out of bounds and move to next chromosome */
	if (indel_pos < outofbounds_start) {
	  /* indel is in eliminated part, so abort */
	  return (T) NULL;
	}
	Univ_IIT_interval_bounds(&chroffset,&chrhigh,&chrlength,chromosome_iit,chrnum,circular_typeint);
	outofbounds_end = 0;
      }
    }

  } else {
    alignstart1 = (left - nindels) + (querylength /*- querystart (0)*/);
    alignend1 = alignstart2 = (left - nindels) + (querylength - indel_pos);
    alignend2 = left /* + (querylength - queryend)*/;

    genomicstart1 = alignstart1;
    genomicend1 = alignend1;
    genomicstart2 = alignstart2;
    genomicend2 = alignend2;

    if (genomicend2 < chroffset && genomicstart1 > chrhigh) {
      /* Out of bounds on both sides */
      return (T) NULL;

    } else if (genomicend2 < chroffset) {
      outofbounds_end = chroffset - genomicend2;
      outofbounds_start = genomicstart1 - chroffset;
      debug0(printf("Out of bounds left (high) %d, out of bounds right (low) %d\n",outofbounds_start,outofbounds_end));
      if (outofbounds_end > outofbounds_start) {
	/* Consider high part to be out of bounds and keep existing chromosome */
	if (indel_pos < outofbounds_start) {
	  /* indel is in eliminated part, so abort */
	  return (T) NULL;
	}
	outofbounds_end = 0;
      } else {
	/* Consider low part to be out of bounds and stay in this chromosome */
	if (querylength - indel_pos - nindels < outofbounds_end) {
	  /* indel is in eliminated part, so abort */
	  return (T) NULL;
	}
	/* Univ_IIT_interval_bounds(&chroffset,&chrhigh,&chrlength,chromosome_iit,chrnum,circular_typeint); */
	outofbounds_start = 0;
      }

    } else if (genomicstart1 > chrhigh) {
      outofbounds_end = chrhigh - genomicend2;
      outofbounds_start = genomicstart1 - chrhigh;
      debug0(printf("Out of bounds left (high) %d, out of bounds right (low) %d\n",outofbounds_start,outofbounds_end));
      if (outofbounds_end > outofbounds_start) {
	/* Consider high part to be out of bounds and keep existing chromosome */
	if (indel_pos < outofbounds_start) {
	  /* indel is in eliminated part, so abort */
	  return (T) NULL;
	}
	outofbounds_end = 0;
      } else if (++chrnum > nchromosomes) {
	return (T) NULL;
      } else {
	/* Consider low part to be out of bounds and move to next chromosome */
	if (querylength - indel_pos - nindels < outofbounds_end) {
	  /* indel is in eliminated part, so abort */
	  return (T) NULL;
	}
	Univ_IIT_interval_bounds(&chroffset,&chrhigh,&chrlength,chromosome_iit,chrnum,circular_typeint);
	outofbounds_start = 0;
      }
    }
  }

  if ((substring1 = Substring_new(nmismatches1_whole,chrnum,chroffset,chrhigh,chrlength,
				  query_compress,/*start_endtype*/END,/*end_endtype*/INS,
				  querystart1,queryend1,querylength,alignstart1,alignend1,genomiclength,
				  /*exactp*/false,plusp,genestrand,first_read_p,
				  /*trim_left_p (previously was end1_indel_p ? false : true)*/true,
				  /*trim_right_p*/false,outofbounds_start,/*outofbounds_end*/0,/*minlength*/0)) == NULL) {
    return (T) NULL;

  } else if ((substring2 = Substring_new(nmismatches2_whole,chrnum,chroffset,chrhigh,chrlength,
					 query_compress,/*start_endtype*/INS,/*end_endtype*/END,
					 querystart2,queryend2,querylength,alignstart2,alignend2,genomiclength,
					 /*exactp*/false,plusp,genestrand,first_read_p,
					 /*trim_left_p*/false,/*trim_right_p (previously was end2_indel_p ? false : true)*/true,
					 /*outofbounds_start*/0,outofbounds_end,/*minlength*/0)) == NULL) {
    Substring_free(&substring1);
    return (T) NULL;

  } else {
    new = (T) MALLOC_OUT(sizeof(*new));
    debug0(printf("Stage3end_new_insertion %p: left %llu, chrnum %d, nmismatches %d+%d, indel_pos %d, nindels %d, sarrayp %d\n",
		  new,(unsigned long long) left,chrnum,nmismatches1_whole,nmismatches2_whole,indel_pos,nindels,sarrayp));

    new->substrings_1toN = List_push(NULL,substring2);
    new->substrings_1toN = List_push(new->substrings_1toN,substring1);

    new->substrings_Nto1 = List_push(NULL,substring1);
    new->substrings_Nto1 = List_push(new->substrings_Nto1,substring2);

    if (plusp == true) {
      new->substrings_LtoH = List_push(NULL,substring2);
      new->substrings_LtoH = List_push(new->substrings_LtoH,substring1);
    } else {
      new->substrings_LtoH = List_push(NULL,substring1);
      new->substrings_LtoH = List_push(new->substrings_LtoH,substring2);
    }
    junction = Junction_new_insertion(nindels);
    new->junctions_LtoH = List_push(NULL,junction);
    new->junctions_1toN = List_push(NULL,junction);
    new->junctions_Nto1 = List_push(NULL,junction);

    new->pairarray = (struct Pair_T *) NULL;
    new->cigar_tokens = (List_T) NULL;
    new->gmap_intronp = false;

    new->querylength = querylength;
    new->querylength_adj =  querylength - nindels;
    new->genomicstart = genomicstart1;
    new->genomicend = genomicend2;

    if (genomicstart1 < genomicend2) {
      new->low = genomicstart1;
      new->high = genomicend2;
    } else {
      new->low = genomicend2;
      new->high = genomicstart1;
    }
    new->genomiclength = new->high - new->low;
    new->guided_insertlength = 0U;

    new->hittype = INSERTION;
    new->genestrand = genestrand;
    new->sarrayp = sarrayp;
    new->gmap_source = GMAP_NOT_APPLICABLE;
    new->improved_by_gmap_p = false;

    new->chrnum = new->effective_chrnum = chrnum;
    new->other_chrnum = 0;
    new->chroffset = chroffset;
    new->chrhigh = chrhigh;
    new->chrlength = chrlength;
    new->plusp = plusp;
    new->sensedir = SENSE_NULL;

#if 0
    new->mapq_loglik = Substring_mapq_loglik(substring1) + Substring_mapq_loglik(substring2) + 
      MAPQ_loglik_exact(quality_string,queryend1,querystart2);
    new->mapq_score = 0;
    new->absmq_score = 0;
#endif

    new->nindels = nindels;
    new->nmismatches_whole = nmismatches1_whole + nmismatches2_whole;
    new->ntscore = indel_penalty + nmismatches1_whole + nmismatches2_whole;
    new->score = new->ntscore;
    new->nsegments = 2;

    new->nmismatches_bothdiff = Substring_nmismatches_bothdiff(substring1) + Substring_nmismatches_bothdiff(substring2);
    /* new->nmismatches_refdiff = Substring_nmismatches_refdiff(new->substring1) + Substring_nmismatches_refdiff(new->substring2); */

#if 0
    /* This method is correct for SNP-tolerant alignment */
    new->nmatches = Substring_match_length(substring1) + Substring_match_length(substring2) - new->total_nmismatches;
#else
    /* This method is now correct for SNP-tolerant alignment */
    new->nmatches = Substring_nmatches(substring1) + Substring_nmatches(substring2);
    new->nmatches_posttrim = Substring_nmatches_posttrim(substring1) + Substring_nmatches_posttrim(substring2);
    new->nmatches_posttrim += nindels; /* for use in goodness_cmp procedures */
    new->nmatches_posttrim -= indel_penalty; /* for use in goodness_cmp procedures */
#endif

    new->trim_left = Substring_trim_left(substring1);
    new->trim_right = Substring_trim_right(substring2);
    new->trim_left_splicep = Substring_trim_left_splicep(substring1);
    new->trim_right_splicep = Substring_trim_right_splicep(substring2);

#if 0
#ifdef SCORE_INDELS
    /* indel_penalty will be counted later */
    new->penalties = 0;
#else
    new->penalties = indel_penalty;
#endif
#endif

    /* new->gene_overlap = NO_KNOWN_GENE; -- initialized later when resolving multimappers */
    new->tally = -1L;

    if (new->score < *found_score) {
      *found_score = new->score;
    }

    new->nsplices = 0;

    new->distance = 0U;
    new->shortexonA_distance = new->shortexonD_distance = 0U;

    new->paired_usedp = false;
    new->paired_seenp = false;
    new->concordantp = false;

    new->circularpos = compute_circularpos(&new->alias,new);
    if (new->alias == +2 || new->alias == -2) {
      Stage3end_free(&new);
      return (T) NULL;
    } else {
      return new;
    }
  }
}


T
Stage3end_new_deletion (int *found_score, int nindels, int indel_pos, int nmismatches1_whole, int nmismatches2_whole,
			Univcoord_T left, int genomiclength, Compress_T query_compress,
			int querylength, bool plusp, int genestrand, bool first_read_p,
			Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh, Chrpos_T chrlength,
			int indel_penalty, bool sarrayp) {
  T new;
  Substring_T substring1, substring2;
  Junction_T junction;
  int querystart1, queryend1, querystart2, queryend2;
  Univcoord_T genomicstart1, genomicend1, genomicstart2, genomicend2;
  Univcoord_T genomicstart_adj_2, genomicend_adj_2;
  Univcoord_T alignstart1, alignend1, alignstart2, alignend2;
  int outofbounds_start = 0, outofbounds_end = 0;


  debug3(printf("Entered with left %llu, querylength %d, genomiclength %d, indel_pos %d\n",
		(unsigned long long) left,querylength,genomiclength,indel_pos));
#if 0
  debug3(printf("q: %s\n",query));
  debug3(printf("g: %s\n",genomicseg));
#endif

  assert(nindels > 0);

  querystart1 = 0;
  queryend1 = indel_pos;
  querystart2 = indel_pos;  /* Do not add nindels */
  queryend2 = querylength;


  if (plusp == true) {
    alignstart1 = left /*+ querystart1 (0)*/;
    alignend1 = left + indel_pos;
    alignstart2 = (left + nindels) + indel_pos;
    alignend2 = (left + nindels) + querylength;

    genomicstart1 = alignstart1;
    genomicend1 = alignend1;
    genomicstart2 = alignstart2;
    genomicend2 = alignend2;

    debug3(printf("plusp is true.  genomicstart %llu, genomicend %llu, alignstart1 %llu, alignend1 %llu, alignstart2 %llu, alignend2 %llu, left1 %llu\n",
		  (unsigned long long) genomicstart,(unsigned long long) genomicend,
		  (unsigned long long) alignstart1,(unsigned long long) alignend1,(unsigned long long) alignstart2,
		  (unsigned long long) alignend2,(unsigned long long) left));


    if (genomicstart1 < chroffset && genomicend2 > chrhigh) {
      /* Out of bounds on both sides */
      return (T) NULL;

    } else if (genomicstart1 < chroffset) {
      outofbounds_start = chroffset - genomicstart1;
      outofbounds_end = genomicend2 - chroffset;
      debug0(printf("Out of bounds left (low) %d, out of bounds right (high) %d\n",outofbounds_start,outofbounds_end));
      if (outofbounds_start > outofbounds_end) {
	/* Consider high part to be out of bounds and keep existing chromosome */
	if (querylength - indel_pos - nindels < outofbounds_end) {
	  /* indel is in eliminated part, so abort */
	  return (T) NULL;
	}
	outofbounds_start = 0;
      } else {
	/* Consider low part to be out of bounds and stay in this chromosome */
	if (indel_pos < outofbounds_start) {
	  /* indel is in eliminated part, so abort */
	  return (T) NULL;
	}
	/* Univ_IIT_interval_bounds(&chroffset,&chrhigh,&chrlength,chromosome_iit,chrnum,circular_typeint); */
	outofbounds_end = 0;
      }

    } else if (genomicend2 > chrhigh) {
      outofbounds_start = chrhigh - genomicstart1;
      outofbounds_end = genomicend2 - chrhigh;
      debug0(printf("Out of bounds left (low) %d, out of bounds right (high) %d\n",outofbounds_start,outofbounds_end));
      if (outofbounds_start > outofbounds_end) {
	/* Consider high part to be out of bounds and keep existing chromosome */
	if (querylength - indel_pos - nindels < outofbounds_end) {
	  /* indel is in eliminated part, so abort */
	  return (T) NULL;
	}
	outofbounds_start = 0;
      } else if (++chrnum > nchromosomes) {
	return (T) NULL;
      } else {
	/* Consider low part to be out of bounds and move to next chromosome */
	if (indel_pos < outofbounds_start) {
	  /* indel is in eliminated part, so abort */
	  return (T) NULL;
	}
	/* Univ_IIT_interval_bounds(&chroffset,&chrhigh,&chrlength,chromosome_iit,chrnum,circular_typeint); */
	outofbounds_end = 0;
      }
    }

  } else {
    alignstart1 = left + (querylength /*- querystart (0)*/);
    alignend1 = left + (querylength - indel_pos);
    alignstart2 = (left - nindels) + (querylength - indel_pos);
    if (left < nindels) {
      /* Handles deletion at start of genome, and avoids left alignend2 becoming negative */
      alignend2 = 0;
    } else {
      alignend2 = (left - nindels) /*+ querylength - queryend (querylength)*/;
    }

    genomicstart1 = alignstart1;
    genomicend1 = alignend1;
    genomicstart2 = alignstart2;
    genomicend2 = alignend2;

    debug3(printf("plusp is false.  genomicstart %llu, genomicend %llu, alignstart1 %llu, alignend1 %llu, alignstart2 %llu, alignend2 %llu, left1 %llu\n",
		  (unsigned long long) genomicstart,(unsigned long long) genomicend,
		  (unsigned long long) alignstart1,(unsigned long long) alignend1,(unsigned long long) alignstart2,
		  (unsigned long long) alignend2,(unsigned long long) left));

    if (genomicend2 < chroffset && genomicstart1 > chrhigh) {
      /* Out of bounds on both sides */
      return (T) NULL;

    } else if (genomicend2 < chroffset) {
      outofbounds_end = chroffset - genomicend2;
      outofbounds_start = genomicstart1 - chroffset;
      debug0(printf("Out of bounds left (high) %d, out of bounds right (low) %d\n",outofbounds_start,outofbounds_end));
      if (outofbounds_end > outofbounds_start) {
	/* Consider high part to be out of bounds and keep existing chromosome */
	if (indel_pos < outofbounds_start) {
	  /* indel is in eliminated part, so abort */
	  return (T) NULL;
	}
	outofbounds_end = 0;
      } else {
	/* Consider low part to be out of bounds.  Stay in this chromosome */
	if (querylength - indel_pos - nindels < outofbounds_end) {
	  /* indel is in eliminated part, so abort */
	  return (T) NULL;
	}
	/* Univ_IIT_interval_bounds(&chroffset,&chrhigh,&chrlength,chromosome_iit,chrnum,circular_typeint); */
	outofbounds_start = 0;
      }

    } else if (genomicstart1 > chrhigh) {
      outofbounds_end = chrhigh - genomicend2;
      outofbounds_start = genomicstart1 - chrhigh;
      debug0(printf("Out of bounds left (high) %d, out of bounds right (low) %d\n",outofbounds_start,outofbounds_end));
      if (outofbounds_end > outofbounds_start) {
	/* Consider high part to be out of bounds and keep existing chromosome */
	if (indel_pos < outofbounds_start) {
	  /* indel is in eliminated part, so abort */
	  return (T) NULL;
	}
	outofbounds_end = 0;
      } else if (++chrnum > nchromosomes) {
	return (T) NULL;
      } else {
	/* Consider low part to be out of bounds and move to next chromosome */
	if (querylength - indel_pos - nindels < outofbounds_end) {
	  /* indel is in eliminated part, so abort */
	  return (T) NULL;
	}
	Univ_IIT_interval_bounds(&chroffset,&chrhigh,&chrlength,chromosome_iit,chrnum,circular_typeint);
	outofbounds_start = 0;
      }
    }
  }

  if ((substring1 = Substring_new(nmismatches1_whole,chrnum,chroffset,chrhigh,chrlength,
				  query_compress,/*start_endtype*/END,/*end_endtype*/DEL,
				  querystart1,queryend1,querylength,alignstart1,alignend1,genomiclength,
				  /*exactp*/false,plusp,genestrand,first_read_p,
				  /*trim_left_p (previously was end1_indel_p ? false : true)*/true,
				  /*trim_right_p*/false,outofbounds_start,/*outofbounds_end*/0,/*minlength*/0)) == NULL) {
    return (T) NULL;

  } else if ((substring2 = Substring_new(nmismatches2_whole,chrnum,chroffset,chrhigh,chrlength,
					 query_compress,/*start_endtype*/DEL,/*end_endtype*/END,
					 querystart2,queryend2,querylength,alignstart2,alignend2,genomiclength,
					 /*exactp*/false,plusp,genestrand,first_read_p,
					 /*trim_left_p*/false,/*trim_right_p (previously was end2_indel_p ? false : true) */true,
					 /*outofbounds_start*/0,outofbounds_end,/*minlength*/0)) == NULL) {
    Substring_free(&substring1);
    return (T) NULL;

  } else {
    new = (T) MALLOC_OUT(sizeof(*new));
    debug0(printf("Stage3end_new_deletion %p: left %llu, chrnum %d, nmismatches %d+%d, indel_pos %d, nindels %d, sarrayp %d\n",
		  new,(unsigned long long) left,chrnum,nmismatches1_whole,nmismatches2_whole,indel_pos,nindels,sarrayp));

    new->pairarray = (struct Pair_T *) NULL;
    new->cigar_tokens = (List_T) NULL;
    new->gmap_intronp = false;

    /* Deletion contents are always from plus genomic strand */
    if (plusp == true) {
      junction = Junction_new_deletion(nindels,/*deletionpos*/left + indel_pos);
    } else {
      junction = Junction_new_deletion(nindels,/*deletionpos*/left + (querylength - indel_pos));
    }
    new->junctions_LtoH = List_push(NULL,junction);
    new->junctions_1toN = List_push(NULL,junction);
    new->junctions_Nto1 = List_push(NULL,junction);


    /* Initialize so Substring_free will not try to free */
    /* Filled in by Stage3end_display_prep */
    new->substrings_1toN = List_push(NULL,substring2);
    new->substrings_1toN = List_push(new->substrings_1toN,substring1);

    new->substrings_Nto1 = List_push(NULL,substring1);
    new->substrings_Nto1 = List_push(new->substrings_Nto1,substring2);

    if (plusp == true) {
      new->substrings_LtoH = List_push(NULL,substring2);
      new->substrings_LtoH = List_push(new->substrings_LtoH,substring1);
    } else {
      new->substrings_LtoH = List_push(NULL,substring1);
      new->substrings_LtoH = List_push(new->substrings_LtoH,substring2);
    }

    new->querylength = querylength;
    new->querylength_adj = querylength + nindels;
    new->genomicstart = genomicstart1;
    new->genomicend = genomicend2;

    if (genomicstart1 < genomicend2) {
      new->low = genomicstart1;
      new->high = genomicend2;
    } else {
      new->low = genomicend2;
      new->high = genomicstart1;
    }
    new->genomiclength = new->high - new->low;
    new->guided_insertlength = 0U;

    new->hittype = DELETION;
    new->genestrand = genestrand;
    new->sarrayp = sarrayp;
    new->gmap_source = GMAP_NOT_APPLICABLE;
    new->improved_by_gmap_p = false;

    new->chrnum = new->effective_chrnum = chrnum;
    new->other_chrnum = 0;
    new->chroffset = chroffset;
    new->chrhigh = chrhigh;
    new->chrlength = chrlength;
    new->plusp = plusp;
    new->sensedir = SENSE_NULL;

#if 0
    new->mapq_loglik = Substring_mapq_loglik(substring1) + Substring_mapq_loglik(substring2);
    new->mapq_score = 0;
    new->absmq_score = 0;
#endif

    new->nindels = nindels;
    new->nmismatches_whole = nmismatches1_whole + nmismatches2_whole;
    new->ntscore = indel_penalty + nmismatches1_whole + nmismatches2_whole;
    new->score = new->ntscore;
    new->nsegments = 2;

    new->nmismatches_bothdiff = Substring_nmismatches_bothdiff(substring1) + Substring_nmismatches_bothdiff(substring2);
    /* new->nmismatches_refdiff = Substring_nmismatches_refdiff(substring1) + Substring_nmismatches_refdiff(substring2); */

#if 0
    /* This method is correct for SNP-tolerant alignment */
    new->nmatches = Substring_match_length(substring1) + Substring_match_length(substring2) - new->total_nmismatches;
#else
    /* This method is now correct for SNP-tolerant alignment */
    new->nmatches = Substring_nmatches(substring1) + Substring_nmatches(substring2);
    new->nmatches_posttrim = Substring_nmatches_posttrim(substring1) + Substring_nmatches_posttrim(substring2);
    new->nmatches_posttrim -= indel_penalty; /* for use in goodness_cmp procedures */
#endif

    new->trim_left = Substring_trim_left(substring1);
    new->trim_right = Substring_trim_right(substring2);
    new->trim_left_splicep = Substring_trim_left_splicep(substring1);
    new->trim_right_splicep = Substring_trim_right_splicep(substring2);

#if 0
#ifdef SCORE_INDELS
    /* indel_penalty will be counted later */
    new->penalties = 0;
#else
    new->penalties = indel_penalty;
#endif
#endif

    /* new->gene_overlap = NO_KNOWN_GENE; -- initialized later when resolving multimappers */
    new->tally = -1L;

    if (new->score < *found_score) {
      *found_score = new->score;
    }

    new->nsplices = 0;

    new->distance = 0U;
    new->shortexonA_distance = new->shortexonD_distance = 0U;

    new->paired_usedp = false;
    new->paired_seenp = false;
    new->concordantp = false;

    new->circularpos = compute_circularpos(&new->alias,new);
    if (new->alias == +2 || new->alias == -2) {
      Stage3end_free(&new);
      return (T) NULL;
    } else {
      return new;
    }
  }
}


/* Never returns NULL */
/* Previously new->substring1 was donor and new->substring2 was acceptor */
T
Stage3end_new_splice (int *found_score, int nmismatches_donor, int nmismatches_acceptor,
		      Substring_T donor, Substring_T acceptor,
		      double donor_prob, double acceptor_prob, Chrpos_T distance,
		      bool shortdistancep, int splicing_penalty, int querylength, int amb_length, double amb_prob,
#ifdef LARGE_GENOMES
		      Uint8list_T ambcoords_donor, Uint8list_T ambcoords_acceptor,
#else
		      Uintlist_T ambcoords_donor, Uintlist_T ambcoords_acceptor,
#endif
		      Intlist_T amb_knowni_donor, Intlist_T amb_knowni_acceptor,
		      Intlist_T amb_nmismatches_donor, Intlist_T amb_nmismatches_acceptor,
		      Doublelist_T amb_probs_donor, Doublelist_T amb_probs_acceptor,
		      bool copy_donor_p, bool copy_acceptor_p, bool first_read_p, int sensedir,
		      bool sarrayp) {
  T new;
  Substring_T substring_for_concordance; /* always the inner substring */
  Substring_T substring_other;		 /* the outer substring */
  Substring_T substring;
  Junction_T junction;
#ifdef DEBUG0
  int i;
#endif

  new = (T) MALLOC_OUT(sizeof(*new));
  debug0(printf("Stage3end_new_splice %p with sensedir %d, donor substring %p and acceptor substring %p, and amb_length %d, sarrayp %d\n",
		new,sensedir,donor,acceptor,amb_length,sarrayp));

#if 0
  assert(Substring_match_length_orig(donor) + Substring_match_length_orig(acceptor) + amb_length == querylength);
#endif

  new->querylength_adj = new->querylength = querylength;

#if 0
  if (donor == NULL) {
    /* new->substring1 = copy_acceptor_p ? Substring_copy(acceptor) : acceptor; */
    /* new->substring2 = (Substring_T) NULL; */
    new->substring_donor = (Substring_T) NULL;
    new->substring_acceptor = new->substring1;

  } else if (acceptor == NULL) {
    /* new->substring1 = copy_donor_p ? Substring_copy(donor) : donor; */
    /* new->substring2 = (Substring_T) NULL; */
    new->substring_donor = new->substring1;
    new->substring_acceptor = (Substring_T) NULL;

  } else {
    if (sensedir != SENSE_ANTI) {
      /* SENSE_FORWARD or SENSE_NULL */
      /* new->substring1 = copy_donor_p ? Substring_copy(donor) : donor; */
      /* new->substring2 = copy_acceptor_p ? Substring_copy(acceptor) : acceptor; */
      new->substring_donor = new->substring1;
      new->substring_acceptor = new->substring2;
    } else {
      /* new->substring1 = copy_acceptor_p ? Substring_copy(acceptor) : acceptor; */
      /* new->substring2 = copy_donor_p ? Substring_copy(donor) : donor; */
      new->substring_donor = new->substring2;
      new->substring_acceptor = new->substring1;
    }
  }
#endif

  new->nindels = 0;

  new->pairarray = (struct Pair_T *) NULL;
  new->cigar_tokens = (List_T) NULL;
  new->gmap_intronp = false;

  new->sarrayp = sarrayp;
  new->gmap_source = GMAP_NOT_APPLICABLE;
  new->improved_by_gmap_p = false;

  if (donor == NULL) {
    new->hittype = SPLICE;
    new->genestrand = Substring_genestrand(acceptor);
    new->chrnum = Substring_chrnum(acceptor);
    new->chroffset = Substring_chroffset(acceptor);
    new->chrhigh = Substring_chrhigh(acceptor);
    new->chrlength = Substring_chrlength(acceptor);
    new->plusp = Substring_plusp(acceptor);

  } else if (acceptor == NULL) {
    new->hittype = SPLICE;
    new->genestrand = Substring_genestrand(donor);
    new->chrnum = Substring_chrnum(donor);
    new->chroffset = Substring_chroffset(donor);
    new->chrhigh = Substring_chrhigh(donor);
    new->chrlength = Substring_chrlength(donor);
    new->plusp = Substring_plusp(donor);

  } else if (shortdistancep == true) {
    new->hittype = SPLICE;
    new->genestrand = Substring_genestrand(donor);
    new->chrnum = Substring_chrnum(donor);
    new->chroffset = Substring_chroffset(donor);
    new->chrhigh = Substring_chrhigh(donor);
    new->chrlength = Substring_chrlength(donor);

#if 0
    if (sensedir == SENSE_FORWARD) {
      if (first_read_p) {
	new->plusp = Substring_plusp(acceptor);
      } else {
	new->plusp = Substring_plusp(donor);
      }

    } else if (sensedir == SENSE_ANTI) {
      if (first_read_p) {
	new->plusp = Substring_plusp(donor);
      } else {
	new->plusp = Substring_plusp(acceptor);
      }

    } else {
      /* No good selection for SENSE_NULL */
      new->plusp = Substring_plusp(donor);
    }
#else
    assert(Substring_plusp(donor) == Substring_plusp(acceptor));
    assert(Substring_chimera_sensedir(donor) == Substring_chimera_sensedir(acceptor));
    new->plusp = Substring_plusp(donor);
#endif


#if 0
  } else if (merge_samechr_p == false) {
    new->hittype = DISTANT_SPLICE;
    new->sarrayp = sarrayp;
    new->gmap_source = GMAP_NOT_APPLICABLE;
    new->improved_by_gmap_p = false;
    new->chrnum = 0;
    new->chroffset = 0;
    new->chrhigh = 0;
    new->chrlength = 0;
#endif

  } else {
    new->sarrayp = sarrayp;
    new->gmap_source = GMAP_NOT_APPLICABLE;
    new->improved_by_gmap_p = false;
    if (Substring_chrnum(donor) == Substring_chrnum(acceptor)) {
      new->hittype = SAMECHR_SPLICE;
      new->genestrand = Substring_genestrand(donor);
      new->chrnum = Substring_chrnum(donor);
      new->chroffset = Substring_chroffset(donor);
      new->chrhigh = Substring_chrhigh(donor);
      new->chrlength = Substring_chrlength(donor);
      new->plusp = Substring_plusp(donor); /* default value, used if merge_samechr_p is true */

    } else {
      new->hittype = TRANSLOC_SPLICE;
      new->genestrand = 0;
      new->chrnum = 0;
      new->chroffset = 0;
      new->chrhigh = 0;
      new->chrlength = 0;
    }

    /* new->plusp assigned below */

#if 0
    if (Substring_plusp(donor) == Substring_plusp(acceptor)) {
      new->plusp = Substring_plusp(donor);
    } else {
      /* Not sure what to do here.  Probably need to have substring->dir rather than substring->plusp. */
      /* Look at ss.samechr for an example.  plusp true => pair_inversion, plusp false => pair_scramble. */
      new->plusp = true;
    }
#endif
  }

  /* printf("Making splice with shortdistancep = %d, donor chrnum %d, and acceptor chrnum %d => chrnum %d\n",
     shortdistancep,Substring_chrnum(donor),Substring_chrnum(acceptor),new->chrnum); */

  donor = copy_donor_p ? Substring_copy(donor) : donor;
  acceptor = copy_acceptor_p ? Substring_copy(acceptor) : acceptor;

  if (sensedir != SENSE_ANTI) {
    /* SENSE_FORWARD or SENSE_NULL */
    if (donor == NULL) {
      new->genomicstart = Substring_genomicstart(acceptor);
      new->genomicend = Substring_genomicend(acceptor);
      
      donor = Substring_new_ambig(/*querystart*/0,/*queryend*/Substring_querystart(acceptor),
				  /*splice_pos*/Substring_querystart(acceptor),querylength,
				  new->chrnum,new->chroffset,new->chrhigh,new->chrlength,
				  /*genomiclength*/querylength,new->plusp,new->genestrand,first_read_p,
				  ambcoords_donor,amb_knowni_donor,amb_nmismatches_donor,amb_probs_donor,
				  /*amb_common_prob*/acceptor_prob,/*amb_donor_common_p*/false,
				  /*substring1p*/true);
      debug0(printf("Making sense ambiguous donor at %d..%d with %d matches\n",
		    0,Substring_querystart(acceptor),Substring_nmatches(donor)));
      donor_prob = Doublelist_max(amb_probs_donor);

    } else if (acceptor == NULL) {
      new->genomicstart = Substring_genomicstart(donor);
      new->genomicend = Substring_genomicend(donor);

      acceptor = Substring_new_ambig(/*querystart*/Substring_queryend(donor),/*queryend*/querylength,
				     /*splice_pos*/Substring_queryend(donor),querylength,
				     new->chrnum,new->chroffset,new->chrhigh,new->chrlength,
				     /*genomiclength*/querylength,new->plusp,new->genestrand,first_read_p,
				     ambcoords_acceptor,amb_knowni_acceptor,amb_nmismatches_acceptor,amb_probs_acceptor,
				     /*amb_common_prob*/donor_prob,/*amb_donor_common_p*/true,
				     /*substring1p*/false);
      debug0(printf("Making sense ambiguous donor at %d..%d with %d matches\n",
		    Substring_queryend(donor),querylength,Substring_nmatches(acceptor)));
      acceptor_prob = Doublelist_max(amb_probs_acceptor);

    } else {
      new->genomicstart = Substring_genomicstart(donor);
      new->genomicend = Substring_genomicend(acceptor);
    }

  } else {
    /* SENSE_ANTI */
    if (donor == NULL) {
      new->genomicstart = Substring_genomicstart(acceptor);
      new->genomicend = Substring_genomicend(acceptor);

      donor = Substring_new_ambig(/*querystart*/Substring_queryend(acceptor),/*queryend*/querylength,
				  /*splice_pos*/Substring_queryend(acceptor),querylength,
				  new->chrnum,new->chroffset,new->chrhigh,new->chrlength,
				  /*genomiclength*/querylength,new->plusp,new->genestrand,first_read_p,
				  ambcoords_donor,amb_knowni_donor,amb_nmismatches_donor,amb_probs_donor,
				  /*amb_common_prob*/acceptor_prob,/*amb_donor_common_p*/false,
				  /*substring1p*/false);
      debug0(printf("Making antisense ambiguous donor at %d..%d with %d matches\n",
		    Substring_queryend(acceptor),querylength,Substring_nmatches(donor)));
      donor_prob = Doublelist_max(amb_probs_donor);

    } else if (acceptor == NULL) {
      new->genomicstart = Substring_genomicstart(donor);
      new->genomicend = Substring_genomicend(donor);

      acceptor = Substring_new_ambig(/*querystart*/0,/*queryend*/Substring_querystart(donor),
				     /*splice_pos*/Substring_querystart(donor),querylength,
				     new->chrnum,new->chroffset,new->chrhigh,new->chrlength,
				     /*genomiclength*/querylength,new->plusp,new->genestrand,first_read_p,
				     ambcoords_acceptor,amb_knowni_acceptor,amb_nmismatches_acceptor,amb_probs_acceptor,
				     /*amb_common_prob*/donor_prob,/*amb_donor_common_p*/true,
				     /*substring1p*/true);
      debug0(printf("Making antisense ambiguous acceptor at %d..%d with %d matches\n",
		    0,Substring_querystart(donor),Substring_nmatches(acceptor)));
      acceptor_prob = Doublelist_max(amb_probs_acceptor);

    } else {
      new->genomicstart = Substring_genomicstart(acceptor);
      new->genomicend = Substring_genomicend(donor);
    }
  }


  if (new->genomicstart < new->genomicend) {
    new->low = new->genomicstart;
    new->high = new->genomicend;

  } else {
    new->low = new->genomicend;
    new->high = new->genomicstart;
  }

  debug0(printf("  hittype is %s, plusp %d, genomicpos %u..%u\n",
		hittype_string(new->hittype),new->plusp,new->genomicstart - new->chroffset,new->genomicend - new->chroffset));
#if 0
  debug0(printf("start_ambiguous_p %d (%d starts), end_ambiguous_p %d (%d ends)\n",
		new->start_ambiguous_p,new->start_nambcoords,new->end_ambiguous_p,new->end_nambcoords));
#endif
#if 0
  for (i = 0; i < new->start_nambcoords; i++) {
    printf("amb start %u\n",new->start_ambcoords[i]);
  }
  for (i = 0; i < new->end_nambcoords; i++) {
    printf("amb end %u\n",new->end_ambcoords[i]);
  }
#endif

  new->genomiclength = new->high - new->low;
  new->guided_insertlength = 0U;

  new->nsplices = 1;

  if (new->chrnum == 0) {
    /* Previously also did this for (donor != NULL && acceptor != NULL && shortdistancep == false), but this led to the wrong chrpos for SAM output */
    /* Checking for merge_samechr_p leads to wrong mappingstart and mappingend for running GMAP */
    /* Always want the original query end */
    junction = Junction_new_chimera(sensedir,donor_prob,acceptor_prob);
    new->junctions_LtoH = List_push(NULL,junction);
    new->junctions_1toN = List_push(NULL,junction);
    new->junctions_Nto1 = List_push(NULL,junction);

    /* For translocations, LtoH makes no sense, so we rely upon 1toN and force LtoH to be the same as 1toN */
    debug0(printf("donor querypos %d..%d\n",Substring_querystart(donor),Substring_queryend(donor)));
    debug0(printf("acceptor querypos %d..%d\n",Substring_querystart(acceptor),Substring_queryend(acceptor)));

    if (Substring_querystart(donor) < Substring_querystart(acceptor)) {
      new->sensedir = SENSE_FORWARD;
      new->substrings_LtoH = List_push(NULL,(void *) acceptor);
      new->substrings_LtoH = List_push(new->substrings_LtoH,(void *) donor);
    } else {
      new->sensedir = SENSE_ANTI;
      new->substrings_LtoH = List_push(NULL,(void *) donor);
      new->substrings_LtoH = List_push(new->substrings_LtoH,(void *) acceptor);
    }
    debug0(printf("sensedir %d\n",new->sensedir));


    new->substrings_1toN = List_copy(new->substrings_LtoH);
    new->substrings_Nto1 = List_reverse(List_copy(new->substrings_LtoH));
    assert(Substring_querystart(List_head(new->substrings_1toN)) < Substring_querystart(List_head(new->substrings_Nto1)));


    if (first_read_p == true) {
      substring_for_concordance = (Substring_T) List_head(new->substrings_Nto1);
      substring_other = (Substring_T) List_head(new->substrings_1toN);
    } else {
      substring_for_concordance = (Substring_T) List_head(new->substrings_1toN);
      substring_other = (Substring_T) List_head(new->substrings_Nto1);
    }
    new->effective_chrnum = Substring_chrnum(substring_for_concordance);
    new->other_chrnum = Substring_chrnum(substring_other);


    /* Redefine based on inner substring */
    new->genomicstart = Substring_genomicstart(substring_for_concordance);
    new->genomicend = Substring_genomicend(substring_for_concordance);
    new->plusp = Substring_plusp(substring_for_concordance);
    
  } else {
    /* Ordinary splice */
    new->effective_chrnum = new->chrnum;
    new->other_chrnum = 0;

    new->substrings_LtoH = (List_T) NULL;
    new->junctions_LtoH = (List_T) NULL;
    new->sensedir = sensedir;

    if (sensedir != SENSE_ANTI) {
      /* SENSE_FORWARD or SENSE_NULL */
      if (new->plusp == true) {
	/* Order is donor, acceptor.  Same as substring1, substring2, as expected */
	new->substrings_LtoH = List_push(new->substrings_LtoH,(void *) acceptor);
	junction = Junction_new_splice(distance,sensedir,donor_prob,acceptor_prob);
	new->junctions_LtoH = List_push(new->junctions_LtoH,(void *) junction);
	new->substrings_LtoH = List_push(new->substrings_LtoH,(void *) donor);
	  
      } else {
	/* Order is acceptor, donor.  Same as substring2, substring1, as expected */
	new->substrings_LtoH = List_push(new->substrings_LtoH,(void *) donor);
	junction = Junction_new_splice(distance,sensedir,donor_prob,acceptor_prob);
	new->junctions_LtoH = List_push(new->junctions_LtoH,(void *) junction);
	new->substrings_LtoH = List_push(new->substrings_LtoH,(void *) acceptor);
      }

    } else {
      /* SENSE_ANTI */
      if (new->plusp == true) {
	/* Order is acceptor, donor.  Same as substring1, substring2, as expected */
	new->substrings_LtoH = List_push(new->substrings_LtoH,(void *) donor);
	junction = Junction_new_splice(distance,sensedir,donor_prob,acceptor_prob);
	new->junctions_LtoH = List_push(new->junctions_LtoH,(void *) junction);
	new->substrings_LtoH = List_push(new->substrings_LtoH,(void *) acceptor);
      } else {
	/* Order is donor, acceptor.  Same as substring2, substring1, as expected */
	new->substrings_LtoH = List_push(new->substrings_LtoH,(void *) acceptor);
	junction = Junction_new_splice(distance,sensedir,donor_prob,acceptor_prob);
	new->junctions_LtoH = List_push(new->junctions_LtoH,(void *) junction);
	new->substrings_LtoH = List_push(new->substrings_LtoH,(void *) donor);
      }
    }

    if (new->plusp == true) {
      new->substrings_1toN = List_copy(new->substrings_LtoH);
      new->substrings_Nto1 = List_reverse(List_copy(new->substrings_LtoH));
      new->junctions_1toN = List_copy(new->junctions_LtoH);
      new->junctions_Nto1 = List_reverse(List_copy(new->junctions_LtoH));
    } else {
      new->substrings_1toN = List_reverse(List_copy(new->substrings_LtoH));
      new->substrings_Nto1 = List_copy(new->substrings_LtoH);
      new->junctions_1toN = List_reverse(List_copy(new->junctions_LtoH));
      new->junctions_Nto1 = List_copy(new->junctions_LtoH);
    }
  }


  new->nmismatches_whole = nmismatches_donor + nmismatches_acceptor;
  new->score = new->ntscore = splicing_penalty + new->nmismatches_whole;
  new->nsegments = 2;
#if 0
  if (sensedir == SENSE_FORWARD) {
    new->score += antistranded_penalty;
  }
#endif

#if 0
  if (donor == NULL) {
    /* new->mapq_loglik = Substring_mapq_loglik(acceptor); */
    new->nmismatches_bothdiff = Substring_nmismatches_bothdiff(acceptor) + nmismatches_donor;
    /* new->nmismatches_refdiff = Substring_nmismatches_refdiff(acceptor) + nmismatches_donor; */
    new->nmatches = Substring_nmatches(acceptor);
    new->nmatches_posttrim = Substring_nmatches_posttrim(acceptor);
    if (favor_ambiguous_p == true) {
      new->nmatches += amb_length;
    }
    new->sensedir_nonamb = SENSE_NULL;	/* Ignore sense based on ambiguous end */
    debug0(printf("New splice has acceptor %d + amb %d matches, sensedir nonamb %d\n",
		  Substring_nmatches(acceptor),amb_length,new->sensedir_nonamb));
  } else if (acceptor == NULL) {
    /* new->mapq_loglik = Substring_mapq_loglik(donor); */
    new->nmismatches_bothdiff = Substring_nmismatches_bothdiff(donor) + nmismatches_acceptor;
    /* new->nmismatches_refdiff = Substring_nmismatches_refdiff(donor) + nmismatches_acceptor; */
    new->nmatches = Substring_nmatches(donor);
    new->nmatches_posttrim = Substring_nmatches_posttrim(donor);
    if (favor_ambiguous_p == true) {
      new->nmatches += amb_length;
    }
    new->sensedir_nonamb = SENSE_NULL;	/* Ignore sense based on ambiguous end */
    debug0(printf("New splice has donor %d + amb %d matches, sensedir nonamb %d\n",
		  Substring_nmatches(donor),amb_length,new->sensedir_nonamb));
  }
#endif

  /* new->mapq_loglik = Substring_mapq_loglik(donor) + Substring_mapq_loglik(acceptor); */
  new->nmismatches_bothdiff = Substring_nmismatches_bothdiff(donor) + Substring_nmismatches_bothdiff(acceptor);
  /* new->nmismatches_refdiff = Substring_nmismatches_refdiff(donor) + Substring_nmismatches_refdiff(acceptor); */
  new->nmatches = Substring_nmatches(donor) + Substring_nmatches(acceptor);
  new->nmatches_posttrim = Substring_nmatches_posttrim(donor) + Substring_nmatches_posttrim(acceptor);
  debug0(printf("New splice has donor %d + acceptor %d matches, sensedir %d\n",
		Substring_nmatches(donor),Substring_nmatches(acceptor),new->sensedir));

  substring = (Substring_T) List_head(new->substrings_1toN);
  new->trim_left = Substring_trim_left(substring);
  new->trim_left_splicep = Substring_trim_left_splicep(substring);

  substring = (Substring_T) List_head(new->substrings_Nto1);
  new->trim_right = Substring_trim_right(substring);
  new->trim_right_splicep = Substring_trim_right_splicep(substring);
  
  /* new->penalties = splicing_penalty; */

  /* new->gene_overlap = NO_KNOWN_GENE; -- initialized later when resolving multimappers */
  new->tally = -1L;

#if 0
  new->mapq_score = 0;
  new->absmq_score = 0;
#endif

  if (new->score < *found_score) {
    *found_score = new->score;
  }

  new->distance = distance;
  new->shortexonA_distance = new->shortexonD_distance = 0U;

  new->paired_usedp = false;
  new->paired_seenp = false;
  new->concordantp = false;

  new->circularpos = compute_circularpos(&new->alias,new);
  if (new->alias == +2 || new->alias == -2) {
    Stage3end_free(&new);
    return (T) NULL;
  } else {
    debug0(printf("Returning new splice %p at genomic %u..%u, donor %p (%u => %u), acceptor %p (%u => %u)\n",
		  new,new->genomicstart - new->chroffset,new->genomicend - new->chroffset,donor,
		  donor == NULL ? 0 : Substring_left_genomicseg(donor),
		  donor == NULL ? 0 : Substring_splicecoord(donor),
		  acceptor,acceptor == NULL ? 0 : Substring_left_genomicseg(acceptor),
		  acceptor == NULL ? 0 : Substring_splicecoord(acceptor)));
    debug0(printf("sensedir %d\n",new->sensedir));
    return new;
  }
}



/* Never returns NULL.  Never copies substrings.  Always shortdistance. */
/* [E Donor D] ----(junction)---- [A Shortexon D] ----(junction)---- [A Acceptor E] */
T
Stage3end_new_shortexon (int *found_score, Substring_T donor, Substring_T acceptor, Substring_T shortexon,
			 double donor_prob, double shortexonA_prob, double shortexonD_prob, double acceptor_prob,
			 int amb_length_donor, int amb_length_acceptor, double amb_prob_donor, double amb_prob_acceptor,
#ifdef LARGE_GENOMES
			 Uint8list_T ambcoords_donor, Uint8list_T ambcoords_acceptor,
#else
			 Uintlist_T ambcoords_donor, Uintlist_T ambcoords_acceptor,
#endif
			 Intlist_T amb_knowni_donor, Intlist_T amb_knowni_acceptor,
			 Intlist_T amb_nmismatches_donor, Intlist_T amb_nmismatches_acceptor,
			 Doublelist_T amb_probs_donor, Doublelist_T amb_probs_acceptor,
			 bool copy_donor_p, bool copy_acceptor_p, bool copy_shortexon_p,
			 int splicing_penalty, int querylength, bool first_read_p, int sensedir, bool sarrayp) {
  T new;
  int ignore;
  Substring_T substring, substring0, substring1, substring2;
  Chrpos_T distance;
  Junction_T junction0, junction2;
  
  new = (T) MALLOC_OUT(sizeof(*new));
  debug0(printf("Stage3end_new_shortexon %p, amb_donor %d, amb_acceptor %d, sensedir %d, sarrayp %d\n",
		new,amb_length_donor,amb_length_acceptor,sensedir,sarrayp));
  assert(Substring_match_length_orig(donor) + Substring_match_length_orig(shortexon) + Substring_match_length_orig(acceptor) +
	 amb_length_donor + amb_length_acceptor == querylength);

  new->querylength_adj = new->querylength = querylength;

  new->sarrayp = sarrayp;
  new->gmap_source = GMAP_NOT_APPLICABLE;
  new->improved_by_gmap_p = false;

#if 0
  if (donor == NULL && acceptor == NULL) {
    new->hittype = ONE_THIRD_SHORTEXON;

    new->shortexonA_distance = 0;
    new->shortexonD_distance = 0;

  }
#endif

#if 0
  if (donor == NULL) {
    new->hittype = TWO_THIRDS_SHORTEXON;
  } else if (acceptor == NULL) {
    new->hittype = TWO_THIRDS_SHORTEXON;
  } else {
    new->hittype = SHORTEXON;
  }
#else
  new->hittype = SUBSTRINGS;
#endif


  new->pairarray = (struct Pair_T *) NULL;
  new->cigar_tokens = (List_T) NULL;
  new->gmap_intronp = false;

  new->nindels = 0;

  new->chrnum = Substring_chrnum(shortexon);
  new->chroffset = Substring_chroffset(shortexon);
  new->chrhigh = Substring_chrhigh(shortexon);
  new->chrlength = Substring_chrlength(shortexon);
  new->plusp = Substring_plusp(shortexon);
  new->genestrand = Substring_genestrand(shortexon);


  /* Compute distances */
  if (donor == NULL) {
    new->shortexonA_distance = 0;
  } else if (Substring_splicecoord_A(shortexon) > Substring_splicecoord(donor)) {
    new->shortexonA_distance = Substring_splicecoord_A(shortexon) - Substring_splicecoord(donor);
  } else {
    new->shortexonA_distance = Substring_splicecoord(donor) - Substring_splicecoord_A(shortexon);
  }

  if (acceptor == NULL) {
    new->shortexonD_distance = 0;
  } else if (Substring_splicecoord_D(shortexon) > Substring_splicecoord(acceptor)) {
    new->shortexonD_distance = Substring_splicecoord_D(shortexon) - Substring_splicecoord(acceptor);
  } else {
    new->shortexonD_distance = Substring_splicecoord(acceptor) - Substring_splicecoord_D(shortexon);
  }
  new->distance = new->shortexonA_distance + new->shortexonD_distance;

  if (sensedir == SENSE_FORWARD) {
    new->genomicstart = (donor != NULL ? Substring_genomicstart(donor) : Substring_genomicstart(shortexon));
    new->genomicend = (acceptor != NULL ? Substring_genomicend(acceptor) : Substring_genomicend(shortexon));

  } else if (sensedir == SENSE_ANTI) {
    new->genomicstart = (acceptor != NULL ? Substring_genomicstart(acceptor) : Substring_genomicstart(shortexon));
    new->genomicend = (donor != NULL ? Substring_genomicend(donor) : Substring_genomicend(shortexon));

  } else {
    abort();
  }

  substring1 = copy_shortexon_p ? Substring_copy(shortexon) : shortexon;
  if (sensedir == SENSE_FORWARD) {
    substring0 = copy_donor_p ? Substring_copy(donor) : donor;
    if (donor == NULL) {
      donor = substring0 = Substring_new_ambig(/*querystart*/0,/*queryend*/Substring_querystart(shortexon),
					       /*splice_pos*/Substring_querystart(shortexon),querylength,
					       new->chrnum,new->chroffset,new->chrhigh,new->chrlength,
					       /*genomiclength*/querylength,new->plusp,new->genestrand,first_read_p,
					       ambcoords_donor,amb_knowni_donor,amb_nmismatches_donor,amb_probs_donor,
					       /*amb_common_prob*/acceptor_prob,/*amb_donor_common_p*/false,
					       /*substring1p*/true);
      /* new->start_amb_prob = Doublelist_max(amb_probs_donor); */
      /* new->start_amb_length = amb_length_donor; */
      junction0 = Junction_new_splice(/*distance*/0,sensedir,Doublelist_max(amb_probs_donor),shortexonA_prob);
    } else if (Substring_splicecoord_A(shortexon) > Substring_splicecoord(donor)) {
      distance = Substring_splicecoord_A(shortexon) - Substring_splicecoord(donor);
      junction0 = Junction_new_splice(distance,sensedir,donor_prob,shortexonA_prob);
    } else {
      distance = Substring_splicecoord(donor) - Substring_splicecoord_A(shortexon);
      junction0 = Junction_new_splice(distance,sensedir,donor_prob,shortexonA_prob);
    }

    substring2 = copy_acceptor_p ? Substring_copy(acceptor) : acceptor;
    if (acceptor == NULL) {
      acceptor = substring2 = Substring_new_ambig(/*querystart*/Substring_queryend(shortexon),/*queryend*/querylength,
						  /*splice_pos*/Substring_queryend(shortexon),querylength,
						  new->chrnum,new->chroffset,new->chrhigh,new->chrlength,
						  /*genomiclength*/querylength,new->plusp,new->genestrand,first_read_p,
						  ambcoords_acceptor,amb_knowni_acceptor,amb_nmismatches_acceptor,amb_probs_acceptor,
						  /*amb_common_prob*/donor_prob,/*amb_donor_common_p*/true,
						  /*substring1p*/false);
      /* new->end_amb_prob = Doublelist_max(amb_probs_acceptor); */
      /* new->end_amb_length = amb_length_acceptor; */
      junction2 = Junction_new_splice(/*distance*/0,sensedir,shortexonD_prob,Doublelist_max(amb_probs_acceptor));
    } else if (Substring_splicecoord_D(shortexon) > Substring_splicecoord(acceptor)) {
      distance = Substring_splicecoord_D(shortexon) - Substring_splicecoord(acceptor);
      junction2 = Junction_new_splice(distance,sensedir,shortexonD_prob,acceptor_prob);
    } else {
      distance = Substring_splicecoord(acceptor) - Substring_splicecoord_D(shortexon);
      junction2 = Junction_new_splice(distance,sensedir,shortexonD_prob,acceptor_prob);
    }

  } else if (sensedir == SENSE_ANTI) {
    substring0 = copy_acceptor_p ? Substring_copy(acceptor) : acceptor;
    if (acceptor == NULL) {
      acceptor = substring0 = Substring_new_ambig(/*querystart*/0,/*queryend*/Substring_querystart(shortexon),
						  /*splice_pos*/Substring_querystart(shortexon),querylength,
						  new->chrnum,new->chroffset,new->chrhigh,new->chrlength,
						  /*genomiclength*/querylength,new->plusp,new->genestrand,first_read_p,
						  ambcoords_acceptor,amb_knowni_acceptor,amb_nmismatches_acceptor,amb_probs_acceptor,
						  /*amb_common_prob*/donor_prob,/*amb_donor_common_p*/true,
						  /*substring1p*/true);
      /* new->start_amb_prob = Doublelist_max(amb_probs_acceptor); */
      /* new->start_amb_length = amb_length_acceptor; */
      junction0 = Junction_new_splice(/*distance*/0,sensedir,shortexonD_prob,Doublelist_max(amb_probs_acceptor));
    } else if (Substring_splicecoord_D(shortexon) > Substring_splicecoord(acceptor)) {
      distance = Substring_splicecoord_D(shortexon) - Substring_splicecoord(acceptor);
      junction0 = Junction_new_splice(distance,sensedir,shortexonD_prob,acceptor_prob);
    } else {
      distance = Substring_splicecoord(acceptor) - Substring_splicecoord_D(shortexon);
      junction0 = Junction_new_splice(distance,sensedir,shortexonD_prob,acceptor_prob);
    }

    substring2 = copy_donor_p ? Substring_copy(donor) : donor;
    if (donor == NULL) {
      donor = substring2 = Substring_new_ambig(/*querystart*/Substring_queryend(shortexon),/*queryend*/querylength,
					       /*splice_pos*/Substring_queryend(shortexon),querylength,
					       new->chrnum,new->chroffset,new->chrhigh,new->chrlength,
					       /*genomiclength*/querylength,new->plusp,new->genestrand,first_read_p,
					       ambcoords_donor,amb_knowni_donor,amb_nmismatches_donor,amb_probs_donor,
					       /*amb_common_prob*/acceptor_prob,/*amb_donor_common_p*/false,
					       /*substring1p*/false);
      /* new->end_amb_prob = Doublelist_max(amb_probs_donor); */
      /* new->end_amb_length = amb_length_donor; */
      junction2 = Junction_new_splice(/*distance*/0,sensedir,Doublelist_max(amb_probs_donor),shortexonA_prob);
    } else if (Substring_splicecoord_A(shortexon) > Substring_splicecoord(donor)) {
      distance = Substring_splicecoord_A(shortexon) - Substring_splicecoord(donor);
      junction2 = Junction_new_splice(distance,sensedir,donor_prob,shortexonA_prob);
    } else {
      distance = Substring_splicecoord(donor) - Substring_splicecoord_A(shortexon);
      junction2 = Junction_new_splice(distance,sensedir,donor_prob,shortexonA_prob);
    }

  } else {
    abort();
  }
  new->sensedir = sensedir;

  /* printf("Making splice with shortdistancep = %d, donor chrnum %d, and acceptor chrnum %d => chrnum %d\n",
     shortdistancep,Substring_chrnum(donor),Substring_chrnum(acceptor),new->chrnum); */


  if (new->genomicstart < new->genomicend) {
    debug0(printf("plus %s\n",print_sense(sensedir)));
    new->low = new->genomicstart;
    new->high = new->genomicend;

  } else {
    debug0(printf("minus %s\n",print_sense(sensedir)));
    new->low = new->genomicend;
    new->high = new->genomicstart;
  }

  debug0(printf("  hittype is %s, genomicpos %u..%u\n",
		hittype_string(new->hittype),new->genomicstart - new->chroffset,new->genomicend - new->chroffset));
  /* debug0(printf("start_ambiguous_p %d, end_ambiguous_p %d\n",new->start_ambiguous_p,new->end_ambiguous_p)); */

  new->genomiclength = new->high - new->low;
  new->guided_insertlength = 0U;

  new->nsplices = 2;

  new->effective_chrnum = new->chrnum;
  new->other_chrnum = 0;

  /* Currently not allowing translocations on shortexons */
  /* substring_for_concordance = (Substring_T) NULL; */

  new->substrings_LtoH = (List_T) NULL;
  new->junctions_LtoH = (List_T) NULL;
  if (new->plusp == true) {
    if (substring2 != NULL) {
      new->substrings_LtoH = List_push(new->substrings_LtoH,(void *) substring2);
      new->junctions_LtoH = List_push(new->junctions_LtoH,(void *) junction2);
    }
    new->substrings_LtoH = List_push(new->substrings_LtoH,(void *) substring1);
    if (substring0 != NULL) {
      new->substrings_LtoH = List_push(new->substrings_LtoH,(void *) substring0);
      new->junctions_LtoH = List_push(new->junctions_LtoH,(void *) junction0);
    }

  } else {
    if (substring0 != NULL) {
      new->substrings_LtoH = List_push(new->substrings_LtoH,(void *) substring0);
      new->junctions_LtoH = List_push(new->junctions_LtoH,(void *) junction0);
    }
    new->substrings_LtoH = List_push(new->substrings_LtoH,(void *) substring1);
    if (substring2 != NULL) {
      new->substrings_LtoH = List_push(new->substrings_LtoH,(void *) substring2);
      new->junctions_LtoH = List_push(new->junctions_LtoH,(void *) junction2);
    }
  }

  if (new->plusp == true) {
    new->substrings_1toN = List_copy(new->substrings_LtoH);
    new->substrings_Nto1 = List_reverse(List_copy(new->substrings_LtoH));
    new->junctions_1toN = List_copy(new->junctions_LtoH);
    new->junctions_Nto1 = List_reverse(List_copy(new->junctions_LtoH));
  } else {
    new->substrings_1toN = List_reverse(List_copy(new->substrings_LtoH));
    new->substrings_Nto1 = List_copy(new->substrings_LtoH);
    new->junctions_1toN = List_reverse(List_copy(new->junctions_LtoH));
    new->junctions_Nto1 = List_copy(new->junctions_LtoH);
  }


#if 0
  new->mapq_loglik = Substring_mapq_loglik(donor) + Substring_mapq_loglik(acceptor) + Substring_mapq_loglik(shortexon);
  new->mapq_score = 0;
  new->absmq_score = 0;
#endif

  new->nmismatches_whole = Substring_nmismatches_whole(shortexon);
  new->nmismatches_whole += Substring_nmismatches_whole(donor);
  new->nmismatches_whole += Substring_nmismatches_whole(acceptor);

  new->ntscore = splicing_penalty + splicing_penalty + new->nmismatches_whole;
  new->score = new->ntscore;
  new->nsegments = 3;
#if 0
  if (sensedir == SENSE_FORWARD) {
    new->score += antistranded_penalty;
  }
#endif

  new->nmismatches_bothdiff = Substring_nmismatches_bothdiff(shortexon);
  new->nmismatches_bothdiff += Substring_nmismatches_bothdiff(donor);
  new->nmismatches_bothdiff += Substring_nmismatches_bothdiff(acceptor);
  /* new->nmismatches_refdiff = Substring_nmismatches_refdiff(donor) + Substring_nmismatches_refdiff(acceptor) + Substring_nmismatches_refdiff(shortexon); */

  new->nmatches = Substring_nmatches(shortexon);
  new->nmatches_posttrim = Substring_nmatches_posttrim(shortexon);
  if (donor == NULL) {
    if (favor_ambiguous_p == true) {
      new->nmatches += amb_length_donor;
    }
  } else {
    /* assert(amb_length_donor == 0); */
    new->nmatches += Substring_nmatches(donor);
  }
  if (acceptor == NULL) {
    if (favor_ambiguous_p == true) {
      new->nmatches += amb_length_acceptor;
    }
  } else {
    /* assert(amb_length_acceptor == 0); */
    new->nmatches += Substring_nmatches(acceptor);
  }

  substring = (Substring_T) List_head(new->substrings_1toN);
  new->trim_left = Substring_trim_left(substring);
  new->trim_left_splicep = Substring_trim_left_splicep(substring);

  substring = (Substring_T) List_head(new->substrings_Nto1);
  new->trim_right = Substring_trim_right(substring);
  new->trim_right_splicep = Substring_trim_right_splicep(substring);

  /* new->penalties = splicing_penalty + splicing_penalty; */

  /* new->gene_overlap = NO_KNOWN_GENE; -- initialized later when resolving multimappers */
  new->tally = -1L;

  if (new->score < *found_score) {
    *found_score = new->score;
  }

  new->paired_usedp = false;
  new->paired_seenp = false;
  new->concordantp = false;

  new->circularpos = compute_circularpos(&new->alias,new);
  if (new->alias == +2 || new->alias == -2) {
    Stage3end_free(&new);
    return (T) NULL;
  } else {
    return new;
  }
}


T
Stage3end_new_terminal (int querystart, int queryend, Univcoord_T left, Compress_T query_compress,
			int querylength, bool plusp, int genestrand, bool first_read_p,
			Endtype_T start_endtype, Endtype_T end_endtype,
			Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh, Chrpos_T chrlength,
			int max_mismatches_allowed, bool sarrayp) {
  T new;
  Substring_T substring;
  Univcoord_T genomicstart, genomicend, alignstart, alignend, alignstart_trim, alignend_trim;
  int nmismatches_whole, minlength;
  bool trim_left_p, trim_right_p;
  int outofbounds_start = 0, outofbounds_end = 0;

  int *mismatch_positions_left, *mismatch_positions_right, nmismatches_left, nmismatches_right;
  int length_left, length_right, pos5, pos3;
  int i;


  debug0(printf("\nStage3end_new_terminal possible: endtypes %s and %s, left %llu, querystart %d, queryend %d, sarrayp %d\n",
		Endtype_string(start_endtype),Endtype_string(end_endtype),(unsigned long long) left,querystart,queryend,sarrayp));

  if (plusp == true) {
    genomicstart = left;
    genomicend = left + querylength;

    alignstart = genomicstart + querystart;
    alignend = genomicstart + queryend;

    if (genomicstart < chroffset && genomicend > chrhigh) {
      /* Out of bounds on both sides */
      return (T) NULL;

    } else if (genomicstart < chroffset) {
      outofbounds_start = chroffset - genomicstart;
      outofbounds_end = genomicend - chroffset;
      debug0(printf("Out of bounds left (low) %d, out of bounds right (high) %d\n",outofbounds_start,outofbounds_end));
      if (outofbounds_start > outofbounds_end) {
	/* Consider high part to be out of bounds and keep existing chromosome */
	outofbounds_start = 0;
      } else {
	/* Consider low part to be out of bounds and stay in this chromosome */
	/* Univ_IIT_interval_bounds(&chroffset,&chrhigh,&chrlength,chromosome_iit,chrnum,circular_typeint); */
	outofbounds_end = 0;
      }

    } else if (genomicend > chrhigh) {
      outofbounds_start = chrhigh - genomicstart;
      outofbounds_end = genomicend - chrhigh;
      debug0(printf("Out of bounds left (low) %d, out of bounds right (high) %d\n",outofbounds_start,outofbounds_end));
      if (outofbounds_start > outofbounds_end) {
	/* Consider high part to be out of bounds and keep existing chromosome */
	outofbounds_start = 0;
      } else if (++chrnum > nchromosomes) {
	return (T) NULL;
      } else {
	/* Consider low part to be out of bounds and move to next chromosome */
	Univ_IIT_interval_bounds(&chroffset,&chrhigh,&chrlength,chromosome_iit,chrnum,circular_typeint);
	outofbounds_end = 0;
      }
    }

  } else {
    genomicend = left;
    genomicstart = left + querylength;

    alignstart = genomicstart - querystart;
    alignend = genomicstart - queryend;

    if (genomicend < chroffset && genomicstart > chrhigh) {
      /* Out of bounds on both sides */
      return (T) NULL;

    } else if (genomicend < chroffset) {
      outofbounds_end = chroffset - genomicend;
      outofbounds_start = genomicstart - chroffset;
      debug0(printf("Out of bounds left (high) %d, out of bounds right (low) %d\n",outofbounds_start,outofbounds_end));
      if (outofbounds_end > outofbounds_start) {
	/* Consider high part to be out of bounds and keep existing chromosome */
	outofbounds_end = 0;
      } else {
	/* Consider low part to be out of bounds and stay in this chromosome */
	/* Univ_IIT_interval_bounds(&chroffset,&chrhigh,&chrlength,chromosome_iit,chrnum,circular_typeint); */
	outofbounds_start = 0;
      }

    } else if (genomicstart > chrhigh) {
      outofbounds_end = chrhigh - genomicend;
      outofbounds_start = genomicstart - chrhigh;
      debug0(printf("Out of bounds left (high) %d, out of bounds right (low) %d\n",outofbounds_start,outofbounds_end));
      if (outofbounds_end > outofbounds_start) {
	/* Consider high part to be out of bounds and keep existing chromosome */
	outofbounds_end = 0;
      } else if (++chrnum > nchromosomes) {
	return (T) NULL;
      } else {
	/* Consider low part to be out of bounds and move to next chromosome */
	Univ_IIT_interval_bounds(&chroffset,&chrhigh,&chrlength,chromosome_iit,chrnum,circular_typeint);
	outofbounds_start = 0;
      }
    }
  }

  if (start_endtype == TERM) {
    trim_left_p = true;
  } else {
    trim_left_p = false;
  }

  if (end_endtype == TERM) {
    trim_right_p = true;
  } else {
    trim_right_p = false;
  }

  /* Note: Changing querylength/3 to querylength/2 loses about 1% of concordant reads */
  minlength = querylength/3;
  if (minlength > TERMINAL_COMPUTE_MINLENGTH) {
    minlength = TERMINAL_COMPUTE_MINLENGTH;
  }

  if ((substring = Substring_new(/*nmismatches_whole*/0,chrnum,chroffset,chrhigh,chrlength,
				 query_compress,start_endtype,end_endtype,querystart,queryend,querylength,
				 alignstart,alignend,/*genomiclength*/querylength,
				 /*exactp*/false,plusp,genestrand,first_read_p,
				 trim_left_p,trim_right_p,outofbounds_start,outofbounds_end,minlength)) == NULL) {
    debug0(printf("returning NULL\n"));
    return (T) NULL;
    
#if 0
  } else if (start_endtype == TERM) {
    if (Substring_trim_left(substring) == 0) {
      /* Not really a terminal */
      Substring_free(&substring);
      return (T) NULL;
    }

  } else if (end_endtype == TERM) {
    if (Substring_trim_right(substring) == 0) {
      /* Not really a terminal */
      Substring_free(&substring);
      return (T) NULL;
    }
#endif
  }

  /* Re-compute nmismatches_whole and nmatches for terminal alignments */
  alignstart_trim = Substring_alignstart_trim(substring);
  alignend_trim = Substring_alignend_trim(substring);

  debug0(printf("alignstart_trim = %llu, alignend_trim = %llu\n",(unsigned long long) alignstart_trim,(unsigned long long) alignend_trim));
  if (plusp == true) {
    debug0(printf("plus: pos5 = %d, pos3 = %d\n",(int) (alignstart_trim-left),(int) (alignend_trim-left)));
    nmismatches_whole =
      Genome_count_mismatches_substring(query_compress,left,/*pos5*/alignstart_trim-left,
					/*pos3*/alignend_trim-left,/*plusp*/true,genestrand,first_read_p);
    debug0(printf("Recomputing nmismatches_whole as %d\n",nmismatches_whole));
    if (nmismatches_whole <= max_mismatches_allowed) {
      /* Code in substring.c suggests that nmismatches_bothdiff would be the same value as nmismatches_whole */
      Substring_set_nmismatches_terminal(substring,nmismatches_whole,/*nmismatches_bothdiff*/nmismatches_whole);

    } else {
      debug0(printf("returning NULL\n"));
      Substring_free(&substring);
      return (T) NULL;
    }

    if (0) {
      /* This may be dependent on the trimming algorithm, but is needed to avoid bad terminal alignments */
      Substring_free(&substring);
      debug0(printf("returning NULL\n"));
      return (T) NULL;
    }

  } else {
    debug0(printf("minus: pos5 = %d, pos3 = %d\n",(int) (alignend_trim-left),(int) (alignstart_trim-left)));
    nmismatches_whole =
      Genome_count_mismatches_substring(query_compress,left,/*pos5*/alignend_trim-left,
					/*pos3*/alignstart_trim-left,/*plusp*/false,genestrand,first_read_p);
    debug0(printf("Recomputing nmismatches_whole as %d\n",nmismatches_whole));
    if (nmismatches_whole <= max_mismatches_allowed) {
      /* Code in substring.c suggests that nmismatches_bothdiff would be the same value as nmismatches_whole */
      Substring_set_nmismatches_terminal(substring,nmismatches_whole,/*nmismatches_bothdiff*/nmismatches_whole);

    } else {
      debug0(printf("returning NULL\n"));
      Substring_free(&substring);
      return (T) NULL;
    }
  }

  new = (T) MALLOC_OUT(sizeof(*new));
  debug0(printf("Stage3end_new_terminal %p: endtypes %s and %s, left %llu, genomicstart/end %llu..%llu, chrhigh %llu, chrnum %d, querystart %d, queryend %d\n",
		new,Endtype_string(start_endtype),Endtype_string(end_endtype),
		(unsigned long long) left,(unsigned long long) genomicstart,(unsigned long long) genomicend,
		(unsigned long long) chrhigh,chrnum,querystart,queryend));

  new->substrings_LtoH = List_push(NULL,(void *) substring);
  new->substrings_1toN = List_push(NULL,(void *) substring);
  new->substrings_Nto1 = List_push(NULL,(void *) substring);

  new->junctions_LtoH = (List_T) NULL;
  new->junctions_1toN = (List_T) NULL;
  new->junctions_Nto1 = (List_T) NULL;


  new->pairarray = (struct Pair_T *) NULL;
  new->cigar_tokens = (List_T) NULL;
  new->gmap_intronp = false;

  new->querylength_adj = new->querylength = querylength;
  new->genomicstart = genomicstart;
  new->genomicend = genomicend;

  if (genomicstart < genomicend) {
    new->low = genomicstart;
    new->high = genomicend;
  } else {
    new->low = genomicend;
    new->high = genomicstart;
  }
  new->genomiclength = new->high - new->low;
  new->guided_insertlength = 0U;


  new->hittype = TERMINAL;
  new->genestrand = genestrand;
  new->sarrayp = sarrayp;
  new->gmap_source = GMAP_NOT_APPLICABLE;
  new->improved_by_gmap_p = false;

  new->chrnum = new->effective_chrnum = chrnum;
  new->other_chrnum = 0;
  new->chroffset = chroffset;
  new->chrhigh = chrhigh;
  new->chrlength = chrlength;
  new->plusp = plusp;

#if 0
  new->mapq_loglik = Substring_mapq_loglik(substring);
  new->mapq_score = 0;
  new->absmq_score = 0;
#endif

  new->nindels = 0;
  new->nmismatches_whole = Substring_nmismatches_whole(substring); /* This value was recomputed to include non-terminal end */
  new->ntscore = /* terminal_penalty + */ nmismatches_whole;

  new->nmismatches_bothdiff = Substring_nmismatches_bothdiff(substring);
  /* new->nmismatches_refdiff = Substring_nmismatches_refdiff(substring); */

#if 0
  new->score = terminal_penalty + nmismatches_whole;
#else
  new->score = /* terminal_penalty + */ Substring_nmismatches_whole(substring);
#endif
  new->nsegments = 1;

#if 0
  new->nmatches = Substring_match_length(substring) - nmismatches;
#else
  new->nmatches = Substring_nmatches(substring);
  new->nmatches_posttrim = Substring_nmatches_posttrim(substring);
#endif

  new->trim_left = Substring_trim_left(substring);
  new->trim_right = Substring_trim_right(substring);
  new->trim_left_splicep = Substring_trim_left_splicep(substring);
  new->trim_right_splicep = Substring_trim_right_splicep(substring);

  /* new->penalties = 0; */

  /* new->gene_overlap = NO_KNOWN_GENE; -- initialized later when resolving multimappers */
  new->tally = -1L;

  new->nsplices = 0;

  new->distance = 0U;
  new->shortexonA_distance = new->shortexonD_distance = 0U;
  new->sensedir = SENSE_NULL;

  new->paired_usedp = false;
  new->paired_seenp = false;
  new->concordantp = false;

  new->circularpos = compute_circularpos(&new->alias,new);
  if (new->alias == +2 || new->alias == -2) {
    Stage3end_free(&new);
    return (T) NULL;
  } else {
    return new;
  }
}


T
Stage3end_new_gmap (int nmismatches_whole, int nmatches_posttrim, int max_match_length,
		    int ambig_end_length_5, int ambig_end_length_3,
		    Splicetype_T ambig_splicetype_5, Splicetype_T ambig_splicetype_3,
		    double ambig_prob_5, double ambig_prob_3, double min_splice_prob,
		    struct Pair_T *pairarray, int npairs, int nsegments, int nintrons, int nindelbreaks,
		    Univcoord_T left, int genomiclength, bool plusp, int genestrand, bool first_read_p,
		    char *accession, int querylength, Chrnum_T chrnum,
		    Univcoord_T chroffset, Univcoord_T chrhigh, Chrpos_T chrlength,
		    int cdna_direction, int sensedir, GMAP_source_T gmap_source) {
  T new;
  Univcoord_T genomicstart, genomicend, genomepos;
  double prob1, prob2;

  Pair_T start, end;
  List_T cigar_tokens;
  bool intronp;
  int hardclip_start, hardclip_end;

  /* In 2012-12-20, removed statements to return NULL, because GMAP alignments seem
     to be okay, at least when starting before coordinate 0 */
  /* Example (when aligned to chrM at beginning of genome) (actually aligns circularly):
     GGATGAGGCAGGAATCAAAGACAGATACTGCGACATAGGGTGCTCCGGCTCCAGCGTCTCGCAATGCTATCGCGTG
     ATAGCCCACACGTTCCCCTTAAATAAGACATCACGATGGATCACAGGTCTATCACCCTATTAACCACTCACGGGAG
  */
  /* However, this leads to fatal bugs later, so restored these statements */

  debug0(printf("Entered Stage3end_new_gmap with first_read_p %d and sensedir %d\n",first_read_p,sensedir));

  start = &(pairarray[0]);
  end = &(pairarray[npairs-1]);
  hardclip_start = start->querypos;
  hardclip_end = (querylength - 1) - end->querypos;

  cigar_tokens = Pair_compute_cigar(&intronp,&hardclip_start,&hardclip_end,pairarray,npairs,querylength,
				    /*watsonp*/plusp,sensedir,/*chimera_part*/0);
  if (Pair_tokens_cigarlength(cigar_tokens) + hardclip_start + hardclip_end != querylength) {
    fprintf(stderr,"Could not compute a valid cigar for %s: %d + %d + %d != %d\n",
	    accession,Pair_tokens_cigarlength(cigar_tokens),hardclip_start,hardclip_end,querylength);
    Pair_dump_array_stderr(pairarray,npairs,/*zerobasedp*/true);
    Pair_tokens_free(&cigar_tokens);
#ifdef CHECK_ASSERTIONS
    abort();
#endif
    return (T) NULL;

  } else if (Stage3_bad_stretch_p(pairarray,npairs,/*pos5*/0,/*pos3*/querylength) == true) {
    debug0(printf("Bad GMAP: bad stretch\n"));
    Pair_tokens_free(&cigar_tokens);
    return (T) NULL;

  } else if (plusp == true) {
    genomicstart = left;
    if ((genomicend = left + genomiclength) > chrhigh) {
      Pair_tokens_free(&cigar_tokens);
      return (T) NULL;
    }
    if (genomicstart > genomicend) {
      /* Must have started before coordinate 0 */
      debug0(printf("plusp and genomicstart %llu > genomicend %llu => started before coordinate 0\n",
		    (unsigned long long) genomicstart,(unsigned long long) genomicend));
      Pair_tokens_free(&cigar_tokens);
      return (T) NULL;
    }
  } else {
    if ((genomicstart = left + genomiclength) > chrhigh) {
      Pair_tokens_free(&cigar_tokens);
      return (T) NULL;
    }
    genomicend = left;
    if (genomicend > genomicstart) {
      /* Must have started before coordinate 0 */
      debug0(printf("minusp and genomicend %llu > genomicstart %llu => started before coordinate 0\n",
		    (unsigned long long) genomicend,(unsigned long long) genomicstart));
      Pair_tokens_free(&cigar_tokens);
      return (T) NULL;
    }
  }

  new = (T) MALLOC_OUT(sizeof(*new));

  debug0(printf("Stage3end_new_gmap %p: left %llu, genomicstart/end %u..%u, chrhigh %llu, chrnum %d, nmismatches %d, cdna_direction %d, sensedir %d, max_match_length %d, gmap_source %d\n",
		new,(unsigned long long) left,(unsigned int) (genomicstart - chroffset),(unsigned int) (genomicend - chroffset),
		(unsigned long long) chrhigh,chrnum,nmismatches_whole,cdna_direction,sensedir,max_match_length,gmap_source));
  debug0(printf("  ambig_end_length_5 %d (prob %f), ambig_end_length_3 %d (prob %f)\n",ambig_end_length_5,ambig_prob_5,ambig_end_length_3,ambig_prob_3));

  new->substrings_LtoH = (List_T) NULL;
  new->substrings_1toN = (List_T) NULL;
  new->substrings_Nto1 = (List_T) NULL;

  new->junctions_LtoH = (List_T) NULL;
  new->junctions_1toN = (List_T) NULL;
  new->junctions_Nto1 = (List_T) NULL;


  new->pairarray = pairarray;
  new->npairs = npairs;
  new->cigar_tokens = cigar_tokens;
  new->gmap_intronp = intronp;
  new->nsegments = nsegments;

  new->querylength_adj = new->querylength = querylength /* - nindels */;
  new->genomicstart = genomicstart;
  new->genomicend = genomicend;

  if (genomicstart < genomicend) {
    new->low = genomicstart;
    new->high = genomicend;
  } else {
    new->low = genomicend;
    new->high = genomicstart;
  }
  new->genomiclength = new->high - new->low;
  new->guided_insertlength = 0U;


  new->hittype = GMAP;
  new->genestrand = genestrand;
  new->sarrayp = false;
  new->gmap_source = gmap_source;
  new->improved_by_gmap_p = false;

  new->chrnum = new->effective_chrnum = chrnum;
  new->other_chrnum = 0;
  new->chroffset = chroffset;
  new->chrhigh = chrhigh;
  new->chrlength = chrlength;
  new->plusp = plusp;
  new->gmap_nindelbreaks = nindelbreaks;
  new->gmap_cdna_direction = cdna_direction;
  new->gmap_nintrons = nintrons;
  new->sensedir = sensedir;

#if 0
  new->mapq_loglik = Substring_mapq_loglik(substring);
  new->mapq_score = 0;
  new->absmq_score = 0;
#endif

  new->nindels = 0;
  new->nmismatches_whole = nmismatches_whole;
  new->ntscore = nmismatches_whole;

#if 1
  /* This favors the trimmed results */
  new->score = nmismatches_whole;
  new->score += localsplicing_penalty * nintrons;
  new->score += indel_penalty_middle * nindelbreaks;
  debug0(printf("gmap score = %d = %d + %d*%d + %d*%d\n",
		new->score,nmismatches_whole,localsplicing_penalty,nintrons,indel_penalty_middle,nindelbreaks));
#else
  /* This is a better way to score GMAP.  Using nmatches_pretrim puts all GMAP entries on an even level. */
  new->score = querylength - nmatches_posttrim;
  if (nindels > 0) {
    /* Account for the fact that a query insertion reduces number of possible posttrim matches */
    new->score -= nindels;
  }
  new->score += localsplicing_penalty * nintrons;
  new->score += indel_penalty_middle * nindelbreaks;
  debug0(printf("gmap score = %d = querylength %d (nindels %d, subtract if pos) - posttrim %d + %d*%d + %d*%d\n",
		new->score,querylength,nindels,nmatches_posttrim,
		localsplicing_penalty,nintrons,indel_penalty_middle,nindelbreaks));
#endif


  new->nmismatches_bothdiff = nmismatches_whole;
  /* new->nmismatches_bothdiff = Substring_nmismatches_bothdiff(new->substring1); */
  /* new->nmismatches_refdiff = Substring_nmismatches_refdiff(new->substring1); */

  /* Adding ambig_end_lengths to nmatches_posttrim would unnecessarily favor long ambig ends when comparing GMAP results */
  new->nmatches = nmatches_posttrim; /* To make addition of ambiguous lengths work, we need to use posttrim, not pretrim */
  new->nmatches_posttrim = nmatches_posttrim;
  if (favor_ambiguous_p == true) {
    new->nmatches += ambig_end_length_5 + ambig_end_length_3;
  }
  debug0(printf("  nmatches %d = posttrim %d + ambig_end_length_5 %d + ambig_end_length_3 %d\n",
		new->nmatches,nmatches_posttrim,ambig_end_length_5,ambig_end_length_3));

  new->nmatches_posttrim -= localsplicing_penalty * nintrons; /* for use in goodness_cmp procedures */
  new->nmatches_posttrim -= indel_penalty_middle * nindelbreaks; /* for use in goodness_cmp procedures */

  if (new->nmatches_posttrim < querylength/2) {
    debug0(printf("  nmatches %d < querylength %d/2, so returning NULL\n",
		  new->nmatches_posttrim,querylength));
    Pair_tokens_free(&cigar_tokens);
    FREE_OUT(new);
    return (T) NULL;
  } else if (max_match_length < gmap_min_nconsecutive) {
    debug0(printf("  max_match_length %d < %d, so returning NULL\n",max_match_length,gmap_min_nconsecutive));
    Pair_tokens_free(&cigar_tokens);
    FREE_OUT(new);
    return (T) NULL;
  }

  new->gmap_max_match_length = max_match_length;
  new->gmap_min_splice_prob = min_splice_prob;


  new->trim_left = Pair_querypos(&(pairarray[0])) - ambig_end_length_5;
  if ((new->gmap_start_amb_length = ambig_end_length_5) > 0) {
    new->trim_left_splicep = true;
  } else if (novelsplicingp == false) {
    new->trim_left_splicep = false;
  } else {
    genomepos = chroffset + Pair_genomepos(&(pairarray[0])) + 1U;
    if (plusp == true) {
      prob1 = Maxent_hr_acceptor_prob(genomepos,chroffset);
      prob2 = Maxent_hr_antidonor_prob(genomepos,chroffset);
      /* printf("At %llu, acceptor prob %f, antidonor prob %f\n",(unsigned long long) genomepos,prob1,prob2); */
    } else {
      prob1 = Maxent_hr_donor_prob(genomepos,chroffset);
      prob2 = Maxent_hr_antiacceptor_prob(genomepos,chroffset);
      /* printf("At %llu, donor prob %f, antiacceptor prob %f\n",(unsigned long long) genomepos,prob1,prob2); */
    }
    if (prob1 > 0.90 || prob2 > 0.90) {
      new->trim_left_splicep = true;
    } else {
      new->trim_left_splicep = false;
    }
  }

  new->trim_right = (querylength - 1) - Pair_querypos(&(pairarray[npairs-1])) - ambig_end_length_3;
  if ((new->gmap_end_amb_length = ambig_end_length_3) > 0) {
    new->trim_right_splicep = true;
  } else if (novelsplicingp == false) {
    new->trim_right_splicep = false;
  } else {
    genomepos = chroffset + Pair_genomepos(&(pairarray[npairs-1])) + 1U;
    if (plusp == true) {
      prob1 = Maxent_hr_donor_prob(genomepos,chroffset);
      prob2 = Maxent_hr_antiacceptor_prob(genomepos,chroffset);
      /* printf("At %llu, donor prob %f, antiacceptor prob %f\n",(unsigned long long) genomepos,prob1,prob2); */
    } else {
      prob1 = Maxent_hr_acceptor_prob(genomepos,chroffset);
      prob2 = Maxent_hr_antidonor_prob(genomepos,chroffset);
      /* printf("At %llu, acceptor prob %f, antidonor prob %f\n",(unsigned long long) genomepos,prob1,prob2); */
    }
    if (prob1 > 0.90 || prob2 > 0.90) {
      new->trim_right_splicep = true;
    } else {
      new->trim_right_splicep = false;
    }
  }

#if 0
  /* new->penalties not used anyway for GMAP alignments */
#ifdef SCORE_INDELS
  /* indel_penalty will be counted later */
  new->penalties = localsplicing_penalty * nintrons;
#else
  new->penalties = localsplicing_penalty * nintrons + indel_penalty_middle * nindelbreaks;
#endif
  /* new->penalties += ambig_end_length_5/ambig_end_interval; */
  /* new->penalties += ambig_end_length_3/ambig_end_interval; */
#endif

  /* new->gene_overlap = NO_KNOWN_GENE; -- initialized later when resolving multimappers */
  new->tally = -1L;

  if (ambig_end_length_5 == 0) {
    new->gmap_start_endtype = END;
  } else if (ambig_splicetype_5 == DONOR || ambig_splicetype_5 == ANTIDONOR) {
    new->gmap_start_endtype = AMB_DON;
  } else if (ambig_splicetype_5 == ACCEPTOR || ambig_splicetype_5 == ANTIACCEPTOR) {
    new->gmap_start_endtype = AMB_ACC;
  } else {
    fprintf(stderr,"Do not recognize splicetype %d for ambig_end_length_5 %d\n",
	    ambig_splicetype_5,ambig_end_length_5);
    abort();
  }
  /* new->start_amb_prob = ambig_prob_5; */
    
  if (ambig_end_length_3 == 0) {
    new->gmap_end_endtype = END;
  } else if (ambig_splicetype_3 == DONOR || ambig_splicetype_3 == ANTIDONOR) {
    new->gmap_end_endtype = AMB_DON;
  } else if (ambig_splicetype_3 == ACCEPTOR || ambig_splicetype_3 == ANTIACCEPTOR) {
    new->gmap_end_endtype = AMB_ACC;
  } else {
    fprintf(stderr,"Do not recognize splicetype %d for ambig_end_length_3 %d\n",
	    ambig_splicetype_3,ambig_end_length_3);
    abort();
  }
  /* new->end_amb_prob = ambig_prob_3; */

  new->nsplices = nintrons;
  if (ambig_end_length_5 > 0) {
    new->nsplices += 1;
  }
  if (ambig_end_length_3 > 0) {
    new->nsplices += 1;
  }

  new->distance = 0U;
  new->shortexonA_distance = new->shortexonD_distance = 0;

  new->paired_usedp = false;
  new->paired_seenp = false;
  new->concordantp = false;

  new->circularpos = compute_circularpos(&new->alias,new);
  if (new->alias == +2 || new->alias == -2) {
    /* Stage3end_free(&new); -- Cannot use, because it frees pairarray */
    Pair_tokens_free(&new->cigar_tokens);
    /* No substrings or junctions */
    FREE(new);
    return (T) NULL;
  } else {
    return new;
  }
}


static int
Stage3end_output_cmp (const void *a, const void *b) {
  T x = * (T *) a;
  T y = * (T *) b;

  if (x->nmatches > y->nmatches) {
    return -1;
  } else if (y->nmatches > x->nmatches) {
    return +1;
  } else if (x->mapq_loglik > y->mapq_loglik) {
    return -1;
  } else if (y->mapq_loglik > x->mapq_loglik) {
    return +1;
  } else if (x->guided_insertlength > 0 && y->guided_insertlength == 0) {
    return -1;
  } else if (y->guided_insertlength > 0 && x->guided_insertlength == 0) {
    return +1;
  } else if (x->guided_insertlength < y->guided_insertlength) {
    return -1;
  } else if (y->guided_insertlength < x->guided_insertlength) {
    return +1;
  } else if (x->nmatches_posttrim > y->nmatches_posttrim) {
    return -1;
  } else if (y->nmatches_posttrim > x->nmatches_posttrim) {
    return +1;
  } else if (x->score < y->score) {
    return -1;
  } else if (y->score < x->score) {
    return +1;
#if 0
  } else if (x->hittype < y->hittype) {
    return -1;
  } else if (y->hittype < x->hittype) {
    return +1;
#endif

    /* This genomic ordering will be undone if want_random_p is true */
  } else if (x->genomicstart < y->genomicstart) {
    return -1;
  } else if (y->genomicstart < x->genomicstart) {
    return +1;
  } else if (x->plusp == true && y->plusp == false) {
    return -1;
  } else if (x->plusp == false && y->plusp == true) {
    return +1;

  } else {
    return 0;
  }
}


static int
Stage3pair_output_cmp (const void *a, const void *b) {
  Stage3pair_T x = * (Stage3pair_T *) a;
  Stage3pair_T y = * (Stage3pair_T *) b;

#ifdef USE_BINGO
  if (x->absdifflength_bingo_p == true && y->absdifflength_bingo_p == false) {
    return -1;
  } else if (y->absdifflength_bingo_p == true && x->absdifflength_bingo_p == false) {
    return +1;
  }
#endif

  if (x->nmatches > y->nmatches) {
    return -1;
  } else if (y->nmatches > x->nmatches) {
    return +1;
  } else if (x->mapq_loglik > y->mapq_loglik) {
    return -1;
  } else if (y->mapq_loglik > x->mapq_loglik) {
    return +1;
  } else if (x->insertlength > 0 && y->insertlength == 0) {
    return -1;
  } else if (y->insertlength > 0 && x->insertlength == 0) {
    return +1;
  } else if (x->insertlength < y->insertlength) {
    return -1;
  } else if (y->insertlength < x->insertlength) {
    return +1;
  } else if (x->nmatches_posttrim > y->nmatches_posttrim) {
    return -1;
  } else if (y->nmatches_posttrim > x->nmatches_posttrim) {
    return +1;
  } else if (x->score < y->score) {
    return -1;
  } else if (y->score < x->score) {
    return +1;

    /* This genomic ordering will be undone if want_random_p is true */
  } else if (x->low < y->low) {
    return -1;
  } else if (y->low < x->low) {
    return +1;

  } else {
    return 0;
  }
}



static float
Stage3end_compute_mapq (Stage3end_T this, Compress_T query_compress_fwd, Compress_T query_compress_rev,
			char *quality_string, bool trim_terminals_p) {
  List_T p;
  Substring_T substring;

  if (this == NULL) {
    return 0.0;

  } else if (this->hittype == GMAP) {
    this->mapq_loglik = Pair_compute_mapq(this->pairarray,this->npairs,
					  this->trim_left,this->trim_right,this->querylength,
					  quality_string,trim_terminals_p);

  } else if (this->plusp == true) {
    this->mapq_loglik = 0.0;
    for (p = this->substrings_LtoH; p != NULL; p = List_next(p)) {
      substring = (Substring_T) List_head(p);
      this->mapq_loglik +=
	Substring_compute_mapq(substring,query_compress_fwd,quality_string,trim_terminals_p);
    }

  } else {
    this->mapq_loglik = 0.0;
    for (p = this->substrings_LtoH; p != NULL; p = List_next(p)) {
      substring = (Substring_T) List_head(p);
      this->mapq_loglik +=
	Substring_compute_mapq(substring,query_compress_rev,quality_string,trim_terminals_p);
    }
  }

  return this->mapq_loglik;
}



static void
Stage3end_display_prep (Stage3end_T this, char *queryuc_ptr, char *queryrc,
			Compress_T query_compress_fwd, Compress_T query_compress_rev,
			int amb_resolve, bool first_read_p) {
  List_T p, q;
  Substring_T substring, anchor;
  Junction_T pre_junction, post_junction, junction;
  int extraleft, extraright, extralow, extrahigh;
  Univcoord_T left, ignore;
  double donor_prob, acceptor_prob;
  int type;


  if (this != NULL) {
    debug0(printf("Doing a display prep of end %p\n",this));
    if (this->hittype == GMAP) {
      this->nmismatches_refdiff = this->nmismatches_bothdiff;

    } else {
      /* Resolve ambiguous end */
      if (amb_resolve >= 0) {
	if (first_read_p == true) {
	  substring = (Substring_T) List_head(this->substrings_Nto1);
	  anchor = (Substring_T) List_head(List_next(this->substrings_Nto1));
	  junction = (Junction_T) List_head(this->junctions_Nto1);
	  left = Substring_set_unambiguous(&donor_prob,&acceptor_prob,&ignore,&this->genomicend,substring,amb_resolve);
	  if (this->plusp == true) {
	    Junction_set_unambiguous(junction,left - Substring_left(anchor),donor_prob,acceptor_prob);
	  } else {
	    Junction_set_unambiguous(junction,Substring_left(anchor) - left,donor_prob,acceptor_prob);
	  }

	} else {
	  substring = (Substring_T) List_head(this->substrings_1toN);
	  anchor = (Substring_T) List_head(List_next(this->substrings_1toN));
	  junction = (Junction_T) List_head(this->junctions_1toN);
	  left = Substring_set_unambiguous(&donor_prob,&acceptor_prob,&this->genomicstart,&ignore,substring,amb_resolve);
	  if (this->plusp == true) {
	    Junction_set_unambiguous(junction,Substring_left(anchor) - left,donor_prob,acceptor_prob);
	  } else {
	    Junction_set_unambiguous(junction,left - Substring_left(anchor),donor_prob,acceptor_prob);
	  }
	}
      }

      this->nmismatches_refdiff = 0;

      /* First segments */
      /* For operations on substrings, proceed in 1toN order, not LtoH order */
      substring = (Substring_T) List_head(this->substrings_1toN);
      extraleft = Substring_querystart(substring); /* terminal start */

      if (List_length(this->substrings_1toN) == 1) {
	post_junction = (Junction_T) NULL;
	extraright = this->querylength - Substring_queryend(substring); /* terminal end */
      } else {
	post_junction = (Junction_T) List_head(this->junctions_1toN);
	if (Junction_type(post_junction) == SPLICE_JUNCTION) {
	  extraright = 2;
	} else {
	  extraright = 0;
	}
      }
      
      if (Substring_ambiguous_p(substring) == true) {
      } else {
	this->nmismatches_refdiff += 
	  Substring_display_prep(substring,queryuc_ptr,this->querylength,
				 extraleft,extraright,query_compress_fwd,query_compress_rev,
				 genome);
      }

      if ((p = List_next(this->substrings_1toN)) == NULL) {
	/* No middle segments */
      } else {
	for (q = List_next(this->junctions_1toN); q != NULL; p = List_next(p), q = List_next(q)) {
	  /* Middle segments */
	  pre_junction = post_junction;
	  post_junction = List_head(q);
#if 0
	  extraleft = 0;
	  if ((type = Junction_type(pre_junction)) == INS_JUNCTION) {
	    ninsertions += Junction_nindels(pre_junction);
	  } else if (type == DEL_JUNCTION) {
	    ndeletions += Junction_nindels(pre_junction);
	  } else if (type == SPLICE_JUNCTION) {
	    extraleft = 2;
	  }
#else
	  if (Junction_type(pre_junction) == SPLICE_JUNCTION) {
	    extraleft = 2;
	  } else {
	    extraleft = 0;
	  }
#endif
	  if (Junction_type(post_junction) == SPLICE_JUNCTION) {
	    extraright = 2;
	  } else {
	    extraright = 0;
	  }

	  substring = (Substring_T) List_head(p);
	  if (Substring_ambiguous_p(substring) == true) {
	    /* Skip */
	  } else {
	    this->nmismatches_refdiff += 
	      Substring_display_prep(substring,queryuc_ptr,this->querylength,
				     extraleft,extraright,query_compress_fwd,query_compress_rev,
				     genome);
	  }
	}

	/* Last segment */
	pre_junction = post_junction;
#if 0
	extraleft = 0;
	if ((type = Junction_type(pre_junction)) == INS_JUNCTION) {
	  ninsertions += Junction_nindels(pre_junction);
	} else if (type == DEL_JUNCTION) {
	  ndeletions += Junction_nindels(pre_junction);
	} else if (type == SPLICE_JUNCTION) {
	  extraleft = 2;
	}
#else
	if (Junction_type(pre_junction) == SPLICE_JUNCTION) {
	  extraleft = 2;
	} else {
	  extraleft = 0;
	}
#endif
	substring = (Substring_T) List_head(p);
	extraright = this->querylength - Substring_queryend(substring);
	
	if (Substring_ambiguous_p(substring) == true) {
	  /* Skip */
	} else {
	  this->nmismatches_refdiff += 
	    Substring_display_prep(substring,queryuc_ptr,this->querylength,
				   extraleft,extraright,query_compress_fwd,query_compress_rev,
				   genome);
	}
      }
    }
  }
  return;
}


static int
end_matches_cmp (const void *a, const void *b) {
  T x = * (T *) a;
  T y = * (T *) b;

  if (x->nmatches > y->nmatches) {
    return -1;
  } else if (x->nmatches < y->nmatches) {
    return +1;
  } else {
    return 0;
  }
}

List_T
Stage3end_sort_bymatches (List_T hits) {
  List_T sorted = NULL;
  T *array;
  int n, i;

  if ((n = List_length(hits)) == 0) {
    return (List_T) NULL;
  } else {
#ifdef USE_ALLOCA_FOR_HITS
    array = (T *) MALLOCA(n * sizeof(T));
    List_fill_array_and_free((void **) array,&hits);
#else
    array = (T *) List_to_array(hits,NULL);
    List_free(&hits);
#endif

    qsort(array,n,sizeof(T),end_matches_cmp);
    for (i = n-1; i >= 0; i--) {
      sorted = List_push(sorted,(void *) array[i]);
    }
#ifdef USE_ALLOCA_FOR_HITS
    FREEA(array);
#else
    FREE(array);
#endif

    return sorted;
  }
}


/* Need to include criteria from end_matches_cmp to work on most likely terminals */
static int
paired_seenp_cmp (const void *a, const void *b) {
  T x = * (T *) a;
  T y = * (T *) b;

  if (x->paired_seenp > y->paired_seenp) {
    return -1;
  } else if (y->paired_seenp > x->paired_seenp) {
    return +1;
  } else if (x->nmatches > y->nmatches) {
    return -1;
  } else if (x->nmatches < y->nmatches) {
    return +1;
  } else {
    return 0;
  }
}

List_T
Stage3end_sort_by_paired_seenp (List_T hits) {
  List_T sorted = NULL;
  T *array;
  int n, i;

  if ((n = List_length(hits)) == 0) {
    return (List_T) NULL;
  } else {
#ifdef USE_ALLOCA_FOR_HITS
    array = (T *) MALLOCA(n * sizeof(T));
    List_fill_array_and_free((void **) array,&hits);
#else
    array = (T *) List_to_array(hits,NULL);
    List_free(&hits);
#endif

    qsort(array,n,sizeof(T),paired_seenp_cmp);
    for (i = n-1; i >= 0; i--) {
      sorted = List_push(sorted,(void *) array[i]);
    }
#ifdef USE_ALLOCA_FOR_HITS
    FREEA(array);
#else
    FREE(array);
#endif

    return sorted;
  }
}



List_T
Stage3end_filter_coverage (List_T hits, int min_coverage) {
  List_T newhits = NULL, p;
  Stage3end_T hit;

  for (p = hits; p != NULL; p = List_next(p)) {
    hit = (Stage3end_T) List_head(p);
    if (hit->querylength - hit->trim_left - hit->trim_right >= min_coverage) {
      newhits = List_push(newhits,(void *) hit);
    } else {
      Stage3end_free(&hit);
    }
  }

  List_free(&hits);
  return newhits;
}




Stage3end_T *
Stage3end_eval_and_sort (int *npaths, int *first_absmq, int *second_absmq,
			 Stage3end_T *stage3array, int maxpaths, Shortread_T queryseq,
			 char *queryuc_ptr, char *queryrc,
			 Compress_T query_compress_fwd, Compress_T query_compress_rev,
			 char *quality_string, bool displayp) {
  float maxlik, loglik;
  float total, q;		/* For Bayesian mapq calculation */
  int compute_npaths;
  bool non_terminal_p;

  int randomi, i;
  Stage3end_T temp;

  if (*npaths == 0) {
    /* Skip */
    *first_absmq = 0;
    *second_absmq = 0;

  } else if (*npaths == 1) {
    stage3array[0]->mapq_loglik = MAPQ_MAXIMUM_SCORE;
    stage3array[0]->mapq_score = 
      MAPQ_max_quality_score(Shortread_quality_string(queryseq),Shortread_fulllength(queryseq));
    stage3array[0]->absmq_score = MAPQ_MAXIMUM_SCORE;

    if (displayp == true) {
      Stage3end_display_prep(stage3array[0],queryuc_ptr,queryrc,query_compress_fwd,query_compress_rev,
			     /*amb_resolve*/-1,/*first_read_p*/true);
    }
    *first_absmq = stage3array[0]->absmq_score;
    *second_absmq = 0;

  } else {
    /* Determine whether to trim terminal ends */
    non_terminal_p = false;
    for (i = 0; i < *npaths; i++) {
      if (stage3array[i]->hittype != TERMINAL) {
	non_terminal_p = true;
      }
    }

    /* Compute mapq_loglik */
    for (i = 0; i < *npaths; i++) {
      Stage3end_compute_mapq(stage3array[i],query_compress_fwd,query_compress_rev,
			     quality_string,/*trim_terminals_p*/non_terminal_p ? false : true);
    }

    /* Sort by nmatches, then mapq */
    qsort(stage3array,*npaths,sizeof(Stage3end_T),Stage3end_output_cmp);

    if (want_random_p) {
      /* Randomize among best alignments */
      i = 1;
      while (i < *npaths && Stage3end_output_cmp(&(stage3array[i]),&(stage3array[0])) == 0) {
	i++;
      }
      if (i > 1) {		/* i is number of ties */
	/* randomi = (int) ((double) i * rand()/((double) RAND_MAX + 1.0)); */
	randomi = (int) (rand() / (((double) RAND_MAX + 1.0) / (double) i));
	/* fprintf(stderr,"%d dups => random %d\n",i,randomi); */
	temp = stage3array[0];
	stage3array[0] = stage3array[randomi];
	stage3array[randomi] = temp;
      }
    }

    /* Enforce monotonicity */
    for (i = *npaths - 1; i > 0; i--) {
      if (stage3array[i-1]->mapq_loglik < stage3array[i]->mapq_loglik) {
	stage3array[i-1]->mapq_loglik = stage3array[i]->mapq_loglik;
      }
    }
    maxlik = stage3array[0]->mapq_loglik;
    
    /* Subtract maxlik to avoid underflow */
    for (i = 0; i < *npaths; i++) {
      stage3array[i]->mapq_loglik -= maxlik;
    }

    /* Save on computation if possible */
    if (*npaths < maxpaths) {
      compute_npaths = *npaths;
    } else {
      compute_npaths = maxpaths;
    }
    if (compute_npaths < 2) {
      compute_npaths = 2;
    }

    /* Compute absolute mapq */
    for (i = 0; i < compute_npaths; i++) {
      loglik = stage3array[i]->mapq_loglik + MAPQ_MAXIMUM_SCORE;
      if (loglik < 0.0) {
	loglik = 0.0;
      }
      stage3array[i]->absmq_score = rint(loglik);
    }
    *first_absmq = stage3array[0]->absmq_score;
    *second_absmq = stage3array[1]->absmq_score;


    /* Compute Bayesian mapq */
    total = 0.0;
    for (i = 0; i < *npaths; i++) {
      total += (stage3array[i]->mapq_loglik = fasterexp(stage3array[i]->mapq_loglik));
    }

    /* Obtain posterior probabilities of being true */
    for (i = 0; i < compute_npaths; i++) {
      stage3array[i]->mapq_loglik /= total;
    }

    /* Convert to Phred scores */
    for (i = 0; i < compute_npaths; i++) {
      if ((q = 1.0 - stage3array[i]->mapq_loglik) < 2.5e-10 /* 10^-9.6 */) {
	stage3array[i]->mapq_score = 96;
      } else {
	stage3array[i]->mapq_score = rint(-10.0 * log10(q));
      }
    }

    if (displayp == true) {
      /* Prepare for display */
      for (i = 0; i < compute_npaths; i++) {
	Stage3end_display_prep(stage3array[i],queryuc_ptr,queryrc,query_compress_fwd,query_compress_rev,
			       /*amb_resolve*/-1,/*first_read_p*/true);
      }
    }

#if 0
    /* Apply filtering for mapq unique -- currently not used since mapq_unique_score is high */
    if (stage3array[0]->mapq_score >= mapq_unique_score &&
	stage3array[1]->mapq_score < mapq_unique_score) {
      for (i = 1; i < *npaths; i++) {
	Stage3end_free(&(stage3array[i]));
      }
      *npaths = 1;
    }
#endif
  }

  return stage3array;
}


static int
insertlength_expected (int insertlength) {
  if (insertlength < expected_pairlength_low) {
    return -1;
  } else if (insertlength > expected_pairlength_very_high) {
    return -1;
  } else if (insertlength > expected_pairlength_high) {
    return 0;
  } else {
    return +1;
  }
}


/* For concordant ends */
static Chrpos_T
pair_insert_length (Stage3end_T hit5, Stage3end_T hit3) {
  List_T p, q;
  Substring_T substring5, substring3;

  if (hit5->plusp != hit3->plusp) {
    debug10(printf("pair_insert_length: hit5->plusp %d != hit3->plusp %d, so returning 0\n",
		   hit5->plusp,hit3->plusp));
    return 0;
  }

  if (hit5->chrnum != 0 && hit3->chrnum != 0) {
    for (q = hit3->substrings_1toN; q != NULL; q = List_next(q)) {
      substring3 = (Substring_T) List_head(q);
      for (p = hit5->substrings_1toN; p != NULL; p = List_next(p)) {
	substring5 = (Substring_T) List_head(p);
	if (Substring_overlap_p(substring5,substring3)) {
	  return Substring_insert_length(substring5,substring3);
	}
      }
    }
  }

  /* No overlap found between any combination of substrings */
  if (hit5->plusp == true) {
    if (hit5->genomicend > hit3->genomicstart + hit5->querylength + hit3->querylength) {
      debug10(printf("pair_insert_length: no overlap found, and %u - %u + %d + %d < 0, so returning 0\n",
		     hit3->genomicstart - hit3->chroffset,hit5->genomicend - hit5->chroffset,
		     hit5->querylength,hit3->querylength));
      return 0;
    } else {
      debug10(printf("pair_insert_length: no overlap found, so returning %u - %u + %d + %d\n",
		     hit3->genomicstart - hit3->chroffset,hit5->genomicend - hit5->chroffset,
		     hit5->querylength,hit3->querylength));
    }
    return hit3->genomicstart - hit5->genomicend + hit5->querylength + hit3->querylength;
  } else {
    if (hit3->genomicstart > hit5->genomicend + hit5->querylength + hit3->querylength) {
      debug10(printf("pair_insert_length: no overlap found, and %u - %u + %d + %d < 0, so returning 0\n",
		     hit5->genomicend - hit5->chroffset,hit3->genomicstart - hit3->chroffset,
		     hit5->querylength,hit3->querylength));
      return 0;
    } else {
      debug10(printf("pair_insert_length: no overlap found, so returning %u - %u + %d + %d\n",
		     hit5->genomicend - hit5->chroffset,hit3->genomicstart - hit3->chroffset,
		     hit5->querylength,hit3->querylength));
      return hit5->genomicend - hit3->genomicstart + hit5->querylength + hit3->querylength;
    }
  }
}



/* For unpaired ends */
static Chrpos_T
pair_insert_length_unpaired (Stage3end_T hit5, Stage3end_T hit3) {

  if (hit5->chrnum != hit3->chrnum) {
    debug10(printf("pair_insert_length: hit5->plusp %d != hit3->plusp %d, so returning 0\n",
		   hit5->plusp,hit3->plusp));
    return 0;
  } else if (hit5->high < hit3->low) {
    return hit3->low - hit5->high + hit5->querylength + hit3->querylength;
  } else if (hit3->high < hit5->low) {
    return hit5->low - hit3->high + hit5->querylength + hit3->querylength;
  } else {
    return hit5->querylength + hit3->querylength;
  }
}


Stage3end_T *
Stage3end_eval_and_sort_guided (int *npaths, int *first_absmq, int *second_absmq, Stage3end_T guide,
				Stage3end_T *stage3array, int maxpaths, Shortread_T queryseq,
				char *queryuc_ptr, char *queryrc,
				Compress_T query_compress_fwd, Compress_T query_compress_rev,
				char *quality_string, bool displayp) {
  float maxlik, loglik;
  float total, q;		/* For Bayesian mapq calculation */
  int compute_npaths;
  bool non_terminal_p;

  int randomi, i;
  Stage3end_T temp;

  if (*npaths == 0) {
    /* Skip */
    *first_absmq = 0;
    *second_absmq = 0;

  } else if (*npaths == 1) {
    stage3array[0]->mapq_loglik = MAPQ_MAXIMUM_SCORE;
    stage3array[0]->mapq_score = 
      MAPQ_max_quality_score(Shortread_quality_string(queryseq),Shortread_fulllength(queryseq));
    stage3array[0]->absmq_score = MAPQ_MAXIMUM_SCORE;

    if (displayp == true) {
      Stage3end_display_prep(stage3array[0],queryuc_ptr,queryrc,query_compress_fwd,query_compress_rev,
			     /*amb_resolve*/-1,/*first_read_p*/true);
    }
    *first_absmq = stage3array[0]->absmq_score;
    *second_absmq = 0;

  } else {
    /* Determine whether to trim terminal ends */
    non_terminal_p = false;
    for (i = 0; i < *npaths; i++) {
      if (stage3array[i]->hittype != TERMINAL) {
	non_terminal_p = true;
      }
    }

    /* Compute mapq_loglik */
    for (i = 0; i < *npaths; i++) {
      Stage3end_compute_mapq(stage3array[i],query_compress_fwd,query_compress_rev,
			     quality_string,/*trim_terminals_p*/non_terminal_p ? false : true);
    }

    /* Compute insert_length relative to guide.  This is the only change from the unguided procedure. */
    for (i = 0; i < *npaths; i++) {
      stage3array[i]->guided_insertlength = pair_insert_length_unpaired(stage3array[i],guide);
    }

    /* Sort by nmatches, then mapq */
    qsort(stage3array,*npaths,sizeof(Stage3end_T),Stage3end_output_cmp);

    if (want_random_p) {
      /* Randomize among best alignments */
      i = 1;
      while (i < *npaths && Stage3end_output_cmp(&(stage3array[i]),&(stage3array[0])) == 0) {
	i++;
      }
      if (i > 1) {		/* i is number of ties */
	/* randomi = (int) ((double) i * rand()/((double) RAND_MAX + 1.0)); */
	randomi = (int) (rand() / (((double) RAND_MAX + 1.0) / (double) i));
	/* fprintf(stderr,"%d dups => random %d\n",i,randomi); */
	temp = stage3array[0];
	stage3array[0] = stage3array[randomi];
	stage3array[randomi] = temp;
      }
    }

    /* Enforce monotonicity */
    for (i = *npaths - 1; i > 0; i--) {
      if (stage3array[i-1]->mapq_loglik < stage3array[i]->mapq_loglik) {
	stage3array[i-1]->mapq_loglik = stage3array[i]->mapq_loglik;
      }
    }
    maxlik = stage3array[0]->mapq_loglik;
    
    /* Subtract maxlik to avoid underflow */
    for (i = 0; i < *npaths; i++) {
      stage3array[i]->mapq_loglik -= maxlik;
    }

    /* Save on computation if possible */
    if (*npaths < maxpaths) {
      compute_npaths = *npaths;
    } else {
      compute_npaths = maxpaths;
    }
    if (compute_npaths < 2) {
      compute_npaths = 2;
    }

    /* Compute absolute mapq */
    for (i = 0; i < compute_npaths; i++) {
      loglik = stage3array[i]->mapq_loglik + MAPQ_MAXIMUM_SCORE;
      if (loglik < 0.0) {
	loglik = 0.0;
      }
      stage3array[i]->absmq_score = rint(loglik);
    }
    *first_absmq = stage3array[0]->absmq_score;
    *second_absmq = stage3array[1]->absmq_score;


    /* Compute Bayesian mapq */
    total = 0.0;
    for (i = 0; i < *npaths; i++) {
      total += (stage3array[i]->mapq_loglik = fasterexp(stage3array[i]->mapq_loglik));
    }

    /* Obtain posterior probabilities of being true */
    for (i = 0; i < compute_npaths; i++) {
      stage3array[i]->mapq_loglik /= total;
    }

    /* Convert to Phred scores */
    for (i = 0; i < compute_npaths; i++) {
      if ((q = 1.0 - stage3array[i]->mapq_loglik) < 2.5e-10 /* 10^-9.6 */) {
	stage3array[i]->mapq_score = 96;
      } else {
	stage3array[i]->mapq_score = rint(-10.0 * log10(q));
      }
    }

    if (displayp == true) {
      /* Prepare for display */
      for (i = 0; i < compute_npaths; i++) {
	Stage3end_display_prep(stage3array[i],queryuc_ptr,queryrc,query_compress_fwd,query_compress_rev,
			       /*amb_resolve*/-1,/*first_read_p*/true);
      }
    }

#if 0
    /* Apply filtering for mapq unique -- currently not used since mapq_unique_score is high */
    if (stage3array[0]->mapq_score >= mapq_unique_score &&
	stage3array[1]->mapq_score < mapq_unique_score) {
      for (i = 1; i < *npaths; i++) {
	Stage3end_free(&(stage3array[i]));
      }
      *npaths = 1;
    }
#endif
  }

  return stage3array;
}


/* Note: single-end terminals can be present with non-terminals when
   paired-end reads are searched for concordance, which can accumulate
   terminal alignments */

/* Pre-final: max (max-terminal, min-other)
   Final: max (min-terminal, max-GMAP, min-other) */


static List_T
Stage3end_optimal_score_aux (bool *eliminatedp, List_T hitlist, int cutoff_level, int suboptimal_mismatches,
			     Compress_T query_compress_fwd, Compress_T query_compress_rev,
			     int querylength, bool keep_gmap_p, bool finalp) {
  List_T optimal = NULL, p, q;
  T hit;
  Substring_T substring;
  int n;
  int minscore = querylength;
  int max_nmatches = 0, max_nmatches_posttrim = 0;
  int trim_left = querylength, trim_right = querylength;
  int nindelbreaks;
  int best_nsegments;

#ifdef TRANSLOC_SPECIAL
  bool non_translocation_p = false;
#endif


  *eliminatedp = false;
  n = List_length(hitlist);
  debug4(printf("\nEntered Stage3end_optimal_score with %d hits: %s\n",
		n,finalp == true ? "FINAL" : "not final"));

  if (n <= 1) {
    return hitlist;
  }

  for (p = hitlist; p != NULL; p = p->rest) {
    hit = (T) p->first;
#ifdef TRANSLOC_SPECIAL
    if (hit->chrnum != 0) {
      non_translocation_p = true;
    }
#endif
#if 0
    if (hit->hittype == GMAP) {
      debug4(printf("Found gmap/terminal\n"));
      gmap_terminal_p = true;
    } else if (hit->hittype == TERMINAL) {
      debug4(printf("Found gmap/terminal\n"));
      gmap_terminal_p = true;
    } else {
      debug4(printf("Found a non-gmap/terminal\n"));
      non_gmap_terminal_p = true;
    }
#endif
  }


  /* Use eventrim for comparing alignments */
  for (p = hitlist; p != NULL; p = p->rest) {
    hit = (T) p->first;

    debug4(printf("hittype: %s, trim_left: %d, trim_right %d\n",
		  hittype_string(hit->hittype),hit->trim_left,hit->trim_right));
    if (hit->hittype == TERMINAL) {
      /* Don't allow terminals to set trims */

#if 0
    } else if ((hit->hittype == INSERTION || hit->hittype == DELETION) &&
	       (hit->indel_pos < 15 || hit->indel_pos > hit->querylength - 15)) {
      /* Don't allow end indels to set trims */
#endif

    } else {
      if (hit->trim_left_splicep == true) {
	/* Skip */
      } else if (hit->trim_left < trim_left) {
	trim_left = hit->trim_left;
      }
      if (hit->trim_right_splicep == true) {
	/* Skip */
      } else if (hit->trim_right < trim_right) {
	trim_right = hit->trim_right;
      }
    }
  }

  if (trim_left == querylength) {
    trim_left = 0;
  }
  if (trim_right == querylength) {
    trim_right = 0;
  }

  debug4(printf("trim_left: %d, trim_right %d\n",trim_left,trim_right));

  for (p = hitlist; p != NULL; p = p->rest) {
    hit = (T) p->first;

    if (hit->hittype == GMAP) {
      debug4(printf("score GMAP:"));
      hit->score_eventrim = Pair_nmismatches_region(&nindelbreaks,hit->pairarray,hit->npairs,
						    trim_left,trim_right,start_amb_length(hit),end_amb_length(hit),
						    hit->querylength);
      debug4(printf("  add mismatches %d.",hit->score_eventrim));
      if (start_amb_length(hit) > 0) {
	debug4(printf("  add penalty for start amb %d.",amb_penalty));
	hit->score_eventrim += amb_penalty;
      }
      if (end_amb_length(hit) > 0) {
	debug4(printf("  add penalty for end amb %d.",amb_penalty));
	hit->score_eventrim += amb_penalty;
      }
      
#ifdef SCORE_INDELS
      hit->score_eventrim += indel_penalty_middle * nindelbreaks;
#endif
      debug4(printf("  RESULT: %d\n",hit->score_eventrim));

    } else {
      hit->score_eventrim = 0; /* was hit->penalties */
      debug4(printf("score OTHER:"));

      for (q = hit->substrings_1toN; q != NULL; q = List_next(q)) {
	substring = (Substring_T) List_head(q);
	hit->score_eventrim += Substring_count_mismatches_region(substring,trim_left,trim_right,
								 query_compress_fwd,query_compress_rev);
	debug4(printf("  substring %d.",Substring_count_mismatches_region(substring,trim_left,trim_right,
									  query_compress_fwd,query_compress_rev)));
      }

#ifdef SCORE_INDELS
      /* Needs to match GMAP scoring */
      if (hit->hittype == INSERTION || hit->hittype == DELETION) {
	debug4(printf("  indel at %d",hit->indel_pos));
	if (hit->indel_pos > trim_left && hit->indel_pos < querylength - trim_right) {
	  hit->score_eventrim += indel_penalty_middle;
	  debug4(printf(" => add %d.",indel_penalty_middle));
	}
      }
#endif
      debug4(printf("  RESULT: %d\n",hit->score_eventrim));
    }
  }

  /* Compute minscore */
  for (p = hitlist; p != NULL; p = p->rest) {
    hit = (T) p->first;
    if (hit->hittype == TERMINAL && finalp == false) {
      /* Don't use to determine minscore */
#ifdef TRANSLOC_SPECIAL
    } else if (hit->chrnum == 0 && non_translocation_p == true) {
      /* Skip, since we will eliminate */
#endif
    } else {
      if (hit->nmatches > max_nmatches) {
	max_nmatches = hit->nmatches;
	max_nmatches_posttrim = hit->nmatches_posttrim;
      }
#ifdef TERMINAL_SECOND_CLASS
      if (non_gmap_terminal_p == true && (hit->hittype == TERMINAL || hit->hittype == GMAP)) {
	/* Skip from setting minscore */
      }
#endif
      if (hit->score_eventrim < minscore) {
	minscore = hit->score_eventrim;
      }
    }
  }

  debug4(printf("Stage3end_optimal_score over %d hits: minscore = %d + subopt:%d\n",
		n,minscore,suboptimal_mismatches));
  minscore += suboptimal_mismatches;
  max_nmatches -= suboptimal_mismatches;
  max_nmatches_posttrim -= suboptimal_mismatches;

#if 0
  if (non_gmap_terminal_p == false && minscore > cutoff_level) {
    /* If we are down to GMAP or terminal hits, keep at least one hit */
    cutoff_level = minscore;
  }
#else
  cutoff_level = minscore;
#endif

  for (p = hitlist; p != NULL; p = p->rest) {
    hit = (T) p->first;

    if (hit->hittype == TERMINAL && finalp == false) {
      debug4(printf("Keeping a hit of type TERMINAL\n"));
      optimal = List_push(optimal,hit);
      
    } else if (keep_gmap_p == true && hit->hittype == GMAP) {
      /* GMAP hits already found to be better than their corresponding terminals */
      debug4(printf("Keeping a hit of type GMAP\n"));
      optimal = List_push(optimal,hit);

#ifdef TRANSLOC_SPECIAL
    } else if (hit->chrnum == 0 && non_translocation_p == true) {
      debug4(printf("Eliminating a hit with splice translocation\n"));
      *eliminatedp = true;
      Stage3end_free(&hit);
#endif

#ifdef TERMINAL_SECOND_CLASS
    } else if ((hit->hittype == TERMINAL || hit->hittype == GMAP) &&
	       non_gmap_terminal_p == true) {
      if (hit->nmatches >= max_nmatches) {
	debug4(printf("Keeping a terminal with nmatches %d\n",hit->nmatches));
	optimal = List_push(optimal,(void *) hit);
      } else {
	debug4(printf("Eliminating a terminal where non-terminals are present\n"));
	*eliminatedp = true;
	Stage3end_free(&hit);
      }
#endif

#if 0
    } else if (hit->score_eventrim > cutoff_level) {
      /* Turning off, because we turned it off for Stage3pair_optimal_score_aux */
      /* For dibasep were previously using hit->ntscore, but gives false positives */
      debug4(printf("Eliminating a hit of type %s with score_eventrim %d > cutoff_level %d\n",
		    hittype_string(hit->hittype),hit->score_eventrim,cutoff_level));
      *eliminatedp = true;
      Stage3end_free(&hit);
#endif

    } else if (hit->score_eventrim > minscore /* && hit->nmatches_posttrim < max_nmatches_posttrim */) {
      debug4(printf("Eliminating a hit with score_eventrim %d and type %s\n",
		    hit->score_eventrim,hittype_string(hit->hittype)));
      *eliminatedp = true;
      Stage3end_free(&hit);

    } else {
      debug4(printf("Keeping a hit with score_eventrim %d and type %s\n",
		    hit->score_eventrim,hittype_string(hit->hittype)));
      optimal = List_push(optimal,hit);
    }
  }
  
  List_free(&hitlist);


  /* Filter on nsegments */
  if (finalp == true && optimal != NULL) {
    hitlist = optimal;
    optimal = (List_T) NULL;

    hit = (T) hitlist->first;
    best_nsegments = hit->nsegments;

    for (p = hitlist; p != NULL; p = p->rest) {
      hit = (T) p->first;
      if (hit->nsegments < best_nsegments) {
	best_nsegments = hit->nsegments;
      }
    }

    for (p = hitlist; p != NULL; p = p->rest) {
      hit = (T) p->first;
      if (hit->nsegments > best_nsegments + 2) {
	debug4(printf("Eliminating a hit with nsegments %d\n",hit->nsegments));
	*eliminatedp = true;
	Stage3end_free(&hit);
      } else {
	debug4(printf("Keeping a hit with nsegments %d\n",hit->nsegments));
	optimal = List_push(optimal,hit);
      }
    }

    List_free(&hitlist);
  }


  debug4(printf("hitlist now has %d entries\n",List_length(optimal)));
  return optimal;
}


List_T
Stage3end_optimal_score (List_T hitlist, int cutoff_level, int suboptimal_mismatches,
			 Compress_T query_compress_fwd, Compress_T query_compress_rev,
			 int querylength, bool keep_gmap_p, bool finalp) {
  List_T optimal;
  bool eliminatedp;


  optimal = Stage3end_optimal_score_aux(&eliminatedp,hitlist,cutoff_level,suboptimal_mismatches,
					query_compress_fwd,query_compress_rev,
					querylength,keep_gmap_p,finalp);
  while (eliminatedp == true) {
    optimal = Stage3end_optimal_score_aux(&eliminatedp,optimal,cutoff_level,suboptimal_mismatches,
					  query_compress_fwd,query_compress_rev,
					  querylength,keep_gmap_p,finalp);
  }

  return optimal;
}


static void
unalias_circular (T hit) {
  Chrpos_T chrlength = hit->chrlength;
  List_T p;
  Substring_T substring;

  debug12(printf("Calling unalias_circular\n"));
  assert(hit->alias == +1);
  if (hit->hittype == GMAP) {
    Pair_unalias_circular(hit->pairarray,hit->npairs,chrlength);

  } else {
    for (p = hit->substrings_1toN; p != NULL; p = List_next(p)) {
      substring = (Substring_T) List_head(p);
      Substring_unalias_circular(substring);
    }
  }

  /* Doesn't fix hitpair->low and hitpair->high */
  hit->genomicstart -= chrlength;
  hit->genomicend -= chrlength;
  hit->low -= chrlength;
  hit->high -= chrlength;

  hit->alias = -1;

  return;
}


#if 0
List_T
Stage3end_unalias_circular (List_T hitlist) {
  List_T p;
  T hit;

  for (p = hitlist; p != NULL; p = p->rest) {
    hit = (T) p->first;
    if (hit->alias == +1) {
      unalias_circular(hit);
    }
  }

  return hitlist;
}
#endif

List_T
Stage3end_remove_circular_alias (List_T hitlist) {
  List_T newlist = NULL, p;
  T hit;
#ifdef SOFT_CLIPS_AVOID_CIRCULARIZATION
  int trim;
#endif

  debug12(printf("Calling Stage3end_remove_circular_alias on %d hits\n",List_length(hitlist)));

  for (p = hitlist; p != NULL; p = p->rest) {
    hit = (T) p->first;

    if (hit->alias == +1) {
      /* First, try to salvage alias +1 */
      unalias_circular(hit);
    }

    if (hit->chrnum == 0) {
      /* Distant splice */
      newlist = List_push(newlist,(void *) hit);

    } else {
#ifdef SOFT_CLIPS_AVOID_CIRCULARIZATION
      /* This allows soft clips to avoid circularization */
      if (hit->plusp == true) {
	trim = hit->trim_left;
      } else {
	trim = hit->trim_right;
      }
#endif

      if (
#ifdef SOFT_CLIPS_AVOID_CIRCULARIZATION
	  hit->low + trim >= hit->chroffset + hit->chrlength
#else
	  hit->low >= hit->chroffset + hit->chrlength
#endif
	  ) {
	/* All in circular alias */
	debug12(printf("Freeing hit because all is in circular alias\n"));
	Stage3end_free(&hit);

      } else {
	newlist = List_push(newlist,(void *) hit);
      }
    }
  }

  List_free(&hitlist);
  return newlist;
}


#if 0
int
Stage3end_noptimal (List_T hitlist, int querylength) {
  int noptimal;
  List_T p;
  T hit;
  int minscore = querylength;

  noptimal = 0;
  for (p = hitlist; p != NULL; p = p->rest) {
    hit = (T) p->first;
    if (hit->score < minscore) {
      minscore = hit->score;
      noptimal = 0;
    }
    if (hit->score == minscore) {
      noptimal++;
    }
  }

  return noptimal;
}
#endif


static Univcoord_T
normalize_coord (Univcoord_T orig, int alias, Chrpos_T chrlength) {
  if (alias == +1) {
    return orig - chrlength;
  } else {
    return orig;
  }
}



static int
duplicate_sort_cmp (const void *a, const void *b) {
  int cmp;
  T x = * (T *) a;
  T y = * (T *) b;
  Univcoord_T x_genomicstart, x_genomicend, y_genomicstart, y_genomicend;
  List_T p, q;
  Substring_T x_substring, y_substring;

  x_genomicstart = normalize_coord(x->genomicstart,x->alias,x->chrlength);
  x_genomicend = normalize_coord(x->genomicend,x->alias,x->chrlength);

  y_genomicstart = normalize_coord(y->genomicstart,y->alias,y->chrlength);
  y_genomicend = normalize_coord(y->genomicend,y->alias,y->chrlength);

  
  if (x_genomicstart < y_genomicstart) {
    return -1;
  } else if (x_genomicstart > y_genomicstart) {
    return +1;
  } else if (x->hittype < y->hittype) {
    return -1;
  } else if (x->hittype > y->hittype) {
    return +1;
  } else if (x->genomicend < y->genomicend) {
    return -1;
  } else if (x->genomicend > y->genomicend) {
    return +1;

  } else {
    for (p = x->substrings_1toN, q = y->substrings_1toN; p != NULL && q != NULL; p = List_next(p), q = List_next(q)) {
      x_substring = (Substring_T) List_head(p);
      y_substring = (Substring_T) List_head(q);
      if ((cmp = Substring_compare(x_substring,y_substring,x->alias,y->alias,x->chrlength,y->chrlength)) != 0) {
	return cmp;
      }
    }
    if (p == NULL && q != NULL) {
      return -1;
    } else if (p != NULL && q == NULL) {
      return +1;
    }

#if 0
    /* Need to change to search on junctions */
    if (x->indel_low < y->indel_low) {
      return -1;
    } else if (y->indel_low < x->indel_low) {
      return +1;
    }
#endif

    if (x->sarrayp == true && y->sarrayp == false) {
      return -1;
    } else if (x->sarrayp == false && y->sarrayp == true) {
      return +1;
    } else {
      return 0;
    }
  }
}

/* Same as duplicate_sort_cmp, except for indel_low */
static int
duplicate_equiv_cmp (const void *a, const void *b) {
  int cmp;
  T x = * (T *) a;
  T y = * (T *) b;
  List_T p, q;
  Substring_T x_substring, y_substring;

  Univcoord_T x_genomicstart, x_genomicend, y_genomicstart, y_genomicend;

  x_genomicstart = normalize_coord(x->genomicstart,x->alias,x->chrlength);
  x_genomicend = normalize_coord(x->genomicend,x->alias,x->chrlength);

  y_genomicstart = normalize_coord(y->genomicstart,y->alias,y->chrlength);
  y_genomicend = normalize_coord(y->genomicend,y->alias,y->chrlength);

  if (x_genomicstart < y_genomicstart) {
    return -1;
  } else if (x_genomicstart > y_genomicstart) {
    return +1;
#if 0
  } else if (x->hittype < y->hittype) {
    return -1;
  } else if (x->hittype > y->hittype) {
    return +1;
#endif
  } else if (x_genomicend < y_genomicend) {
    return -1;
  } else if (x_genomicend > y_genomicend) {
    return +1;
  } else {
    for (p = x->substrings_1toN, q = y->substrings_1toN; p != NULL && q != NULL; p = List_next(p), q = List_next(q)) {
      x_substring = (Substring_T) List_head(p);
      y_substring = (Substring_T) List_head(q);
      if ((cmp = Substring_compare(x_substring,y_substring,x->alias,y->alias,x->chrlength,y->chrlength)) != 0) {
	return cmp;
      }
    }
    if (p == NULL && q != NULL) {
      return -1;
    } else if (p != NULL && q == NULL) {
      return +1;
    } else {
      return 0;
    }
  }
}


static int
genomicstart_cmp (const void *a, const void *b) {
  T x = * (T *) a;
  T y = * (T *) b;

  if (x->genomicstart < y->genomicstart) {
    return -1;
  } else if (x->genomicstart > y->genomicstart) {
    return +1;
  } else if (x->genomicend < y->genomicend) {
    return -1;
  } else if (x->genomicend > y->genomicend) {
    return +1;
  } else {
    return 0;
  }
}

static int
genomicend_cmp (const void *a, const void *b) {
  T x = * (T *) a;
  T y = * (T *) b;

  if (x->genomicend < y->genomicend) {
    return -1;
  } else if (x->genomicend > y->genomicend) {
    return +1;
  } else if (x->genomicstart < y->genomicstart) {
    return -1;
  } else if (x->genomicstart > y->genomicstart) {
    return +1;
  } else {
    return 0;
  }
}


#if 0

List_T
Stage3end_mark_ambiguous_splices (bool *ambiguousp, List_T hitlist) {
  T x, y, *hits;
  int n, i, j;
  Chrpos_T splice_distance_1, splice_distance_2;

#ifndef CLEAN_SINGLE_END_AMBIGUOUS
  return hitlist;
#endif

  *ambiguousp = false;
  
  debug9(printf("Entered Stage3end_mark_ambiguous_splices with %d pairs\n",n));
  if ((n = List_length(hitlist)) == 0) {
    return NULL;
  } else {
#ifdef USE_ALLOCA_FOR_HITS
    hits = (T *) MALLOCA(n * sizeof(T));
    List_fill_array((void **) hits,hitlist); /* hitlist is a return value */
#else
    hits = (T *) List_to_array(hitlist,NULL);
#endif
  }

  /* By genomicstart */
  debug9(printf("Stage3end_mark_ambiguous_splices: checking %d hits by genomicstart\n",n));
  qsort(hits,n,sizeof(T),genomicstart_cmp);
  for (i = 0; i < n; i++) {
    x = hits[i];
    if (x->hittype == SPLICE) {
      j = i+1;
      while (j < n && hits[j]->genomicstart == x->genomicstart) {
	y = hits[j];
	if (y->hittype == SPLICE) {
	  debug9(printf("  #%d has common start with #%d\n",i,j));
	  if (SENSE_INCONSISTENT_P(x->sensedir,y->sensedir)) {
	    debug9(printf("  #%d and #%d have different sense, so skipping\n",i,j));
	  } else if (x->score != y->score) {
	    debug9(printf("  #%d overlaps #%d and score %d != %d, so skipping\n",
			  i,j,x->score,y->score));
	  } else {
	    debug9(printf("  #%d overlaps #%d and score %d == %d, so both are ambiguous\n",
			  i,j,x->score,y->score));
	    if (x->genomicstart < x->genomicend) {
	      splice_distance_1 = x->genomicend - x->genomicstart;
	    } else {
	      splice_distance_1 = x->genomicstart - x->genomicend;
	    }

	    if (y->genomicstart < y->genomicend) {
	      splice_distance_2 = y->genomicend - y->genomicstart;
	    } else {
	      splice_distance_2 = y->genomicstart - y->genomicend;
	    }

	    debug9(printf("    splice distances: %u and %u\n",splice_distance_1,splice_distance_2));
	    if (splice_distance_1 > splice_distance_2) {
	      debug9(printf("    first is ambiguous\n"));
	      x->chimera_ambiguous_p = true;
	      *ambiguousp = true;
	    } else if (splice_distance_2 > splice_distance_1) {
	      debug9(printf("    second is ambiguous\n"));
	      y->chimera_ambiguous_p = true;
	      *ambiguousp = true;
	    }
	  }
	}
	j++;
      }
      i = j;
    }
  }
    
  /* By genomicend */
  debug9(printf("Stage3end_mark_ambiguous_splices: checking %d hits by genomicend\n",n));
  qsort(hits,n,sizeof(T),genomicend_cmp);
  for (i = 0; i < n; i++) {
    x = hits[i];
    if (x->hittype == SPLICE) {
      j = i+1;
      while (j < n && hits[j]->genomicend == x->genomicend) {
	y = hits[j];
	if (y->hittype == SPLICE) {
	  debug9(printf("  #%d has common end with #%d\n",i,j));
	  if (SENSE_INCONSISTENT_P(x->sensedir,y->sensedir)) {
	    debug9(printf("  #%d and #%d have different sense, so skipping\n",i,j));
	  } else if (x->score != y->score) {
	    debug9(printf("  #%d overlaps #%d and score %d != %d, so skipping\n",
			  i,j,x->score,y->score));
	  } else {
	    debug9(printf("  #%d overlaps #%d and score %d == %d, so both are ambiguous\n",
			  i,j,x->score,y->score));
	    if (x->genomicstart < x->genomicend) {
	      splice_distance_1 = x->genomicend - x->genomicstart;
	    } else {
	      splice_distance_1 = x->genomicstart - x->genomicend;
	    }

	    if (y->genomicstart < y->genomicend) {
	      splice_distance_2 = y->genomicend - y->genomicstart;
	    } else {
	      splice_distance_2 = y->genomicstart - y->genomicend;
	    }

	    debug9(printf("    splice distances: %u and %u\n",splice_distance_1,splice_distance_2));
	    if (splice_distance_1 > splice_distance_2) {
	      debug9(printf("    first is ambiguous\n"));
	      x->chimera_ambiguous_p = true;
	      *ambiguousp = true;
	    } else if (splice_distance_2 > splice_distance_1) {
	      debug9(printf("    second is ambiguous\n"));
	      y->chimera_ambiguous_p = true;
	      *ambiguousp = true;
	    }
	  }
	}
	j++;
      }
      i = j;
    }
  }

  debug9(
	 for (i = 0; i < n; i++) {
	   x = hits[i];
	   if (x->hittype == SPLICE) {
	     printf("  %d: %u..%u (plusp = %d) ambiguousp:%d\n",
		    i,x->genomicstart - x->chroffset,x->genomicend - x->chroffset,x->plusp,x->chimera_ambiguous_p);
	   }
	 }
	 );

#ifdef USE_ALLOCA_FOR_HITS
  FREEA(hits);
#else
  FREE(hits);
#endif

  return hitlist;
}

#endif

#ifdef DEBUG7
static void
Stage3end_print_substrings (Stage3end_T hit) {
  List_T p;
  Substring_T substring;

  for (p = hit->substrings_1toN; p != NULL; p = List_next(p)) {
    substring = (Substring_T) List_head(p);
    Substring_print_ends(substring,hit->chrnum);
  }
  return;
}
#endif


const Except_T Duplicate_Pairing = { "Duplicates both seen in pairing" };

List_T
Stage3end_remove_duplicates (List_T hitlist) {
#ifdef DEBUG7
  List_T p;
#endif
  T x, y, *hits;
  int n, usedi, i, j, k;
  bool *eliminate, eliminatep;

  
  debug7(printf("Entered Stage3end_remove_duplicates with %d hits\n",List_length(hitlist)));
  if ((n = List_length(hitlist)) == 0) {
    return NULL;
  } else {
#ifdef USE_ALLOCA_FOR_HITS
    eliminate = (bool *) CALLOCA(n,sizeof(bool));
    hits = (T *) MALLOCA(n * sizeof(T));    
    List_fill_array((void **) hits,hitlist); /* hitlist is a return value */
#else
    eliminate = (bool *) CALLOC(n,sizeof(bool));
    hits = (T *) List_to_array(hitlist,NULL);
#endif
  }


  /* By equivalence */
  debug7(printf("Stage3end_remove_duplicates: checking %d hits by equivalence class\n",n));
  qsort(hits,n,sizeof(T),duplicate_sort_cmp);

  debug7(
	 for (i = 0; i < n; i++) {
	   x = hits[i];
	   printf("  Initial %d (%s): %p #%d:%u..%u, alias %d, nmatches %d, score %d ",
		  i,hittype_string(x->hittype),x,x->chrnum,x->genomicstart - x->chroffset,x->genomicend - x->chroffset,
		  x->alias,x->nmatches,x->score);
	   Stage3end_print_substrings(x);
	   printf("\n");
	 }
	 );

  eliminatep = false;
  i = 0;
  while (i < n) {
    j = i+1;
    while (j < n && duplicate_equiv_cmp(&(hits[j]),&(hits[i])) == 0) {
      j++;
    }

    if (j > i+1) {
      debug7(printf("Equivalence class #%d through #%d.  ",i,j-1));

      x = hits[i];
      if (x->paired_usedp == true) {
	usedi = i;
      } else {
	usedi = -1;
      }
      
      for (k = i+1; k < j; k++) {
	y = hits[k];
	if (y->paired_usedp == true) {
	  if (usedi >= 0) {
	    debug7(printf("  #%d equivalent to #%d and both used (%p and %p)\n",k,usedi,hits[k],hits[usedi]));
#if 0
	    /* This doesn't matter anymore.  Example from NM_001033853:
	       TTGCCCTTGGTCACCCCGATGACGTCGATCATCTCATCCTGCCCAAACACTTGGTTCACAGGTACCTGCTGCTCA
	       AGTGATGAATCCAAGAGGCGTTTCTATAAGAATTGGCATAAATCTAAGAAGAAGGCCCACCTGATGGAGATCCAG */
	    fprintf(stderr,"Duplicates of Stage3end_T both seen\n");
#if 0
	    /* No longer providing queryseq1 and queryseq2 */
	    Shortread_print_query_pairedend_fasta(stderr,queryseq1,queryseq2,
						  /*invert_first_p*/false,/*invert_second_p*/true);
#endif
	    Except_raise(&Duplicate_Pairing, __FILE__, __LINE__);
#endif
	  } else {
	    usedi = k;
	  }
	}
      }

      if (usedi < 0) {
	debug7(printf("None used yet so eliminating #%d through #%d\n",i+1,j-1));
	for (k = i+1; k < j; k++) {
	  eliminate[k] = true;
	  eliminatep = true;
	}
      } else {
	debug7(printf("One used already so eliminating all but #%d\n",usedi));
	for (k = i; k < j; k++) {
	  if (k != usedi) {
	    eliminate[k] = true;
	    eliminatep = true;
	  }
	}
      }
    }

    i = j;
  }
    

#if 0
  nkept = 0;
  for (i = 0; i < n; i++) {
    if (eliminate[i] == false) {
      nkept++;
    }
  }
  if (nkept == 0) {
    /* All entries eliminated one another, so keep the first one */
    eliminate[0] = false;
  }
#endif

  if (eliminatep == false) {
    debug7(printf("No eliminations, so hitlist is unchanged\n"));
  } else {
    List_free(&hitlist);
    hitlist = (List_T) NULL;
    for (i = n-1; i >= 0; i--) {
      x = hits[i];
      if (eliminate[i] == false) {
	debug7(printf("  Keeping #%d:%u..%u, nmatches %d (nindels %d, distance %u, chrnum %d) (plusp = %d)\n",
		      x->chrnum,x->genomicstart - x->chroffset,x->genomicend - x->chroffset,
		      x->nmatches,x->nindels,x->distance,x->chrnum,x->plusp));
	hitlist = List_push(hitlist,x);
      } else {
	debug7(printf("  Eliminating #%d:%u..%u, nmatches %d (nindels %d, distance %u, chrnum %d) (plusp = %d)\n",
		      x->chrnum,x->genomicstart - x->chroffset,x->genomicend - x->chroffset,
		      x->nmatches,x->nindels,x->distance,x->chrnum,x->plusp));
	Stage3end_free(&x);
      }
    }
  }

#ifdef USE_ALLOCA_FOR_HITS
  FREEA(hits);
  FREEA(eliminate);
#else
  FREE(hits);
  FREE(eliminate);
#endif

  debug7(
	 for (p = hitlist, i = 0; p != NULL; p = p->rest, i++) {
	   x = (T) p->first;
	   printf("  Final %d: #%d:%u..%u (plusp = %d)\n",
		  i,x->chrnum,x->genomicstart - x->chroffset,x->genomicend - x->chroffset,x->plusp);
	 }
	 );

  debug7(printf("Exited Stage3end_remove_duplicates with %d hits\n",List_length(hitlist)));
  return hitlist;
}


#if 0
List_T
Stage3end_reject_trimlengths (List_T hits) {
  List_T filtered = NULL, p;
  T hit;

  for (p = hits; p != NULL; p = p->rest) {
    hit = (T) p->first;
    if (hit->trim_left + hit->trim_right >= reject_trimlength) {
      Stage3end_free(&hit);
    } else {
      filtered = List_push(filtered,(void *) hit);
    }
  }

  List_free(&hits);
  return filtered;
}
#endif


/* Used for eliminating exact duplicates.  Also sorts secondarily by hittype. */
static int
hit_sort_cmp (const void *a, const void *b) {
  Stage3end_T x = * (Stage3end_T *) a;
  Stage3end_T y = * (Stage3end_T *) b;
  
  debug7(printf("Comparing %s: #%d:%u..%u, alias %d, nmatches %d, score %d with %s: #%d:%u..%u, alias %d, nmatches %d, score %d\n",
		hittype_string(x->hittype),x->chrnum,x->genomicstart-x->chroffset,x->genomicend-x->chroffset,
		x->alias,x->nmatches,x->score,
		hittype_string(y->hittype),y->chrnum,y->genomicstart-y->chroffset,y->genomicend-y->chroffset,
		y->alias,y->nmatches,y->score));

  if (x->plusp > y->plusp) {
    return -1;
  } else if (y->plusp > x->plusp) {
    return +1;

#if 0
  } else if (x->high < y->low) {
    return -1;
  } else if (y->high < x->low) {
    return +1;
#else
  } else if (x->low < y->low) {
    debug7(printf("Returning -1 for low\n"));
    return -1;
  } else if (y->low < x->low) {
    debug7(printf("Returning +1 for low\n"));
    return +1;

  } else if (x->high < y->high) {
    debug7(printf("Returning -1 for high\n"));
    return -1;
  } else if (y->high < x->high) {
    debug7(printf("Returning +1 for high\n"));
    return +1;
#endif

#if 1
    /* Rank terminals last, so terminals cannot win */
  } else if (x->hittype != TERMINAL && y->hittype == TERMINAL) {
    return -1;
  } else if (x->hittype == TERMINAL && y->hittype != TERMINAL) {
    return +1;
#endif

  } else if (x->score < y->score) {
    return -1;
  } else if (y->score < x->score) {
    return +1;
  } else if (x->nmatches > y->nmatches) {
    return -1;
  } else if (y->nmatches > x->nmatches) {
    return +1;
#if 0
  } else if (x->nmatches_posttrim > y->nmatches_posttrim) {
    return -1;
  } else if (y->nmatches_posttrim > x->nmatches_posttrim) {
    return +1;
#endif

#if 0
  } else if (x->nsplices < y->nsplices) {
    return -1;
  } else if (y->nsplices < x->nsplices) {
    return +1;
#endif

  } else if (x->hittype < y->hittype) {
    return -1;
  } else if (y->hittype < x->hittype) {
    return +1;

#if 0
  } else if (y->start_amb_length + y->end_amb_length == 0 &&
	     x->start_amb_length + x->end_amb_length > 0) {
    return -1;
  } else if (x->start_amb_length + x->end_amb_length == 0 &&
	     y->start_amb_length + y->end_amb_length > 0) {
    return +1;
#endif

#if 0
  } else if (x->indel_low < y->indel_low) {
    return -1;
  } else if (y->indel_low < x->indel_low) {
    return +1;
#endif

  } else if (x->sensedir != 0 && y->sensedir == 0) {
    return -1;
  } else if (y->sensedir != 0 && x->sensedir == 0) {
    return +1;

  } else if (x->sarrayp == true && y->sarrayp == false) {
    return -1;
  } else if (x->sarrayp == false && y->sarrayp == true) {
    return +1;
  } else {
    return 0;
  }
}

/* Same as hit_sort_cmp, except for hittype, nmatches_posttrim, and indel_low */
static int
hit_equiv_cmp (Stage3end_T x, Stage3end_T y) {

  if (x->plusp > y->plusp) {
    return -1;
  } else if (y->plusp > x->plusp) {
    return +1;
  } else if (x->low < y->low) {
    return -1;
  } else if (y->low < x->low) {
    return +1;
  } else if (x->high < y->high) {
    return +1;
  } else if (y->high < x->high) {
    return -1;

#if 0
  } else if (x->hittype != TERMINAL && y->hittype == TERMINAL) {
    return -1;
  } else if (x->hittype == TERMINAL && y->hittype != TERMINAL) {
    return +1;
#endif

  } else if (x->score < y->score) {
    return -1;
  } else if (y->score < x->score) {
    return +1;
  } else if (x->nmatches > y->nmatches) {
    return -1;
  } else if (y->nmatches > x->nmatches) {
    return +1;
#if 0
  } else if (x->nmatches_posttrim > y->nmatches_posttrim) {
    return -1;
  } else if (y->nmatches_posttrim > x->nmatches_posttrim) {
    return +1;
#endif

#if 0
    /* Causes hits to not be recognized as equivalent */
  } else if (x->nsplices < y->nsplices) {
    return -1;
  } else if (y->nsplices < x->nsplices) {
    return +1;
#endif

#if 0
  } else if (y->start_amb_length + y->end_amb_length == 0 &&
	     x->start_amb_length + x->end_amb_length > 0) {
    return -1;
  } else if (x->start_amb_length + x->end_amb_length == 0 &&
	     y->start_amb_length + y->end_amb_length > 0) {
    return +1;
#endif

#if 0
  } else if (x->indel_low < y->indel_low) {
    return -1;
  } else if (y->indel_low < x->indel_low) {
    return +1;
#endif

#if 0
  } else if (x->sarrayp == true && y->sarrayp == false) {
    return -1;
  } else if (x->sarrayp == false && y->sarrayp == true) {
    return +1;
#endif

#if 0
    /* Used for sorting but not equiv */
  } else if (x->sensedir != 0 && y->sensedir == 0) {
    return -1;
  } else if (y->sensedir != 0 && x->sensedir == 0) {
    return +1;
#endif

  } else {
    return 0;
  }
}


static int
hit_goodness_cmp (bool *equalp, Stage3end_T hit,
#ifdef DEBUG7
		  int k,
#endif
		  Stage3end_T best_hit, bool finalp) {
  double prob1, prob2;

#ifdef PRE_RESOLVE_MULTIMAPPING
  if (Stage3end_tally(x) > TALLY_RATIO*Stage3end_tally(y)) {
    debug7(printf("  #%d overlaps #%d and tally %ld > %f*%ld, so marking %d for elimination\n",
		  i,j,x->tally,TALLY_RATIO,y->tally,j));
    eliminate[j] = true;
  } else if (Stage3end_tally(y) > TALLY_RATIO*Stage3end_tally(x)) {
    debug7(printf("  #%d overlaps #%d and tally %f*%ld < %ld, so marking %d for elimination\n",
		  i,j,TALLY_RATIO,x->tally,y->tally,i));
    eliminate[i] = true;
  }
#endif

  *equalp = false;

#if 1
  if (finalp == true) {
    /* Skip */
  } else if (hit->hittype == TERMINAL || best_hit->hittype == TERMINAL) {
    /* Do not allow terminal to win or lose in pre-final stages */
    debug7(printf(" => %d ties by terminal\n",k));
    return 0;
  }
#endif

#if 0
  if (hit->hittype == TERMINAL || best_hit->hittype == TERMINAL) {
    /* Skip: Don't use scores if terminal is involved */
  } else if (hit->nindels == 0 && best_hit->nindels == 0) {
    /* Skip: Use scores only if indel is involved */
  } else if (hit->score > best_hit->score) {
    debug7(printf("  => %d loses by score\n",k));
    return -1;
  } else if (hit->score < best_hit->score) {
    debug7(printf("  => %d wins by score\n",k));
    return +1;
  }
#endif

  if (hit->nmatches < best_hit->nmatches) {
    debug7(printf("  => %d loses by nmatches\n",k));
    return -1;
  } else if (hit->nmatches > best_hit->nmatches) {
    debug7(printf("  => %d wins by nmatches\n",k));
    return +1;

#if 0
    /* Causes ambiguous splices to lose to a definitive splice, which could be wrong */
  } else if (hit->nmatches_posttrim < best_hit->nmatches_posttrim) {
    debug7(printf("  => %d loses by nmatches (post-trim)\n",k));
    return -1;
  } else if (hit->nmatches_posttrim > best_hit->nmatches_posttrim) {
    debug7(printf("  => %d wins by nmatches (post-trim)\n",k));
    return +1;
#endif

#if 0
  } else if (hit->nsplices > best_hit->nsplices) {
    debug7(printf("  => %d loses by nsplices: %d > %d in best\n",k,hit->nsplices,best_hit->nsplices));
    return -1;
  } else if (hit->nsplices < best_hit->nsplices) {
    debug7(printf("  => %d wins by nsplices: %d < %d in best\n",k,hit->nsplices,best_hit->nsplices));
    return +1;
#endif

  } else if (hit->hittype > best_hit->hittype) {
    debug7(printf("  => %d loses by hittype\n",k));
    return -1;
  } else if (hit->hittype < best_hit->hittype) {
    debug7(printf("  => %d wins by hittype\n",k));
    return +1;

#if 0
  } else if (hit->start_amb_length + hit->end_amb_length == 0 &&
	     best_hit->start_amb_length + best_hit->end_amb_length > 0) {
    debug7(printf("  => %d loses by ambiguity\n",k));
    return -1;
  } else if (hit->start_amb_length + hit->end_amb_length > 0 &&
	     best_hit->start_amb_length + best_hit->end_amb_length == 0) {
    debug7(printf("  => %d wins by ambiguity\n",k));
    return +1;
#endif

  } else if (hit->nindels > best_hit->nindels) {
    debug7(printf("  => %d loses by nindels\n",k));
    return -1;
  } else if (hit->nindels < best_hit->nindels) {
    debug7(printf("  => %d wins by nindels\n",k));
    return +1;

  } else if (hit->chrnum == 0 && best_hit->chrnum != 0) {
    debug7(printf("  => %d loses because distant splice\n",k));
    return -1;
  } else if (hit->chrnum > 0 && best_hit->chrnum == 0) {
    debug7(printf("  => %d wins because not distant splice\n",k));
    return +1;

  } else if (finalp == false) {
    debug7(printf("  => indistinguishable\n"));
    return 0;

  } else {
    prob1 = Stage3end_prob(hit);
    prob2 = Stage3end_prob(best_hit);
    if (prob1 < prob2) {
      debug7(printf("  => %d loses by splice prob %f vs %f\n",k,prob1,prob2));
      return -1;
    } else if (prob1 > prob2) {
      debug7(printf("  => %d wins by splice prob %f vs %f\n",k,prob1,prob2));
      return +1;
    }

    if (hit->genomiclength > best_hit->genomiclength) {
      debug7(printf("  => %d loses by genomiclength\n",k));
      return -1;
    } else if (hit->genomiclength < best_hit->genomiclength) {
      debug7(printf("  => %d wins by genomiclength\n",k));
      return +1;

    } else {
      debug7(printf("  => equal\n"));
      *equalp = true;
      return 0;
    }
  }
}


static bool
hit_subsumption (Stage3end_T x, Stage3end_T y) {
  if (x->plusp != y->plusp) {
    return false;		/* Different strands */
  } else if (x->low <= y->low && x->high >= y->high) {
    return true;
  } else if (y->low <= x->low && y->high >= x->high) {
    return true;
  } else {
    return false;
  }
}

static bool
hit_bad_superstretch_p (Stage3end_T hit_k, Stage3end_T *hits, int k, int j, bool finalp) {
  int a;
  bool equalp;

  for (a = k+1; a <= j; a++) {
    if (hit_subsumption(hit_k,hits[a]) == true) {
      debug7(printf("Testing %d because stretches over %d",k,a));
      if (hit_goodness_cmp(&equalp,hits[a],
#ifdef DEBUG7
			   a,
#endif
			   hit_k,finalp) > 0 || equalp == true) {
	debug7(printf(" => eliminating\n"));
	return true;
      }
      debug7(printf("\n"));
    }
  }
  return false;
}



List_T
Stage3end_remove_overlaps (List_T hitlist, bool finalp) {
  List_T unique = NULL;
  T best_hit, hit, *hits, *prev;
  int cmp;
  int nkept, n, i, j, k, besti;
  bool *eliminate, equalp;
#ifdef PRE_RESOLVE_MULTIMAPPING
  long int best_tally;
#endif

  
  debug7(printf("Entered Stage3end_remove_overlaps with %d hits: %s\n",
		List_length(hitlist),finalp == true ? "FINAL" : "not final"));
  if ((n = List_length(hitlist)) == 0) {
    return NULL;
  } else {
#ifdef USE_ALLOCA_FOR_HITS
    eliminate = (bool *) CALLOCA(n,sizeof(bool));
    hits = (T *) MALLOCA(n * sizeof(T));
    List_fill_array_and_free((void **) hits,&hitlist);
#else
    eliminate = (bool *) CALLOC(n,sizeof(bool));
    hits = (T *) List_to_array(hitlist,NULL);
    List_free(&hitlist);
#endif
  }


  /* Step 1.  Check for exact duplicates */
  /* Probably don't want to eliminate aliases at this point */
  debug7(printf("Step 1.  Checking for exact duplicates\n"));
  qsort(hits,n,sizeof(Stage3end_T),hit_sort_cmp);

  debug7(
	 for (i = 0; i < n; i++) {
	   hit = hits[i];
	   printf("  Initial %d (%s): %p #%d:%u..%u, alias %d, nmatches %d, score %d\n",
		  i,hittype_string(hit->hittype),hit,hit->chrnum,hit->genomicstart-hit->chroffset,hit->genomicend-hit->chroffset,
		  hit->alias,hit->nmatches,hit->score);
	 }
	 );

  i = 0;
  while (i < n) {
    j = i+1;
    debug7(printf(" %d,%d",i,j));
    while (j < n && hit_equiv_cmp(hits[j],hits[i]) == 0) {
      debug7(printf("  %d is identical to %d => eliminating\n",j,i));
      eliminate[j] = true;
      j++;
    }
    i = j;
  }
  debug7(printf("\n"));


  nkept = 0;
  for (i = 0; i < n; i++) {
    if (eliminate[i] == false) {
      nkept++;
    } else if (hits[i]->paired_usedp == true) {
      nkept++;
    }
  }
  if (nkept == 0) {
    /* All entries eliminated one another, so keep the first one */
    eliminate[0] = false;
    nkept = 1;
  }

  prev = hits;
#ifdef USE_ALLOCA_FOR_HITS
  hits = (Stage3end_T *) MALLOCA(nkept * sizeof(Stage3end_T));
#else
  hits = (Stage3end_T *) MALLOC(nkept * sizeof(Stage3end_T));
#endif

  for (i = 0, j = 0; i < n; i++) {
    hit = prev[i];
    if (eliminate[i] == false) {
      debug7(printf("  Keeping %u..%u, nmatches (trimmed) %d (plusp = %d)\n",
		    hit->low - hit->chroffset,hit->high - hit->chroffset,hit->nmatches,hit->plusp));
      hits[j++] = hit;
    } else if (hit->paired_usedp == true) {
      debug7(printf("  Already paired %u..%u, nmatches (trimmed) %d (plusp = %d)\n",
		    hit->low - hit->chroffset,hit->high - hit->chroffset,hit->nmatches,hit->plusp));
      hits[j++] = hit;
    } else {
      debug7(printf("  Eliminating %u..%u, nmatches (trimmed) %d (plusp = %d)\n",
		    hit->low - hit->chroffset,hit->high - hit->chroffset,hit->nmatches,hit->plusp));
      Stage3end_free(&hit);
    }
  }

#ifdef USE_ALLOCA_FOR_HITS
  FREEA(prev);
#else
  FREE(prev);
#endif


  /* Step 2: Check for superstretches */
  n = nkept;
  debug7(printf("Step 2.  Checking for superstretches among %d hits within subsumption clusters\n",n));

  for (i = 0; i < n; i++) {
    eliminate[i] = false;
  }

  debug7(
	 for (i = 0; i < n; i++) {
	   hit = hits[i];
	   printf("  Initial %d (%s): %p #%d:%u..%u, nmatches (trimmed) %d, nmatches (post-trim) %d, score %d\n",
		  i,hittype_string(hit->hittype),hit,hit->chrnum,hit->genomicstart-hit->chroffset,hit->genomicend-hit->chroffset,
		  hit->nmatches,hit->nmatches_posttrim,hit->score);
	 }
	 );

  /* Find clusters */
  i = 0;
  while (i < n) {
    j = i;
    if (hits[i]->chrnum != 0) {
      while (j+1 < n && hit_subsumption(hits[i],hits[j+1]) == true) {
	j = j+1;
      }

      if (j > i) {
	debug7(printf("Cluster from %d up through %d\n",i,j));

	/* Find bad superstretches */
	for (k = i; k <= j; k++) {
	  if (hits[k]->chrnum == 0) {
	    /* Skip */
	  } else if (hit_bad_superstretch_p(hits[k],hits,k,j,finalp) == true) {
	    eliminate[k] = true;
	  }
	}
      }
    }

    i = j+1;
  }

  nkept = 0;
  for (i = 0; i < n; i++) {
    if (eliminate[i] == false) {
      nkept++;
    } else if (hits[i]->paired_usedp == true) {
      nkept++;
    }
  }
  if (nkept == 0) {
    /* All entries eliminated one another, so keep the first one */
    eliminate[0] = false;
    nkept = 1;
  }

  prev = hits;
#ifdef USE_ALLOCA_FOR_HITS
  hits = (Stage3end_T *) MALLOCA(nkept * sizeof(Stage3end_T));
#else
  hits = (Stage3end_T *) MALLOC(nkept * sizeof(Stage3end_T));
#endif

  for (i = 0, j = 0; i < n; i++) {
    hit = prev[i];
    if (eliminate[i] == false) {
      debug7(printf("  Keeping %u..%u, nmatches (trimmed) %d (plusp = %d)\n",
		    hit->low - hit->chroffset,hit->high - hit->chroffset,hit->nmatches,hit->plusp));
      hits[j++] = hit;
    } else if (hit->paired_usedp == true) {
      debug7(printf("  Already paired %u..%u, nmatches (trimmed) %d (plusp = %d)\n",
		    hit->low - hit->chroffset,hit->high - hit->chroffset,hit->nmatches,hit->plusp));
      hits[j++] = hit;
    } else {
      debug7(printf("  Eliminating %u..%u, nmatches (trimmed) %d (plusp = %d)\n",
		    hit->low - hit->chroffset,hit->high - hit->chroffset,hit->nmatches,hit->plusp));
      Stage3end_free(&hit);
    }
  }

#ifdef USE_ALLOCA_FOR_HITS
  FREEA(prev);
#else
  FREE(prev);
#endif


  /* Step 3: Check for best within subsumption clusters */
  n = nkept;
  debug7(printf("Checking for best among %d hits within subsumption clusters\n",n));

  for (i = 0; i < n; i++) {
    eliminate[i] = false;
  }
  /* qsort(hits,n,sizeof(Stage3end_T),hit_sort_cmp); -- No need since original order was kept */
  
  debug7(
	 for (i = 0; i < n; i++) {
	   hit = hits[i];
	   printf("  Initial %d (%s): %p #%d:%u..%u, nmatches %d, score %d\n",
		  i,hittype_string(hit->hittype),hit,hit->chrnum,hit->genomicstart-hit->chroffset,hit->genomicend-hit->chroffset,
		  hit->nmatches,hit->score);
	 }
	 );

  /* Find clusters from left */
  i = 0;
  while (i < n) {
    j = i;
    if (hits[i]->chrnum != 0) {
      while (j+1 < n && hit_subsumption(hits[i],hits[j+1]) == true) {
	j = j+1;
      }

      if (j > i) {
	debug7(printf("Cluster from %d up through %d\n",i,j));

	best_hit = hits[i];
	besti = i;
	debug7(printf("Assume best is %d\n",besti));

	for (k = i+1; k <= j; k++) {
	  if (hits[k]->chrnum == 0) {
	    /* Skip */
	  } else {
	    cmp = hit_goodness_cmp(&equalp,hits[k],
#ifdef DEBUG7
				   k,
#endif
				   best_hit,finalp);
	    debug7(printf("Comparison of %d with best %d yields %d\n",k,besti,cmp));
	    if (cmp > 0) {
	      best_hit = hits[k];
	      besti = k;
	      debug7(printf("Best is now %d\n",besti));
	    }
	  }
	}

	for (k = i; k <= j; k++) {
	  if (k == besti) {
	    /* Skip */
	  } else if (hits[k]->chrnum == 0) {
	    /* Skip */
	  } else if (hit_goodness_cmp(&equalp,hits[k],
#ifdef DEBUG7
				      k,
#endif
				      best_hit,finalp) < 0 || equalp == true) {
	    debug7(printf("  Eliminating hit %d from left, because beaten by %d\n",k,besti));
	    eliminate[k] = true;
	  }
	}
      }
    }
      
    i = j+1;
  }


  /* Find clusters starting from right */
  j = n - 1;
  while (j >= 0) {
    i = j;
    if (hits[j]->chrnum != 0) {
      while (i-1 >= 0 && hit_subsumption(hits[j],hits[i-1]) == true) {
	i = i-1;
      }

      if (i < j) {
	debug7(printf("Cluster from %d down through %d\n",j,i));
	best_hit = hits[i];
	besti = i;
	debug7(printf("Assume best is %d\n",besti));

	for (k = i+1; k <= j; k++) {
	  if (hits[k]->chrnum == 0) {
	    /* Skip */
	  } else {
	    cmp = hit_goodness_cmp(&equalp,hits[k],
#ifdef DEBUG7
				   k,
#endif
				   best_hit,finalp);
	    debug7(printf("Comparison of %d with best %d yields %d\n",k,besti,cmp));
	    if (cmp > 0) {
	      best_hit = hits[k];
	      besti = k;
	      debug7(printf("Best is now %d\n",besti));
	    }
	  }
	}

	for (k = i; k <= j; k++) {
	  if (k == besti) {
	    /* Skip */
	  } else if (hits[k]->chrnum == 0) {
	    /* Skip */
	  } else if (hit_goodness_cmp(&equalp,hits[k],
#ifdef DEBUG7
				      k,
#endif
				      best_hit,finalp) < 0 || equalp == true) {
	    debug7(printf("  Eliminating hit %d from right, because beaten by %d\n",k,besti));
	    eliminate[k] = true;
	  }
	}
      }
    }
      
    j = i-1;
  }


  nkept = 0;
  for (i = 0; i < n; i++) {
    if (eliminate[i] == false) {
      nkept++;
    } else if (hits[i]->paired_usedp == true) {
      nkept++;
    }
  }
  if (nkept == 0) {
    eliminate[0] = false;
    nkept = 1;
  }

  prev = hits;
#ifdef USE_ALLOCA_FOR_HITS
  hits = (Stage3end_T *) MALLOCA(nkept * sizeof(Stage3end_T));
#else
  hits = (Stage3end_T *) MALLOC(nkept * sizeof(Stage3end_T));
#endif

  for (i = 0, j = 0; i < n; i++) {
    hit = prev[i];
    if (eliminate[i] == false) {
      debug7(printf("  Keeping %u..%u, nmatches (trimmed) %d (plusp = %d)\n",
		    hit->low - hit->chroffset,hit->high - hit->chroffset,hit->nmatches,hit->plusp));
      hits[j++] = hit;
    } else if (hit->paired_usedp == true) {
      debug7(printf("  Already paired %u..%u, nmatches (trimmed) %d (plusp = %d)\n",
		    hit->low - hit->chroffset,hit->high - hit->chroffset,hit->nmatches,hit->plusp));
      hits[j++] = hit;
    } else {
      debug7(printf("  Eliminating %u..%u, nmatches (trimmed) %d (plusp = %d)\n",
		    hit->low - hit->chroffset,hit->high - hit->chroffset,hit->nmatches,hit->plusp));
      Stage3end_free(&hit);
    }
  }

#ifdef USE_ALLOCA_FOR_HITS
  FREEA(prev);
#else
  FREE(prev);
#endif


  /* Step 4: Check for identity */
  n = nkept;
  debug7(printf("Checking for duplicates among %d hits by identity\n",n));

  for (i = 0; i < n; i++) {
    eliminate[i] = false;
  }
  /* qsort(hits,n,sizeof(Stage3end_T),hit_sort_cmp); -- No need since original order was kept */

  debug7(
	 for (i = 0; i < n; i++) {
	   hit = hits[i];
	   printf("  Initial %d (%s): %p #%d:%u..%u, nmatches %d, score %d\n",
		  i,hittype_string(hit->hittype),hit,hit->chrnum,hit->genomicstart-hit->chroffset,hit->genomicend-hit->chroffset,
		  hit->nmatches,hit->score);
	 }
	 );

  i = 0;
  while (i < n) {
    debug7(printf("Looking at %d with score %d\n",i,hits[i]->score));
    j = i+1;
    while (j < n && hit_equiv_cmp(hits[j],hits[i]) == 0) {
      debug7(printf("  %d equal to %d\n",j,i));
      eliminate[j] = true;
      j++;
    }

    i = j;
  }

  for (i = n-1; i >= 0; i--) {
    hit = hits[i];
    if (eliminate[i] == false) {
      unique = List_push(unique,hit);
    } else if (hit->paired_usedp == true) {
      unique = List_push(unique,hit);
    } else {
      Stage3end_free(&hit);
    }
  }

#ifdef USE_ALLOCA_FOR_HITS
  FREEA(hits);
  FREEA(eliminate);
#else
  FREE(hits);
  FREE(eliminate);
#endif


#ifdef PRE_RESOLVE_MULTIMAPPING
  if (use_tally_p == true && tally_iit != NULL) {
    if ((n = List_length(unique)) > 1) {
#ifdef USE_ALLOCA_FOR_HITS
      hits = (T *) MALLOCA(n * sizeof(T));
      List_fill_array_and_free((void **) hits,&unique);
#else
      hits = (T *) List_to_array(unique,NULL);
      List_free(&unique);
#endif

      best_tally = 0;
      for (i = 0; i < n; i++) {
	if (hits[i]->tally < 0) {
	  hits[i]->tally = Stage3end_compute_tally(hits[i]);
	}
	if (hits[i]->tally > best_tally) {
	  best_tally = hits[i]->tally;
	}
      }

      unique = (List_T) NULL;
      for (i = 0; i < n; i++) {
	if (hits[i]->tally < best_tally) {
	  /* Stage3end_free(&(hits[i])); */
	} else {
	  unique = List_push(unique,(void *) hits[i]);
	}
      }

#ifdef USE_ALLOCA_FOR_HITS
      FREEA(hits);
#else
      FREE(hits);
#endif
    }
  }
#endif

  debug7(printf("Exited Stage3end_remove_overlaps with %d hits\n",List_length(unique)));
  return unique;
}


List_T
Stage3end_resolve_multimapping (List_T hits) {
  List_T resolve1, resolve2, resolve3, p;
  Stage3end_T hit;

  Overlap_T best_overlap;
  long int best_tally;
  double tally_threshold;
  bool runlengthp;

  if (List_length(hits) <= 1) {
    return hits;
  }

  if (genes_iit == NULL) {
    resolve1 = hits;
  } else {
    best_overlap = NO_KNOWN_GENE;
    for (p = hits; p != NULL; p = p->rest) {
      hit = (Stage3end_T) p->first;
      if ((hit->gene_overlap = Stage3end_gene_overlap(hit)) > best_overlap) {
	best_overlap = hit->gene_overlap;
      }
    }
    if (best_overlap == NO_KNOWN_GENE) {
      resolve1 = hits;
    } else {
      resolve1 = (List_T) NULL;
      for (p = hits; p != NULL; p = p->rest) {
	hit = (Stage3end_T) p->first;
	if (hit->gene_overlap < best_overlap) {
	  Stage3end_free(&hit);
	} else {
	  resolve1 = List_push(resolve1,(void *) hit);
	}
      }
      List_free(&hits);
    }
  }
      
  if (List_length(resolve1) <= 1) {
    return resolve1;
  }

  if (tally_iit == NULL) {
    resolve2 = resolve1;
  } else {
    best_tally = 0L;
    for (p = resolve1; p != NULL; p = p->rest) {
      hit = (Stage3end_T) p->first;
      if ((hit->tally = Stage3end_compute_tally(hit)) > best_tally) {
	best_tally = hit->tally;
      }
    }
    if (best_tally == 0L) {
      resolve2 = resolve1;
    } else {
      resolve2 = (List_T) NULL;
#ifdef USE_TALLY_RATIO
      tally_threshold = (double) best_tally / TALLY_RATIO;
#else
      tally_threshold = 1.0;
#endif
      for (p = resolve1; p != NULL; p = p->rest) {
	hit = (Stage3end_T) p->first;
	if ((double) hit->tally < tally_threshold) {
	  Stage3end_free(&hit);
	} else {
	  resolve2 = List_push(resolve2,(void *) hit);
	}
      }
      List_free(&resolve1);
    }
  }


  if (List_length(resolve2) <= 1) {
    return resolve2;
  }

  if (runlength_iit == NULL) {
    resolve3 = resolve2;
  } else {
    runlengthp = false;
    for (p = resolve2; p != NULL; p = p->rest) {
      hit = (Stage3end_T) p->first;
      if (Stage3end_runlength_p(hit) == true) {
	runlengthp = true;
      }
    }
    if (runlengthp == false) {
      resolve3 = resolve2;
    } else {
      resolve3 = (List_T) NULL;
      for (p = resolve2; p != NULL; p = p->rest) {
	hit = (Stage3end_T) p->first;
	if (Stage3end_runlength_p(hit) == false) {
	  Stage3end_free(&hit);
	} else {
	  resolve3 = List_push(resolve3,(void *) hit);
	}
      }
      List_free(&resolve2);
    }
  }


  return resolve3;
}


static void
print_alignment_info (Filestring_T fp, int nblocks, int score, int mapq_score, bool sarrayp) {
  FPRINTF(fp,"segs:%d,align_score:%d,mapq:%d",nblocks,score,mapq_score);
  if (sarrayp == true) {
    FPRINTF(fp,",method:sa");
  }
  return;
}


Pairtype_T
Stage3_determine_pairtype (T hit5, T hit3) {
  debug14(printf("Entered Stage3_determine_pairtype\n"));
  if (hit5->effective_chrnum != hit3->effective_chrnum) {
    debug14(printf("Returning unpaired\n"));
    return UNPAIRED;
  } else if (hit5->plusp != hit3->plusp) {
    debug14(printf("Returning paired_inversion\n"));
    return PAIRED_INVERSION;
  } else if (hit5->plusp == true) {
    if (hit3->genomicend < hit5->genomicstart) {
      debug14(printf("Returning paired_scramble\n"));
      return PAIRED_SCRAMBLE;
    } else if (hit3->genomicstart > hit5->genomicend + pairmax) {
      debug14(printf("Returning paired_toolong\n"));
      return PAIRED_TOOLONG;
    } else {
      debug14(printf("Returning concordant\n"));
      return CONCORDANT;
    }
  } else {
    if (hit3->genomicend > hit5->genomicstart) {
      debug14(printf("Returning paired_scramble\n"));
      return PAIRED_SCRAMBLE;
    } else if (hit3->genomicstart + pairmax < hit5->genomicend) {
      debug14(printf("Returning paired_toolong\n"));
      return PAIRED_TOOLONG;
    } else {
      debug14(printf("Returning concordant\n"));
      return CONCORDANT;
    }
  }
}


Pairtype_T
Stage3pair_pairtype (Stage3pair_T this) {
  return this->pairtype;
}

bool
Stage3pair_circularp (Stage3pair_T this) {
  return this->circularp;
}



#if 0
static char *
unpaired_type_text (T hit5, T hit3) {
  if (hit5->chrnum != hit3->chrnum) {
    return UNPAIRED_INTERCHROM_TEXT;
  } else if (hit5->plusp != hit3->plusp) {
    return PAIRED_INVERSION_TEXT;
  } else if (hit5->plusp == true) {
    if (hit3->genomicstart < hit5->genomicstart) {
      return PAIRED_SCRAMBLE_TEXT;
    } else {
      return UNPAIRED_TOOLONG_TEXT;
    }
  } else {
    if (hit5->genomicstart < hit3->genomicstart) {
      return PAIRED_SCRAMBLE_TEXT;
    } else {
      return UNPAIRED_TOOLONG_TEXT;
    }
  }
}
#endif


/* Has a copy in pair.c */
static void
print_pair_info (Filestring_T fp, T hit5, T hit3, int insertlength, int pairscore,
		 Pairtype_T pairtype) {

  assert(hit5->effective_chrnum == hit3->effective_chrnum); /* Same chromosomes */

#if 0
  /* Doesn't hold for paired (inversion) */
  assert(hit5->plusp == hit3->plusp);	/* Same direction */
#endif

  FPRINTF(fp,"pair_score:%d",pairscore);
  FPRINTF(fp,",insert_length:%d",insertlength);

  switch (pairtype) {
  case CONCORDANT: break;
  case PAIRED_SCRAMBLE: FPRINTF(fp,",pairtype:scramble"); break;
  case PAIRED_INVERSION: FPRINTF(fp,",pairtype:inversion"); break;
  case PAIRED_TOOLONG: FPRINTF(fp,",pairtype:toolong"); break;
  case CONCORDANT_TRANSLOCATIONS: break;
  case CONCORDANT_TERMINAL: break;
  case PAIRED_UNSPECIFIED: abort();
  case UNPAIRED: abort();
  case UNSPECIFIED: abort();
  }

  return;
}


static void
print_substrings (Filestring_T fp, T this, int score, Univ_IIT_T chromosome_iit, Shortread_T queryseq,
		  Shortread_T headerseq, char *acc_suffix, bool invertp, T hit5, T hit3, int insertlength,
		  int pairscore, Pairtype_T pairtype, int mapq_score) {
  char *single_chr, *chr;
  bool allocp, alloc1p, pairinfo_printed_p = false;
  List_T substrings, junctions, p, q;
  Substring_T substring;
  Junction_T pre_junction, post_junction;
  int nblocks;

  if (this->chrnum == 0) {
    single_chr = (char *) NULL;
  } else {
    single_chr = Univ_IIT_label(chromosome_iit,this->chrnum,&alloc1p);
  }
  if (invertp == true) {
    substrings = this->substrings_Nto1;
    junctions = this->junctions_Nto1;
  } else {
    substrings = this->substrings_1toN;
    junctions = this->junctions_1toN;
  }

  if (print_m8_p) {
    for (p = substrings; p != NULL; p = List_next(p)) {
      substring = (Substring_T) List_head(p);
      if (Substring_ambiguous_p(substring) == true) {
	/* Skip */
      } else {
	if ((chr = single_chr) == NULL) {
	  chr = Univ_IIT_label(chromosome_iit,Substring_chrnum(substring),&allocp);
	}
	Substring_print_m8(fp,substring,headerseq,acc_suffix,chr,invertp);
	if (single_chr == NULL && allocp == true) {
	  FREE(chr);
	}
      }
    }

  } else {
    if ((nblocks = List_length(substrings)) == 1) {
      post_junction = (Junction_T) NULL;
    } else {
      post_junction = (Junction_T) List_head(junctions);
    }
    substring = (Substring_T) List_head(substrings);
    if (Substring_ambiguous_p(substring) == true) {
      nblocks -= 1;
    }
    substring = (Substring_T) List_last_value(substrings);
    if (Substring_ambiguous_p(substring) == true) {
      nblocks -= 1;
    }


    /* First line */
    substring = (Substring_T) List_head(substrings);
    if (Substring_ambiguous_p(substring) == true) {
      /* Skip */
    } else {
      if ((chr = single_chr) == NULL) {
	chr = Univ_IIT_label(chromosome_iit,Substring_chrnum(substring),&allocp);
      }
      FPRINTF(fp," ");
      Substring_print_alignment(fp,/*pre_junction*/NULL,substring,post_junction,queryseq,genome,chr,invertp);
      if (single_chr == NULL && allocp == true) {
	FREE(chr);
      }

      /* Alignment info */
      FPRINTF(fp,"\t");
      print_alignment_info(fp,nblocks,score,mapq_score,this->sarrayp);
    
      /* Pairing info */
      if (hit5 != NULL && hit3 != NULL) {
	FPRINTF(fp,"\t");
	print_pair_info(fp,hit5,hit3,insertlength,pairscore,pairtype);
      }
      pairinfo_printed_p = true;

      FPRINTF(fp,"\n");
    }

    if ((p = List_next(substrings)) == NULL) {
      /* Done */
    } else {
      /* Middle lines */
      for (q = List_next(junctions); q != NULL; p = List_next(p), q = List_next(q)) {
	pre_junction = post_junction;
	post_junction = List_head(q);

	substring = (Substring_T) List_head(p);
	if (Substring_ambiguous_p(substring) == true) {
	  /* Skip */
	} else {
	  if (pairinfo_printed_p == true) {
	    FPRINTF(fp,",");
	  } else {
	    FPRINTF(fp," ");
	  }
	  if ((chr = single_chr) == NULL) {
	    chr = Univ_IIT_label(chromosome_iit,Substring_chrnum(substring),&allocp);
	  }
	  Substring_print_alignment(fp,pre_junction,substring,post_junction,queryseq,genome,chr,invertp);
	  if (single_chr == NULL && allocp == true) {
	    FREE(chr);
	  }

	  if (pairinfo_printed_p == false) {
	    /* Alignment info */
	    FPRINTF(fp,"\t");
	    print_alignment_info(fp,nblocks,score,mapq_score,this->sarrayp);
    
	    /* Pairing info */
	    if (hit5 != NULL && hit3 != NULL) {
	      FPRINTF(fp,"\t");
	      print_pair_info(fp,hit5,hit3,insertlength,pairscore,pairtype);
	    }
	    pairinfo_printed_p = true;
	  }
	  FPRINTF(fp,"\n");
	}
      }

      /* Last line */
      pre_junction = post_junction;

      substring = (Substring_T) List_head(p);
      if (Substring_ambiguous_p(substring) == true) {
	/* Skip */
      } else {
	if (pairinfo_printed_p == true) {
	  FPRINTF(fp,",");
	} else {
	  FPRINTF(fp," ");
	}
	if ((chr = single_chr) == NULL) {
	  chr = Univ_IIT_label(chromosome_iit,Substring_chrnum(substring),&allocp);
	}
	Substring_print_alignment(fp,pre_junction,substring,/*post_junction*/NULL,queryseq,genome,chr,invertp);
	if (single_chr == NULL && allocp == true) {
	  FREE(chr);
	}

	if (pairinfo_printed_p == false) {
	  /* Alignment info */
	  FPRINTF(fp,"\t");
	  print_alignment_info(fp,nblocks,score,mapq_score,this->sarrayp);
	  
	  /* Pairing info */
	  if (hit5 != NULL && hit3 != NULL) {
	    FPRINTF(fp,"\t");
	    print_pair_info(fp,hit5,hit3,insertlength,pairscore,pairtype);
	  }
	  /* pairinfo_printed_p = true; */
	}
	FPRINTF(fp,"\n");
      }
    }
  }

  if (alloc1p == true) {
    FREE(single_chr);
  }
}



/* May substitute paired-end loglik for single-end loglik */
void
Stage3end_print (Filestring_T fp, T this, int score,
		 Univ_IIT_T chromosome_iit, Shortread_T queryseq, Shortread_T headerseq,
		 char *acc_suffix, bool invertp, T hit5, T hit3, int insertlength,
		 int pairscore, Pairtype_T pairtype, int mapq_score) {
  bool pairedp;

  if (this->hittype == GMAP) {
    if (print_m8_p) {
      Pair_print_m8(fp,this->pairarray,this->npairs,/*invertedp*/false,
		    this->chrnum,queryseq,headerseq,acc_suffix,chromosome_iit);
    } else if (Shortread_invertedp(queryseq) == false) {
      if (pairtype == UNPAIRED) {
	pairedp = false;
      } else {
	pairedp = true;
      }
      Substring_print_gmap(fp,this->pairarray,this->npairs,this->nsegments,/*invertedp*/false,
			   this->gmap_start_endtype,this->gmap_end_endtype,
			   this->chrnum,this->chroffset,this->chrhigh,Shortread_fulllength(queryseq),
			   this->plusp,this->gmap_cdna_direction,this->score,insertlength,pairscore,mapq_score,
			   chromosome_iit,pairedp,this->gmap_source);
    } else {
      if (pairtype == UNPAIRED) {
	pairedp = false;
      } else {
	pairedp = true;
      }
      Substring_print_gmap(fp,this->pairarray,this->npairs,this->nsegments,/*invertedp*/true,
			   this->gmap_end_endtype,this->gmap_start_endtype,
			   this->chrnum,this->chroffset,this->chrhigh,Shortread_fulllength(queryseq),
			   this->plusp,this->gmap_cdna_direction,this->score,insertlength,pairscore,mapq_score,
			   chromosome_iit,pairedp,this->gmap_source);
    }

  } else {
    print_substrings(fp,this,score,chromosome_iit,queryseq,headerseq,acc_suffix,invertp,
		     hit5,hit3,insertlength,pairscore,pairtype,mapq_score);
  }

  return;
}


static void
print_query_header (Filestring_T fp, char initchar, Shortread_T queryseq, bool invertp) {
  FPRINTF(fp,"%c",initchar);
  if (invertp == false) {
    Shortread_print_oneline(fp,queryseq);
  } else {
    Shortread_print_oneline_revcomp(fp,queryseq);
  }

  return;
}



static void
print_barcode_and_quality (Filestring_T fp, Shortread_T queryseq, bool invertp, int quality_shift) {
  char *barcode;

  if ((barcode = Shortread_barcode(queryseq)) != NULL) {
    FPRINTF(fp,"\tbarcode:%s",barcode);
  }

  if (Shortread_quality_string(queryseq) != NULL) {
    FPRINTF(fp,"\t");
    if (invertp == false) {
      Shortread_print_quality(fp,queryseq,/*hardclip_low*/0,/*hardclip_high*/0,
			      quality_shift,/*show_chopped_p*/true);
    } else {
      Shortread_print_quality_revcomp(fp,queryseq,/*hardclip_low*/0,/*hardclip_high*/0,
				      quality_shift,/*show_chopped_p*/true);
    }
  }

  return;
}


void
Stage3pair_print_end (Filestring_T fp, Filestring_T fp_failedinput,
		      Result_T result, Resulttype_T resulttype,
		      char initchar, bool firstp, Univ_IIT_T chromosome_iit,
		      Shortread_T queryseq, Shortread_T headerseq1, Shortread_T headerseq2,
		      int maxpaths, bool quiet_if_excessive_p,
		      bool invertp, int quality_shift) {
  Stage3pair_T *stage3pairarray, stage3pair;
  T *stage3array, *stage3array_mate, this, hit5, hit3;
  int npaths, npaths_mate, pathnum, first_absmq, second_absmq;
  bool excessivep, translocationp;


  if (resulttype == PAIREDEND_NOMAPPING) {
    if (print_m8_p == false) {
      Filestring_set_split_output(fp,OUTPUT_NM);
      print_query_header(fp,initchar,queryseq,invertp);
      FPRINTF(fp,"\t0 %s",UNPAIRED_TEXT);

      print_barcode_and_quality(fp,queryseq,invertp,quality_shift);
    
      FPRINTF(fp,"\t");
      Shortread_print_header(fp,headerseq1,headerseq2);
      FPRINTF(fp,"\n");
    }
    /* If failedinput_root != NULL, then this case is handled by calling procedure */

  } else if (resulttype == CONCORDANT_UNIQ) {
    stage3pairarray = (Stage3pair_T *) Result_array(&npaths,&first_absmq,&second_absmq,result);
    stage3pair = stage3pairarray[0];
    hit5 = stage3pair->hit5;
    hit3 = stage3pair->hit3;

    if (stage3pair->circularp == true) {
      Filestring_set_split_output(fp,OUTPUT_CC);
    } else {
      Filestring_set_split_output(fp,OUTPUT_CU);
    }

    if (print_m8_p == false) {
      print_query_header(fp,initchar,queryseq,invertp);
      FPRINTF(fp,"\t1 %s",CONCORDANT_TEXT);
    
      print_barcode_and_quality(fp,queryseq,invertp,quality_shift);

      FPRINTF(fp,"\t");
      Shortread_print_header(fp,headerseq1,headerseq2);
    }
    
    if (firstp == true) {
      Stage3end_print(fp,hit5,hit5->score,
		      chromosome_iit,queryseq,headerseq1,/*acc_suffix*/"/1",
		      invertp,hit5,hit3,stage3pair->insertlength,
		      stage3pair->score,stage3pair->pairtype,stage3pair->mapq_score);
    } else {
      Stage3end_print(fp,hit3,hit3->score,
		      chromosome_iit,queryseq,headerseq1,/*acc_suffix*/"/2",
		      invertp,hit5,hit3,stage3pair->insertlength,
		      stage3pair->score,stage3pair->pairtype,stage3pair->mapq_score);
    }

    if (print_m8_p == false) {
      FPRINTF(fp,"\n");
    }

  } else if (resulttype == CONCORDANT_TRANSLOC) {
    Filestring_set_split_output(fp,OUTPUT_CT);
    stage3pairarray = (Stage3pair_T *) Result_array(&npaths,&first_absmq,&second_absmq,result);

    if (quiet_if_excessive_p && npaths > maxpaths) {
      if (print_m8_p == false) {
	/* No xs category for transloc, so ignore quiet-if-excessive_p */
	print_query_header(fp,initchar,queryseq,invertp);
	FPRINTF(fp,"\t%d %s",npaths,CONCORDANT_TEXT);
	FPRINTF(fp," (transloc)");
	
	print_barcode_and_quality(fp,queryseq,invertp,quality_shift);
      
	FPRINTF(fp,"\t");
	Shortread_print_header(fp,headerseq1,headerseq2);

	/* No further output */
	FPRINTF(fp,"\n");
      }

      if (failedinput_root != NULL) {
	Shortread_print_query_singleend(fp_failedinput,queryseq,headerseq1);
      }

    } else {
      if (print_m8_p == false) {
	print_query_header(fp,initchar,queryseq,invertp);
	FPRINTF(fp,"\t%d %s",npaths,CONCORDANT_TEXT);
	FPRINTF(fp," (transloc)");

	print_barcode_and_quality(fp,queryseq,invertp,quality_shift);
      
	FPRINTF(fp,"\t");
	Shortread_print_header(fp,headerseq1,headerseq2);
      }

      for (pathnum = 1; pathnum <= npaths && pathnum <= maxpaths; pathnum++) {
	stage3pair = stage3pairarray[pathnum-1];
	hit5 = stage3pair->hit5;
	hit3 = stage3pair->hit3;

	if (firstp == true) {
	  Stage3end_print(fp,hit5,hit5->score,
			  chromosome_iit,queryseq,headerseq1,/*acc_suffix*/"/1",
			  invertp,hit5,hit3,stage3pair->insertlength,
			  stage3pair->score,stage3pair->pairtype,stage3pair->mapq_score);
	} else {
	  Stage3end_print(fp,hit3,hit3->score,
			  chromosome_iit,queryseq,headerseq1,/*acc_suffix*/"/2",
			  invertp,hit5,hit3,stage3pair->insertlength,
			  stage3pair->score,stage3pair->pairtype,stage3pair->mapq_score);
	}
      }

      if (print_m8_p == false) {
	FPRINTF(fp,"\n");
      }
    }


  } else if (resulttype == CONCORDANT_MULT) {
    stage3pairarray = (Stage3pair_T *) Result_array(&npaths,&first_absmq,&second_absmq,result);

    if (quiet_if_excessive_p && npaths > maxpaths) {
      Filestring_set_split_output(fp,OUTPUT_CX);
      if (print_m8_p == false) {
	print_query_header(fp,initchar,queryseq,invertp);
	FPRINTF(fp,"\t%d %s",npaths,CONCORDANT_TEXT);

	print_barcode_and_quality(fp,queryseq,invertp,quality_shift);

	FPRINTF(fp,"\t");
	Shortread_print_header(fp,headerseq1,headerseq2);

	/* No further output */
	FPRINTF(fp,"\n");

	if (failedinput_root != NULL) {
	  Shortread_print_query_singleend(fp_failedinput,queryseq,headerseq1);
	}
      }

    } else {
      Filestring_set_split_output(fp,OUTPUT_CM);
      if (print_m8_p == false) {
	print_query_header(fp,initchar,queryseq,invertp);
	FPRINTF(fp,"\t%d %s",npaths,CONCORDANT_TEXT);
	
	print_barcode_and_quality(fp,queryseq,invertp,quality_shift);

	FPRINTF(fp,"\t");
	Shortread_print_header(fp,headerseq1,headerseq2);
      }

      for (pathnum = 1; pathnum <= npaths && pathnum <= maxpaths; pathnum++) {
	stage3pair = stage3pairarray[pathnum-1];
	hit5 = stage3pair->hit5;
	hit3 = stage3pair->hit3;

	if (firstp == true) {
	  Stage3end_print(fp,hit5,hit5->score,
			  chromosome_iit,queryseq,headerseq1,/*acc_suffix*/"/1",
			  invertp,hit5,hit3,stage3pair->insertlength,
			  stage3pair->score,stage3pair->pairtype,stage3pair->mapq_score);
	} else {
	  Stage3end_print(fp,hit3,hit3->score,
			  chromosome_iit,queryseq,headerseq1,/*acc_suffix*/"/2",
			  invertp,hit5,hit3,stage3pair->insertlength,
			  stage3pair->score,stage3pair->pairtype,stage3pair->mapq_score);
	}
      }

      if (print_m8_p == false) {
	FPRINTF(fp,"\n");
      }
    }


  } else if (resulttype == PAIRED_UNIQ) {
    stage3pairarray = (Stage3pair_T *) Result_array(&npaths,&first_absmq,&second_absmq,result);
    stage3pair = stage3pairarray[0];

    if (stage3pair->circularp == true) {
      Filestring_set_split_output(fp,OUTPUT_PC);
    } else if (stage3pair->pairtype == PAIRED_INVERSION) {
      Filestring_set_split_output(fp,OUTPUT_PI);
    } else if (stage3pair->pairtype == PAIRED_SCRAMBLE) {
      Filestring_set_split_output(fp,OUTPUT_PS);
    } else if (stage3pair->pairtype == PAIRED_TOOLONG) {
      Filestring_set_split_output(fp,OUTPUT_PL);
    } else {
      fprintf(stderr,"Unexpected pairtype %d\n",stage3pair->pairtype);
      abort();
    }
    
    if (print_m8_p == false) {
      print_query_header(fp,initchar,queryseq,invertp);
      FPRINTF(fp,"\t1 %s",PAIRED_TEXT);

      print_barcode_and_quality(fp,queryseq,invertp,quality_shift);

      FPRINTF(fp,"\t");
      Shortread_print_header(fp,headerseq1,headerseq2);
    }

    hit5 = stage3pair->hit5;
    hit3 = stage3pair->hit3;
    
    if (firstp == true) {
      Stage3end_print(fp,hit5,hit5->score,
		      chromosome_iit,queryseq,headerseq1,/*acc_suffix*/"/1",
		      invertp,hit5,hit3,stage3pair->insertlength,
		      stage3pair->score,stage3pair->pairtype,stage3pair->mapq_score);
    } else {
      Stage3end_print(fp,hit3,hit3->score,
		      chromosome_iit,queryseq,headerseq1,/*acc_suffix*/"/2",
		      invertp,hit5,hit3,stage3pair->insertlength,
		      stage3pair->score,stage3pair->pairtype,stage3pair->mapq_score);
    }

    if (print_m8_p == false) {
      FPRINTF(fp,"\n");
    }

  } else if (resulttype == PAIRED_MULT) {
    stage3pairarray = (Stage3pair_T *) Result_array(&npaths,&first_absmq,&second_absmq,result);

    if (quiet_if_excessive_p && npaths > maxpaths) {
      Filestring_set_split_output(fp,OUTPUT_PX);
      if (print_m8_p == false) {
	print_query_header(fp,initchar,queryseq,invertp);
	FPRINTF(fp,"\t%d %s",npaths,PAIRED_TEXT);
	
	print_barcode_and_quality(fp,queryseq,invertp,quality_shift);

	FPRINTF(fp,"\t");
	Shortread_print_header(fp,headerseq1,headerseq2);

	/* No further output */
	FPRINTF(fp,"\n");

	if (failedinput_root != NULL) {
	  Shortread_print_query_singleend(fp_failedinput,queryseq,headerseq1);
	}
      }

    } else {
      Filestring_set_split_output(fp,OUTPUT_PM);
      if (print_m8_p == false) {
	print_query_header(fp,initchar,queryseq,invertp);
	FPRINTF(fp,"\t%d %s",npaths,PAIRED_TEXT);

	print_barcode_and_quality(fp,queryseq,invertp,quality_shift);

	FPRINTF(fp,"\t");
	Shortread_print_header(fp,headerseq1,headerseq2);
      }

      for (pathnum = 1; pathnum <= npaths && pathnum <= maxpaths; pathnum++) {
	stage3pair = stage3pairarray[pathnum-1];
	hit5 = stage3pair->hit5;
	hit3 = stage3pair->hit3;

	if (firstp == true) {
	  Stage3end_print(fp,hit5,hit5->score,
			  chromosome_iit,queryseq,headerseq1,/*acc_suffix*/"/1",
			  invertp,hit5,hit3,stage3pair->insertlength,
			  stage3pair->score,stage3pair->pairtype,stage3pair->mapq_score);
	} else {
	  Stage3end_print(fp,hit3,hit3->score,
			  chromosome_iit,queryseq,headerseq1,/*acc_suffix*/"/2",
			  invertp,hit5,hit3,stage3pair->insertlength,
			  stage3pair->score,stage3pair->pairtype,stage3pair->mapq_score);
	}
      }

      if (print_m8_p == false) {
	FPRINTF(fp,"\n");
      }
    }


  } else {
    /* Print as singles */
    if (firstp == true) {
      /* Get stage3array_mate first to avoid incorrect values for npaths */
      stage3array_mate = (T *) Result_array2(&npaths_mate,&first_absmq,&second_absmq,result);
      stage3array = (T *) Result_array(&npaths,&first_absmq,&second_absmq,result);
    } else {
      /* Get stage3array_mate first to avoid incorrect values for npaths */
      stage3array_mate = (T *) Result_array(&npaths_mate,&first_absmq,&second_absmq,result);
      stage3array = (T *) Result_array2(&npaths,&first_absmq,&second_absmq,result);
    }

    excessivep = false;
    translocationp = false;
    if (resulttype == HALFMAPPING_UNIQ) {
      if (npaths > 0 && Stage3end_circularpos(stage3array[0]) > 0) {
	Filestring_set_split_output(fp,OUTPUT_HC);
      } else if (npaths_mate > 0 && Stage3end_circularpos(stage3array_mate[0]) > 0) {
	Filestring_set_split_output(fp,OUTPUT_HC);
      } else {
	Filestring_set_split_output(fp,OUTPUT_HU);
      }

    } else if (resulttype == HALFMAPPING_TRANSLOC) {
      Filestring_set_split_output(fp,OUTPUT_HT);
      translocationp = true;

    } else if (resulttype == HALFMAPPING_MULT) {
      if (quiet_if_excessive_p && npaths > maxpaths) {
	Filestring_set_split_output(fp,OUTPUT_HX);
	excessivep = true;
      } else {
	Filestring_set_split_output(fp,OUTPUT_HM);
      }

    } else if (resulttype == UNPAIRED_UNIQ) {
      if (npaths > 0 && Stage3end_circularpos(stage3array[0]) > 0) {
	Filestring_set_split_output(fp,OUTPUT_UC);
      } else if (npaths_mate > 0 && Stage3end_circularpos(stage3array_mate[0]) > 0) {
	Filestring_set_split_output(fp,OUTPUT_UC);
      } else {
	Filestring_set_split_output(fp,OUTPUT_UU);
      }

    } else if (resulttype == UNPAIRED_TRANSLOC) {
      Filestring_set_split_output(fp,OUTPUT_UT);
      translocationp = true;

    } else if (resulttype == UNPAIRED_MULT) {
      if (quiet_if_excessive_p && npaths > maxpaths) {
	Filestring_set_split_output(fp,OUTPUT_UX);
	excessivep = true;
      } else {
	Filestring_set_split_output(fp,OUTPUT_UM);
      }

    } else {
      fprintf(stderr,"Resulttype is %s\n",Resulttype_string(resulttype));
      abort();
    }

    if (print_m8_p == false) {
      print_query_header(fp,initchar,queryseq,invertp);
      FPRINTF(fp,"\t%d %s",npaths,UNPAIRED_TEXT);
      if (translocationp == true) {
	FPRINTF(fp," (transloc)");
      }

      print_barcode_and_quality(fp,queryseq,invertp,quality_shift);

      FPRINTF(fp,"\t");
      Shortread_print_header(fp,headerseq1,headerseq2);
    }

    if (excessivep == true) {
      /* No output */
      if (failedinput_root != NULL) {
	Shortread_print_query_singleend(fp_failedinput,queryseq,headerseq1);
      }
					      
    } else {
      if (firstp == true) {
	for (pathnum = 1; pathnum <= npaths && pathnum <= maxpaths; pathnum++) {
	  this = stage3array[pathnum-1];
	  Stage3end_print(fp,this,this->score,
			  chromosome_iit,queryseq,headerseq1,/*acc_suffix*/"/1",
			  invertp,/*hit5*/(T) NULL,/*hit3*/(T) NULL,
			  /*insertlength*/0,/*pairscore*/0,/*pairtype*/UNPAIRED,this->mapq_score);
	}
      } else {
	for (pathnum = 1; pathnum <= npaths && pathnum <= maxpaths; pathnum++) {
	  this = stage3array[pathnum-1];
	  Stage3end_print(fp,this,this->score,
			  chromosome_iit,queryseq,headerseq1,/*acc_suffix*/"/2",
			  invertp,/*hit5*/(T) NULL,/*hit3*/(T) NULL,
			  /*insertlength*/0,/*pairscore*/0,/*pairtype*/UNPAIRED,this->mapq_score);
	}
      }
    }

    if (print_m8_p == false) {
      FPRINTF(fp,"\n");
    }
  }


  return;
}


static List_T
strip_gaps_at_head (List_T pairs) {
  Pair_T pair;

  while (pairs != NULL && (pair = pairs->first) != NULL &&
	 (pair->gapp == true || pair->cdna == ' ' || pair->genome == ' ')) {
    pairs = List_pop_out(pairs,(void **) &pair);
    Pair_free_out(&pair);
  }

  return pairs;
}

static List_T
strip_gaps_at_tail (List_T pairs) {
  Pair_T pair;

  if (pairs != NULL) {
    pairs = List_reverse(pairs);

    while (pairs != NULL && (pair = pairs->first) != NULL &&
	   (pair->gapp == true || pair->cdna == ' ' || pair->genome == ' ')) {
      pairs = List_pop_out(pairs,(void **) &pair);
      Pair_free_out(&pair);
    }

    pairs = List_reverse(pairs);
  }

  return pairs;
}


/* If use querylength_adj, ss.bug.4 fails.  If use querylength, ss.bug.3 fails */
static List_T
Stage3end_convert_to_pairs (List_T pairs, T hit, Shortread_T queryseq,
			    int hardclip_low, int hardclip_high, int queryseq_offset) {
  Pair_T pair;
  List_T newpairs = NULL, p, q;
  Chrpos_T genomicpos1, genomicpos2;
  Substring_T substring, prev_substring;
  Junction_T junction;
  Junctiontype_T type;
  char *deletion_string;

  if (hit->hittype == TRANSLOC_SPLICE) {
    /* Cannot handle translocations within a single GMAP alignment */
    abort();
    return NULL;
    
  } else if (hit->hittype == GMAP) {
    debug15(printf("Converting gmap to pairs\n"));
    /* Use querylength here, but adj elsewhere */
    return Pair_convert_array_to_pairs(pairs,hit->pairarray,hit->npairs,hit->plusp,hit->querylength,
				       hardclip_low,hardclip_high,queryseq_offset);
  } else {
    p = hit->substrings_1toN;
    prev_substring = (Substring_T) List_head(p);
    pairs = Substring_convert_to_pairs(pairs,prev_substring,hit->querylength,
				       queryseq,hardclip_low,hardclip_high,queryseq_offset);

    for (q = hit->junctions_1toN, p = List_next(p); p != NULL; q = List_next(q), p = List_next(p)) {
      junction = (Junction_T) List_head(q);
      substring = (Substring_T) List_head(p);
    
      if ((type = Junction_type(junction)) == INS_JUNCTION) {
	pairs = Substring_add_insertion(pairs,prev_substring,substring,hit->querylength,
					/*insertionlength*/Junction_nindels(junction),queryseq,
					hardclip_low,hardclip_high,queryseq_offset);
      } else if (type == DEL_JUNCTION) {
	deletion_string = Junction_deletion_string(junction,genome,hit->plusp);
	pairs = Substring_add_deletion(pairs,prev_substring,substring,hit->querylength,
				       deletion_string,/*deletionlength*/Junction_nindels(junction),
				       hardclip_low,hardclip_high,queryseq_offset);
      } else if (type == SPLICE_JUNCTION) {
	pairs = Substring_add_intron(pairs,prev_substring,substring,hit->querylength,
				     hardclip_low,hardclip_high,queryseq_offset);
	
      } else {
	abort();
      }
    
      pairs = Substring_convert_to_pairs(pairs,substring,hit->querylength,
					 queryseq,hardclip_low,hardclip_high,queryseq_offset);
      prev_substring = substring;
    }

    debug15(Pair_dump_list(pairs,true));
    return pairs;
  }
}


/* Don't want querylength_adj */
struct Pair_T *
Stage3pair_merge (int *npairs, int *querylength_merged, char **queryseq_merged, char **quality_merged,
		  Stage3pair_T this, Shortread_T queryseq5, Shortread_T queryseq3,
		  int querylength5, int querylength3, int clipdir,
		  int hardclip5_low, int hardclip5_high, int hardclip3_low, int hardclip3_high) {
  struct Pair_T *pairarray, *newpair;
  Pair_T oldpair;
  List_T pairs, pairs5, pairs3, p;
  T hit5, hit3;
  int querylengthA, querylengthB;
  char *queryseq_ptr_5, *queryseq_ptr_3, *quality_ptr_5, *quality_ptr_3;
#ifdef CHECK_ASSERTIONS
  Chrpos_T genomicpos1, genomicpos2;
#endif

  hit5 = this->hit5;
  hit3 = this->hit3;
  queryseq_ptr_5 = Shortread_fullpointer_uc(queryseq5);
  queryseq_ptr_3 = Shortread_fullpointer_uc(queryseq3);
  quality_ptr_5 = Shortread_quality_string(queryseq5);
  quality_ptr_3 = Shortread_quality_string(queryseq3);

  if (hit5->plusp == true) {
    if (clipdir > 0) {
      pairs5 = Stage3end_convert_to_pairs(NULL,hit5,queryseq5,hardclip5_low,hardclip5_high,/*queryseq_offset*/0);
      pairs5 = strip_gaps_at_head(pairs5);

      pairs3 = Stage3end_convert_to_pairs(NULL,hit3,queryseq3,hardclip3_low,hardclip3_high,
					  /*queryseq_offset*/querylength5-hardclip5_low-hardclip5_high-hardclip3_low-hardclip3_high);
      pairs3 = strip_gaps_at_tail(pairs3);

#ifdef CHECK_ASSERTIONS
      genomicpos1 = ((Pair_T) List_head(pairs5))->genomepos;
      genomicpos2 = ((Pair_T) List_last_value(pairs3))->genomepos;
      if (genomicpos2 != genomicpos1 + 1U) {
	printf("Accession %s, plus\n",Shortread_accession(queryseq5));
	printf("Expected genomicpos2 %u == genomicpos1 %u + 1\n",genomicpos2,genomicpos1);
	Pair_dump_list(pairs5,true);
	Pair_dump_list(pairs3,true);
	abort();
      }
#endif      

      pairs = List_append(pairs3,pairs5);

      querylengthA = querylength5 - hardclip5_low - hardclip5_high;
      querylengthB = querylength3 - hardclip3_low - hardclip3_high;
      *querylength_merged = querylengthA + querylengthB;

      *queryseq_merged = (char *) MALLOC_OUT((querylengthA+querylengthB+1) * sizeof(char));
      strncpy(*queryseq_merged,queryseq_ptr_5,querylengthA);
      strncpy(&((*queryseq_merged)[querylengthA]),&(queryseq_ptr_3[querylength3 - querylengthB]),querylengthB);

      if (quality_ptr_5 == NULL || quality_ptr_3 == NULL) {
	*quality_merged = (char *) NULL;
      } else {
	*quality_merged = (char *) MALLOC_OUT((querylengthA+querylengthB+1) * sizeof(char));
	strncpy(*quality_merged,quality_ptr_5,querylengthA);
	strncpy(&((*quality_merged)[querylengthA]),&(quality_ptr_3[querylength3 - querylengthB]),querylengthB);
      }

    } else if (clipdir < 0) {
      pairs3 = Stage3end_convert_to_pairs(NULL,hit3,queryseq3,hardclip3_low,hardclip3_high,/*queryseq_offset*/0);
      pairs3 = strip_gaps_at_head(pairs3);

      pairs5 = Stage3end_convert_to_pairs(NULL,hit5,queryseq5,hardclip5_low,hardclip5_high,
					  /*queryseq_offset*/querylength3-hardclip3_low-hardclip3_high-hardclip5_low-hardclip5_high);
      pairs5 = strip_gaps_at_tail(pairs5);

#ifdef CHECK_ASSERTIONS
      genomicpos1 = ((Pair_T) List_head(pairs3))->genomepos;
      genomicpos2 = ((Pair_T) List_last_value(pairs5))->genomepos;
      if (genomicpos2 != genomicpos1 + 1U) {
	printf("Accession %s, plus, clipdir %d\n",Shortread_accession(queryseq5),clipdir);
	printf("Expected genomicpos2 %u == genomicpos1 %u + 1\n",genomicpos2,genomicpos1);
	printf("Begin of pairs3\n");
	Pair_dump_list(pairs3,true);
	printf("Begin of pairs5\n");
	Pair_dump_list(pairs5,true);
	abort();
      }
#endif      

      pairs = List_append(pairs5,pairs3);

      querylengthA = querylength3 - hardclip3_low - hardclip3_high;
      querylengthB = querylength5 - hardclip5_low - hardclip5_high;
      *querylength_merged = querylengthA + querylengthB;

      *queryseq_merged = (char *) MALLOC_OUT((querylengthA+querylengthB+1) * sizeof(char));
      strncpy(*queryseq_merged,queryseq_ptr_3,querylengthA);
      strncpy(&((*queryseq_merged)[querylengthA]),&(queryseq_ptr_5[querylength5 - querylengthB]),querylengthB);

      if (quality_ptr_5 == NULL || quality_ptr_3 == NULL) {
	*quality_merged = (char *) NULL;
      } else {
	*quality_merged = (char *) MALLOC_OUT((querylengthA+querylengthB+1) * sizeof(char));
	strncpy(*quality_merged,quality_ptr_3,querylengthA);
	strncpy(&((*quality_merged)[querylengthA]),&(quality_ptr_5[querylength5 - querylengthB]),querylengthB);
      }

    } else {
      abort();
    }

  } else {
    if (clipdir > 0) {
      pairs3 = Stage3end_convert_to_pairs(NULL,hit3,queryseq3,hardclip3_low,hardclip3_high,/*queryseq_offset*/0);
      pairs3 = strip_gaps_at_head(pairs3);

      pairs5 = Stage3end_convert_to_pairs(NULL,hit5,queryseq5,hardclip5_low,hardclip5_high,
					  /*queryseq_offset*/querylength3-hardclip3_low-hardclip3_high-hardclip5_low-hardclip5_high);
      pairs5 = strip_gaps_at_tail(pairs5);

#ifdef CHECK_ASSERTIONS
      genomicpos1 = ((Pair_T) List_head(pairs3))->genomepos;
      genomicpos2 = ((Pair_T) List_last_value(pairs5))->genomepos;
      if (genomicpos2 != genomicpos1 - 1U) {
	printf("Accession %s, minus\n",Shortread_accession(queryseq5));
	printf("Expected genomicpos2 %u == genomicpos1 %u - 1\n",genomicpos2,genomicpos1);
	Pair_dump_list(pairs3,true);
	Pair_dump_list(pairs5,true);
	abort();
      }
#endif      

      pairs = List_append(pairs5,pairs3);

      querylengthA = querylength3 - hardclip3_low - hardclip3_high;
      querylengthB = querylength5 - hardclip5_low - hardclip5_high;
      *querylength_merged = querylengthA + querylengthB;

      *queryseq_merged = (char *) MALLOC_OUT((querylengthA+querylengthB+1) * sizeof(char));
      strncpy(*queryseq_merged,queryseq_ptr_3,querylengthA);
      strncpy(&((*queryseq_merged)[querylengthA]),&(queryseq_ptr_5[querylength5 - querylengthB]),querylengthB);

      if (quality_ptr_5 == NULL || quality_ptr_3 == NULL) {
	*quality_merged = (char *) NULL;
      } else {
	*quality_merged = (char *) MALLOC_OUT((querylengthA+querylengthB+1) * sizeof(char));
	strncpy(*quality_merged,quality_ptr_3,querylengthA);
	strncpy(&((*quality_merged)[querylengthA]),&(quality_ptr_5[querylength5 - querylengthB]),querylengthB);
      }

    } else if (clipdir < 0) {
      pairs5 = Stage3end_convert_to_pairs(NULL,hit5,queryseq5,hardclip5_low,hardclip5_high,/*queryseq_offset*/0);
      pairs5 = strip_gaps_at_head(pairs5);

      pairs3 = Stage3end_convert_to_pairs(NULL,hit3,queryseq3,hardclip3_low,hardclip3_high,
					  /*queryseq_offset*/querylength5-hardclip5_low-hardclip5_high-hardclip3_low-hardclip3_high);
      pairs3 = strip_gaps_at_tail(pairs3);

#ifdef CHECK_ASSERTIONS
      genomicpos1 = ((Pair_T) List_head(pairs5))->genomepos;
      genomicpos2 = ((Pair_T) List_last_value(pairs3))->genomepos;
      if (genomicpos2 != genomicpos1 - 1U) {
	printf("Accession %s, minus\n",Shortread_accession(queryseq5));
	printf("Expected genomicpos2 %u == genomicpos1 %u - 1\n",genomicpos2,genomicpos1);
	Pair_dump_list(pairs5,true);
	Pair_dump_list(pairs3,true);
	abort();
      }
#endif      

      pairs = List_append(pairs3,pairs5);

      querylengthA = querylength5 - hardclip5_low - hardclip5_high;
      querylengthB = querylength3 - hardclip3_low - hardclip3_high;
      *querylength_merged = querylengthA + querylengthB;

      *queryseq_merged = (char *) MALLOC_OUT((querylengthA+querylengthB+1) * sizeof(char));
      strncpy(*queryseq_merged,queryseq_ptr_5,querylengthA);
      strncpy(&((*queryseq_merged)[querylengthA]),&(queryseq_ptr_3[querylength3 - querylengthB]),querylengthB);

      if (quality_ptr_5 == NULL || quality_ptr_3 == NULL) {
	*quality_merged = (char *) NULL;
      } else {
	*quality_merged = (char *) MALLOC_OUT((querylengthA+querylengthB+1) * sizeof(char));
	strncpy(*quality_merged,quality_ptr_5,querylengthA);
	strncpy(&((*quality_merged)[querylengthA]),&(quality_ptr_3[querylength3 - querylengthB]),querylengthB);
      }

    } else {
      abort();
    }
  }

  pairs = List_reverse(pairs);
  /* Pair_dump_list(pairs,true); */

  *npairs = List_length(pairs);
  newpair = pairarray = (struct Pair_T *) MALLOC_OUT((*npairs)*sizeof(struct Pair_T));
  for (p = pairs; p != NULL; p = p->rest) {
    oldpair = (Pair_T) p->first;
    memcpy(newpair++,oldpair,sizeof(struct Pair_T));
    Pair_free_out(&oldpair);
  }
  List_free_out(&pairs);

  return pairarray;
}



#if 0
List_T
Stage3end_filter_bymatch (List_T hitlist) {
  List_T filtered = NULL, p;
  T hit;
  int min_nmismatches_whole = 1000;

  for (p = hitlist; p != NULL; p = p->rest) {
    hit = (T) p->first;
    if (hit->nmismatches_whole < min_nmismatches_whole) {
      min_nmismatches_whole = hit->nmismatches_whole;
    }
  }

  for (p = hitlist; p != NULL; p = p->rest) {
    hit = (T) p->first;
    if (hit->nmismatches_whole == min_nmismatches_whole) {
      filtered = List_push(filtered,hit);
    } else {
      Stage3end_free(&hit);
    }
  }
  List_free(&hitlist);

  return filtered;
}
#endif



static Chrpos_T
overlap5_gmap_plus (int *querypos, Chrpos_T *genomicstart, Chrpos_T *genomicend,
		    Stage3end_T hit5, Stage3end_T gmap) {
  Chrpos_T chrpos;
  Substring_T substring;
  List_T p;

  debug10(printf("Entered overlap5_gmap_plus with gmap %d..%d\n",
		 gmap->pairarray[0].querypos,gmap->pairarray[gmap->npairs - 1].querypos));
  for (p = hit5->substrings_LtoH; p != NULL; p = List_next(p)) {
    substring = (Substring_T) List_head(p);
    if (Substring_ambiguous_p(substring) == false) {
      *genomicstart = Substring_alignstart_chr(substring);
      *genomicend = Substring_alignend_chr(substring);
      if ((chrpos = Pair_binary_search_ascending(&(*querypos),/*lowi*/0,/*highi*/gmap->npairs,gmap->pairarray,
						 *genomicstart,*genomicend)) > 0) {
	return chrpos;
      }
    }
  }

  return 0;
}

static Chrpos_T
overlap3_gmap_plus (int *querypos, Chrpos_T *genomicstart, Chrpos_T *genomicend,
		    Stage3end_T hit3, Stage3end_T gmap) {
  Chrpos_T chrpos;
  Substring_T substring;
  List_T p;

  debug10(printf("Entered overlap3_gmap_plus with gmap %d..%d\n",
		 gmap->pairarray[0].querypos,gmap->pairarray[gmap->npairs - 1].querypos));
  for (p = hit3->substrings_LtoH; p != NULL; p = List_next(p)) {
    substring = (Substring_T) List_head(p);
    if (Substring_ambiguous_p(substring) == false) {
      *genomicstart = Substring_alignstart_chr(substring);
      *genomicend = Substring_alignend_chr(substring);
      if ((chrpos = Pair_binary_search_ascending(&(*querypos),/*lowi*/0,/*highi*/gmap->npairs,gmap->pairarray,
						 *genomicstart,*genomicend)) > 0) {
	return chrpos;
      }
    }
  }

  return 0;
}


static Chrpos_T
overlap5_gmap_minus (int *querypos, Chrpos_T *genomicstart, Chrpos_T *genomicend,
		     Stage3end_T hit5, Stage3end_T gmap) {
  Chrpos_T chrpos;
  Substring_T substring;
  List_T p;

  debug10(printf("Entered overlap5_gmap_minus with gmap %d..%d\n",
		 gmap->pairarray[0].querypos,gmap->pairarray[gmap->npairs - 1].querypos));
  for (p = hit5->substrings_LtoH; p != NULL; p = List_next(p)) {
    substring = (Substring_T) List_head(p);
    if (Substring_ambiguous_p(substring) == false) {
      *genomicstart = Substring_alignstart_chr(substring);
      *genomicend = Substring_alignend_chr(substring);
      if ((chrpos = Pair_binary_search_descending(&(*querypos),/*lowi*/0,/*highi*/gmap->npairs,gmap->pairarray,
						  *genomicstart,*genomicend)) > 0) {
	return chrpos;
      }
    }
  }

  return 0;
}

static Chrpos_T
overlap3_gmap_minus (int *querypos, Chrpos_T *genomicstart, Chrpos_T *genomicend,
		     Stage3end_T hit3, Stage3end_T gmap) {
  Chrpos_T chrpos;
  Substring_T substring;
  List_T p;

  debug10(printf("Entered overlap3_gmap_minus with gmap %d..%d\n",
		 gmap->pairarray[0].querypos,gmap->pairarray[gmap->npairs - 1].querypos));
  for (p = hit3->substrings_LtoH; p != NULL; p = List_next(p)) {
    substring = (Substring_T) List_head(p);
    if (Substring_ambiguous_p(substring) == false) {
      *genomicstart = Substring_alignstart_chr(substring);
      *genomicend = Substring_alignend_chr(substring);
      if ((chrpos = Pair_binary_search_descending(&(*querypos),/*lowi*/0,/*highi*/gmap->npairs,gmap->pairarray,
						  *genomicstart,*genomicend)) > 0) {
	return chrpos;
      }
    }
  }

  return 0;
}


/* Should not set ambiguous flag in substrings, because resolution of
   an ambiguity depends on a particular pair of ends */

static void
resolve_inside_ambiguous_splice_plus (int *unresolved_amb_length, int *amb_resolve_5, int *amb_resolve_3,
				      int *amb_status_inside, T hit5, T hit3, Univcoord_T *splicesites,
				      Compress_T query5_compress_fwd, Compress_T query3_compress_fwd,
				      int localsplicing_penalty, int querylength5, int querylength3,
				      int genestrand) {
  int insertlength;
  Univcoord_T genomicstart, genomicend;
  int nbingo, bingoi5, bingoi3;
  int nbest, besti5, besti3, i, j;
  int best_nmismatches, nmismatches;

  Substring_T substring5, substring3;
#ifdef LARGE_GENOMES
  Uint8list_T ambcoords;
#else
  Uintlist_T ambcoords;
#endif
  Univcoord_T *end_ambcoords, *start_ambcoords;
  int *end_amb_nmismatches, *start_amb_nmismatches;
  int end_amb_length_5, start_amb_length_3;


  *unresolved_amb_length = 0;

  if (hit5->hittype == GMAP) {
    substring5 = (Substring_T) NULL;
  } else {
    substring5 = (Substring_T) List_head(hit5->substrings_Nto1);
  }
  if (hit3->hittype == GMAP) {
    substring3 = (Substring_T) NULL;
  } else {
    substring3 = (Substring_T) List_head(hit3->substrings_1toN);
  }
  debug9(printf("resolve plus: hit5 %s and hit3 %s\n",
		hittype_string(hit5->hittype),hittype_string(hit3->hittype)));

  if (substring5 != NULL && Substring_ambiguous_p(substring5) == true && 
      substring3 != NULL && Substring_ambiguous_p(substring3) == true) {
    debug9(printf("Resolve plus: Got ambiguous at 5' and ambiguous at 3':"));
    end_ambcoords = Substring_ambcoords(substring5);
    end_amb_nmismatches = Substring_amb_nmismatches(substring5);
    start_ambcoords = Substring_ambcoords(substring3);
    start_amb_nmismatches = Substring_amb_nmismatches(substring3);
    end_amb_length_5 = end_amb_length(hit5);
    start_amb_length_3 = start_amb_length(hit3);
    
    nbingo = nbest = 0;
    best_nmismatches = querylength5 + querylength3;
    for (i = 0; i < Substring_nambcoords(substring5); i++) {
      genomicend = end_ambcoords[i] + end_amb_length_5;
      for (j = 0; j < Substring_nambcoords(substring3); j++) {
	genomicstart = start_ambcoords[j] - start_amb_length_3;
	debug9(printf(" %u,%u",(Chrpos_T) (genomicend - hit5->chroffset),(Chrpos_T) (genomicstart - hit3->chroffset)));
	if (genomicend < genomicstart) {
	  /* Look for valid insertlength */
	  insertlength = genomicstart - genomicend + querylength5 + querylength3;
	  debug9(printf(" (%u)",insertlength));
	  if (insertlength >= expected_pairlength_low && insertlength <= expected_pairlength_high) {
	    nbingo++;
	    bingoi5 = i;
	    bingoi3 = j;
	    debug9(printf("*"));
	  }

	  if ((nmismatches = end_amb_nmismatches[i] + start_amb_nmismatches[j]) < best_nmismatches) {
	    best_nmismatches = nmismatches;
	    besti5 = i;
	    besti3 = j;
	    nbest = 1;
	  } else if (nmismatches == best_nmismatches) {
	    nbest++;
	  }
	}
      }
    }

    if (nbingo == 1) {
      *amb_resolve_5 = bingoi5;
      *amb_resolve_3 = bingoi3;
      *amb_status_inside = AMB_RESOLVED_BYLENGTH;

    } else if (nbest == 0) {
      debug9(printf("\nnbest is zero: amb_length = %d...%d",
		    end_amb_length(hit5),start_amb_length(hit3)));
      *unresolved_amb_length = end_amb_length_5 + start_amb_length_3;
      *amb_status_inside = AMB_UNRESOLVED_TOOCLOSE;

    } else if (nbest == 1) {
      debug9(printf("\nnbest is 1, with nmismatches %d\n",best_nmismatches));
      *amb_resolve_5 = besti5;
      *amb_resolve_3 = besti3;
      *amb_status_inside = AMB_RESOLVED_BYMATCHES;

    } else {
      *amb_resolve_5 = -1;	/* Signifies cannot resolve */
      *amb_resolve_3 = -1;	/* Signifies cannot resolve */
      *amb_status_inside = AMB_UNRESOLVED_MULTIPLE;
    }
    debug9(printf("\n"));

  } else if (substring5 != NULL && Substring_ambiguous_p(substring5) == true) {
    debug9(printf("hit3 %u..%u\n",hit3->genomicstart - hit3->chroffset,hit3->genomicend - hit3->chroffset));
    debug9(printf("Resolve plus: Got ambiguous at 5' (%s):",hittype_string(hit5->hittype)));
    debug9(printf(" (?< %u):",hit3->genomicstart + querylength3 - hit3->chroffset));
    end_ambcoords = Substring_ambcoords(substring5);
    end_amb_nmismatches = Substring_amb_nmismatches(substring5);
    end_amb_length_5 = end_amb_length(hit5);

    nbingo = nbest = 0;
    best_nmismatches = querylength5;
    for (i = 0; i < Substring_nambcoords(substring5); i++) {
      genomicend = end_ambcoords[i] + end_amb_length_5;
      debug9(printf(" %u (%d mismatches)",(Chrpos_T) (genomicend - hit5->chroffset),end_amb_nmismatches[i]));
      if (genomicend < hit3->genomicstart /*allow overlap*/+ querylength3) {
	/* Look for valid insertlength */
	insertlength = hit3->genomicstart - genomicend + querylength5 + querylength3;
	debug9(printf(" (%u)",insertlength));
	if (insertlength >= expected_pairlength_low && insertlength <= expected_pairlength_high) {
	  nbingo++;
	  bingoi5 = i;
	  debug9(printf("*"));
	}

	if ((nmismatches = end_amb_nmismatches[i]) < best_nmismatches) {
	  best_nmismatches = nmismatches;
	  besti5 = i;
	  nbest = 1;
	} else if (nmismatches == best_nmismatches) {
	  nbest++;
	}
      }
    }

    if (nbingo == 1) {
      debug9(printf("\nnbingo is 1\n"));
      *amb_resolve_5 = bingoi5;
      *amb_status_inside = AMB_RESOLVED_BYLENGTH;

    } else if (nbest == 0) {
      debug9(printf("\nnbest is zero: amb_length = %d...%d",
		    end_amb_length(hit5),start_amb_length(hit3)));
      *unresolved_amb_length = end_amb_length_5;
      *amb_status_inside = AMB_UNRESOLVED_TOOCLOSE;

    } else if (nbest == 1) {
      debug9(printf("\nnbest is 1, with nmismatches %d\n",best_nmismatches));
      *amb_resolve_5 = besti5;
      *amb_status_inside = AMB_RESOLVED_BYMATCHES;

    } else {
      *amb_resolve_5 = -1;
      *amb_status_inside = AMB_UNRESOLVED_MULTIPLE;
    }
    debug9(printf("\n"));

  } else if (substring3 != NULL && Substring_ambiguous_p(substring3) == true) {
    debug9(printf("Resolve plus: Got ambiguous at 3':"));
    start_ambcoords = Substring_ambcoords(substring3);
    start_amb_nmismatches = Substring_amb_nmismatches(substring3);
    start_amb_length_3 = start_amb_length(hit3);
    
    nbingo = nbest = 0;
    best_nmismatches = querylength3;
    for (j = 0; j < Substring_nambcoords(substring3); j++) {
      genomicstart = start_ambcoords[j] - start_amb_length_3;
      debug9(printf(" %u",(Chrpos_T) (genomicstart - hit3->chroffset)));
      if (hit5->genomicend < genomicstart /*allow overlap*/+ querylength5) {
	/* Look for valid insertlength */
	insertlength = genomicstart - hit5->genomicend + querylength5 + querylength3;
	debug9(printf(" (%u)",insertlength));
	if (insertlength >= expected_pairlength_low && insertlength <= expected_pairlength_high) {
	  nbingo++;
	  bingoi3 = j;
	  debug9(printf("*"));
	}

	if ((nmismatches = start_amb_nmismatches[j]) < best_nmismatches) {
	  best_nmismatches = nmismatches;
	  besti3 = j;
	  nbest = 1;
	} else if (nmismatches == best_nmismatches) {
	  nbest++;
	}
      }
    }

    if (nbingo == 1) {
      debug9(printf("\nnbingo is 1\n"));
      *amb_resolve_3 = bingoi3;
      *amb_status_inside = AMB_RESOLVED_BYLENGTH;

    } else if (nbest == 0) {
      debug9(printf("\nnbest is zero: amb_length = %d...%d",
		    end_amb_length(hit5),start_amb_length(hit3)));
      *unresolved_amb_length = start_amb_length_3;
      *amb_status_inside = AMB_UNRESOLVED_TOOCLOSE;

    } else if (nbest == 1) {
      debug9(printf("\nnbest is 1, with nmismatches %d\n",best_nmismatches));
      *amb_resolve_3 = besti3;
      *amb_status_inside = AMB_RESOLVED_BYMATCHES;

    } else {
      *amb_resolve_3 = -1;
      *amb_status_inside = AMB_UNRESOLVED_MULTIPLE;
    }
     
    debug9(printf("\n"));
  }

  return;
}


static void
resolve_inside_ambiguous_splice_minus (int *unresolved_amb_length, int *amb_resolve_5, int *amb_resolve_3,
				       int *amb_status_inside, T hit5, T hit3, Univcoord_T *splicesites,
				       Compress_T query5_compress_rev, Compress_T query3_compress_rev,
				       int localsplicing_penalty, int querylength5, int querylength3,
				       int genestrand) {
  int insertlength;
  Univcoord_T genomicstart, genomicend;
  int nbingo, bingoi5, bingoi3;
  int nbest, besti5, besti3, i, j;
  int best_nmismatches, nmismatches;
  bool new5p = false, new3p = false;

  Substring_T substring5, substring3;
#ifdef LARGE_GENOMES
  Uint8list_T ambcoords;
#else
  Uintlist_T ambcoords;
#endif
  Univcoord_T *end_ambcoords, *start_ambcoords;
  int *end_amb_nmismatches, *start_amb_nmismatches;
  int end_amb_length_5, start_amb_length_3;

  *unresolved_amb_length = 0;

  debug9(printf("resolve minus: hit5 %s and hit3 %s\n",
		hittype_string(hit5->hittype),hittype_string(hit3->hittype)));
  if (hit5->hittype == GMAP) {
    substring5 = (Substring_T) NULL;
  } else {
    substring5 = (Substring_T) List_head(hit5->substrings_Nto1);
    debug9(printf("hit5 ambiguous_p %d\n",Substring_ambiguous_p(substring5)));
  }
  if (hit3->hittype == GMAP) {
    substring3 = (Substring_T) NULL;
  } else {
    substring3 = (Substring_T) List_head(hit3->substrings_1toN);
    debug9(printf("hit3 ambiguous_p %d\n",Substring_ambiguous_p(substring3)));
  }

  if (substring5 != NULL && Substring_ambiguous_p(substring5) == true &&
      substring3 != NULL && Substring_ambiguous_p(substring3) == true) {
    debug9(printf("Resolve minus: Got ambiguous at 5' and ambiguous at 3':"));
    end_ambcoords = Substring_ambcoords(substring5);
    end_amb_nmismatches = Substring_amb_nmismatches(substring5);
    start_ambcoords = Substring_ambcoords(substring3);
    start_amb_nmismatches = Substring_amb_nmismatches(substring3);
    end_amb_length_5 = end_amb_length(hit5);
    start_amb_length_3 = start_amb_length(hit3);

    nbingo = nbest = 0;
    best_nmismatches = querylength5 + querylength3;
    for (i = 0; i < Substring_nambcoords(substring5); i++) {
      genomicend = end_ambcoords[i] - end_amb_length_5;
      for (j = 0; j < Substring_nambcoords(substring3); j++) {
	genomicstart = start_ambcoords[j] + start_amb_length_3;
	debug9(printf(" %u,%u",(Chrpos_T) (genomicend - hit5->chroffset),(Chrpos_T) (genomicstart - hit3->chroffset)));
	if (genomicstart < genomicend) {
	  /* Look for valid insertlength */
	  insertlength = genomicend - genomicstart + querylength5 + querylength3;
	  debug9(printf(" (%u)",insertlength));
	  if (insertlength >= expected_pairlength_low && insertlength <= expected_pairlength_high) {
	    nbingo++;
	    bingoi5 = i;
	    bingoi3 = j;
	    debug9(printf("*"));
	  }

	  if ((nmismatches = end_amb_nmismatches[i] + start_amb_nmismatches[j]) < best_nmismatches) {
	    best_nmismatches = nmismatches;
	    besti5 = i;
	    besti3 = j;
	    nbest = 1;
	  } else if (nmismatches == best_nmismatches) {
	    nbest++;
	  }
	}
      }
    }

    if (nbingo == 1) {
      debug9(printf("\nnbingo is 1\n"));
      *amb_resolve_5 = bingoi5;
      *amb_resolve_3 = bingoi3;
      *amb_status_inside = AMB_RESOLVED_BYLENGTH;

    } else if (nbest == 0) {
      debug9(printf("\nnbest is zero: amb_length = %d...%d",
		    end_amb_length(hit5),start_amb_length(hit3)));
      *unresolved_amb_length = end_amb_length_5 + start_amb_length_3;
      *amb_status_inside = AMB_UNRESOLVED_TOOCLOSE;

    } else if (nbest == 1) {
      debug9(printf("\nnbest is 1, with nmismatches %d\n",best_nmismatches));
      *amb_resolve_5 = besti5;
      *amb_resolve_3 = besti3;
      *amb_status_inside = AMB_RESOLVED_BYMATCHES;

    } else {
      *amb_resolve_5 = -1;	/* Signifies cannot resolve */
      *amb_resolve_3 = -1;
      *amb_status_inside = AMB_UNRESOLVED_MULTIPLE;

    }
    debug9(printf("\n"));

  } else if (substring5 != NULL && Substring_ambiguous_p(substring5) == true) {
    debug9(printf("Resolve minus: Got ambiguous at 5':"));
    end_ambcoords = Substring_ambcoords(substring5);
    end_amb_nmismatches = Substring_amb_nmismatches(substring5);
    end_amb_length_5 = end_amb_length(hit5);

    nbingo = nbest = 0;
    best_nmismatches = querylength5;
    for (i = 0; i < Substring_nambcoords(substring5); i++) {
      genomicend = end_ambcoords[i] - end_amb_length_5;
      debug9(printf(" %u",(Chrpos_T) (genomicend - hit5->chroffset)));
      if (hit3->genomicstart < genomicend /*allow overlap*/+ querylength3) {
	/* Look for valid insertlength */
	insertlength = genomicend - hit3->genomicstart + querylength5 + querylength3;
	debug9(printf(" (%u)",insertlength));
	if (insertlength >= expected_pairlength_low && insertlength <= expected_pairlength_high) {
	  nbingo++;
	  bingoi5 = i;
	  debug9(printf("*"));
	}

	if ((nmismatches = end_amb_nmismatches[i]) < best_nmismatches) {
	  best_nmismatches = nmismatches;
	  besti5 = i;
	  nbest = 1;
	} else if (nmismatches == best_nmismatches) {
	  nbest++;
	}
      }
    }

    if (nbingo == 1) {
      debug9(printf("\nnbingo is 1\n"));
      *amb_resolve_5 = bingoi5;
      *amb_status_inside = AMB_RESOLVED_BYLENGTH;

    } else if (nbest == 0) {
      debug9(printf("\nnbest is zero: amb_length = %d...%d",
		    end_amb_length(hit5),start_amb_length(hit3)));
      *unresolved_amb_length = end_amb_length_5;
      *amb_status_inside = AMB_UNRESOLVED_TOOCLOSE;

    } else if (nbest == 1) {
      debug9(printf("\nnbest is 1, with nmismatches %d\n",best_nmismatches));
      *amb_resolve_5 = besti5;
      *amb_status_inside = AMB_RESOLVED_BYMATCHES;

    } else {
      *amb_resolve_5 = -1;
      *amb_status_inside = AMB_UNRESOLVED_MULTIPLE;
    }
    debug9(printf("\n"));

  } else if (substring3 != NULL && Substring_ambiguous_p(substring3) == true) {
    debug9(printf("Resolve minus: Got ambiguous at 3':"));
    start_ambcoords = Substring_ambcoords(substring3);
    start_amb_nmismatches = Substring_amb_nmismatches(substring3);
    start_amb_length_3 = start_amb_length(hit3);

    nbingo = nbest = 0;
    best_nmismatches = querylength3;
    for (j = 0; j < Substring_nambcoords(substring3); j++) {
      genomicstart = start_ambcoords[j] + start_amb_length_3;
      debug9(printf(" %u",(Chrpos_T) (genomicstart - hit3->chroffset)));
      if (genomicstart < hit5->genomicend /*allow overlap*/+ querylength5) {
	/* Look for valid insertlength */
	insertlength = hit5->genomicend - genomicstart + querylength5 + querylength3;
	debug9(printf(" (%u)",insertlength));
	if (insertlength >= expected_pairlength_low && insertlength <= expected_pairlength_high) {
          nbingo++;
          bingoi3 = j;
          debug9(printf("*"));
	}

	if ((nmismatches = start_amb_nmismatches[j]) < best_nmismatches) {
	  best_nmismatches = nmismatches;
	  besti3 = j;
	  nbest = 1;
	} else if (nmismatches == best_nmismatches) {
	  nbest++;
	}
      }
    }

    if (nbingo == 1) {
      debug9(printf("\nnbingo is 1\n"));
      *amb_resolve_3 = bingoi3;
      *amb_status_inside = AMB_RESOLVED_BYLENGTH;

    } else if (nbest == 0) {
      debug9(printf("\nnbest is zero: amb_length = %d...%d",
		    end_amb_length(hit5),start_amb_length(hit3)));
      *unresolved_amb_length = start_amb_length_3;
      *amb_status_inside = AMB_UNRESOLVED_TOOCLOSE;

    } else if (nbest == 1) {
      debug9(printf("\nnbest is 1, with nmismatches %d\n",best_nmismatches));
      *amb_resolve_3 = besti3;
      *amb_status_inside = AMB_RESOLVED_BYMATCHES;

    } else {
      *amb_resolve_3 = -1;
      *amb_status_inside = AMB_UNRESOLVED_MULTIPLE;
    }
    debug9(printf("\n"));
  }


  return;
}



static void
alias_circular (T hit) {
  Chrpos_T chrlength = hit->chrlength;
  List_T p;
  Substring_T substring;

  assert(hit->alias == -1);
  if (hit->hittype == GMAP) {
    Pair_alias_circular(hit->pairarray,hit->npairs,chrlength);

  } else {
    for (p = hit->substrings_1toN; p != NULL; p = List_next(p)) {
      substring = (Substring_T) List_head(p);
      Substring_alias_circular(substring);
    }
  }

  /* Doesn't fix hitpair->low and hitpair->high */
  hit->genomicstart += chrlength;
  hit->genomicend += chrlength;
  hit->low += chrlength;
  hit->high += chrlength;

  hit->alias = +1;

  return;
}


static int
compute_insertlength (Stage3pair_T this) {
  T hit5, hit3;
  Chrpos_T chrstart, chrend, chrpos;
  int querypos;
  int querylength5, querylength3;


  hit5 = this->hit5;
  hit3 = this->hit3;
  querylength5 = hit5->querylength;
  querylength3 = hit3->querylength;

  if (hit5->hittype == GMAP && hit3->hittype == GMAP) {
    debug10(printf("Got hit5 and hit3 both of type GMAP\n"));

    /* Do not try to resolve ambiguity on inside of concordant ends */
    if (hit5->plusp == true && hit3->plusp == true) {
      return (hit3->genomicstart - hit5->genomicend) + querylength5 + querylength3;
    } else if (hit5->plusp == false && hit3->plusp == false) {
      return (hit5->genomicend - hit3->genomicstart) + querylength5 + querylength3;
    } else {
      return pair_insert_length_unpaired(hit5,hit3);
    }

  } else if (hit5->hittype == GMAP) {
    debug10(printf("Got hit5 of type GMAP\n"));
    if (hit5->plusp == true && hit3->plusp == true) {
      /* Have 5-start..end and 3-start..end */
      debug10(printf("plus: comparing hit5->genomicend %u <= hit3->genomicstart %u\n",
		     hit5->genomicend - hit5->chroffset,hit3->genomicstart - hit3->chroffset));

      if (hit5->genomicend <= hit3->genomicstart) {
	/* No overlap */
	return (hit3->genomicstart - hit5->genomicend) + querylength5 + querylength3;
      } else if ((chrpos = overlap3_gmap_plus(&querypos,&chrstart,&chrend,/*hit*/hit3,/*gmap*/hit5)) > 0U) {
	return /* end3 */ chrend - /* start5 */ (chrpos - querypos);
      } else {
	/* Still no overlap */
	return (hit3->genomicstart - hit5->genomicend) + querylength5 + querylength3;
      }

    } else if (hit5->plusp == false && hit3->plusp == false) {
      /* Have 3-end..start and 5-end..start */
      debug10(printf("minus: comparing hit3->genomicstart %u <= hit5->genomicend %u\n",
		     hit3->genomicstart - hit3->chroffset,hit5->genomicend - hit5->chroffset));

      if (hit3->genomicstart <= hit5->genomicend) {
	return (hit5->genomicend - hit3->genomicstart) + querylength5 + querylength3;
      } else if ((chrpos = overlap3_gmap_minus(&querypos,&chrstart,&chrend,/*hit*/hit3,/*gmap*/hit5)) > 0U) {
	return /* start5 */ (chrpos + querypos) - /* end3 */ chrend + 1;
      } else {
	/* Still no overlap */
	return (hit5->genomicend - hit3->genomicstart) + querylength5 + querylength3;
      }
    } else {
      return pair_insert_length_unpaired(hit5,hit3);
    }

  } else if (hit3->hittype == GMAP) {
    debug10(printf("Got hit3 of type GMAP\n"));
    if (hit5->plusp == true && hit3->plusp == true) {
      /* Have 5-start..end and 3-start..end */
      debug10(printf("plus: comparing hit5->genomicend %u <= hit3->genomicstart %u\n",
		     hit5->genomicend - hit5->chroffset,hit3->genomicstart - hit3->chroffset));

      if (hit5->genomicend <= hit3->genomicstart) {
	/* No overlap */
	return (hit3->genomicstart - hit5->genomicend) + querylength5 + querylength3;
      } else if ((chrpos = overlap5_gmap_plus(&querypos,&chrstart,&chrend,/*hit*/hit5,/*gmap*/hit3)) > 0U) {
	return /* end3 */ (chrpos - querypos + querylength3) - /* start5 */ chrstart;
      } else {
	/* Still no overlap */
	return (hit3->genomicstart - hit5->genomicend) + querylength5 + querylength3;
      }

    } else if (hit5->plusp == false && hit3->plusp == false) {
      /* Have 3-end..start and 5-end..start */
      debug10(printf("minus: comparing hit3->genomicstart %u <= hit5->genomicend %u\n",
		     hit3->genomicstart - hit3->chroffset,hit5->genomicend - hit5->chroffset));
      if (hit3->genomicstart <= hit5->genomicend) {
	/* No overlap */
	return (hit5->genomicend - hit3->genomicstart) + querylength5 + querylength3;
      } else if ((chrpos = overlap5_gmap_minus(&querypos,&chrstart,&chrend,/*hit*/hit5,/*gmap*/hit3)) > 0U) {
	return /* start5 */ chrstart - /* end3 */ (chrpos + querypos - querylength3) - 1;
      } else {
	/* Still no overlap */
	return (hit5->genomicend - hit3->genomicstart) + querylength5 + querylength3;
      }
    } else {
      return pair_insert_length_unpaired(hit5,hit3);
    }

  } else if (hit5->plusp == true && hit3->plusp == false) {
    /* Have 5-start..end and 3-end..start */
    /*   or 3-end..start and 5-start..end */

    if (hit5->genomicend < hit3->genomicend) {
      return (hit3->genomicend - hit5->genomicend) + querylength5 + querylength3;
    } else if (hit3->genomicstart < hit5->genomicstart) {
      return (hit5->genomicstart - hit3->genomicstart) + querylength5 + querylength3;
    } else {
      return pair_insert_length_unpaired(hit5,hit3);
    }

  } else if (hit5->plusp == false && hit3->plusp == true) {
    /* Have 5-end..start and 3-start..end */
    /*   or 3-start..end and 5-end..start */

    if (hit5->genomicstart < hit3->genomicstart) {
      return (hit3->genomicstart - hit5->genomicstart) + querylength5 + querylength3;
    } else if (hit3->genomicend < hit5->genomicend) {
      return (hit5->genomicend - hit3->genomicend) + querylength5 + querylength3;
    } else {
      return pair_insert_length_unpaired(hit5,hit3);
    }

  } else if (hit5->plusp == true) {
    /* Concordant directions on same chromosome (plus) */
    debug10(printf("Concordant on plus strand\n"));
    /* Have 5-start..end and 3-start..end */
    if (hit5->genomicend < hit3->genomicstart) {
      /* No overlap */
      return (hit3->genomicstart - hit5->genomicend) + querylength5 + querylength3;
    } else {
      return pair_insert_length(hit5,hit3);
    }


  } else {
    /* Concordant directions on same chromosome (minus) */
    debug10(printf("Concordant on minus strand\n"));
    /* Have 3-end..start and 5-end..start */
    if (hit3->genomicstart < hit5->genomicend) {
      /* No overlap */
      return (hit5->genomicend - hit3->genomicstart) + querylength5 + querylength3;
    } else {
      return pair_insert_length(hit5,hit3);
    }
  }
}



Stage3pair_T
Stage3pair_new (T hit5, T hit3,	Univcoord_T *splicesites,
		Compress_T query5_compress_fwd, Compress_T query5_compress_rev,
		Compress_T query3_compress_fwd, Compress_T query3_compress_rev,
		int genestrand,	Pairtype_T pairtype, int localsplicing_penalty,
		bool private5p, bool private3p, bool expect_concordant_p) {
  Stage3pair_T new;
  Stage3end_T copy;
  Substring_T substring1, substringN;
  Chrpos_T chrstart, chrend, chrpos;
  int querypos;
  int unresolved_amb_length = 0;
  int found_score = 0;
  bool overreach5p, overreach3p;

  int querylength5 = hit5->querylength;
  int querylength3 = hit3->querylength;

  debug10(printf("\nStage3pair_new called with pairtype %s and chrnum %d, %d (effective %d, %d)\n",
		 Pairtype_string(pairtype),hit5->chrnum,hit3->chrnum,hit5->effective_chrnum,hit3->effective_chrnum));

#if 0
  if (hit5->hittype == TERMINAL && hit5->trim_left + hit5->trim_right >= reject_trimlength) {
    debug10(printf("5' rejected by trimlength %d + %d\n",hit5->trim_left,hit5->trim_right));
    if (private5p == true) {
      Stage3end_free(&hit5);
    }
    if (private3p == true) {
      Stage3end_free(&hit3);
    }
    debug5(printf("Rejecting terminal as NULL because hit5 trim %d+%d > reject_trimlength %d\n",hit5->trim_left,hit5->trim_right,reject_trimlength));
    return (Stage3pair_T) NULL;

  } else if (hit3->hittype == TERMINAL && hit3->trim_left + hit3->trim_right >= reject_trimlength) {
    debug10(printf("3' rejected by trimlength %d + %d\n",hit3->trim_left,hit3->trim_right));
    if (private5p == true) {
      Stage3end_free(&hit5);
    }
    if (private3p == true) {
      Stage3end_free(&hit3);
    }
    debug5(printf("Rejecting terminal as NULL because hit3 trim %d+%d > reject_trimlength %d\n",hit3->trim_left,hit3->trim_right,reject_trimlength));
    return (Stage3pair_T) NULL;
  } else {
#endif
    new = (Stage3pair_T) MALLOC_OUT(sizeof(*new));
#if 0
  }
#endif

  if (pairtype == PAIRED_UNSPECIFIED || pairtype == UNSPECIFIED) {
    /* Can get here from running GMAP improvement on a paired result */
    pairtype = Stage3_determine_pairtype(hit5,hit3);
    debug10(printf("  Changing pairtype to %s\n",Pairtype_string(pairtype)));
    if (pairtype == CONCORDANT) {
      expect_concordant_p = true;
    }
  }
  new->pairtype = pairtype;
  new->genestrand = genestrand;
  new->amb_resolve_5 = -1;
  new->amb_resolve_3 = -1;
  new->amb_status_inside = AMB_NOT_AMBIGUOUS;


#if 0
  new->mapq_loglik = hit5->mapq_loglik + hit3->mapq_loglik;
  new->mapq_score = 0;
  new->absmq_score = 0;
#endif

  if (hit5->hittype == GMAP && hit3->hittype == GMAP) {
    debug10(printf("Got hit5 and hit3 both of type GMAP\n"));

    /* Do not try to resolve ambiguity on inside of concordant ends */
    if (hit5->plusp == true && hit3->plusp == true) {
      new->dir = +1;
      new->insertlength = (hit3->genomicstart - hit5->genomicend) + querylength5 + querylength3;
      new->insertlength_expected_sign = insertlength_expected(new->insertlength);
      debug10(printf("plus, no overlap: insert length %d = start3 %u - end5 %u + %d + %d\n",
		     new->insertlength,hit3->genomicstart - hit3->chroffset,
		     hit5->genomicend - hit5->chroffset,querylength5,querylength3));
    } else if (hit5->plusp == false && hit3->plusp == false) {
      new->dir = -1;
      new->insertlength = (hit5->genomicend - hit3->genomicstart) + querylength5 + querylength3;
      new->insertlength_expected_sign = insertlength_expected(new->insertlength);
      debug10(printf("minus, no overlap: insert length %d = end5 %u - start3 %u + %d + %d\n",
		     new->insertlength,hit5->genomicend - hit5->chroffset,
		     hit3->genomicstart - hit3->chroffset,querylength5,querylength3));
    } else {
      new->dir = 0;
      new->insertlength = pair_insert_length_unpaired(hit5,hit3);
      new->insertlength_expected_sign = false;
    }

  } else if (hit5->hittype == GMAP) {
    debug10(printf("Got hit5 of type GMAP\n"));
    if (hit5->plusp == true && hit3->plusp == true) {
      new->dir = +1;

      if (expect_concordant_p == true) {
	/* Try to resolve ambiguity on inside of concordant ends */
	resolve_inside_ambiguous_splice_plus(&unresolved_amb_length,&new->amb_resolve_5,&new->amb_resolve_3,
					     &new->amb_status_inside,hit5,hit3,
					     splicesites,query5_compress_fwd,query3_compress_fwd,
					     localsplicing_penalty,querylength5,querylength3,genestrand);
      }

      /* Have 5-start..end and 3-start..end */
      debug10(printf("plus: comparing hit5->genomicend %u <= hit3->genomicstart %u\n",
		     hit5->genomicend - hit5->chroffset,hit3->genomicstart - hit3->chroffset));

      if (hit5->genomicend <= hit3->genomicstart) {
	/* No overlap */
	new->insertlength = (hit3->genomicstart - hit5->genomicend) + querylength5 + querylength3;
	new->insertlength_expected_sign = insertlength_expected(new->insertlength);
	debug10(printf("plus, no overlap: insert length %d = start3 %u - end5 %u + %d + %d\n",
		       new->insertlength,hit3->genomicstart - hit3->chroffset,
		       hit5->genomicend - hit5->chroffset,querylength5,querylength3));
      } else if ((chrpos = overlap3_gmap_plus(&querypos,&chrstart,&chrend,/*hit*/hit3,/*gmap*/hit5)) > 0U) {
	new->insertlength = /* end3 */ chrend - /* start5 */ (chrpos - querypos);
	new->insertlength_expected_sign = insertlength_expected(new->insertlength);
	debug10(printf("plus, overlap: insert length %d = end3 %u - start5 (%u - %d)\n",
		       new->insertlength,chrend,chrpos,querypos));
      } else {
	/* Still no overlap */
	new->insertlength = (hit3->genomicstart - hit5->genomicend) + querylength5 + querylength3;
	new->insertlength_expected_sign = insertlength_expected(new->insertlength);

	if (new->insertlength <= 0) {
	  /* Overreach */
	  debug5(printf("  Returning NULL because of overreach\n"));
	  if (private5p == true) {
	    Stage3end_free(&hit5);
	  }
	  if (private3p == true) {
	    Stage3end_free(&hit3);
	  }
	  FREE_OUT(new);
	  return (Stage3pair_T) NULL;
	}
      }

    } else if (hit5->plusp == false && hit3->plusp == false) {
      new->dir = -1;

      if (expect_concordant_p == true) {
	/* Try to resolve ambiguity on inside of concordant ends */
	resolve_inside_ambiguous_splice_minus(&unresolved_amb_length,&new->amb_resolve_5,&new->amb_resolve_3,
					      &new->amb_status_inside,hit5,hit3,
					      splicesites,query5_compress_rev,query3_compress_rev,
					      localsplicing_penalty,querylength5,querylength3,genestrand);
      }

      /* Have 3-end..start and 5-end..start */
      debug10(printf("minus: comparing hit3->genomicstart %u <= hit5->genomicend %u\n",
		     hit3->genomicstart - hit3->chroffset,hit5->genomicend - hit5->chroffset));

      if (hit3->genomicstart <= hit5->genomicend) {
	/* No overlap */
	new->insertlength = (hit5->genomicend - hit3->genomicstart) + querylength5 + querylength3;
	new->insertlength_expected_sign = insertlength_expected(new->insertlength);
	debug10(printf("minus, no overlap: insert length %d = end5 %u - start3 %u + %d + %d\n",
		       new->insertlength,hit5->genomicend - hit5->chroffset,
		       hit3->genomicstart - hit3->chroffset,querylength5,querylength3));
      } else if ((chrpos = overlap3_gmap_minus(&querypos,&chrstart,&chrend,/*hit*/hit3,/*gmap*/hit5)) > 0U) {
	new->insertlength = /* start5 */ (chrpos + querypos) - /* end3 */ chrend + 1;
	new->insertlength_expected_sign = insertlength_expected(new->insertlength);
	debug10(printf("minus, overlap: insert length %d = start5 (%u + %d) - end3 %u + 1\n",
		       new->insertlength,chrpos,querypos,chrend));
      } else {
	/* Still no overlap */
	new->insertlength = (hit5->genomicend - hit3->genomicstart) + querylength5 + querylength3;
	new->insertlength_expected_sign = insertlength_expected(new->insertlength);

	if (new->insertlength <= 0) {
	  /* Overreach */
	  debug5(printf("  Returning NULL because of overreach\n"));
	  if (private5p == true) {
	    Stage3end_free(&hit5);
	  }
	  if (private3p == true) {
	    Stage3end_free(&hit3);
	  }
	  FREE_OUT(new);
	  return (Stage3pair_T) NULL;
	}
      }

    } else {
      new->dir = 0;
      new->insertlength = pair_insert_length_unpaired(hit5,hit3); /* was 0 */
      new->insertlength_expected_sign = false;
    }

  } else if (hit3->hittype == GMAP) {
    debug10(printf("Got hit3 of type GMAP\n"));
    if (hit5->plusp == true && hit3->plusp == true) {
      new->dir = +1;

      if (expect_concordant_p == true) {
	/* Try to resolve ambiguity on inside of concordant ends */
	resolve_inside_ambiguous_splice_plus(&unresolved_amb_length,&new->amb_resolve_5,&new->amb_resolve_3,
					     &new->amb_status_inside,hit5,hit3,
					     splicesites,query5_compress_fwd,query3_compress_fwd,
					     localsplicing_penalty,querylength5,querylength3,genestrand);
      }

      /* Have 5-start..end and 3-start..end */
      debug10(printf("plus: comparing hit5->genomicend %u <= hit3->genomicstart %u\n",
		     hit5->genomicend - hit5->chroffset,hit3->genomicstart - hit3->chroffset));

      if (hit5->genomicend <= hit3->genomicstart) {
	/* No overlap */
	new->insertlength = (hit3->genomicstart - hit5->genomicend) + querylength5 + querylength3;
	new->insertlength_expected_sign = insertlength_expected(new->insertlength);
	debug10(printf("plus, no overlap: insert length %d = start3 %u - end5 %u + %d + %d\n",
		       new->insertlength,hit3->genomicstart - hit3->chroffset,
		       hit5->genomicend - hit5->chroffset,querylength5,querylength3));
      } else if ((chrpos = overlap5_gmap_plus(&querypos,&chrstart,&chrend,/*hit*/hit5,/*gmap*/hit3)) > 0U) {
	new->insertlength = /* end3 */ (chrpos - querypos + querylength3) - /* start5 */ chrstart;
	new->insertlength_expected_sign = insertlength_expected(new->insertlength);
	debug10(printf("plus, overlap: insert length %d = end3 (%u - %d + %d) - start5 %u\n",
		       new->insertlength,chrpos,querypos,querylength3,chrstart));
      } else {
	/* Still no overlap */
	new->insertlength = (hit3->genomicstart - hit5->genomicend) + querylength5 + querylength3;
	new->insertlength_expected_sign = insertlength_expected(new->insertlength);

	if (new->insertlength <= 0) {
	  /* Overreach */
	  debug5(printf("  Returning NULL because of overreach\n"));
	  if (private5p == true) {
	    Stage3end_free(&hit5);
	  }
	  if (private3p == true) {
	    Stage3end_free(&hit3);
	  }
	  FREE_OUT(new);
	  return (Stage3pair_T) NULL;
	}
      }

    } else if (hit5->plusp == false && hit3->plusp == false) {
      new->dir = -1;

      if (expect_concordant_p == true) {
	/* Try to resolve ambiguity on inside of concordant ends */
	resolve_inside_ambiguous_splice_minus(&unresolved_amb_length,&new->amb_resolve_5,&new->amb_resolve_3,
					      &new->amb_status_inside,hit5,hit3,
					      splicesites,query5_compress_rev,query3_compress_rev,
					      localsplicing_penalty,querylength5,querylength3,genestrand);
      }

      /* Have 3-end..start and 5-end..start */
      debug10(printf("minus: comparing hit3->genomicstart %u <= hit5->genomicend %u\n",
		     hit3->genomicstart - hit3->chroffset,hit5->genomicend - hit5->chroffset));
      if (hit3->genomicstart <= hit5->genomicend) {
	/* No overlap */
	new->insertlength = (hit5->genomicend - hit3->genomicstart) + querylength5 + querylength3;
	new->insertlength_expected_sign = insertlength_expected(new->insertlength);
	debug10(printf("minus, no overlap: insert length %d = end5 %u - start3 %u + %d + %d\n",
		       new->insertlength,hit5->genomicend - hit5->chroffset,
		       hit3->genomicstart - hit3->chroffset,querylength5,querylength3));
      } else if ((chrpos = overlap5_gmap_minus(&querypos,&chrstart,&chrend,/*hit*/hit5,/*gmap*/hit3)) > 0U) {
	new->insertlength = /* start5 */ chrstart - /* end3 */ (chrpos + querypos - querylength3) - 1;
	new->insertlength_expected_sign = insertlength_expected(new->insertlength);
	debug10(printf("minus, overlap: insert length %d = start5 %u - end3 (%u + %d - %d) - 1\n",
		       new->insertlength,chrstart,chrpos,querypos,querylength3));
      } else {
	/* Still no overlap */
	new->insertlength = (hit5->genomicend - hit3->genomicstart) + querylength5 + querylength3;
	new->insertlength_expected_sign = insertlength_expected(new->insertlength);

	if (new->insertlength <= 0) {
	  /* Overreach */
	  debug5(printf("  Returning NULL because of overreach\n"));
	  if (private5p == true) {
	    Stage3end_free(&hit5);
	  }
	  if (private3p == true) {
	    Stage3end_free(&hit3);
	  }
	  FREE_OUT(new);
	  return (Stage3pair_T) NULL;
	}
      }
    } else {
      new->dir = 0;
      new->insertlength = pair_insert_length_unpaired(hit5,hit3); /* was 0 */
      new->insertlength_expected_sign = false;
    }

  } else if (hit5->plusp == true && hit3->plusp == false) {
    new->dir = 0;
    
    /* Have 5-start..end and 3-end..start */
    /*   or 3-end..start and 5-start..end */

    if (hit5->genomicend < hit3->genomicend) {
      new->insertlength = (hit3->genomicend - hit5->genomicend) + querylength5 + querylength3;
      new->insertlength_expected_sign = insertlength_expected(new->insertlength);
    } else if (hit3->genomicstart < hit5->genomicstart) {
      new->insertlength = (hit5->genomicstart - hit3->genomicstart) + querylength5 + querylength3;
      new->insertlength_expected_sign = insertlength_expected(new->insertlength);
    } else {
      new->insertlength = pair_insert_length_unpaired(hit5,hit3); /* was 0 */
      new->insertlength_expected_sign = false;
    }

  } else if (hit5->plusp == false && hit3->plusp == true) {
    new->dir = 0;
    
    /* Have 5-end..start and 3-start..end */
    /*   or 3-start..end and 5-end..start */

    if (hit5->genomicstart < hit3->genomicstart) {
      new->insertlength = (hit3->genomicstart - hit5->genomicstart) + querylength5 + querylength3;
      new->insertlength_expected_sign = insertlength_expected(new->insertlength);
    } else if (hit3->genomicend < hit5->genomicend) {
      new->insertlength = (hit5->genomicend - hit3->genomicend) + querylength5 + querylength3;
      new->insertlength_expected_sign = insertlength_expected(new->insertlength);
    } else {
      new->insertlength = pair_insert_length_unpaired(hit5,hit3); /* was 0 */
      new->insertlength_expected_sign = false;
    }

  } else if (hit5->plusp == true) {
    /* Concordant directions on same chromosome (plus) */
    debug10(printf("Concordant on plus strand\n"));
    new->dir = +1;

    if (expect_concordant_p == true) {
      overreach5p = overreach3p = false;
      if (hit5->hittype == SPLICE) {

	substringN = (Substring_T) List_head(hit5->substrings_Nto1);
	if (Substring_alignstart(substringN) > hit3->genomicend) {
	  substring1 = (Substring_T) List_head(hit5->substrings_1toN);
	  if (Substring_alignend(substring1) < hit3->genomicstart) {
	    overreach5p = true;
	  }
	}
      }
      if (hit3->hittype == SPLICE) {
	substring1 = (Substring_T) List_head(hit3->substrings_1toN);
	if (Substring_alignend(substring1) < hit5->genomicstart) {
	  substringN = (Substring_T) List_head(hit3->substrings_Nto1);
	  if (Substring_alignstart(substringN) > hit5->genomicend) {
	    overreach3p = true;
	  }
	}
      }

      if (overreach5p == true || overreach3p == true) {
	/* Either overreach */
	debug5(printf("  Returning NULL because of dual overreach\n"));
	if (private5p == true) {
	  Stage3end_free(&hit5);
	}
	if (private3p == true) {
	  Stage3end_free(&hit3);
	}
	FREE_OUT(new);
	return (Stage3pair_T) NULL;

#if 0
      } else if (overreach5p == true) {
	/* Overreach of hit5 */
	debug9(printf("Overreach of hit5 of type SPLICE.  Removing substring2\n"));
	if (hit5->sensedir == SENSE_FORWARD) {
	  copy = Stage3end_new_splice(&found_score,/*nmismatches_donor*/Substring_nmismatches_whole(hit5->substring1),
				      /*nmismatches_acceptor*/0,/*donor*/hit5->substring1,/*acceptor*/NULL,/*distance*/0U,
				      /*shortdistancep*/true,localsplicing_penalty,hit5->querylength,/*amb_length*/0,/*amb_prob*/0.0,
				      /*ambcoords_donor*/NULL,/*ambcoords_acceptor*/NULL,
				      /*amb_nmismatches_donor*/NULL,/*amb_nmismatches_acceptor*/NULL,
				      /*amb_probs_donor*/NULL,/*amb_probs_acceptor*/NULL,
				      /*copy_donor_p*/true,/*copy_acceptor_p*/false,/*first_read_p*/true,
				      /*sensedir*/hit5->sensedir,hit5->sarrayp);
	} else if (hit5->sensedir == SENSE_ANTI) {
	  copy = Stage3end_new_splice(&found_score,/*nmismatches_donor*/0,
				      /*nmismatches_acceptor*/Substring_nmismatches_whole(hit5->substring1),/*donor*/NULL,
				      /*acceptor*/hit5->substring1,/*distance*/0U,
				      /*shortdistancep*/true,localsplicing_penalty,hit5->querylength,/*amb_length*/0,/*amb_prob*/0.0,
				      /*ambcoords_donor*/NULL,/*ambcoords_acceptor*/NULL,
				      /*amb_nmismatches_donor*/NULL,/*amb_nmismatches_acceptor*/NULL,
				      /*amb_probs_donor*/NULL,/*amb_probs_acceptor*/NULL,
				      /*copy_donor_p*/false,/*copy_acceptor_p*/true,/*first_read_p*/true,
				      /*sensedir*/hit5->sensedir,hit5->sarrayp);
	} else {
	  abort();
	}
	if (private5p == true) {
	  Stage3end_free(&hit5);
	} else {
	  private5p = true;
	}
	hit5 = copy;

      } else if (overreach3p == true) {
	/* Overreach of hit3 */
	debug9(printf("Overreach of hit3 of type SPLICE.  Removing substring1\n"));
	if (hit3->sensedir == SENSE_FORWARD) {
	  copy = Stage3end_new_splice(&found_score,/*nmismatches_donor*/0,
				      /*nmismatches_acceptor*/Substring_nmismatches_whole(hit3->substring2),/*donor*/NULL,
				      /*acceptor*/hit3->substring2,/*distance*/0U,
				      /*shortdistancep*/true,localsplicing_penalty,hit3->querylength,/*amb_length*/0,/*amb_prob*/0.0,
				      /*ambcoords_donor*/NULL,/*ambcoords_acceptor*/NULL,
				      /*amb_nmismatches_donor*/NULL,/*amb_nmismatches_acceptor*/NULL,
				      /*amb_probs_donor*/NULL,/*amb_probs_acceptor*/NULL,
				      /*copy_donor_p*/false,/*copy_acceptor_p*/true,/*first_read_p*/true,
				      /*sensedir*/hit3->sensedir,hit3->sarrayp);
	} else if (hit3->sensedir == SENSE_ANTI) {
	  copy = Stage3end_new_splice(&found_score,/*nmismatches_donor*/Substring_nmismatches_whole(hit3->substring2),
				      /*nmismatches_acceptor*/0,/*donor*/hit3->substring2,/*acceptor*/NULL,/*distance*/0U,
				      /*shortdistancep*/true,localsplicing_penalty,hit3->querylength,/*amb_length*/0,/*amb_prob*/0.0,
				      /*ambcoords_donor*/NULL,/*ambcoords_acceptor*/NULL,
				      /*amb_nmismatches_donor*/NULL,/*amb_nmismatches_acceptor*/NULL,
				      /*amb_probs_donor*/NULL,/*amb_probs_acceptor*/NULL,
				      /*copy_donor_p*/true,/*copy_acceptor_p*/false,/*first_read_p*/true,
				      /*sensedir*/hit3->sensedir,hit3->sarrayp);
	} else {
	  abort();
	}
	if (private3p == true) {
	  Stage3end_free(&hit3);
	} else {
	  private3p = true;
	}
	hit3 = copy;
#endif
      }

      /* Try to resolve ambiguity on inside of concordant ends */
      resolve_inside_ambiguous_splice_plus(&unresolved_amb_length,&new->amb_resolve_5,&new->amb_resolve_3,
					   &new->amb_status_inside,hit5,hit3,
					   splicesites,query5_compress_fwd,query3_compress_fwd,
					   localsplicing_penalty,querylength5,querylength3,genestrand);
    }

    /* Have 5-start..end and 3-start..end */
    if (hit5->genomicend < hit3->genomicstart) {
      /* No overlap */
      new->insertlength = (hit3->genomicstart - hit5->genomicend) + querylength5 + querylength3;
      new->insertlength_expected_sign = insertlength_expected(new->insertlength);
      debug10(printf("plus, no overlap: insert length %d = start3 %u - end5 %u + %d + %d\n",
		     new->insertlength,hit3->genomicstart - hit3->chroffset,
		     hit5->genomicend - hit5->chroffset,querylength5,querylength3));
#if 0
    } else if (hit5->genomicend > hit3->genomicend + SUBSUMPTION_SLOP) {
      /* hit5 subsumes hit3 */
      debug10(printf("plus, subsumption %u > %u\n",
		     hit5->genomicend - hit5->chroffset,hit3->genomicend - hit3->chroffset));
      new->insertlength = 0;
      new->insertlength_expected_sign = false;
#endif
    } else {
      new->insertlength = pair_insert_length(hit5,hit3);
      new->insertlength_expected_sign = insertlength_expected(new->insertlength);
    }


  } else {
    /* Concordant directions on same chromosome (minus) */
    debug10(printf("Concordant on minus strand\n"));
    new->dir = -1;

    if (expect_concordant_p == true) {
      overreach5p = overreach3p = false;
      if (hit5->hittype == SPLICE) {
	substringN = (Substring_T) List_head(hit5->substrings_Nto1);
	if (Substring_alignstart(substringN) < hit3->genomicend) {
	  substring1 = (Substring_T) List_head(hit5->substrings_1toN);
	  if (Substring_alignend(substring1) > hit3->genomicstart) {
	    overreach5p = true;
	  }
	}
      }
      if (hit3->hittype == SPLICE) {
	substring1 = (Substring_T) List_head(hit3->substrings_1toN);
	if (Substring_alignend(substring1) > hit5->genomicstart) {
	  substringN = (Substring_T) List_head(hit3->substrings_Nto1);
	  if (Substring_alignstart(substringN) < hit5->genomicend) {
	    overreach3p = true;
	  }
	}
      }

      if (overreach5p == true || overreach3p == true) {
	/* Either overreach */
	debug5(printf("  Returning NULL because of dual overreach\n"));
	if (private5p == true) {
	  Stage3end_free(&hit5);
	}
	if (private3p == true) {
	  Stage3end_free(&hit3);
	}
	FREE_OUT(new);
	return (Stage3pair_T) NULL;

#if 0
      } else if (overreach5p == true) {
	/* Overreach of hit5 */
	debug9(printf("Overreach of hit5 of type SPLICE.  Removing substring2\n"));
	if (hit5->sensedir == SENSE_FORWARD) {
	  copy = Stage3end_new_splice(&found_score,/*nmismatches_donor*/Substring_nmismatches_whole(hit5->substring1),
				      /*nmismatches_acceptor*/0,/*donor*/hit5->substring1,/*acceptor*/NULL,/*distance*/0U,
				      /*shortdistancep*/true,localsplicing_penalty,hit5->querylength,/*amb_length*/0,/*amb_prob*/0.0,
				      /*ambcoords_donor*/NULL,/*ambcoords_acceptor*/NULL,
				      /*amb_nmismatches_donor*/NULL,/*amb_nmismatches_acceptor*/NULL,
				      /*amb_probs_donor*/NULL,/*amb_probs_acceptor*/NULL,
				      /*copy_donor_p*/true,/*copy_acceptor_p*/false,/*first_read_p*/false,
				      /*sensedir*/hit5->sensedir,hit5->sarrayp);
	} else if (hit5->sensedir == SENSE_ANTI) {
	  copy = Stage3end_new_splice(&found_score,/*nmismatches_donor*/0,
				      /*nmismatches_acceptor*/Substring_nmismatches_whole(hit5->substring1),/*donor*/NULL,
				      /*acceptor*/hit5->substring1,/*distance*/0U,
				      /*shortdistancep*/true,localsplicing_penalty,hit5->querylength,/*amb_length*/0,/*amb_prob*/0.0,
				      /*ambcoords_donor*/NULL,/*ambcoords_acceptor*/NULL,
				      /*amb_nmismatches_donor*/NULL,/*amb_nmismatches_acceptor*/NULL,
				      /*amb_probs_donor*/NULL,/*amb_probs_acceptor*/NULL,
				      /*copy_donor_p*/false,/*copy_acceptor_p*/true,/*first_read_p*/false,
				      /*sensedir*/hit5->sensedir,hit5->sarrayp);
	} else {
	  abort();
	}
	if (private5p == true) {
	  Stage3end_free(&hit5);
	} else {
	  private5p = true;
	}
	hit5 = copy;

      } else if (overreach3p == true) {
	/* Overreach of hit3 */
	debug9(printf("Overreach of hit3 of type SPLICE.  Removing substring1\n"));
	if (hit3->sensedir == SENSE_FORWARD) {
	  copy = Stage3end_new_splice(&found_score,/*nmismatches_donor*/0,
				      /*nmismatches_acceptor*/Substring_nmismatches_whole(hit3->substring2),/*donor*/NULL,
				      /*acceptor*/hit3->substring2,/*distance*/0U,
				      /*shortdistancep*/true,localsplicing_penalty,hit3->querylength,/*amb_length*/0,/*amb_prob*/0.0,
				      /*ambcoords_donor*/NULL,/*ambcoords_acceptor*/NULL,
				      /*amb_nmismatches_donor*/NULL,/*amb_nmismatches_acceptor*/NULL,
				      /*amb_probs_donor*/NULL,/*amb_probs_acceptor*/NULL,
				      /*copy_donor_p*/false,/*copy_acceptor_p*/true,/*first_read_p*/false,
				      /*sensedir*/hit3->sensedir,hit3->sarrayp);
	} else if (hit3->sensedir == SENSE_ANTI) {
	  copy = Stage3end_new_splice(&found_score,/*nmismatches_donor*/Substring_nmismatches_whole(hit3->substring2),
				      /*nmismatches_acceptor*/0,/*donor*/hit3->substring2,/*acceptor*/NULL,/*distance*/0U,
				      /*shortdistancep*/true,localsplicing_penalty,hit3->querylength,/*amb_length*/0,/*amb_prob*/0.0,
				      /*ambcoords_donor*/NULL,/*ambcoords_acceptor*/NULL,
				      /*amb_nmismatches_donor*/NULL,/*amb_nmismatches_acceptor*/NULL,
				      /*amb_probs_donor*/NULL,/*amb_probs_acceptor*/NULL,
				      /*copy_donor_p*/true,/*copy_acceptor_p*/false,/*first_read_p*/false,
				      /*sensedir*/hit3->sensedir,hit3->sarrayp);
	} else {
	  abort();
	}
	if (private3p == true) {
	  Stage3end_free(&hit3);
	} else {
	  private3p = true;
	}
	hit3 = copy;
#endif
      }

      /* Try to resolve ambiguity on inside of concordant ends */
      resolve_inside_ambiguous_splice_minus(&unresolved_amb_length,&new->amb_resolve_5,&new->amb_resolve_3,
					    &new->amb_status_inside,hit5,hit3,
					    splicesites,query5_compress_rev,query3_compress_rev,
					    localsplicing_penalty,querylength5,querylength3,genestrand);
    }

    /* Have 3-end..start and 5-end..start */
    if (hit3->genomicstart < hit5->genomicend) {
      /* No overlap */
      new->insertlength = (hit5->genomicend - hit3->genomicstart) + querylength5 + querylength3;
      new->insertlength_expected_sign = insertlength_expected(new->insertlength);
      debug10(printf("minus, no overlap: insert length %d = end5 %u - start3 %u + %d + %d\n",
		     new->insertlength,hit5->genomicend - hit5->chroffset,
		     hit3->genomicstart - hit3->chroffset,querylength5,querylength3));
#if 0
    } else if (hit3->genomicstart > hit5->genomicstart + SUBSUMPTION_SLOP) {
      /* hit3 subsumes hit5 */
      debug10(printf("minus, subsumption %u > %u\n",
		     hit3->genomicstart - hit3->chroffset,hit5->genomicstart - hit5->chroffset));
      new->insertlength = 0;
      new->insertlength_expected_sign = false;
#endif
    } else {
      new->insertlength = pair_insert_length(hit5,hit3);
      new->insertlength_expected_sign = insertlength_expected(new->insertlength);
    }
  }

  debug5(printf("\nGot insertlength of %d\n",new->insertlength));
  /* Was new->insertlength <= 0, but this eliminates legitimate overlaps */
  /* Was new->insertlength < -pairmax, but this allows overreach */
  if (new->insertlength <= 0) {
    /* Not concordant */
#ifdef USE_BINGO
    new->absdifflength_bingo_p = false;
#endif
#ifdef USE_ABSDIFFLENGTH
    new->absdifflength = -1U;
#endif

    if (expect_concordant_p == true) {
      debug5(printf("  Returning NULL, because not concordant\n"));
      if (private5p == true) {
	Stage3end_free(&hit5);
      }
      if (private3p == true) {
	Stage3end_free(&hit3);
      }
      FREE_OUT(new);
      return (Stage3pair_T) NULL;
    }

  } else if (new->insertlength > pairmax && expect_concordant_p == true) {
    debug5(printf("  Returning NULL because insertlength %d > pairmax %d\n",new->insertlength,pairmax));
    if (private5p == true) {
      Stage3end_free(&hit5);
    }
    if (private3p == true) {
      Stage3end_free(&hit3);
    }
    FREE_OUT(new);
    return (Stage3pair_T) NULL;

  } else {
#ifdef USE_ABSDIFFLENGTH
    if (new->insertlength < expected_pairlength) {
      new->absdifflength = expected_pairlength - new->insertlength;
    } else {
      new->absdifflength = new->insertlength - expected_pairlength;
    }
#endif
#ifdef USE_BINGO
    if (new->absdifflength <= pairlength_deviation) {
      new->absdifflength_bingo_p = true;
    } else {
      new->absdifflength_bingo_p = false;
    }
#endif
  }

  if (SENSE_CONSISTENT_P(hit5->sensedir,hit3->sensedir)) {
    debug0(printf("senses are consistent\n"));
    new->sense_consistent_p = true;
  } else {
    debug0(printf("senses are inconsistent\n"));
    new->sense_consistent_p = false;
  }

  /* Do not alter score, so the alignmnent terminates at the known splice site  */
  new->score = hit5->score + hit3->score /* + unresolved_amb_length */;

  new->nmatches = hit5->nmatches + hit3->nmatches - unresolved_amb_length;
  new->nmatches_posttrim = hit5->nmatches_posttrim + hit3->nmatches_posttrim;
  /* new->overlap_known_gene_p = false; -- initialized later when resolving multimappers */
  new->tally = -1L;

  new->low = (hit5->low < hit3->low) ? hit5->low : hit3->low;
  new->high = (hit5->high > hit3->high) ? hit5->high : hit3->high;

#if 0
  if (new->low > new->high) {
    fprintf(stderr,"new->low %u > new->high %u, hit5->chrnum %d\n",
	    new->low - new->chroffset,new->high - new->chroffset,hit5->chrnum);
    abort();
  }
#endif

  if (hit5->chrnum == 0 || hit3->chrnum == 0) {
    new->outerlength = querylength5 + querylength3;
  } else {
    new->outerlength = new->high - new->low;
  }

  new->hit5 = hit5;
  new->hit3 = hit3;

  new->private5p = private5p;
  new->private3p = private3p;


  if (expect_concordant_p == true) {
    hit5->paired_usedp = true;
    hit3->paired_usedp = true;
  }

  new->nsplices = hit5->nsplices + hit3->nsplices;

  debug0(printf("Created new pair %p from %p and %p with private %d, %d\n",new,hit5,hit3,private5p,private3p));
  debug0(printf("  hittypes %s and %s\n",hittype_string(hit5->hittype),hittype_string(hit3->hittype)));
  debug0(printf("  sensedirs %d and %d\n",hit5->sensedir,hit3->sensedir));
  debug0(printf("  chrpos %u..%u and %u..%u\n",
		hit5->genomicstart - hit5->chroffset,hit5->genomicend - hit5->chroffset,
		hit3->genomicstart - hit3->chroffset,hit3->genomicend - hit3->chroffset));

  if (hit5->circularpos < 0 && hit3->circularpos < 0) {
    new->circularp = false;
  } else {
    new->circularp = true;
  }

  /* Fixing insertlength for circular pairs */
  if (new->insertlength > hit5->chrlength) {
    new->insertlength -= hit5->chrlength;
  }

  if (hit5->alias == +1) {
    debug0(printf("Unaliasing 5' end\n"));
    if (private5p == false) {
      new->hit5 = Stage3end_copy(hit5);
      new->private5p = true;
    }
    unalias_circular(new->hit5);
  }

  if (hit3->alias == +1) {
    debug0(printf("Unaliasing 3' end\n"));
    if (private3p == false) {
      new->hit3 = Stage3end_copy(hit3);
      new->private3p = true;
    }
    unalias_circular(new->hit3);
  }

  /* assert((int) new->insertlength >= 0); */
  return new;
}


void
Stage3pair_privatize (Stage3pair_T *array, int npairs) {
  Stage3pair_T hitpair;
  int i;

  debug0(printf("Call to Stage3pair_privatize on %d hitpairs\n",npairs));

  for (i = 0; i < npairs; i++) {
    hitpair = array[i];
    debug0(printf("  Pair with hitpairs %p (private %d), %p (private %d)\n",
		  hitpair->hit5,hitpair->private5p,hitpair->hit3,hitpair->private3p));
  
    if (hitpair->private5p == false) {
      hitpair->hit5 = Stage3end_copy(hitpair->hit5);
      hitpair->private5p = true;
    }
    if (hitpair->private3p == false) {
      hitpair->hit3 = Stage3end_copy(hitpair->hit3);
      hitpair->private3p = true;
    }
  }

  return;
}


/* Used for eliminating exact duplicates.  Also sorts secondarily by hittype. */
static int
hitpair_sort_cmp (const void *a, const void *b) {
  Stage3pair_T x = * (Stage3pair_T *) a;
  Stage3pair_T y = * (Stage3pair_T *) b;
  
  Univcoord_T x_hit5_high, x_hit5_low, y_hit5_high, y_hit5_low;
  Univcoord_T x_hit3_high, x_hit3_low, y_hit3_high, y_hit3_low;
  Univcoord_T x_low, x_high, y_low, y_high;
  
  debug8(printf("  Comparing (%s, %s-%s): %p, %u..%u|%u..%u (dir = %d), alias %d|%d, nmatches: %d (%d posttrim), amb_lengths %d and %d, sensedirs %d-%d\n",
		Pairtype_string(x->pairtype),hittype_string(x->hit5->hittype),
		hittype_string(x->hit3->hittype),x,
		x->hit5->low - x->hit5->chroffset,x->hit5->high - x->hit5->chroffset,
		x->hit3->low - x->hit3->chroffset,x->hit3->high - x->hit3->chroffset,
		x->dir,x->hit5->alias,x->hit3->alias,x->nmatches,x->nmatches_posttrim,
		start_amb_length(x->hit5) + end_amb_length(x->hit5),start_amb_length(x->hit3) + end_amb_length(x->hit3),
		x->hit5->sensedir,x->hit3->sensedir));

  debug8(printf("       with (%s, %s-%s): %p, %u..%u|%u..%u (dir = %d), alias %d|%d, nmatches: %d (%d posttrim), amb_lengths %d and %d, sensedirs %d-%d\n",
		Pairtype_string(y->pairtype),hittype_string(y->hit5->hittype),
		hittype_string(y->hit3->hittype),y,
		y->hit5->low - y->hit5->chroffset,y->hit5->high - y->hit5->chroffset,
		y->hit3->low - y->hit3->chroffset,y->hit3->high - y->hit3->chroffset,
		y->dir,y->hit5->alias,y->hit3->alias,y->nmatches,y->nmatches_posttrim,
		start_amb_length(y->hit5) + end_amb_length(y->hit5),start_amb_length(y->hit3) + end_amb_length(y->hit3),
		y->hit5->sensedir,y->hit3->sensedir));


  x_hit5_low = normalize_coord(x->hit5->low,x->hit5->alias,x->hit5->chrlength);
  x_hit5_high = normalize_coord(x->hit5->high,x->hit5->alias,x->hit5->chrlength);

  x_hit3_low = normalize_coord(x->hit3->low,x->hit3->alias,x->hit3->chrlength);
  x_hit3_high = normalize_coord(x->hit3->high,x->hit3->alias,x->hit3->chrlength);

  x_low = (x_hit5_low < x_hit3_low) ? x_hit5_low : x_hit3_low;
  x_high = (x_hit5_high > x_hit3_high) ? x_hit5_high : x_hit3_high;


  y_hit5_low = normalize_coord(y->hit5->low,y->hit5->alias,y->hit5->chrlength);
  y_hit5_high = normalize_coord(y->hit5->high,y->hit5->alias,y->hit5->chrlength);

  y_hit3_low = normalize_coord(y->hit3->low,y->hit3->alias,y->hit3->chrlength);
  y_hit3_high = normalize_coord(y->hit3->high,y->hit3->alias,y->hit3->chrlength);

  y_low = (y_hit5_low < y_hit3_low) ? y_hit5_low : y_hit3_low;
  y_high = (y_hit5_high > y_hit3_high) ? y_hit5_high : y_hit3_high;


  if (x->dir != 0 && y->dir == 0) {
    return -1;
  } else if (x->dir == 0 && y->dir != 0) {
    return +1;
  } else if (x->dir > 0 && y->dir < 0) {
    return -1;
  } else if (x->dir < 0 && y->dir > 0) {
    return +1;

#if 0
  } else if (x->high < y->low) {
    return -1;
  } else if (y->high < x->low) {
    return +1;

  } else if (x->hit5->high < y->hit5->low) {
    return -1;
  } else if (y->hit5->high < x->hit5->low) {
    return +1;

  } else if (x->hit3->high < y->hit3->low) {
    return -1;
  } else if (y->hit3->high < x->hit3->low) {
    return +1;
#else
    /* low to high pattern needed for finding overlaps */
  } else if (x_low < y_low) {
    debug8(printf("Returning -1 for low\n"));
    return -1;
  } else if (y_low < x_low) {
    debug8(printf("Returning +1 for low\n"));
    return +1;

  } else if (x_high > y_high) {
    debug8(printf("Returning -1 for high\n"));
    return -1;
  } else if (y_high > x_high) {
    debug8(printf("Returning +1 for high\n"));
    return +1;

    /* Need to check inside ends to avoid declaring unequal hitpairs equal */
  } else if (x_hit5_low < y_hit5_low) {
    return -1;
  } else if (y_hit5_low < x_hit5_low) {
    return +1;

  } else if (x_hit5_high < y_hit5_high) {
    return -1;
  } else if (y_hit5_high < x_hit5_high) {
    return +1;

  } else if (x_hit3_low < y_hit3_low) {
    return -1;
  } else if (y_hit3_low < x_hit3_low) {
    return +1;

  } else if (x_hit3_high < y_hit3_high) {
    return -1;
  } else if (y_hit3_high < x_hit3_high) {
    return +1;
#endif


#if 1
    /* Rank terminals last, so terminals cannot win */
  } else if (x->hit5->hittype != TERMINAL && x->hit3->hittype != TERMINAL && 
	     (y->hit5->hittype == TERMINAL || y->hit3->hittype == TERMINAL)) {
    return -1;
  } else if ((x->hit5->hittype == TERMINAL || x->hit3->hittype == TERMINAL) &&
	     y->hit5->hittype != TERMINAL && y->hit3->hittype != TERMINAL) {
    return +1;
#endif

  } else if (x->score < y->score) {
    return -1;
  } else if (y->score < x->score) {
    return +1;
  } else if (x->nmatches > y->nmatches) {
    return -1;
  } else if (y->nmatches > x->nmatches) {
    return +1;
#if 0
  } else if (x->nmatches_posttrim > y->nmatches_posttrim) {
    return -1;
  } else if (y->nmatches_posttrim > x->nmatches_posttrim) {
    return +1;
#endif

#if 0
  } else if (x->nsplices < y->nsplices) {
    return -1;
  } else if (y->nsplices < x->nsplices) {
    return +1;
#endif

  } else if (x->amb_status_inside < y->amb_status_inside) {
    return -1;
  } else if (y->amb_status_inside < x->amb_status_inside) {
    return +1;

  } else if (x->hit5->hittype < y->hit5->hittype) {
    return -1;
  } else if (y->hit5->hittype < x->hit5->hittype) {
    return +1;
  } else if (x->hit3->hittype < y->hit3->hittype) {
    return -1;
  } else if (y->hit3->hittype < x->hit3->hittype) {
    return +1;

#if 0
  } else if ((x->amb_resolve_5 != -1 && x->amb_resolve_3 != -1) &&
	     (y->amb_resolve_5 == -1 || y->amb_resolve_3 == -1)) {
    /* x is resolved, y is ambiguous.  x wins */
    return -1;
  } else if ((y->amb_resolve_5 != -1 && y->amb_resolve_3 != -1) &&
	     (x->amb_resolve_5 == -1 || x->amb_resolve_3 == -1)) {
    /* y is resolved, x is ambiguous.  y wins */
    return +1;
#endif

#if 0
  } else if (x->hit5->start_amb_length + x->hit5->end_amb_length +
	     x->hit3->start_amb_length + x->hit3->end_amb_length == 0 &&
	     y->hit5->start_amb_length + y->hit5->end_amb_length +
	     y->hit3->start_amb_length + y->hit3->end_amb_length > 0) {
    /* x is resolved, y is ambiguous.  x wins */
    return -1;
  } else if (y->hit5->start_amb_length + y->hit5->end_amb_length +
	     y->hit3->start_amb_length + y->hit3->end_amb_length == 0 &&
	     x->hit5->start_amb_length + x->hit5->end_amb_length +
	     x->hit3->start_amb_length + x->hit3->end_amb_length > 0) {
    /* y is resolved, x is ambiguous.  y wins */
    return +1;
#endif

  } else if (x->sense_consistent_p == true && y->sense_consistent_p == false) {
    return -1;
  } else if (x->sense_consistent_p == false && y->sense_consistent_p == true) {
    return +1;
#if 0
  } else if (x->indel_low < y->indel_low) {
    return -1;
  } else if (y->indel_low < x->indel_low) {
    return +1;
#endif

  } else if (x->sense_consistent_p == true) {
    if ((x->hit5->sensedir != 0 || x->hit3->sensedir != 0) &&
	(y->hit5->sensedir == 0 && y->hit3->sensedir == 0)) {
      return -1;
    } else if ((y->hit5->sensedir != 0 || y->hit3->sensedir != 0) &&
	       (x->hit5->sensedir == 0 && x->hit3->sensedir == 0)) {
      return +1;
    } else {
      return 0;
    }

  } else {
    return 0;
  }
}


/* Same as hitpair_sort_cmp, except for hittype, nmatches_posttrim, and indel_low */
static int
hitpair_equiv_cmp (Stage3pair_T x, Stage3pair_T y) {
  Univcoord_T x_hit5_high, x_hit5_low, y_hit5_high, y_hit5_low;
  Univcoord_T x_hit3_high, x_hit3_low, y_hit3_high, y_hit3_low;
  Univcoord_T x_low, x_high, y_low, y_high;
  
  x_hit5_low = normalize_coord(x->hit5->low,x->hit5->alias,x->hit5->chrlength);
  x_hit5_high = normalize_coord(x->hit5->high,x->hit5->alias,x->hit5->chrlength);

  x_hit3_low = normalize_coord(x->hit3->low,x->hit3->alias,x->hit3->chrlength);
  x_hit3_high = normalize_coord(x->hit3->high,x->hit3->alias,x->hit3->chrlength);

  x_low = (x_hit5_low < x_hit3_low) ? x_hit5_low : x_hit3_low;
  x_high = (x_hit5_high > x_hit3_high) ? x_hit5_high : x_hit3_high;


  y_hit5_low = normalize_coord(y->hit5->low,y->hit5->alias,y->hit5->chrlength);
  y_hit5_high = normalize_coord(y->hit5->high,y->hit5->alias,y->hit5->chrlength);

  y_hit3_low = normalize_coord(y->hit3->low,y->hit3->alias,y->hit3->chrlength);
  y_hit3_high = normalize_coord(y->hit3->high,y->hit3->alias,y->hit3->chrlength);

  y_low = (y_hit5_low < y_hit3_low) ? y_hit5_low : y_hit3_low;
  y_high = (y_hit5_high > y_hit3_high) ? y_hit5_high : y_hit3_high;


  if (x->dir != 0 && y->dir == 0) {
    return -1;
  } else if (x->dir == 0 && y->dir != 0) {
    return +1;
  } else if (x->dir > 0 && y->dir < 0) {
    return -1;
  } else if (x->dir < 0 && y->dir > 0) {
    return +1;
  } else if (x_low < y_low) {
    return -1;
  } else if (y_low < x_low) {
    return +1;
  } else if (x_high < y_high) {
    return -1;
  } else if (y_high < x_high) {
    return +1;

  } else if (x_hit5_low < y_hit5_low) {
    return -1;
  } else if (y_hit5_low < x_hit5_low) {
    return +1;
  } else if (x_hit5_high < y_hit5_high) {
    return -1;
  } else if (y_hit5_high < x_hit5_high) {
    return +1;

  } else if (x_hit3_low < y_hit3_low) {
    return -1;
  } else if (y_hit3_low < x_hit3_low) {
    return +1;
  } else if (x_hit3_high < y_hit3_high) {
    return -1;
  } else if (y_hit3_high < x_hit3_high) {
    return +1;

#if 0
  } else if (x->hit5->hittype != TERMINAL && x->hit3->hittype != TERMINAL && 
	     (y->hit5->hittype == TERMINAL || y->hit3->hittype == TERMINAL)) {
    return -1;
  } else if ((x->hit5->hittype == TERMINAL || x->hit3->hittype == TERMINAL) &&
	     y->hit5->hittype != TERMINAL && y->hit3->hittype != TERMINAL) {
    return +1;
#endif

#if 0
  } else if (x->score < y->score) {
    return -1;
  } else if (y->score < x->score) {
    return +1;
  } else if (x->nmatches > y->nmatches) {
    return -1;
  } else if (y->nmatches > x->nmatches) {
    return +1;
  } else if (x->nmatches_posttrim > y->nmatches_posttrim) {
    return -1;
  } else if (y->nmatches_posttrim > x->nmatches_posttrim) {
    return +1;
#endif

#if 0
    /* Causes hits to not be recognized as equivalent */
  } else if (x->nsplices < y->nsplices) {
    return -1;
  } else if (y->nsplices < x->nsplices) {
    return +1;
#endif

  } else if (x->amb_status_inside < y->amb_status_inside) {
    return -1;
  } else if (y->amb_status_inside < x->amb_status_inside) {
    return +1;

#if 0
  } else if (x->hit5->start_amb_length + x->hit5->end_amb_length +
	     x->hit3->start_amb_length + x->hit3->end_amb_length > 0 &&
	     y->hit5->start_amb_length + y->hit5->end_amb_length +
	     y->hit3->start_amb_length + y->hit3->end_amb_length == 0) {
    return -1;
  } else if (y->hit5->start_amb_length + y->hit5->end_amb_length +
	     y->hit3->start_amb_length + y->hit3->end_amb_length > 0 &&
	     x->hit5->start_amb_length + x->hit5->end_amb_length +
	     x->hit3->start_amb_length + x->hit3->end_amb_length == 0) {
    return +1;
#endif

  } else if (x->sense_consistent_p == true && y->sense_consistent_p == false) {
    return -1;
  } else if (x->sense_consistent_p == false && y->sense_consistent_p == true) {
    return +1;

#if 0
  } else if (x->indel_low < y->indel_low) {
    return -1;
  } else if (y->indel_low < x->indel_low) {
    return +1;
#endif

#if 0
  } else if (x->sense_consistent_p == true) {
    /* Used for sorting, but not equiv */
    if ((x->hit5->sensedir != 0 || x->hit3->sensedir != 0) &&
	(y->hit5->sensedir == 0 && y->hit3->sensedir == 0)) {
      return -1;
    } else if ((y->hit5->sensedir != 0 || y->hit3->sensedir != 0) &&
	       (x->hit5->sensedir == 0 && x->hit3->sensedir == 0)) {
      return +1;
    } else {
      return 0;
    }
#endif

  } else {
    return 0;
  }
}


static int
hitpair_position_cmp (const void *a, const void *b) {
  Stage3pair_T x = * (Stage3pair_T *) a;
  Stage3pair_T y = * (Stage3pair_T *) b;
  
  if (x->low < y->low) {
    return -1;
  } else if (y->low < x->low) {
    return +1;
  } else if (x->high < y->high) {
    return +1;
  } else if (y->high < x->high) {
    return -1;
  } else {
    return 0;
  }
}


#if 0
static bool
hitpair_equal (Stage3pair_T x, Stage3pair_T y) {
  if (x->dir != y->dir) {
    debug8(printf("=>F "));
    return false;		/* Different strands */
  } else if (x->hit5->low != y->hit5->low) {
    debug8(printf("=>F "));
    return false;
  } else if (x->hit5->high != y->hit5->high) {
    debug8(printf("=>F "));
    return false;
  } else if (x->hit3->low != y->hit3->low) {
    debug8(printf("=>F "));
    return false;
  } else if (x->hit3->high != y->hit3->high) {
    debug8(printf("=>F "));
    return false;
  } else {
    debug8(printf("=>T "));
    return true;
  }
}
#endif


#if 0
static bool
hitpair_overlap (Stage3pair_T x, Stage3pair_T y) {
#if 0
  if ((x->hit5->hittype == SPLICE || x->hit3->hittype == SPLICE) &&
      (y->hit5->hittype == SPLICE || y->hit3->hittype == SPLICE)) {
    /* Special case: pairs involving splices don't overlap */
    return false;
  }
#endif
  if (x->dir != y->dir) {
    return false;		/* Different strands */
  } else if (x->high < y->low) {
    return false;
  } else if (x->low > y->high) {
    return false;
  } else {
    return true;
  }
}
#endif


static bool
hitpair_subsumption (Stage3pair_T x, Stage3pair_T y) {
  if (x->dir != y->dir) {
    return false;		/* Different strands */
  } else if (x->low <= y->low && x->high >= y->high) {
    return true;
  } else if (y->low <= x->low && y->high >= x->high) {
    return true;
    
    /* Test each end of the pair.  Example: 1586..1512 and 1400..1468 should subsume 1586..1512 and 1564..1617 */
  } else if (x->hit5->low <= y->hit5->low && x->hit5->high >= y->hit5->high) {
    return true;
  } else if (y->hit5->low <= x->hit5->low && y->hit5->high >= x->hit5->high) {
    return true;

  } else if (x->hit3->low <= y->hit3->low && x->hit3->high >= y->hit3->high) {
    return true;
  } else if (y->hit3->low <= x->hit3->low && y->hit3->high >= x->hit3->high) {
    return true;

  } else {
    return false;
  }
}


static int
pair_matches_cmp (const void *a, const void *b) {
  Stage3pair_T x = * (Stage3pair_T *) a;
  Stage3pair_T y = * (Stage3pair_T *) b;

  if (x->nmatches > y->nmatches) {
    return -1;
  } else if (x->nmatches < y->nmatches) {
    return +1;
  } else {
    return 0;
  }
}

List_T
Stage3pair_sort_bymatches (List_T hits) {
  List_T sorted = NULL;
  Stage3pair_T *array;
  int n, i;

  
  if ((n = List_length(hits)) == 0) {
    return (List_T) NULL;
  } else {
#ifdef USE_ALLOCA_FOR_HITS
    array = (Stage3pair_T *) MALLOCA(n * sizeof(Stage3pair_T));
    List_fill_array_and_free((void **) array,&hits);
#else
    array = (Stage3pair_T *) List_to_array(hits,NULL);
    List_free(&hits);
#endif

    qsort(array,n,sizeof(Stage3pair_T),pair_matches_cmp);
    for (i = n-1; i >= 0; i--) {
      sorted = List_push(sorted,(void *) array[i]);
    }
#ifdef USE_ALLOCA_FOR_HITS
    FREEA(array);
#else
    FREE(array);
#endif

    return sorted;
  }
}


#if 0
static int
hitpair_distance_cmp (const void *a, const void *b) {
  Stage3pair_T x = * (Stage3pair_T *) a;
  Stage3pair_T y = * (Stage3pair_T *) b;

  debug(printf("Comparing %u with %u\n",x->absdifflength,y->absdifflength));
  if (x->absdifflength < y->absdifflength) {
    return -1;
  } else if (x->absdifflength > y->absdifflength) {
    return +1;
  } else {
    return 0;
  }
}
#endif


#if 0
List_T
Stage3pair_sort_distance (List_T hitpairlist) {
#ifdef DEBUG
  Stage3pair_T hitpair;
  List_T p;
#endif
  List_T sorted = NULL;
  Stage3pair_T *hitpairs;
  int n, i;

  if ((n = List_length(hitpairlist)) == 0) {
    return NULL;
  } else {
#ifdef USE_ALLOCA_FOR_HITS
    hitpairs = (Stage3pair_T *) MALLOCA(n * sizeof(Stage3pair_T));
    List_fill_array_and_free((void **) hitpairs,&hitpairlist);
#else
    hitpairs = (Stage3pair_T *) List_to_array(hitpairlist,NULL);
    List_free(&hitpairlist);
#endif
  }

  qsort(hitpairs,n,sizeof(Stage3pair_T),hitpair_distance_cmp);

  for (i = n-1; i >= 0; i--) {
    sorted = List_push(sorted,hitpairs[i]);
  }
#ifdef USE_ALLOCA_FOR_HITS
  FREEA(hitpairs);
#else
  FREE(hitpairs);
#endif

  debug(
	for (p = sorted, i = 0; p != NULL; p = p->rest, i++) {
	  hitpair = (Stage3pair_T) p->first;
	  printf("  Final %d: %llu-%llu (dir = %d), insert length %u\n",
		 i,(unsigned long long) hitpair->low,(unsigned long long) hitpair->high,
		 hitpair->dir,hitpair->insertlength);
	}
	);

  return sorted;
}
#endif


#if 0
List_T
Stage3pair_remove_duplicates_exact (List_T hitpairlist) {
  List_T unique = NULL;
  Stage3pair_T hitpair, *hitpairs;
  int n, i, j;
  bool *eliminate;

  debug8(printf("Entered Stage3pair_remove_duplicates_exact with %d pairs\n",n));
  if ((n = List_length(hitpairlist)) == 0) {
    return NULL;
  } else {
#ifdef USE_ALLOCA_FOR_HITS
    eliminate = (bool *) CALLOCA(n,sizeof(bool));
    hitpairs = (Stage3pair_T *) MALLOCA(n * sizeof(Stage3pair_T));
    List_fill_array_and_free((void **) hitpairs,&hitpairlist);
#else
    eliminate = (bool *) CALLOC(n,sizeof(bool));
    hitpairs = (Stage3pair_T *) List_to_array(hitpairlist,NULL);
    List_free(&hitpairlist);
#endif
  }

  debug8(printf("Checking for exact duplicates\n"));
  qsort(hitpairs,n,sizeof(Stage3pair_T),hitpair_sort_cmp);

  debug8(
	 for (i = 0; i < n; i++) {
	   hitpair = hitpairs[i];
	   printf("  Initial %d (%s, %s-%s): %p, %u..%u|%u..%u (dir = %d), alias %d|%d, nmatches: %d\n",
		  i,Pairtype_string(hitpair->pairtype),hittype_string(hitpair->hit5->hittype),
		  hittype_string(hitpair->hit3->hittype),hitpair,
		  hitpair->hit5->low - hitpair->hit5->chroffset,hitpair->hit5->high - hitpair->hit5->chroffset,
		  hitpair->hit3->low - hitpair->hit3->chroffset,hitpair->hit3->high - hitpair->hit3->chroffset,
		  hitpair->dir,hitpair->hit5->alias,hitpair->hit3->alias,hitpair->nmatches);
	 }
	 );

  i = 0;
  while (i < n) {
    j = i+1;
    while (j < n && hitpair_equiv_cmp(hitpairs[j],hitpairs[i]) == 0) {
      debug8(printf("  %d is identical to %d => eliminating\n",j,i));
      eliminate[j] = true;
      j++;
    }
    i = j;
  }

  for (i = n-1; i >= 0; i--) {
    hitpair = hitpairs[i];
    if (eliminate[i] == false) {
      unique = List_push(unique,hitpair);
    } else {
      Stage3pair_free(&hitpair);
    }
  }

#ifdef USE_ALLOCA_FOR_HITS
  FREEA(hitpairs);
  FREEA(eliminate);
#else
  FREE(hitpairs);
  FREE(eliminate);
#endif

  debug8(printf("Exited Stage3pair_remove_duplicates_exact with %d pairs\n",List_length(unique)));
  return unique;
}
#endif


static int
hitpair_goodness_cmp (bool *equalp, Stage3pair_T hitpair,
		      Stage3pair_T best_hitpair, bool finalp) {
  double prob1, prob2;

#if 0
  int hitpair_nmatches, best_hitpair_nmatches;
  int max_trim_left, max_trim_right;
  Stage3end_T hit5, besthit5, hit3, besthit3;

  if (hitpair->absdifflength_bingo_p < best_hitpair->absdifflength_bingo_p) {
    /* k is worse */
    debug8(printf(" => loses by absdifflength (bingo)\n"));
    return -1;
  } else if (hitpair->absdifflength_bingo_p > best_hitpair->absdifflength_bingo_p) {
    /* k is better */
    debug8(printf(" => wins by absdifflength (bingo)\n"));
    return +1;
  }
#endif

#ifdef PRE_RESOLVE_MULTIMAPPING
  if (TALLY_RATIO*Stage3pair_tally(hitpair) < Stage3pair_tally(best_hitpair)) {
    /* k is worse */
    debug8(printf(" => loses by tally\n"));
    return -1;
  } else if (Stage3pair_tally(hitpair) > TALLY_RATIO*Stage3pair_tally(best_hitpair)) {
    /* k is better */
    debug8(printf(" => wins by tally\n"));
    return +1;
  }
#endif

  *equalp = false;

#if 1
  if (finalp == true) {
    /* Skip */
  } else if (hitpair->hit5->hittype == TERMINAL || hitpair->hit3->hittype == TERMINAL ||
	     best_hitpair->hit5->hittype == TERMINAL || best_hitpair->hit3->hittype == TERMINAL) {
    /* Do not allow terminal to win or lose in pre-final stages */
    debug8(printf(" => ties by terminal\n"));
    return 0;
  }
#endif

#if 0
  /* hitpair_nmatches = hitpair->nmatches; */
  /* best_hitpair_nmatches = best_hitpair->nmatches; */

  if (hitpair->hit5->hittype == TERMINAL || best_hitpair->hit5->hittype == TERMINAL ||
      hitpair->hit3->hittype == TERMINAL || best_hitpair->hit3->hittype == TERMINAL) {
    /* Skip: Don't use scores if terminal is involved */

#if 0
  } else if (hitpair->hit5->hittype == GMAP && best_hitpair->hit5->hittype == GMAP &&
	     hitpair->hit3->hittype == GMAP && best_hitpair->hit3->hittype == GMAP) {
    /* Dual GMAP alignments: Compare only in trimmed region */
    hit5 = hitpair->hit5;
    besthit5 = best_hitpair->hit5;
    max_trim_left = (hit5->trim_left > besthit5->trim_left) ? hit5->trim_left : besthit5->trim_left;
    max_trim_right = (hit5->trim_right > besthit5->trim_right) ? hit5->trim_right : besthit5->trim_right;
    hitpair_nmatches = Pair_array_nmatches_posttrim(hit5->pairarray,hit5->npairs,
						    /*pos5*/max_trim_left,/*pos3*/hit5->querylength - max_trim_right);
    best_hitpair_nmatches = Pair_array_nmatches_posttrim(besthit5->pairarray,besthit5->npairs,
							 /*pos5*/max_trim_left,/*pos3*/besthit5->querylength - max_trim_right);
    debug8(printf(" gmap/gmap on 5' end with trim %d left, %d right: %d versus %d",
		  max_trim_left,max_trim_right,hitpair_nmatches,best_hitpair_nmatches));

    hit3 = hitpair->hit3;
    besthit3 = best_hitpair->hit3;
    max_trim_left = (hit3->trim_left > besthit3->trim_left) ? hit3->trim_left : besthit3->trim_left;
    max_trim_right = (hit3->trim_right > besthit3->trim_right) ? hit3->trim_right : besthit3->trim_right;
    hitpair_nmatches += Pair_array_nmatches_posttrim(hit3->pairarray,hit3->npairs,
						     /*pos5*/max_trim_left,/*pos3*/hit3->querylength - max_trim_right);
    best_hitpair_nmatches += Pair_array_nmatches_posttrim(besthit3->pairarray,besthit3->npairs,
							  /*pos5*/max_trim_left,/*pos3*/besthit3->querylength - max_trim_right);
    debug8(printf(" gmap/gmap on 3' end with trim %d left, %d right: %d versus %d",
		  max_trim_left,max_trim_right,hitpair_nmatches,best_hitpair_nmatches));

  } else if (hitpair->hit5->hittype == GMAP && best_hitpair->hit5->hittype == GMAP) {
    /* 5' GMAP alignments: Compare only in trimmed region */
    hit5 = hitpair->hit5;
    besthit5 = best_hitpair->hit5;
    max_trim_left = (hit5->trim_left > besthit5->trim_left) ? hit5->trim_left : besthit5->trim_left;
    max_trim_right = (hit5->trim_right > besthit5->trim_right) ? hit5->trim_right : besthit5->trim_right;
    hitpair_nmatches = Pair_array_nmatches_posttrim(hit5->pairarray,hit5->npairs,
						    /*pos5*/max_trim_left,/*pos3*/hit5->querylength - max_trim_right);
    best_hitpair_nmatches = Pair_array_nmatches_posttrim(besthit5->pairarray,besthit5->npairs,
							 /*pos5*/max_trim_left,/*pos3*/besthit5->querylength - max_trim_right);
    debug8(printf(" gmap/gmap on 5' end with trim %d left, %d right: %d versus %d",
		  max_trim_left,max_trim_right,hitpair_nmatches,best_hitpair_nmatches));

    hitpair_nmatches += hitpair->hit3->nmatches;
    best_hitpair_nmatches += best_hitpair->hit3->nmatches;

  } else if (hitpair->hit3->hittype == GMAP && best_hitpair->hit3->hittype == GMAP) {
    /* 3' GMAP alignments: Compare only in trimmed region */
    hit3 = hitpair->hit3;
    besthit3 = best_hitpair->hit3;
    max_trim_left = (hit3->trim_left > besthit3->trim_left) ? hit3->trim_left : besthit3->trim_left;
    max_trim_right = (hit3->trim_right > besthit3->trim_right) ? hit3->trim_right : besthit3->trim_right;
    hitpair_nmatches = Pair_array_nmatches_posttrim(hit3->pairarray,hit3->npairs,
						    /*pos5*/max_trim_left,/*pos3*/hit3->querylength - max_trim_right);
    best_hitpair_nmatches = Pair_array_nmatches_posttrim(besthit3->pairarray,besthit3->npairs,
							 /*pos5*/max_trim_left,/*pos3*/besthit3->querylength - max_trim_right);
    debug8(printf(" gmap/gmap on 3' end with trim %d left, %d right: %d versus %d",
		  max_trim_left,max_trim_right,hitpair_nmatches,best_hitpair_nmatches));

    hitpair_nmatches += hitpair->hit5->nmatches;
    best_hitpair_nmatches += best_hitpair->hit5->nmatches;
#endif

  } else if (hitpair->hit5->nindels == 0 && best_hitpair->hit5->nindels == 0 &&
	     hitpair->hit3->nindels == 0 && best_hitpair->hit3->nindels == 0) {
    /* Skip: Use scores only if indel is involved */
  } else if (hitpair->score > best_hitpair->score) {
    /* k is worse */
    debug8(printf(" => loses by score\n"));
    return -1;
  } else if (hitpair->score < best_hitpair->score) {
    /* k is better */
    debug8(printf(" => wins by score\n"));
    return +1;
  }
#endif


  if (hitpair->nmatches < best_hitpair->nmatches) {
    /* k is worse */
    debug8(printf(" => loses by nmatches\n"));
    return -1;
  } else if (hitpair->nmatches > best_hitpair->nmatches) {
    /* k is better */
    debug8(printf(" => wins by nmatches\n"));
    return +1;

#if 0
  } else if (hitpair->nmatches_posttrim < best_hitpair->nmatches_posttrim) {
    /* k is worse */
    debug8(printf(" => loses by nmatches_posttrim\n"));
    return -1;
  } else if (hitpair->nmatches_posttrim > best_hitpair->nmatches_posttrim) {
    /* k is better */
    debug8(printf(" => wins by nmatches_posttrim\n"));
    return +1;
#endif

#if 0
  } else if (hitpair->nsplices > best_hitpair->nsplices) {
    /* k is worse */
    debug8(printf(" => loses by nsplices: %d > %d in best\n",hitpair->nsplices,best_hitpair->nsplices));
    return -1;
  } else if (hitpair->nsplices < best_hitpair->nsplices) {
    /* k is better */
    debug8(printf(" => wins by nsplices: %d < %d in best\n",hitpair->nsplices,best_hitpair->nsplices));
    return +1;
#endif

  } else if (hitpair->amb_status_inside > best_hitpair->amb_status_inside) {
    /* k is worse */
    debug8(printf(" => loses by amb_status_inside\n"));
    return -1;
  } else if (hitpair->amb_status_inside < best_hitpair->amb_status_inside) {
    /* k is better */
    debug8(printf(" => wins by amb_status_inside\n"));
    return +1;


  } else if (hitpair->hit5->hittype > best_hitpair->hit5->hittype &&
	     hitpair->hit3->hittype >= best_hitpair->hit3->hittype) {
    /* k is worse */
    debug8(printf(" => loses by hittype\n"));
    return -1;

  } else if (hitpair->hit5->hittype >= best_hitpair->hit5->hittype &&
	     hitpair->hit3->hittype > best_hitpair->hit3->hittype) {
    /* k is worse */
    debug8(printf(" => loses by hittype\n"));
    return -1;

  } else if (hitpair->hit5->hittype < best_hitpair->hit5->hittype &&
	     hitpair->hit3->hittype <= best_hitpair->hit3->hittype) {
    /* k is better */
    debug8(printf(" => wins by hittype\n"));
    return +1;

  } else if (hitpair->hit5->hittype <= best_hitpair->hit5->hittype &&
	     hitpair->hit3->hittype < best_hitpair->hit3->hittype) {
    /* k is better */
    debug8(printf(" => wins by hittype\n"));
    return +1;

#if 0
  } else if ((hitpair->amb_resolve_5 == -1 || hitpair->amb_resolve_3 == -1) &&
	     (best_hitpair->amb_resolve_5 != -1 && best_hitpair->amb_resolve_3 != -1)) {
    /* best_hitpair is resolved, hitpair is ambiguous.  best_hitpair wins */
    debug8(printf(" => loses by resolve_inside\n"));
    return -1;

  } else if ((hitpair->amb_resolve_5 != -1 && hitpair->amb_resolve_3 != -1) &&
	     (best_hitpair->amb_resolve_5 == -1 || best_hitpair->amb_resolve_3 == -1)) {
    /* hitpair is resolved, best_hitpair is ambiguous.  hitpair wins */
    debug8(printf(" => wins by resolve_inside: %d, %d, %d, %d\n",
		  hitpair->amb_resolve_5,hitpair->amb_resolve_3,
		  best_hitpair->amb_resolve_5,best_hitpair->amb_resolve_3));
    return +1;
#endif

#if 0
  } else if (hitpair->hit5->start_amb_length + hitpair->hit5->end_amb_length +
	     hitpair->hit3->start_amb_length + hitpair->hit3->end_amb_length > 0 &&
	     best_hitpair->hit5->start_amb_length + best_hitpair->hit5->end_amb_length +
	     best_hitpair->hit3->start_amb_length + best_hitpair->hit3->end_amb_length == 0) {
    debug8(printf(" => loses by ambiguity\n"));
    return -1;

  } else if (hitpair->hit5->start_amb_length + hitpair->hit5->end_amb_length +
	     hitpair->hit3->start_amb_length + hitpair->hit3->end_amb_length == 0 &&
	     best_hitpair->hit5->start_amb_length + best_hitpair->hit5->end_amb_length +
	     best_hitpair->hit3->start_amb_length + best_hitpair->hit3->end_amb_length > 0) {
    debug8(printf(" => wins by ambiguity\n"));
    return +1;
#endif

#if 0
  } else if (hitpair->absdifflength < best_hitpair->absdifflength) {
    /* k is worse */
    debug8(printf(" => loses by absdifflength\n"));
    return -1;
  } else if (hitpair->absdifflength > best_hitpair->absdifflength) {
    /* k is better */
    debug8(printf(" => wins by absdifflength\n"));
    return +1;
#endif

  } else if (finalp == false) {
    debug8(printf("  => indistinguishable\n"));
    return 0;

#ifdef USE_ABSDIFFLENGTH
    /* If insert length is within deviation of expected pairlength, favor it */
  } else if (best_hitpair->absdifflength <= (Chrpos_T) pairlength_deviation &&
	     hitpair->absdifflength > (Chrpos_T) pairlength_deviation) {
    /* k is worse */
    debug8(printf(" => loses by absdifflength within deviation %d\n",pairlength_deviation));
    return -1;
  } else if (hitpair->absdifflength <= (Chrpos_T) pairlength_deviation &&
	     best_hitpair->absdifflength > (Chrpos_T) pairlength_deviation) {
    /* k is better */
    debug8(printf(" => wins by absdifflength within deviation %d\n",pairlength_deviation));
    return +1;
#endif

#if 0
    /* Previously favored longer insert lengths to give more compact
       splices.  However, we now accept splices first that give
       expected pairlength */
  } else if (hitpair->insertlength_expected_sign == -1 && best_hitpair->insertlength_expected_sign == +1) {
    /* k is worse */
    debug8(printf(" => loses by insertlength_expected_sign\n"));
    return -1;
  } else if (hitpair->insertlength_expected_sign == +1 && best_hitpair->insertlength_expected_sign == -1) {
    /* k is better */
    debug8(printf(" => wins by insertlength_expected_sign\n"));
    return +1;
#endif

    /* Next we look at splice probability */
  } else {
    debug8(printf(" => prob"));
    prob1 = Stage3end_prob(hitpair->hit5) + Stage3end_prob(hitpair->hit3);
    prob2 = Stage3end_prob(best_hitpair->hit5) + Stage3end_prob(best_hitpair->hit3);
    if (prob1 + 0.3 < prob2) {
      /* k is worse */
      debug8(printf(" => loses by dual splice prob %f vs %f\n",prob1,prob2));
      return -1;
    } else if (prob1 > prob2 + 0.3) {
      /* k is better */
      debug8(printf(" => wins by dual splice prob %f vs %f\n",prob1,prob2));
      return +1;
    }

    /* Overlapping ends worse than separate ends */
    if (hitpair->insertlength <= hitpair->hit5->querylength + hitpair->hit3->querylength &&
	best_hitpair->insertlength > best_hitpair->hit5->querylength + best_hitpair->hit3->querylength) {
      debug8(printf(" => loses by being overlapping\n"));
      return -1;
    } else if (hitpair->insertlength > hitpair->hit5->querylength + hitpair->hit3->querylength &&
	       best_hitpair->insertlength <= best_hitpair->hit5->querylength + best_hitpair->hit3->querylength) {
      debug8(printf(" => wins by being separate\n"));
      return +1;

      /* Next, favor shorter outerlengths to give more compact splices or closer pairs */
    } else if (hitpair->outerlength > best_hitpair->outerlength + OUTERLENGTH_SLOP) {
      /* k is worse */
      debug8(printf(" => loses by outerlength\n"));
      return -1;
    } else if (hitpair->outerlength + OUTERLENGTH_SLOP < best_hitpair->outerlength) {
      /* k is better */
      debug8(printf(" => wins by outerlength\n"));
      return +1;
      
    } else {
#if 0
      if (hitpair->insertlength_expected_sign >= 0 && best_hitpair->insertlength_expected_sign >= 0) {
	/* Both insert lengths are short, so favor shorter insert length */
	debug8(printf(" => short insertlengths"));
	/* Favor shorter insert lengths */
	if (hitpair->insertlength > best_hitpair->insertlength) {
	  /* k is worse */
	  debug8(printf(" => loses by insertlength\n"));
	  return -1;
	} else if (hitpair->insertlength < best_hitpair->insertlength) {
	  /* k is better */
	  debug8(printf(" => wins by insertlength\n"));
	  return +1;
	}
      }
#endif

      /* Both insert lengths are long, so favor longer insert length to give more compact splices */
      debug8(printf(" => long insertlengths"));
      if (hitpair->insertlength < best_hitpair->insertlength) {
	/* k is worse */
	debug8(printf(" => loses by insertlength\n"));
	return -1;
      } else if (hitpair->insertlength > best_hitpair->insertlength) {
	/* k is better */
	debug8(printf(" => wins by insertlength\n"));
	return +1;
      }

      debug8(printf("  => equal\n"));
      *equalp = true;
      return 0;
    }
  }
}


#if 0
static bool
hitpair_bad_superstretch_p (Stage3pair_T hitpair_k, Stage3pair_T *hitpairs, int k, int j,
			    bool finalp) {
  int a;
  bool equalp;

  for (a = k+1; a <= j; a++) {
    if (hitpair_subsumption(hitpair_k,hitpairs[a]) == true) {
      debug8(printf("Testing %d because stretches over %d",k,a));
      if (hitpair_goodness_cmp(&equalp,hitpairs[a],
			       hitpair_k,finalp) > 0 || equalp == true) {
	debug8(printf(" => eliminating\n"));
	return true;
      }
      debug8(printf("\n"));
    }
  }
  return false;
}
#endif


/* Recursive, list-based approach */
static List_T
pair_remove_bad_superstretches (bool *keep_p, Stage3pair_T superstretch, List_T list, bool finalp) {
  List_T result = NULL, better, equal, p, q, r;
  Stage3pair_T stage3pair, hitpair;
  bool equalp, this_kept_p;
  int cmp;

  *keep_p = true;

  p = list;
  while (p != NULL) {
    stage3pair = (Stage3pair_T) List_head(p);

    q = List_next(p);
    while (q != NULL && hitpair_subsumption(stage3pair,(Stage3pair_T) List_head(q)) == true) {
#ifdef DEBUG8
      printf("  This (%s, %s-%s): %p, %u..%u|%u..%u (dir = %d), nmatches: %d (%d posttrim), insertlength %d, amb_status_inside %d, amb_lengths %d and %d\n",
	     Pairtype_string(stage3pair->pairtype),hittype_string(stage3pair->hit5->hittype),
	     hittype_string(stage3pair->hit3->hittype),stage3pair,
	     stage3pair->hit5->low - stage3pair->hit5->chroffset,stage3pair->hit5->high - stage3pair->hit5->chroffset,
	     stage3pair->hit3->low - stage3pair->hit3->chroffset,stage3pair->hit3->high - stage3pair->hit3->chroffset,
	     stage3pair->dir,stage3pair->nmatches,stage3pair->nmatches_posttrim,
	     stage3pair->insertlength,stage3pair->amb_status_inside,
	     start_amb_length(stage3pair->hit5)+ end_amb_length(stage3pair->hit5),start_amb_length(stage3pair->hit3) + end_amb_length(stage3pair->hit3));
      hitpair = (Stage3pair_T) List_head(q);
      printf("subsumes that (%s, %s-%s): %p, %u..%u|%u..%u (dir = %d), nmatches: %d (%d posttrim), insertlength %d, amb_status_inside %d, amb_lengths %d and %d\n",
	     Pairtype_string(hitpair->pairtype),hittype_string(hitpair->hit5->hittype),
	     hittype_string(hitpair->hit3->hittype),hitpair,
	     hitpair->hit5->low - hitpair->hit5->chroffset,hitpair->hit5->high - hitpair->hit5->chroffset,
	     hitpair->hit3->low - hitpair->hit3->chroffset,hitpair->hit3->high - hitpair->hit3->chroffset,
	     hitpair->dir,hitpair->nmatches,hitpair->nmatches_posttrim,
	     hitpair->insertlength,hitpair->amb_status_inside,
	     start_amb_length(hitpair->hit5) + end_amb_length(hitpair->hit5),start_amb_length(hitpair->hit3) + end_amb_length(hitpair->hit3));
#endif
      q = List_next(q);
    }

    if (q == p) {
      result = List_push(result,(void *) stage3pair);
      if (superstretch != NULL && 
	  (hitpair_goodness_cmp(&equalp,stage3pair,superstretch,finalp) > 0 || equalp == true)) {
	*keep_p = false;
      }
      p = List_next(q);

    } else {
      /* Cluster */
      debug8(printf("Processing cluster of size %d - %d\n",List_length(p),List_length(q)));
      better = equal = (List_T) NULL;
      for (r = List_next(p); r != q; r = List_next(r)) {
	debug8(printf("Calling hitpair_goodness_cmp\n"));
	cmp = hitpair_goodness_cmp(&equalp,(Stage3pair_T) List_head(r),stage3pair,finalp);
	debug8(printf("cmp = %d, equalp = %d\n",cmp,equalp));
	if (cmp > 0) {
	  better = List_push(better,(void *) List_head(r));
	} else if (cmp < 0) {
	  hitpair = (Stage3pair_T) List_head(r);
	  Stage3pair_free(&hitpair);
	} else {
	  equal = List_push(equal,(void *) List_head(r));
	}
      }

      debug8(printf("Found %d better, %d equal\n",List_length(better),List_length(equal)));

      if (better == NULL) {
	/* All children are equal to parent */
	debug8(printf("All children are equivalent, so keeping parent and all (equal) children\n"));
	result = List_push(result,(void *) stage3pair);
	equal = List_reverse(equal); /* Keep original order */
	for (r = equal; r != NULL; r = List_next(r)) {
	  hitpair = (Stage3pair_T) List_head(r);
	  result = List_push(result,(void *) hitpair);
	}
	List_free(&equal);

      } else {
	/* Exists a child better than parent */
	debug8(printf("Exists a child better than parent, so deleting parent and equal and calling recursively among all (better) children\n"));
	Stage3pair_free(&stage3pair);
	for (r = equal; r != NULL; r = List_next(r)) {
	  hitpair = (Stage3pair_T) List_head(r);
	  Stage3pair_free(&hitpair);
	}
	List_free(&equal);

	if (List_length(better) == 1) {
	  hitpair = (Stage3pair_T) List_head(better);
	  result = List_push(result,(void *) hitpair);
	  List_free(&better);

	} else {
	  /* Don't call List_reverse(better) */
	  result = List_append(result,pair_remove_bad_superstretches(&this_kept_p,/*superstretch*/NULL,better,finalp));
#if 0
	  /* Already deleted parent */
	  if (this_kept_p == false) {
	    Stage3pair_free(&stage3pair);
	  } else {
	    /* Compare stage3pair against the current parent */
	    result = List_push(result,(void *) stage3pair);
	    if (superstretch != NULL && 
		(hitpair_goodness_cmp(&equalp,stage3pair,superstretch,finalp) >= 0 || equalp == true)) {
	      *keep_p = false;
	    }
	  }
#endif
	}
      }

      p = q;
    }
  }

  List_free(&list);

  debug8(printf("Returning result of length %d\n",List_length(result)));
  return List_reverse(result);
}


static List_T
pair_remove_overlaps (List_T hitpairlist, bool translocp, bool finalp) {
  List_T unique = NULL;
  Stage3pair_T hitpair, *hitpairs;
  int nkept, n, i, j;
  bool *eliminate;
  bool keep_p;

  n = List_length(hitpairlist);
  debug8(printf("  Entering pair_remove_overlaps with %d pairs: %s\n",
		n,finalp == true ? "FINAL" : "not final"));

  if (n < 2) {
    debug8(printf("  Exiting pair_remove_overlaps with %d < 2 pairs\n",n));
    return hitpairlist;
  } else {
#ifdef USE_ALLOCA_FOR_HITS
    eliminate = (bool *) CALLOCA(n,sizeof(bool));
    hitpairs = (Stage3pair_T *) MALLOCA(n * sizeof(Stage3pair_T));
    List_fill_array_and_free((void **) hitpairs,&hitpairlist);
#else
    eliminate = (bool *) CALLOC(n,sizeof(bool));
    hitpairs = (Stage3pair_T *) List_to_array(hitpairlist,NULL);
    List_free(&hitpairlist);
#endif
  }

  /* Step 1.  Check for exact duplicates */
  debug8(printf("  Step 1.  Checking for exact duplicates\n"));
  qsort(hitpairs,n,sizeof(Stage3pair_T),hitpair_sort_cmp);

  debug8(
	 for (i = 0; i < n; i++) {
	   hitpair = hitpairs[i];
	   printf("  Initial %d (%s, %s-%s): %p, %u..%u|%u..%u (dir = %d), alias %d|%d, nmatches: %d (%d posttrim), amb_lengths %d and %d, sensedirs %d and %d\n",
		  i,Pairtype_string(hitpair->pairtype),hittype_string(hitpair->hit5->hittype),
		  hittype_string(hitpair->hit3->hittype),hitpair,
		  hitpair->hit5->low - hitpair->hit5->chroffset,hitpair->hit5->high - hitpair->hit5->chroffset,
		  hitpair->hit3->low - hitpair->hit3->chroffset,hitpair->hit3->high - hitpair->hit3->chroffset,
		  hitpair->dir,hitpair->hit5->alias,hitpair->hit3->alias,hitpair->nmatches,hitpair->nmatches_posttrim,
		  start_amb_length(hitpair->hit5) + end_amb_length(hitpair->hit5),start_amb_length(hitpair->hit3) + end_amb_length(hitpair->hit3),
		  hitpair->hit5->sensedir,hitpair->hit3->sensedir);
	 }
	 );

  i = 0;
  while (i < n) {
    j = i+1;
    debug8(printf(" %d,%d",i,j));
    while (j < n && hitpair_equiv_cmp(hitpairs[j],hitpairs[i]) == 0) {
      debug8(printf("  %d is identical to %d => eliminating\n",j,i));
      eliminate[j] = true;
      j++;
    }
    i = j;
  }
  debug8(printf("\n"));

  nkept = 0;
  for (i = 0; i < n; i++) {
    if (eliminate[i] == false) {
      nkept++;
    }
  }
  if (nkept == 0) {
    /* All entries eliminated one another, so keep the first one */
    eliminate[0] = false;
    nkept = 1;
  }

  for (i = n - 1; i >= 0; --i) {
    hitpair = hitpairs[i];
    if (eliminate[i] == false) {
      debug8(printf("  Keeping %u..%u|%u..%u, nmatches (trimmed) %d, score %d, (dir = %d)\n",
		    hitpair->hit5->low - hitpair->hit5->chroffset,hitpair->hit5->high - hitpair->hit5->chroffset,
		    hitpair->hit3->low - hitpair->hit3->chroffset,hitpair->hit3->high - hitpair->hit3->chroffset,
		    hitpair->nmatches,hitpair->score,hitpair->dir));
      unique = List_push(unique,(void *) hitpair);
    } else {
      debug8(printf("  Eliminating %u..%u|%u..%u, nmatches (trimmed) %d, score %d, (dir = %d)\n",
		    hitpair->hit5->low - hitpair->hit5->chroffset,hitpair->hit5->high - hitpair->hit5->chroffset,
		    hitpair->hit3->low - hitpair->hit3->chroffset,hitpair->hit3->high - hitpair->hit3->chroffset,
		    hitpair->nmatches,hitpair->score,hitpair->dir));
      Stage3pair_free(&hitpair);
    }
  }

#ifdef USE_ALLOCA_FOR_HITS
  FREEA(hitpairs);
  FREEA(eliminate);
#else
  FREE(hitpairs);
  FREE(eliminate);
#endif

  debug8(printf("  Step 2.  Checking for bad superstretches\n"));
  if (translocp == true) {
    return unique;
  } else {
    return pair_remove_bad_superstretches(&keep_p,/*superstretch*/NULL,unique,finalp);
  }
}


List_T
Stage3pair_remove_overlaps (List_T hitpairlist, bool translocp, bool finalp) {
  List_T unique_separate, unique_overlapping,
    separate = NULL, overlapping = NULL, p;
  Stage3pair_T hitpair;

  List_T indep_overlapping = NULL;
  Stage3pair_T *array_separate, *array_overlapping;
  Stage3pair_T hitpair_overlapping;
  Univcoord_T low, high;
  bool subsumedp, equalp;
  int n_separate, n_overlapping, i, j;


  for (p = hitpairlist; p != NULL; p = List_next(p)) {
    hitpair = (Stage3pair_T) List_head(p);
    if (hitpair->insertlength <= hitpair->hit5->querylength + hitpair->hit3->querylength) {
      overlapping = List_push(overlapping,(void *) hitpair);
    } else {
      separate = List_push(separate,(void *) hitpair);
    }
  }
  List_free(&hitpairlist);

  debug8(printf("Calling Stage3pair_remove_overlaps for separate pair ends\n"));
  unique_separate = pair_remove_overlaps(separate,translocp,finalp);

  debug8(printf("Calling Stage3pair_remove_overlaps for overlapping pair ends\n"));
  unique_overlapping = pair_remove_overlaps(overlapping,translocp,finalp);
  
  if (unique_overlapping == NULL) {
    return unique_separate;
  } else if (unique_separate == NULL) {
    return unique_overlapping;
  } else {
    debug8(printf("Have both overlapping and separate\n"));
    n_overlapping = List_length(unique_overlapping);
#ifdef USE_ALLOCA_FOR_HITS
    array_overlapping = (Stage3pair_T *) MALLOCA(n_overlapping * sizeof(Stage3pair_T));
    List_fill_array_and_free((void **) array_overlapping,&unique_overlapping);
#else
    array_overlapping = (Stage3pair_T *) List_to_array(unique_overlapping,NULL);
    List_free(&unique_overlapping);
#endif

    n_separate = List_length(unique_separate);
#ifdef USE_ALLOCA_FOR_HITS
    array_separate = (Stage3pair_T *) MALLOCA(n_separate * sizeof(Stage3pair_T));
    List_fill_array((void **) array_separate,unique_separate);
#else
    array_separate = (Stage3pair_T *) List_to_array(unique_separate,NULL);
    /* List_free(&unique_separate); -- save for final result */
#endif

    qsort(array_overlapping,n_overlapping,sizeof(Stage3pair_T),hitpair_position_cmp);
    qsort(array_separate,n_separate,sizeof(Stage3pair_T),hitpair_position_cmp);

    i = j = 0;
    for (i = 0; i < n_overlapping; i++) {
      hitpair_overlapping = array_overlapping[i];
      low = hitpair_overlapping->low;
      high = hitpair_overlapping->high;
      while (j >= 0 && array_separate[j]->high >= low) {
	j--;
      }
      j += 1;

      subsumedp = false;
      while (j < n_separate && subsumedp == false && array_separate[j]->low <= high) {
	if (hitpair_goodness_cmp(&equalp,array_separate[j],
				 hitpair_overlapping,finalp) > 0) {
	  debug8(printf("separate pair %d better than overlapping pair %d\n",j,i));
	  subsumedp = hitpair_subsumption(array_separate[j],hitpair_overlapping);
	  debug8(printf("  checking if separate pair %d subsumes overlapping pair %d => %d\n",
			j,i,subsumedp));
	}
	j++;
      }
      j -= 1;

      if (subsumedp == true) {
	Stage3pair_free(&hitpair_overlapping);
      } else {
	indep_overlapping = List_push(indep_overlapping,(void *) hitpair_overlapping);
      }
    }

#ifdef USE_ALLOCA_FOR_HITS
    FREEA(array_separate);
    FREEA(array_overlapping);
#else
    FREE(array_separate);
    FREE(array_overlapping);
#endif

    return List_append(unique_separate,indep_overlapping);
  }
}


List_T
Stage3pair_resolve_multimapping (List_T hitpairs) {
  List_T resolve1, resolve2, resolve3, p;
  Stage3pair_T hitpair;

  Overlap_T best_overlap;
  long int best_tally;
  double tally_threshold;
  bool runlengthp;


  if (List_length(hitpairs) <= 1) {
    return hitpairs;
  }

  if (genes_iit == NULL) {
    resolve1 = hitpairs;
  } else {
    best_overlap = NO_KNOWN_GENE;
    for (p = hitpairs; p != NULL; p = p->rest) {
      hitpair = (Stage3pair_T) p->first;
      if ((hitpair->gene_overlap = Stage3pair_gene_overlap(hitpair)) > best_overlap) {
	best_overlap = hitpair->gene_overlap;
      }
    }
    if (best_overlap == NO_KNOWN_GENE) {
      resolve1 = hitpairs;
    } else {
      resolve1 = (List_T) NULL;
      for (p = hitpairs; p != NULL; p = p->rest) {
	hitpair = (Stage3pair_T) p->first;
	if (hitpair->gene_overlap < best_overlap) {
	  Stage3pair_free(&hitpair);
	} else {
	  resolve1 = List_push(resolve1,(void *) hitpair);
	}
      }
      List_free(&hitpairs);
    }
  }
      
  if (List_length(resolve1) <= 1) {
    return resolve1;
  }

  if (tally_iit == NULL) {
    resolve2 = resolve1;
  } else {
    best_tally = 0L;
    for (p = resolve1; p != NULL; p = p->rest) {
      hitpair = (Stage3pair_T) p->first;
      if ((hitpair->tally = Stage3end_compute_tally(hitpair->hit5) + Stage3end_compute_tally(hitpair->hit3)) > best_tally) {
	best_tally = hitpair->tally;
      }
    }
    if (best_tally == 0L) {
      resolve2 = resolve1;
    } else {
      resolve2 = (List_T) NULL;
#ifdef USE_TALLY_RATIO
      tally_threshold = (double) best_tally / TALLY_RATIO;
#else
      tally_threshold = 1.0;
#endif
      for (p = resolve1; p != NULL; p = p->rest) {
	hitpair = (Stage3pair_T) p->first;
	if ((double) hitpair->tally < tally_threshold) {
	  Stage3pair_free(&hitpair);
	} else {
	  resolve2 = List_push(resolve2,(void *) hitpair);
	}
      }
      List_free(&resolve1);
    }
  }


  if (List_length(resolve2) <= 1) {
    return resolve2;
  }

  if (runlength_iit == NULL) {
    resolve3 = resolve2;
  } else {
    runlengthp = false;
    for (p = resolve2; p != NULL; p = p->rest) {
      hitpair = (Stage3pair_T) p->first;
      if (Stage3end_runlength_p(hitpair->hit5) == true || Stage3end_runlength_p(hitpair->hit3) == true) {
	runlengthp = true;
      }
    }
    if (runlengthp == false) {
      resolve3 = resolve2;
    } else {
      resolve3 = (List_T) NULL;
      for (p = resolve2; p != NULL; p = p->rest) {
	hitpair = (Stage3pair_T) p->first;
	if (Stage3end_runlength_p(hitpair->hit5) == false && Stage3end_runlength_p(hitpair->hit3) == false) {
	  Stage3pair_free(&hitpair);
	} else {
	  resolve3 = List_push(resolve3,(void *) hitpair);
	}
      }
      List_free(&resolve2);
    }
  }


  return resolve3;
}


List_T
Stage3pair_filter_coverage (List_T hits, int min_coverage_5, int min_coverage_3) {
  List_T newhits = NULL, p;
  Stage3end_T hit5, hit3;
  Stage3pair_T hitpair;

  for (p = hits; p != NULL; p = List_next(p)) {
    hitpair = (Stage3pair_T) List_head(p);
    hit5 = hitpair->hit5;
    hit3 = hitpair->hit3;
    if (hit5->querylength - hit5->trim_left - hit5->trim_right >= min_coverage_5 ||
	hit3->querylength - hit3->trim_left - hit3->trim_right >= min_coverage_3) {
      newhits = List_push(newhits,(void *) hitpair);
    } else {
      Stage3pair_free(&hitpair);
    }
  }

  List_free(&hits);
  return newhits;
}


Stage3pair_T *
Stage3pair_eval_and_sort (int *npaths, int *first_absmq, int *second_absmq,
			  Stage3pair_T *stage3pairarray, int maxpaths,
			  Shortread_T queryseq5, Shortread_T queryseq3,
			  char *queryuc_ptr_5, char *queryrc5,
			  char *queryuc_ptr_3, char *queryrc3,
			  Compress_T query5_compress_fwd, Compress_T query5_compress_rev, 
			  Compress_T query3_compress_fwd, Compress_T query3_compress_rev, 
			  char *quality_string_5, char *quality_string_3) {
  float maxlik, loglik;

  float total, q;
  int mapq_score;
  bool non_terminal_5p, non_terminal_3p;

  int compute_npaths;
  int randomi, i;
  Stage3pair_T temp;

  if (*npaths == 0) {
    /* Skip */
    *first_absmq = 0;
    *second_absmq = 0;

  } else if (*npaths == 1) {
    stage3pairarray[0]->mapq_loglik = MAPQ_MAXIMUM_SCORE;
    stage3pairarray[0]->mapq_score = MAPQ_max_quality_score(quality_string_5,Shortread_fulllength(queryseq5));
    if ((mapq_score = MAPQ_max_quality_score(quality_string_3,Shortread_fulllength(queryseq3))) > stage3pairarray[0]->mapq_score) {
      stage3pairarray[0]->mapq_score = mapq_score;
    }
    stage3pairarray[0]->absmq_score = MAPQ_MAXIMUM_SCORE;

    assert(stage3pairarray[0]->private5p == true);
    assert(stage3pairarray[0]->private3p == true);
    Stage3end_display_prep(stage3pairarray[0]->hit5,queryuc_ptr_5,queryrc5,query5_compress_fwd,query5_compress_rev,
			   stage3pairarray[0]->amb_resolve_5,/*first_read_p*/true);
    Stage3end_display_prep(stage3pairarray[0]->hit3,queryuc_ptr_3,queryrc3,query3_compress_fwd,query3_compress_rev,
			   stage3pairarray[0]->amb_resolve_3,/*first_read_p*/false);
    if (stage3pairarray[0]->amb_resolve_5 >= 0 || stage3pairarray[0]->amb_resolve_3 >= 0) {
      stage3pairarray[0]->insertlength = compute_insertlength(stage3pairarray[0]);
      assert((int) stage3pairarray[0]->insertlength > 0);
    }


    *first_absmq = stage3pairarray[0]->absmq_score;
    *second_absmq = 0;

  } else {
    /* Determine whether to trim terminal ends */
    non_terminal_5p = non_terminal_3p = false;
    for (i = 0; i < *npaths; i++) {
      if (stage3pairarray[i]->hit5->hittype != TERMINAL) {
	non_terminal_5p = true;
      }
      if (stage3pairarray[i]->hit3->hittype != TERMINAL) {
	non_terminal_3p = true;
      }
    }

    /* Resolve ambiguities, needed for computing mapq */
    for (i = 0; i < *npaths; i++) {
      assert(stage3pairarray[i]->private5p == true);
      assert(stage3pairarray[i]->private3p == true);
      Stage3end_display_prep(stage3pairarray[i]->hit5,queryuc_ptr_5,queryrc5,query5_compress_fwd,query5_compress_rev,
			     stage3pairarray[i]->amb_resolve_5,/*first_read_p*/true);
      Stage3end_display_prep(stage3pairarray[i]->hit3,queryuc_ptr_3,queryrc3,query3_compress_fwd,query3_compress_rev,
			     stage3pairarray[i]->amb_resolve_3,/*first_read_p*/false);
      if (stage3pairarray[i]->amb_resolve_5 >= 0 || stage3pairarray[i]->amb_resolve_3 >= 0) {
	stage3pairarray[i]->insertlength = compute_insertlength(stage3pairarray[i]);
      }
    }


    /* Compute mapq_loglik */
    for (i = 0; i < *npaths; i++) {
      stage3pairarray[i]->mapq_loglik =
	Stage3end_compute_mapq(stage3pairarray[i]->hit5,query5_compress_fwd,query5_compress_rev,
			       quality_string_5,/*trim_terminals_p*/non_terminal_5p ? false : true);
      stage3pairarray[i]->mapq_loglik +=
	Stage3end_compute_mapq(stage3pairarray[i]->hit3,query3_compress_fwd,query3_compress_rev,
			       quality_string_3,/*trim_terminals_p*/non_terminal_3p ? false : true);
    }

    /* Sort by nmatches, then mapq, and then insert length */
    qsort(stage3pairarray,*npaths,sizeof(Stage3pair_T),Stage3pair_output_cmp);

    if (want_random_p) {
      /* Randomize among best alignments */
      i = 1;
      while (i < *npaths && Stage3pair_output_cmp(&(stage3pairarray[i]),&(stage3pairarray[0])) == 0) {
	i++;
      }
      if (i > 1) {		/* i is number of ties */
	/* randomi = (int) ((double) i * rand()/((double) RAND_MAX + 1.0)); */
	randomi = (int) (rand() / (((double) RAND_MAX + 1.0) / (double) i));
	/* fprintf(stderr,"%d dups => random %d\n",i,randomi); */
	temp = stage3pairarray[0];
	stage3pairarray[0] = stage3pairarray[randomi];
	stage3pairarray[randomi] = temp;
      }
    }

    /* Enforce monotonicity */
    for (i = *npaths - 1; i > 0; i--) {
      if (stage3pairarray[i-1]->mapq_loglik < stage3pairarray[i]->mapq_loglik) {
	stage3pairarray[i-1]->mapq_loglik = stage3pairarray[i]->mapq_loglik;
      }
    }
    maxlik = stage3pairarray[0]->mapq_loglik;

    /* Subtract maxlik to avoid underflow */
    for (i = 0; i < *npaths; i++) {
      stage3pairarray[i]->mapq_loglik -= maxlik;
    }

    /* Save on computation if possible */
    if (*npaths < maxpaths) {
      compute_npaths = *npaths;
    } else {
      compute_npaths = maxpaths;
    }
    if (compute_npaths < 2) {
      compute_npaths = 2;
    }

    /* Compute absolute mapq */
    for (i = 0; i < compute_npaths; i++) {
      loglik = stage3pairarray[i]->mapq_loglik + MAPQ_MAXIMUM_SCORE;
      if (loglik < 0.0) {
	loglik = 0.0;
      }
      stage3pairarray[i]->absmq_score = rint(loglik);
    }
    *first_absmq = stage3pairarray[0]->absmq_score;
    *second_absmq = stage3pairarray[1]->absmq_score;


    /* Compute Bayesian mapq */
    total = 0.0;
    for (i = 0; i < *npaths; i++) {
      total += (stage3pairarray[i]->mapq_loglik = fasterexp(stage3pairarray[i]->mapq_loglik));
    }

    /* Obtain posterior probabilities of being true */
    for (i = 0; i < compute_npaths; i++) {
      stage3pairarray[i]->mapq_loglik /= total;
    }

    /* Convert to Phred scores */
    for (i = 0; i < compute_npaths; i++) {
      if ((q = 1.0 - stage3pairarray[i]->mapq_loglik) < 2.5e-10 /* 10^-9.6 */) {
	stage3pairarray[i]->mapq_score = 96;
      } else {
	stage3pairarray[i]->mapq_score = rint(-10.0 * log10(q));
      }
    }

#if 0
    /* Apply filtering for mapq unique -- currently not used since mapq_unique_score is high */
    if (stage3pairarray[0]->mapq_score >= mapq_unique_score &&
	stage3pairarray[1]->mapq_score < mapq_unique_score) {
      for (i = 1; i < *npaths; i++) {
	Stage3pair_free(&(stage3pairarray[i]));
      }
      *npaths = 1;
    }
#endif
  }

  return stage3pairarray;
}


List_T
Stage3pair_remove_excess_terminals (List_T hitpairlist) {
  List_T cleaned = NULL, p;
  Stage3pair_T hitpair;
  int splice5p = false, splice3p = false, terminal5p = false, terminal3p = false;
  int best_splice5_score = 1000000, best_splice3_score = 1000000;

  for (p = hitpairlist; p != NULL; p = p->rest) {
    hitpair = (Stage3pair_T) p->first;
    if (hitpair->hit5->hittype == SPLICE || hitpair->hit5->hittype == SHORTEXON) {
      splice5p = true;
      if (hitpair->hit5->score < best_splice5_score) {
	best_splice5_score = hitpair->hit5->score;
      }
    } else if (hitpair->hit5->hittype == TERMINAL) {
      terminal5p = true;
    }

    if (hitpair->hit3->hittype == SPLICE || hitpair->hit3->hittype == SHORTEXON) {
      splice3p = true;
      if (hitpair->hit3->score < best_splice3_score) {
	best_splice3_score = hitpair->hit3->score;
      }
    } else if (hitpair->hit3->hittype == TERMINAL) {
      terminal3p = true;
    }
  }

  if ((splice5p == true && terminal5p == true) ||
      (splice3p == true && terminal3p == true)) {
    for (p = hitpairlist; p != NULL; p = p->rest) {
      hitpair = (Stage3pair_T) p->first;
      if (hitpair->hit5->hittype == TERMINAL && splice5p == true && hitpair->hit5->score >= best_splice5_score) {
	Stage3pair_free(&hitpair);
      } else if (hitpair->hit3->hittype == TERMINAL && splice3p == true && hitpair->hit3->score >= best_splice3_score) {
	Stage3pair_free(&hitpair);
      } else {
	cleaned = List_push(cleaned,hitpair);
      }
    }

    List_free(&hitpairlist);
    return cleaned;

  } else {
    return hitpairlist;
  }

}



static List_T
Stage3pair_optimal_score_aux (bool *eliminatedp, List_T hitpairlist, int cutoff_level, int suboptimal_mismatches,
			      Compress_T query5_compress_fwd, Compress_T query5_compress_rev,
			      Compress_T query3_compress_fwd, Compress_T query3_compress_rev,
			      int querylength5, int querylength3, bool keep_gmap_p, bool finalp) {
  List_T optimal = NULL, p, q;
  Stage3pair_T hitpair;
  T hit5, hit3;
  Substring_T substring;
  int cutoff_level_5, cutoff_level_3, score;
  int n;
  int minscore5 = querylength5, minscore3 = querylength3, minscore = querylength5 + querylength3;
  int best_nsegments;
  /* int max_nmatches = 0, max_nmatches_posttrim; */
#ifdef USE_OPTIMAL_SCORE_BINGO
  int minscore_bingo = querylength5 + querylength3;
#endif
  int trim_left_5 = querylength5, trim_right_5 = querylength5,
    trim_left_3 = querylength3, trim_right_3 = querylength3;
  int nindelbreaks;

#ifdef TRANSLOC_SPECIAL
  bool non_translocation_p = false;
#endif


  *eliminatedp = false;

  n = List_length(hitpairlist);
  debug6(printf("\nEntered Stage3pair_optimal_score with %d hitpairs: %s\n",
		n,finalp == true ? "FINAL" : "not final"));
  
  if (n <= 1) {
    return hitpairlist;
  }


  /* Use eventrim for comparing alignments */
  for (p = hitpairlist; p != NULL; p = p->rest) {
    hitpair = (Stage3pair_T) p->first;
    hit5 = hitpair->hit5;
    hit3 = hitpair->hit3;

    debug6(printf("hit5 %u..%u type %s, nsegments %d, trim_left: %d%s, trim_right %d%s, start_ambig %d, end_ambig %d.  hit3 %u..%u type %s, nsegments %d, trim_left %d%s, trim_right %d%s, start_ambig %d, end_ambig %d, sensedirs %d and %d.\n",
		  hit5->genomicstart - hit5->chroffset,hit5->genomicend - hit5->chroffset,hittype_string(hit5->hittype),
		  hit5->nsegments,hit5->trim_left,hit5->trim_left_splicep ? " (splice)" : "",
		  hit5->trim_right,hit5->trim_right_splicep ? " (splice)" : "",
		  start_amb_length(hit5),end_amb_length(hit5),
		  hit3->genomicstart - hit3->chroffset,hit3->genomicend - hit3->chroffset,hittype_string(hit3->hittype),
		  hit3->nsegments,hit3->trim_left,hit3->trim_left_splicep ? " (splice)" : "",
		  hit3->trim_right,hit3->trim_right_splicep ? " (splice)" : "",
		  start_amb_length(hit3),end_amb_length(hit3),hit5->sensedir,hit3->sensedir));

    if (hit5->hittype == TERMINAL) {
      /* Don't allow terminals to set trims */

#if 0
    } else if ((hit5->hittype == INSERTION || hit5->hittype == DELETION) &&
	(hit5->indel_pos < 15 || hit5->indel_pos > hit5->querylength - 15)) {
      /* Don't allow end indels to set trims */
#endif

    } else {
      if (hit5->trim_left_splicep == true) {
	/* Skip */
      } else if (hit5->trim_left < trim_left_5) {
	trim_left_5 = hit5->trim_left;
      }
      if (hit5->trim_right_splicep == true) {
	/* Skip */
      } else if (hit5->trim_right < trim_right_5) {
	trim_right_5 = hit5->trim_right;
      }
    }

    if (hit3->hittype == TERMINAL) {
      /* Don't allow terminals to set trims */

#if 0
    } else if ((hit3->hittype == INSERTION || hit3->hittype == DELETION) &&
	(hit3->indel_pos < 15 || hit3->indel_pos > hit3->querylength - 15)) {
      /* Don't allow end indels to set trims */
#endif

    } else {
      if (hit3->trim_left_splicep == true) {
	/* Skip */
      } else if (hit3->trim_left < trim_left_3) {
	trim_left_3 = hit3->trim_left;
      }
      if (hit3->trim_right_splicep == true) {
	/* Skip */
      } else if (hit3->trim_right < trim_right_3) {
	trim_right_3 = hit3->trim_right;
      }
    }
  }

  if (trim_left_5 == querylength5) {
    trim_left_5 = 0;
  }
  if (trim_right_5 == querylength5) {
    trim_right_5 = 0;
  }
  if (trim_left_3 == querylength3) {
    trim_left_3 = 0;
  }
  if (trim_right_3 == querylength3) {
    trim_right_3 = 0;
  }

  debug6(printf("overall 5': trim_left %d, trim_right %d\n",trim_left_5,trim_right_5));
  debug6(printf("overall 3': trim_left %d, trim_right %d\n",trim_left_3,trim_right_3));


  for (p = hitpairlist; p != NULL; p = p->rest) {
    hitpair = (Stage3pair_T) p->first;
    hit5 = hitpair->hit5;
    hit3 = hitpair->hit3;

    if (hit5->hittype == GMAP) {
      hit5->score_eventrim = 0;  /* was hit5->penalties */
      debug6(printf("score 5' GMAP:"));
#if 0
      if (Stage3end_bad_stretch_p(hit5,query5_compress_fwd,query5_compress_rev) == true) {
	hit5->score_eventrim += 2;
	debug6(printf("  bad stretch 2."));
      }
#endif

#if 0
      if (0 && hit5->trim_left <= 8) {
	/* Ignore small trims */
      } else if (hit5->trim_left > trim_left_5) {
	hit5->score_eventrim += hit5->trim_left - trim_left_5;
	debug6(printf("  add trim left (%d - %d).",hit5->trim_left,trim_left_5));
      }
      if (0 && hit5->trim_right <= 8) {
	/* Ignore small trims */
      } else if (hit5->trim_right > trim_right_5) {
	hit5->score_eventrim += hit5->trim_right - trim_right_5;
	debug6(printf("  add trim right (%d - %d).",hit5->trim_right,trim_right_5));
      }
#endif

      hit5->score_eventrim += Pair_nmismatches_region(&nindelbreaks,hit5->pairarray,hit5->npairs,
						      trim_left_5,trim_right_5,start_amb_length(hit5),end_amb_length(hit5),
						      hit5->querylength);
      debug6(printf("  add nmismatches %d.",Pair_nmismatches_region(&nindelbreaks,hit5->pairarray,hit5->npairs,
								    trim_left_5,trim_right_5,start_amb_length(hit5),end_amb_length(hit5),
								    hit5->querylength)));
      if (start_amb_length(hit5) > 0) {
	debug6(printf("  add penalty for start amb %d.",amb_penalty));
	hit5->score_eventrim += amb_penalty;
      }
      if (end_amb_length(hit5) > 0) {
	debug6(printf("  add penalty for end amb %d.",amb_penalty));
	hit5->score_eventrim += amb_penalty;
      }

#ifdef SCORE_INDELS
      hit5->score_eventrim += indel_penalty_middle * nindelbreaks;
      debug6(printf("  add indelbreaks %d.",indel_penalty_middle * nindelbreaks));
#endif
      debug6(printf("  RESULT: %d\n",hit5->score_eventrim));
      
    } else {
      hit5->score_eventrim = 0;	/* was hit5->penalties */
      debug6(printf("score 5' OTHER:"));

      for (q = hit5->substrings_1toN; q != NULL; q = List_next(q)) {
	substring = (Substring_T) List_head(q);
	hit5->score_eventrim += Substring_count_mismatches_region(substring,trim_left_5,trim_right_5,
								  query5_compress_fwd,query5_compress_rev);
	debug6(printf("  substring (%d..%d) %d.",trim_left_5,trim_right_5,
		      Substring_count_mismatches_region(substring,trim_left_5,trim_right_5,
							query5_compress_fwd,query5_compress_rev)));
      }

#ifdef SCORE_INDELS
      /* Needs to match GMAP scoring */
      if (hit5->hittype == INSERTION || hit5->hittype == DELETION) {
	debug6(printf("  indel at %d",hit5->indel_pos));
	if (hit5->indel_pos > trim_left_5 && hit5->indel_pos < querylength5 - trim_right_5) {
	  hit5->score_eventrim += indel_penalty_middle;
	  debug6(printf(" => add %d.",indel_penalty_middle));
	}
      }
#endif
      debug6(printf("  RESULT: %d\n",hit5->score_eventrim));
    }

    if (hit3->hittype == GMAP) {
      hit3->score_eventrim = 0;  /* was hit3->penalties */
      debug6(printf("score 3' GMAP:"));
#if 0
      if (Stage3end_bad_stretch_p(hit3,query3_compress_fwd,query3_compress_rev) == true) {
	hit3->score_eventrim += 2;
	debug6(printf("  bad stretch 2."));
      }
#endif

#if 0
      if (0 && hit3->trim_left <= 8) {
	/* Ignore small trims */
      } else if (hit3->trim_left > trim_left_3) {
	hit3->score_eventrim += hit3->trim_left - trim_left_3;
	debug6(printf("  add trim left (%d - %d).",hit3->trim_left,trim_left_3));
      }
      if (0 && hit3->trim_right <= 8) {
	/* Ignore small trims */
      } else if (hit3->trim_right > trim_right_3) {
	hit3->score_eventrim += hit3->trim_right - trim_right_3;
	debug6(printf("  add trim right (%d - %d).",hit3->trim_right,trim_right_3));
      }
#endif

      hit3->score_eventrim += Pair_nmismatches_region(&nindelbreaks,hit3->pairarray,hit3->npairs,
						      trim_left_3,trim_right_3,start_amb_length(hit3),end_amb_length(hit3),
						      hit3->querylength);
      debug6(printf("  add nmismatches %d.",Pair_nmismatches_region(&nindelbreaks,hit3->pairarray,hit3->npairs,
								    trim_left_3,trim_right_3,start_amb_length(hit3),end_amb_length(hit3),
								    hit3->querylength)));

      if (start_amb_length(hit3) > 0) {
	debug6(printf("  add penalty for start amb %d.",amb_penalty));
	hit3->score_eventrim += amb_penalty;
      }
      if (end_amb_length(hit3) > 0) {
	debug6(printf("  add penalty for end amb %d.",amb_penalty));
	hit3->score_eventrim += amb_penalty;
      }

#ifdef SCORE_INDELS
      hit3->score_eventrim += indel_penalty_middle * nindelbreaks;
      debug6(printf("  add indelbreaks %d.",indel_penalty_middle * nindelbreaks));
#endif
#if 0
      if (hit3->start_amb_prob < 0.9) {
	hit3->score_eventrim += hit3->start_amb_length / ambig_end_interval;
	debug6(printf("  add amb start %d/%d (prob %f).",hit3->start_amb_length,ambig_end_interval,hit3->start_amb_prob));
      }
      if (hit3->end_amb_prob < 0.9) {
	hit3->score_eventrim += hit3->end_amb_length / ambig_end_interval;
	debug6(printf("  add amb end %d/%d (prob %f).",hit3->end_amb_length,ambig_end_interval,hit3->end_amb_prob));
      }
#endif
      debug6(printf("  RESULT: %d\n",hit3->score_eventrim));

    } else {
      hit3->score_eventrim = 0;  /* was hit3->penalties */
      debug6(printf("score 3' OTHER:"));

      for (q = hit3->substrings_1toN; q != NULL; q = List_next(q)) {
	substring = (Substring_T) List_head(q);
	hit3->score_eventrim += Substring_count_mismatches_region(substring,trim_left_3,trim_right_3,
								  query3_compress_fwd,query3_compress_rev);
	debug6(printf("  substring (%d..%d) %d.",trim_left_3,trim_right_3,
		      Substring_count_mismatches_region(substring,trim_left_3,trim_right_3,
							query3_compress_fwd,query3_compress_rev)));
      }

#ifdef SCORE_INDELS
      /* Needs to match GMAP scoring */
      if (hit3->hittype == INSERTION || hit3->hittype == DELETION) {
	debug6(printf("  indel at %d",hit3->indel_pos));
	if (hit3->indel_pos > trim_left_3 && hit3->indel_pos < querylength3 - trim_right_3) {
	  hit3->score_eventrim += indel_penalty_middle;
	  debug6(printf(" => add %d.",indel_penalty_middle));
	}
      }
#endif
      debug6(printf("  RESULT: %d\n",hit3->score_eventrim));
    }

    hitpair->score_eventrim = hit5->score_eventrim + hit3->score_eventrim;
    if (hitpair->score_eventrim < minscore) {
      minscore = hitpair->score_eventrim;
    }
  }


  for (p = hitpairlist; p != NULL; p = p->rest) {
    hitpair = (Stage3pair_T) p->first;
    debug6(printf("%u..%u|%u..%u types %s and %s, score_eventrim %d+%d, pairlength %d, outerlength %u\n",
		  hitpair->hit5->low - hitpair->hit5->chroffset,hitpair->hit5->high - hitpair->hit5->chroffset,
		  hitpair->hit3->low - hitpair->hit3->chroffset,hitpair->hit3->high - hitpair->hit3->chroffset,
		  hittype_string(hitpair->hit5->hittype),hittype_string(hitpair->hit3->hittype),
		  hitpair->hit5->score_eventrim,hitpair->hit3->score_eventrim,
		  hitpair->insertlength,hitpair->outerlength));

    if (hitpair->hit5->score_eventrim < minscore5) {
      minscore5 = hitpair->hit5->score_eventrim;
    }
    if (hitpair->hit3->score_eventrim < minscore3) {
      minscore3 = hitpair->hit3->score_eventrim;
    }
  }
  debug6(printf("Stage3pair_optimal_score over %d pairs: minscore = %d and %d + subopt:%d\n",
		n,minscore5,minscore3,suboptimal_mismatches));

  if (finalp == false) {
    /* finalp == false.  Add suboptimal_mismatches to each end. */
    minscore5 += suboptimal_mismatches;
    minscore3 += suboptimal_mismatches;
    cutoff_level_5 = minscore5;
    cutoff_level_3 = minscore3;

    for (p = hitpairlist; p != NULL; p = p->rest) {
      hitpair = (Stage3pair_T) p->first;

      if (keep_gmap_p == true && (hitpair->hit5->hittype == GMAP || hitpair->hit3->hittype == GMAP)) {
	/* GMAP hits already found to be better than their corresponding terminals */
	debug6(printf("Prefinal: Keeping a hit pair of type %s-%s with score_eventrim %d and %d, because keep_gmap_p is true\n",
		      hittype_string(hitpair->hit5->hittype),hittype_string(hitpair->hit3->hittype),
		      hitpair->hit5->score_eventrim,hitpair->hit3->score_eventrim));
	optimal = List_push(optimal,hitpair);

      } else if (hitpair->hit5->score_eventrim > cutoff_level_5 && hitpair->hit3->score_eventrim > cutoff_level_3) {
	debug6(printf("Prefinal: Eliminating a hit pair at %u..%u|%u..%u with score_eventrim_5 %d > cutoff_level_5 %d and score_eventrim_3 %d > cutoff_level_3 %d (finalp %d)\n",
		      hitpair->hit5->low - hitpair->hit5->chroffset,hitpair->hit5->high - hitpair->hit5->chroffset,
		      hitpair->hit3->low - hitpair->hit3->chroffset,hitpair->hit3->high - hitpair->hit3->chroffset,
		      hitpair->hit5->score_eventrim,cutoff_level_5,hitpair->hit3->score_eventrim,cutoff_level_3,finalp));
	*eliminatedp = true;
	Stage3pair_free(&hitpair);

      } else {
	debug6(printf("Prefinal: Keeping a hit pair with score_eventrim %d and %d (cutoff_level %d and %d)\n",
		      hitpair->hit5->score_eventrim,hitpair->hit3->score_eventrim,cutoff_level_5,cutoff_level_3));
	optimal = List_push(optimal,hitpair);
      }
    }

  } else {
    /* finalp == true.  Add suboptimal_mismatches to overall score. */
#if 0
    if (minscore5 + minscore3 < minscore) {
      cutoff_level = minscore + suboptimal_mismatches;
      debug6(printf("cutoff level %d = minscore %d + subopt %d\n",cutoff_level,minscore,suboptimal_mismatches));
    } else {
      cutoff_level = minscore5 + minscore3 + suboptimal_mismatches;
      debug6(printf("cutoff level %d = minscore5 %d + minscore3 %d + subopt %d\n",cutoff_level,minscore5,minscore3,suboptimal_mismatches));
    }
#else
    cutoff_level = minscore + suboptimal_mismatches;
    debug6(printf("cutoff level %d = minscore %d + subopt %d\n",cutoff_level,minscore,suboptimal_mismatches));
#endif


    for (p = hitpairlist; p != NULL; p = p->rest) {
      hitpair = (Stage3pair_T) p->first;
      debug6(printf("Final: %u..%u|%u..%u types %s and %s, score_eventrim %d (%d+%d), pairlength %d, outerlength %u\n",
		    hitpair->hit5->low - hitpair->hit5->chroffset,hitpair->hit5->high - hitpair->hit5->chroffset,
		    hitpair->hit3->low - hitpair->hit3->chroffset,hitpair->hit3->high - hitpair->hit3->chroffset,
		    hittype_string(hitpair->hit5->hittype),hittype_string(hitpair->hit3->hittype),
		    hitpair->score_eventrim,hitpair->hit5->score_eventrim,hitpair->hit3->score_eventrim,
		    hitpair->insertlength,hitpair->outerlength));

      score = hitpair->score_eventrim;

      if (keep_gmap_p == true && (hitpair->hit5->hittype == GMAP || hitpair->hit3->hittype == GMAP)) {
	/* GMAP hits already found to be better than their corresponding terminals */
	debug6(printf("Final: Keeping a hit pair of type %s-%s with score_eventrim %d, because keep_gmap_p is true\n",
		      hittype_string(hitpair->hit5->hittype),hittype_string(hitpair->hit3->hittype),hitpair->score_eventrim));
	optimal = List_push(optimal,hitpair);

#if 0
      } else if (score > cutoff_level) {
	/* Turning off, because it eliminates some shorter, but good, alignments unnecessarily */
	debug6(printf("Final: Eliminating a hit pair at %u..%u|%u..%u with score %d > cutoff_level %d (finalp %d)\n",
		      hitpair->hit5->low - hitpair->hit5->chroffset,hitpair->hit5->high - hitpair->hit5->chroffset,
		      hitpair->hit3->low - hitpair->hit3->chroffset,hitpair->hit3->high - hitpair->hit3->chroffset,
		      score,cutoff_level,finalp));
	*eliminatedp = true;
	Stage3pair_free(&hitpair);
#endif

      } else {
	debug6(printf("Final: Keeping a hit pair with score_eventrim %d (cutoff_level %d)\n",
		      hitpair->score_eventrim,cutoff_level));
	optimal = List_push(optimal,hitpair);
      }
    }
  }


  List_free(&hitpairlist);

  /* Filter on nsegments */
  if (finalp == true && optimal != NULL) {
    hitpairlist = optimal;
    optimal = (List_T) NULL;

    hitpair = (Stage3pair_T) hitpairlist->first;
    best_nsegments = hitpair->hit5->nsegments + hitpair->hit3->nsegments;

    for (p = hitpairlist; p != NULL; p = p->rest) {
      hitpair = (Stage3pair_T) p->first;
      if (hitpair->hit5->nsegments + hitpair->hit3->nsegments < best_nsegments) {
	best_nsegments = hitpair->hit5->nsegments + hitpair->hit3->nsegments;
      }
    }

    for (p = hitpairlist; p != NULL; p = p->rest) {
      hitpair = (Stage3pair_T) p->first;
      if (hitpair->hit5->nsegments + hitpair->hit3->nsegments > best_nsegments + 2) {
	debug6(printf("Eliminating a hit pair with nsegments %d+%d\n",hitpair->hit5->nsegments,hitpair->hit3->nsegments));
	*eliminatedp = true;
	Stage3pair_free(&hitpair);
      } else {
	debug6(printf("Keeping a hit pair with nsegments %d+%d\n",hitpair->hit5->nsegments,hitpair->hit3->nsegments));
	optimal = List_push(optimal,hitpair);
      }
    }

    List_free(&hitpairlist);
  }


#if 0
  /* Filter on pairlength */
  if (optimal != NULL) {
    hitpairlist = optimal;
    optimal = (List_T) NULL;

    hitpair = (Stage3pair_T) hitpairlist->first;
    best_absdifflength = hitpair->absdifflength;
    best_outerlength = hitpair->outerlength;

    for (p = hitpairlist; p != NULL; p = p->rest) {
      hitpair = (Stage3pair_T) p->first;
      if (hitpair->absdifflength < best_absdifflength) {
	best_absdifflength = hitpair->absdifflength;
	best_outerlength = hitpair->outerlength;
      } else if (hitpair->absdifflength > best_absdifflength) {
	/* Skip */
      } else if (hitpair->outerlength < best_outerlength) {
	best_outerlength = hitpair->outerlength;
      }
    }

    for (p = hitpairlist; p != NULL; p = p->rest) {
      hitpair = (Stage3pair_T) p->first;
      if (hitpair->absdifflength > best_absdifflength) {
	debug6(printf("Eliminating a hit pair with absdifflength %d\n",hitpair->absdifflength));
	*eliminatedp = true;
	Stage3pair_free(&hitpair);
      } else if (hitpair->outerlength > best_outerlength + OUTERLENGTH_SLOP) {
	debug6(printf("Eliminating a hit pair with outerlength %u\n",hitpair->outerlength));
	*eliminatedp = true;
	Stage3pair_free(&hitpair);
      } else {
	debug6(printf("Keeping a hit pair with absdifflength %d and outerlength %d\n",
		      hitpair->absdifflength,hitpair->outerlength));
	optimal = List_push(optimal,hitpair);
      }
    }

    List_free(&hitpairlist);
  }
#endif

  return optimal;
}


List_T
Stage3pair_optimal_score (List_T hitpairlist, int cutoff_level, int suboptimal_mismatches,
			  Compress_T query5_compress_fwd, Compress_T query5_compress_rev,
			  Compress_T query3_compress_fwd, Compress_T query3_compress_rev,
			  int querylength5, int querylength3, bool keep_gmap_p, bool finalp) {
  List_T optimal;
  bool eliminatedp;

  optimal = Stage3pair_optimal_score_aux(&eliminatedp,hitpairlist,cutoff_level,suboptimal_mismatches,
					 query5_compress_fwd,query5_compress_rev,
					 query3_compress_fwd,query3_compress_rev,
					 querylength5,querylength3,keep_gmap_p,finalp);
  while (eliminatedp == true) {
    optimal = Stage3pair_optimal_score_aux(&eliminatedp,optimal,cutoff_level,suboptimal_mismatches,
					   query5_compress_fwd,query5_compress_rev,
					   query3_compress_fwd,query3_compress_rev,
					   querylength5,querylength3,keep_gmap_p,finalp);
  }

  return optimal;
}



bool
Stage3pair_sense_consistent_p (List_T hitpairlist) {
  Stage3pair_T hitpair;
  T hit5, hit3;
  List_T p;

  for (p = hitpairlist; p != NULL; p = List_next(p)) {
    hitpair = (Stage3pair_T) List_head(p);
    hit5 = hitpair->hit5;
    hit3 = hitpair->hit3;
    if (hit5->sensedir == hit3->sensedir) {
      return true;
    } else if (hit5->hittype == GMAP && hit5->gmap_nintrons > 0 && hit5->sensedir == 0) {
      /* false */
    } else if (hit3->hittype == GMAP && hit3->gmap_nintrons > 0 && hit3->sensedir == 0) {
      /* false */
    } else {
      return true;
    }
  }
  return false;
}


/* Want to unalias plus and alias minus */
List_T
Stage3end_linearize_5 (List_T hitlist) {
  Chrpos_T chrlength;
  T hit;
  List_T p;

  for (p = hitlist; p != NULL; p = List_next(p)) {
    hit = (T) List_head(p);
    debug12(chrlength = hit->chrlength);
    debug12(printf("Looking at 5' end %u..%u against chrlength %u\n",
		   hit->genomicstart - hit->chroffset,hit->genomicend - hit->chroffset,chrlength));

    if (hit->alias == 0) {
      /* Skip */

    } else if (hit->alias == +1) {
      if (hit->plusp == true) {
	unalias_circular(hit);
      }

    } else if (hit->alias == -1) {
      if (hit->plusp == false) {
	alias_circular(hit);
      }
    }
  }

  return hitlist;
}


/* Want to alias plus and unalias minus */
List_T
Stage3end_linearize_3 (List_T hitlist) {
  Chrpos_T chrlength;
  T hit;
  List_T p;

  for (p = hitlist; p != NULL; p = List_next(p)) {
    hit = (T) List_head(p);
    debug12(chrlength = hit->chrlength);
    debug12(printf("Looking at 3' end %u..%u against chrlength %u\n",
		   hit->genomicstart - hit->chroffset,hit->genomicend - hit->chroffset,chrlength));

    if (hit->alias == 0) {
      /* Skip */

    } else if (hit->alias == -1) {
      if (hit->plusp == true) {
	alias_circular(hit);
      }

    } else if (hit->alias == +1) {
      if (hit->plusp == false) {
	unalias_circular(hit);
      }
    }
  }

  return hitlist;
}



List_T
Stage3pair_remove_circular_alias (List_T hitpairlist) {
  List_T newlist = NULL, p;
  Stage3pair_T hitpair;
  int trim;

  debug12(printf("Stage3pair_remove_circular_alias called with %d hitpairs\n",
		 List_length(hitpairlist)));
  for (p = hitpairlist; p != NULL; p = p->rest) {
    hitpair = (Stage3pair_T) p->first;

#if 0
    /* Not sure if this is necessary */
    if (hitpair->hit5->alias == +1 && hitpair->hit3->alias == +1) {
      /* First, try to salvage alias +1 */
      unalias_circular(hitpair->hit5);
      unalias_circular(hitpair->hit3);
    }
#endif

    if (hitpair->hit5->plusp == true) {
      trim = hitpair->hit5->trim_left;
    } else {
      trim = hitpair->hit3->trim_right;
    }

    if (hitpair->low + trim >= hitpair->hit5->chroffset + hitpair->hit5->chrlength) {
      /* Both ends in circular alias */
      debug12(printf("Both ends in circular alias\n"));
      Stage3pair_free(&hitpair);

    } else {
      newlist = List_push(newlist,(void *) hitpair);
    }
  }

  List_free(&hitpairlist);
  return newlist;
}


static List_T
pair_up_concordant_aux (bool *abort_pairing_p, int *found_score, int *nconcordant, int *nsamechr,
			List_T *samechr, List_T *conc_transloc, List_T hitpairs,
			T **hits5_plus, int *nhits5_plus, T **hits5_minus, int *nhits5_minus,
			T **hits3_plus, int *nhits3_plus, T **hits3_minus, int *nhits3_minus,
			bool *sorted5p, bool *sorted3p,
			int cutoff_level_5, int cutoff_level_3, int subopt_levels,

			Univcoord_T *splicesites,
			Compress_T query5_compress_fwd, Compress_T query5_compress_rev,
			Compress_T query3_compress_fwd, Compress_T query3_compress_rev,
			int querylength5, int querylength3, int maxpairedpaths,
			int splicing_penalty, int genestrand) {
  int new_found_score = *found_score;
  int frontier_score, score5_start, score5_end, score5, score3, i, j;
  List_T q, prev_start;
  Stage3pair_T stage3pair;
  T *hits5, *hits3, hit5, hit3;
  int nhits5, nhits3;
  Univcoord_T insert_start;


  prev_start = hitpairs;
  frontier_score = 0;
  while (*abort_pairing_p == false && frontier_score <= *found_score + subopt_levels &&
	 frontier_score <= cutoff_level_5 + cutoff_level_3) {
    debug5a(printf("frontier_score = %d\n",frontier_score));
    if ((score5_start = frontier_score - cutoff_level_3) < 0) {
      score5_start = 0;
    }
    score5_end = (cutoff_level_5 < frontier_score) ? cutoff_level_5 : frontier_score;
    for (score5 = score5_start; score5 <= score5_end; score5++) {
      debug5a(printf("score5 = %d (cutoff %d), score3 = %d (cutoff %d)\n",
		     score5,cutoff_level_5,frontier_score-score5,cutoff_level_3));
      score3 = frontier_score - score5;
      assert(score3 <= cutoff_level_3);
      if (1 || (score5 <= cutoff_level_5 && ((score3 = frontier_score - score5) <= cutoff_level_3))) {
	/* Sort this level if necessary: 5' by genomicend, 3' by genomicstart */
	if (sorted5p[score5] == false) {
	  if (nhits5_plus[score5] > 0) {
	    qsort(hits5_plus[score5],nhits5_plus[score5],sizeof(T),genomicend_cmp);
	  }
	  if (nhits5_minus[score5] > 0) {
	    qsort(hits5_minus[score5],nhits5_minus[score5],sizeof(T),genomicend_cmp);
	  }
	  sorted5p[score5] = true;
	}
	if (sorted3p[score3] == false) {
	  if (nhits3_plus[score3] > 0) {
	    qsort(hits3_plus[score3],nhits3_plus[score3],sizeof(T),genomicstart_cmp);
	  }
	  if (nhits3_minus[score3] > 0) {
	    qsort(hits3_minus[score3],nhits3_minus[score3],sizeof(T),genomicstart_cmp);
	  }
	  sorted3p[score3] = true;
	}

	/* plus/plus: hits5_plus against hits3_plus (really on minus) */
	hits5 = hits5_plus[score5];
	hits3 = hits3_plus[score3];
	nhits5 = nhits5_plus[score5];
	nhits3 = nhits3_plus[score3];
	debug5(printf("at score %d, nhits5_plus = %d; at score %d, nhits3_plus = %d\n",
		      score5,nhits5,score3,nhits3));

	if (nhits5 > 0 && nhits3 > 0) {
	  i = j = 0;
	  while (*abort_pairing_p == false && i < nhits5) {
	    hit5 = hits5[i];
	    insert_start = hit5->genomicend - querylength5;
	    debug5(printf("plus/plus: i=%d/%d %u..%u %s %s %p\n",
			  i,nhits5,hit5->genomicstart - hit5->chroffset,hit5->genomicend - hit5->chroffset,
			  print_sense(hit5->sensedir),hittype_string(hit5->hittype),hit5));

	    while (j >= 0 && 
		   hits3[j]->genomicstart + querylength3 /* for scramble: */ + pairmax > insert_start) {
	      debug5(printf("  backup: j=%d/%d %u..%u %s %s %p\n",
			    j,nhits3,hits3[j]->genomicstart - hits3[j]->chroffset,hits3[j]->genomicend - hits3[j]->chroffset,
			    print_sense(hits3[j]->sensedir),hittype_string(hits3[j]->hittype),hits3[j]));
	      j--;
	    }
	    j++;		/* Finish backup */

	    while (j < nhits3 && 
		   hits3[j]->genomicstart + querylength3 /* for scramble: */ + pairmax <= insert_start) {
	      debug5(printf("  advance: j=%d/%d %u..%u %s %s %p\n",
			    j,nhits3,hits3[j]->genomicstart - hits3[j]->chroffset,hits3[j]->genomicend - hits3[j]->chroffset,
			    print_sense(hits3[j]->sensedir),hittype_string(hits3[j]->hittype),hits3[j]));
	      j++;
	    }

	    while (j < nhits3 && hits3[j]->genomicstart + querylength3 <= pairmax + insert_start) {
	      debug5(printf("  overlap: j=%d/%d %u..%u %s %s %p",
			    j,nhits3,hits3[j]->genomicstart - hits3[j]->chroffset,hits3[j]->genomicend - hits3[j]->chroffset,
			    print_sense(hits3[j]->sensedir),hittype_string(hits3[j]->hittype),hits3[j]));
	      hit3 = hits3[j];
		
	      /* Want only pairs not previously seen */
	      if (hit5->paired_seenp == false || hit3->paired_seenp == false) {
		if (hit5->effective_chrnum != hit3->effective_chrnum) {
		  debug5(printf(" => diff chrs %d and %d",hit5->effective_chrnum,hit3->effective_chrnum));
		} else if (hit5->chrnum == 0 && hit3->chrnum == 0 /* && hit5->other_chrnum != hit3->other_chrnum */) {
		  /* Could potentially miss an alignment if the two ends overlap */
		  debug5(printf(" => double splice translocations"));
		} else if (hit5->chrnum == 0 || hit3->chrnum == 0) {
		  debug5(printf(" => conc_transloc effchr %d (chrnum5 %d, chrnum3 %d)",
				hit5->effective_chrnum,hit5->chrnum,hit3->chrnum));
		  if ((stage3pair = Stage3pair_new(hit5,hit3,splicesites,
						   query5_compress_fwd,query5_compress_rev,
						   query3_compress_fwd,query3_compress_rev,genestrand,
						   /*pairtype*/CONCORDANT_TRANSLOCATIONS,splicing_penalty,
						   /*private5p*/false,/*private3p*/false,/*expect_concordant_p*/true)) != NULL) {
		    *conc_transloc = List_push(*conc_transloc,(void *) stage3pair);
		  }

		} else if (SENSE_INCONSISTENT_P(hit5->sensedir,hit3->sensedir)) {
		  debug5(printf(" => sense inconsistent: %d | %d = %d",hit5->sensedir,hit3->sensedir,hit5->sensedir|hit3->sensedir));
		} else if (hit3->genomicend < hit5->genomicstart) {
		  debug5(printf(" => scramble because end3 %llu < start5 %llu\n",
				(unsigned long long) hit3->genomicend,(unsigned long long) hit5->genomicstart));
		  if (*nsamechr <= maxpairedpaths &&
		      (stage3pair = Stage3pair_new(Stage3end_copy(hit5),Stage3end_copy(hit3),splicesites,
						   query5_compress_fwd,query5_compress_rev,
						   query3_compress_fwd,query3_compress_rev,genestrand,
						   /*pairtype*/PAIRED_SCRAMBLE,splicing_penalty,
						   /*private5p*/true,/*private3p*/true,/*expect_concordant_p*/false)) != NULL) {
		    *samechr = List_push(*samechr,(void *) stage3pair);
		    (*nsamechr)++;
		  }
		} else {
		  debug5(printf(" => concordant effchr %d (chrnum5 %d, chrnum3 %d)",
				hit5->effective_chrnum,hit5->chrnum,hit3->chrnum));
		  if ((stage3pair = Stage3pair_new(hit5,hit3,splicesites,
						   query5_compress_fwd,query5_compress_rev,
						   query3_compress_fwd,query3_compress_rev,genestrand,
						   /*pairtype*/CONCORDANT,splicing_penalty,
						   /*private5p*/false,/*private3p*/false,/*expect_concordant_p*/true)) != NULL) {

		    debug5(printf("Have new pair with scores %d + %d, compared with new_found_score %d\n",hit5->score,hit3->score,new_found_score));
#if 0
		    /* Doesn't work with new substrings */
		    if (start_amb_length(hit5) > 0 || end_amb_length(hit5) > 0 || 
			start_amb_length(hit3) > 0 || end_amb_length(hit3) > 0) {
		      /* Don't use ambiguous splices to update found_score*/
		      hitpairs = List_push(hitpairs,(void *) stage3pair);
		      (*nconcordant)++;

		    } /*else*/
#endif
		    if (hit5->score + hit3->score < new_found_score) {
		      /* Don't use frontier_score here, which is the trimmed_score.  Use the full score, to motivate stage1hr to find longer alignments */
		      new_found_score = hit5->score + hit3->score;
		      debug5(printf(" => tentatively updating found_score to be %d = %d + %d\n",new_found_score,hit5->score,hit3->score));
		      hitpairs = List_push(hitpairs,(void *) stage3pair);
		      (*nconcordant)++;

		    } else {
		      hitpairs = List_push(hitpairs,(void *) stage3pair);
		      (*nconcordant)++;
		    }

		    if (0 && *nconcordant > maxpairedpaths) {
		      debug(printf(" -- %d concordant paths exceeds %d",*nconcordant,maxpairedpaths));
		      *abort_pairing_p = true;
		    }
		  }
		}
	      }
	      debug5(printf("\n"));

	      j++;
	    }
	    j--;		/* Finish advance */

	    i++;
	  }
	}

	/* minus/minus: hits3_minus (really on plus) against hits5_minus */
	hits3 = hits3_minus[score3];
	hits5 = hits5_minus[score5];
	nhits3 = nhits3_minus[score3];
	nhits5 = nhits5_minus[score5];
	debug5(printf("at score %d, nhits5_minus = %d; at score %d, nhits3_minus = %d\n",
		      score5,nhits5,score3,nhits3));

	if (nhits3 > 0 && nhits5 > 0) {
	  i = j = 0;
	  while (*abort_pairing_p == false && i < nhits3) {
	    hit3 = hits3[i];
	    insert_start = hit3->genomicstart - querylength3;
	    debug5(printf("minus/minus: i=%d/%d %u..%u %s %s %p\n",
			  i,nhits3,hit3->genomicstart - hit3->chroffset,hit3->genomicend - hit3->chroffset,
			  print_sense(hit3->sensedir),hittype_string(hit3->hittype),hit3));

	    while (j >= 0 && 
		   hits5[j]->genomicend + querylength5 /* for scramble: */ + pairmax > insert_start) {
	      debug5(printf("  backup: j=%d/%d %u..%u %s %s %p\n",
			    j,nhits5,hits5[j]->genomicstart - hits5[j]->chroffset,hits5[j]->genomicend - hits5[j]->chroffset,
			    print_sense(hits5[j]->sensedir),hittype_string(hits5[j]->hittype),hits5[j]));
	      j--;
	    }
	    j++;			/* Finish backup */

	    while (j < nhits5 && 
		   hits5[j]->genomicend + querylength5 /* for scramble: */ + pairmax <= insert_start) {
	      debug5(printf("  advance: j=%d/%d %u..%u %s %s %p\n",
			    j,nhits5,hits5[j]->genomicstart - hits5[j]->chroffset,hits5[j]->genomicend - hits5[j]->chroffset,
			    print_sense(hits5[j]->sensedir),hittype_string(hits5[j]->hittype),hits5[j]));
	      j++;
	    }

	    while (j < nhits5 && hits5[j]->genomicend + querylength5 <= pairmax + insert_start) {
	      debug5(printf("  overlap: j=%d/%d %u..%u %s %s %p",
			    j,nhits5,hits5[j]->genomicstart - hits5[j]->chroffset,hits5[j]->genomicend - hits5[j]->chroffset,
			    print_sense(hits5[j]->sensedir),hittype_string(hits5[j]->hittype),hits5[j]));
	      hit5 = hits5[j];

	      /* Want only pairs not previously seen */
	      if (hit3->paired_seenp == false || hit5->paired_seenp == false) {
		if (hit3->effective_chrnum != hit5->effective_chrnum) {
		  debug5(printf(" => diff chrs %d and %d",hit5->effective_chrnum,hit3->effective_chrnum));
		} else if (hit5->chrnum == 0 && hit3->chrnum == 0 /* && hit5->other_chrnum != hit3->other_chrnum */) {
		  /* Could potentially miss an alignment if the two ends overlap */
		  debug5(printf(" => double splice translocations"));

		} else if (hit5->chrnum == 0 || hit3->chrnum == 0) {
		  debug5(printf(" => conc_transloc effchr %d (chrnum5 %d, chrnum3 %d)",
				hit3->effective_chrnum,hit5->chrnum,hit3->chrnum));
		  if ((stage3pair = Stage3pair_new(hit5,hit3,splicesites,
						   query5_compress_fwd,query5_compress_rev,
						   query3_compress_fwd,query3_compress_rev,genestrand,
						   /*pairtype*/CONCORDANT_TRANSLOCATIONS,splicing_penalty,
						   /*private5p*/false,/*private3p*/false,/*expect_concordant_p*/true)) != NULL) {
		    *conc_transloc = List_push(*conc_transloc,(void *) stage3pair);
		  }

		} else if (SENSE_INCONSISTENT_P(hit3->sensedir,hit5->sensedir)) {
		  debug5(printf(" => sense inconsistent: %d | %d = %d",hit5->sensedir,hit3->sensedir,hit5->sensedir|hit3->sensedir));
		} else if (hit5->genomicstart < hit3->genomicend) {
		  debug5(printf(" => scramble because start5 %llu < end3 %llu\n",
				(unsigned long long) hit5->genomicstart,(unsigned long long) hit3->genomicend));
		  if (*nsamechr <= maxpairedpaths &&
		      (stage3pair = Stage3pair_new(Stage3end_copy(hit5),Stage3end_copy(hit3),splicesites,
						   query5_compress_fwd,query5_compress_rev,
						   query3_compress_fwd,query3_compress_rev,genestrand,
						   /*pairtype*/PAIRED_SCRAMBLE,splicing_penalty,
						   /*private5p*/true,/*private3p*/true,/*expect_concordant_p*/false)) != NULL) {
		    *samechr = List_push(*samechr,(void *) stage3pair);
		    (*nsamechr)++;
		  }
		} else {
		  debug5(printf(" => concordant effchr %d (chrnum5 %d, chrnum3 %d)",
				hit3->effective_chrnum,hit5->chrnum,hit3->chrnum));
		  if ((stage3pair = Stage3pair_new(hit5,hit3,splicesites,
						   query5_compress_fwd,query5_compress_rev,
						   query3_compress_fwd,query3_compress_rev,genestrand,
						   /*pairtype*/CONCORDANT,splicing_penalty,
						   /*private5p*/false,/*private3p*/false,/*expect_concordant_p*/true)) != NULL) {

		    debug5(printf("Have new pair with scores %d + %d, compared with new_found_score %d\n",hit5->score,hit3->score,new_found_score));
#if 0
		    /* Doesn't work with new substrings */
		    if (start_amb_length(hit5) > 0 || end_amb_length(hit5) > 0 || 
			start_amb_length(hit3) > 0 || end_amb_length(hit3) > 0) {
		      /* Don't use ambiguous splices to update found_score*/
		      hitpairs = List_push(hitpairs,(void *) stage3pair);
		      (*nconcordant)++;

		    } /*else*/
#endif
		    if (hit5->score + hit3->score < new_found_score) {
		      /* Don't use frontier_score here, which is the trimmed_score.  Use the full score, to motivate stage1hr to find longer alignments */
		      new_found_score = hit5->score + hit3->score;
		      debug5(printf(" => tentatively updating found_score to be %d = %d + %d\n",new_found_score,hit5->score,hit3->score));
		      hitpairs = List_push(hitpairs,(void *) stage3pair);
		      (*nconcordant)++;

		    } else {
		      hitpairs = List_push(hitpairs,(void *) stage3pair);
		      (*nconcordant)++;
		    }

		    if (0 && *nconcordant > maxpairedpaths) {
		      debug(printf(" -- %d concordant paths exceeds %d",*nconcordant,maxpairedpaths));
		      *abort_pairing_p = true;
		    }
		  }
		}
	      }
	      debug5(printf("\n"));

	      j++;
	    }
	    j--;		/* Finish advance */

	    i++;
	  }
	}

	/* plus/minus (inversions): hits5_plus against hits3_minus */
	hits5 = hits5_plus[score5];
	hits3 = hits3_minus[score3];
	nhits5 = nhits5_plus[score5];
	nhits3 = nhits3_minus[score3];
	debug5(printf("at score %d, nhits5_plus = %d; at score %d, nhits3_minus = %d\n",
		      score5,nhits5,score3,nhits3));

	if (nhits5 > 0 && nhits3 > 0) {
	  i = j = 0;
	  while (*abort_pairing_p == false && i < nhits5) {
	    hit5 = hits5[i];
	    insert_start = hit5->genomicend - querylength5;
	    debug5(printf("plus/minus: i=%d/%d %u..%u %s %s %p\n",
			  i,nhits5,hit5->genomicstart - hit5->chroffset,hit5->genomicend - hit5->chroffset,
			  print_sense(hit5->sensedir),hittype_string(hit5->hittype),hit5));

	    while (j >= 0 && 
		   hits3[j]->genomicstart + querylength3 /* for scramble: */ + pairmax > insert_start) {
	      debug5(printf("  backup: j=%d/%d %u..%u %s %s %p\n",
			    j,nhits3,hits3[j]->genomicstart - hits3[j]->chroffset,hits3[j]->genomicend - hits3[j]->chroffset,
			    print_sense(hits3[j]->sensedir),hittype_string(hits3[j]->hittype),hits3[j]));
	      j--;
	    }
	    j++;		/* Finish backup */

	    while (j < nhits3 && 
		   hits3[j]->genomicstart + querylength3 /* for scramble: */ + pairmax <= insert_start) {
	      debug5(printf("  advance: j=%d/%d %u..%u %s %s %p\n",
			    j,nhits3,hits3[j]->genomicstart - hits3[j]->chroffset,hits3[j]->genomicend - hits3[j]->chroffset,
			    print_sense(hits3[j]->sensedir),hittype_string(hits3[j]->hittype),hits3[j]));
	      j++;
	    }

	    while (j < nhits3 && hits3[j]->genomicstart + querylength3 <= pairmax + insert_start) {
	      debug5(printf("  overlap: j=%d/%d %u..%u %s %s %p",
			    j,nhits3,hits3[j]->genomicstart - hits3[j]->chroffset,hits3[j]->genomicend - hits3[j]->chroffset,
			    print_sense(hits3[j]->sensedir),hittype_string(hits3[j]->hittype),hits3[j]));
	      hit3 = hits3[j];
		
	      /* Want only pairs not previously seen */
	      if (hit5->paired_seenp == false || hit3->paired_seenp == false) {
		if (hit5->effective_chrnum != hit3->effective_chrnum) {
		  debug5(printf(" => diff chrs %d and %d",hit5->effective_chrnum,hit3->effective_chrnum));
		} else if (hit5->chrnum == 0 && hit3->chrnum == 0 /* && hit5->other_chrnum != hit3->other_chrnum */) {
		  debug5(printf(" => double splice translocations"));
		} else if (SENSE_INCONSISTENT_FOR_INVERSION_P(hit5->sensedir,hit3->sensedir)) {
		  debug5(printf(" => sense inconsistent for inversion"));
#if 0
		} else if (hits3[j]->genomicstart + querylength3 <= insert_start) {
		  debug5(printf(" => scramble"));
		  if (*nsamechr <= maxpairedpaths &&
		      (stage3pair = Stage3pair_new(Stage3end_copy(hit5),Stage3end_copy(hit3),splicesites,
						   query5_compress_fwd,query5_compress_rev,
						   query3_compress_fwd,query3_compress_rev,genestrand,
						   /*pairtype*/PAIRED_SCRAMBLE,splicing_penalty,
						   /*private5p*/true,/*private3p*/true,/*expect_concordant_p*/false)) != NULL) {
		    *samechr = List_push(*samechr,(void *) stage3pair);
		    (*nsamechr)++;
		  }
#endif
		} else {
		  debug5(printf(" => inversion effchr %d (chrnum5 %d, chrnum3 %d)",
				hit5->effective_chrnum,hit5->chrnum,hit3->chrnum));
		  if (*nsamechr <= maxpairedpaths &&
		      (stage3pair = Stage3pair_new(Stage3end_copy(hit5),Stage3end_copy(hit3),splicesites,
						   query5_compress_fwd,query5_compress_rev,
						   query3_compress_fwd,query3_compress_rev,genestrand,
						   /*pairtype*/PAIRED_INVERSION,splicing_penalty,
						   /*private5p*/true,/*private3p*/true,/*expect_concordant_p*/false)) != NULL) {
		    *samechr = List_push(*samechr,(void *) stage3pair);
		    (*nsamechr)++;
		  }
		}
	      }
	      debug5(printf("\n"));

	      j++;
	    }
	    j--;		/* Finish advance */

	    i++;
	  }
	}

	/* minus/plus: hits3_plus against hits5_minus */
	hits3 = hits3_plus[score3];
	hits5 = hits5_minus[score5];
	nhits3 = nhits3_plus[score3];
	nhits5 = nhits5_minus[score5];
	debug5(printf("at score %d, nhits5_minus = %d; at score %d, nhits3_plus = %d\n",
		      score5,nhits5,score3,nhits3));

	if (nhits3 > 0 && nhits5 > 0) {
	  i = j = 0;
	  while (*abort_pairing_p == false && i < nhits3) {
	    hit3 = hits3[i];
	    insert_start = hit3->genomicstart - querylength3;
	    debug5(printf("minus/plus: i=%d/%d %u..%u %s %s %p\n",
			  i,nhits3,hit3->genomicstart - hit3->chroffset,hit3->genomicend - hit3->chroffset,
			  print_sense(hit3->sensedir),hittype_string(hit3->hittype),hit3));

	    while (j >= 0 && 
		   hits5[j]->genomicend + querylength5 /* for scramble: */ + pairmax > insert_start) {
	      debug5(printf("  backup: j=%d/%d %u..%u %s %s %p\n",
			    j,nhits5,hits5[j]->genomicstart - hits5[j]->chroffset,hits5[j]->genomicend - hits5[j]->chroffset,
			    print_sense(hits5[j]->sensedir),hittype_string(hits5[j]->hittype),hits5[j]));
	      j--;
	    }
	    j++;			/* Finish backup */

	    while (j < nhits5 && 
		   hits5[j]->genomicend + querylength5 /* for scramble: */ + pairmax <= insert_start) {
	      debug5(printf("  advance: j=%d/%d %u..%u %s %s %p\n",
			    j,nhits5,hits5[j]->genomicstart - hits5[j]->chroffset,hits5[j]->genomicend - hits5[j]->chroffset,
			    print_sense(hits5[j]->sensedir),hittype_string(hits5[j]->hittype),hits5[j]));
	      j++;
	    }

	    while (j < nhits5 && hits5[j]->genomicend + querylength5 <= pairmax + insert_start) {
	      debug5(printf("  overlap: j=%d/%d %u..%u %s %s %p",
			    j,nhits5,hits5[j]->genomicstart - hits5[j]->chroffset,hits5[j]->genomicend - hits5[j]->chroffset,
			    print_sense(hits5[j]->sensedir),hittype_string(hits5[j]->hittype),hits5[j]));
	      hit5 = hits5[j];

	      /* Want only pairs not previously seen */
	      if (hit3->paired_seenp == false || hit5->paired_seenp == false) {
		if (hit3->effective_chrnum != hit5->effective_chrnum) {
		  debug5(printf(" => diff chrs %d and %d",hit5->effective_chrnum,hit3->effective_chrnum));
		} else if (hit5->chrnum == 0 && hit3->chrnum == 0 /* && hit5->other_chrnum != hit3->other_chrnum */) {
		  debug5(printf(" => double splice translocations"));
		} else if (SENSE_INCONSISTENT_FOR_INVERSION_P(hit3->sensedir,hit5->sensedir)) {
		  debug5(printf(" => sense inconsistent for inversion"));
#if 0
		} else if (hits5[j]->genomicend + querylength5 <= insert_start) {
		  debug5(printf(" => scramble"));
		  if (*nsamechr <= maxpairedpaths &&
		      (stage3pair = Stage3pair_new(Stage3end_copy(hit5),Stage3end_copy(hit3),splicesites,
						   query5_compress_fwd,query5_compress_rev,
						   query3_compress_fwd,query3_compress_rev,genestrand,
						   /*pairtype*/PAIRED_SCRAMBLE,splicing_penalty,
						   /*private5p*/true,/*private3p*/true,/*expect_concordant_p*/false)) != NULL) {
		    *samechr = List_push(*samechr,(void *) stage3pair);
		    (*nsamechr)++;
		  }
#endif
		} else {
		  debug5(printf(" => inversion effchr %d (chrnum5 %d, chrnum3 %d)",
				hit3->effective_chrnum,hit5->chrnum,hit3->chrnum));
		  if (*nsamechr <= maxpairedpaths &&
		      (stage3pair = Stage3pair_new(Stage3end_copy(hit5),Stage3end_copy(hit3),splicesites,
						   query5_compress_fwd,query5_compress_rev,
						   query3_compress_fwd,query3_compress_rev,genestrand,
						   /*pairtype*/PAIRED_INVERSION,splicing_penalty,
						   /*private5p*/true,/*private3p*/true,/*expect_concordant_p*/false)) != NULL) {
		    *samechr = List_push(*samechr,(void *) stage3pair);
		    (*nsamechr)++;
		  }
		}
	      }
	      debug5(printf("\n"));

	      j++;
	    }
	    j--;		/* Finish advance */

	    i++;
	  }
	}
      }
    }

    /* Mark all concordant pairs (and terminal and double_terminal) found at this pairscore level */
    for (q = hitpairs; q != prev_start; q = List_next(q)) {
      stage3pair = (Stage3pair_T) List_head(q);
      stage3pair->hit5->paired_seenp = true;
      stage3pair->hit3->paired_seenp = true;
    }
    prev_start = hitpairs;

    if (*abort_pairing_p == false) {
      *found_score = new_found_score;
    }

    frontier_score++;
  }

  return hitpairs;
}


/* Use nmismatches_bothdiff (which is a trimmed_score), instead of stage3end->score */
static int
sort_hits_by_trimmed_score (T **hits_plus, T **hits_minus, int *nhits_plus, int *nhits_minus,
			    List_T *hitarray, int narray, int cutoff_level) {
  int score;
  int nhits, i;
  List_T q;
  T hit;

  /* Find sizes for allocating memory */
  debug5(printf("Sizes of pieces by score level:\n"));
  for (i = 0; i < narray; i++) {
    debug5(printf("  array score level %d with %d hits\n",i,List_length(hitarray[i])));
    for (q = hitarray[i]; q != NULL; q = q->rest) {
      hit = (T) q->first;
      debug5(printf(" : %p score %d, type %s\n",hit,hit->score,hittype_string(hit->hittype)));
      assert(hit->nmismatches_bothdiff >= 0);
      if (hit->nmismatches_bothdiff > cutoff_level) {
	debug5(printf("Skipping hit with trimmed score %d > cutoff level %d\n",hit->nmismatches_bothdiff,cutoff_level));
      } else if (hit->plusp == true) {
	nhits_plus[hit->nmismatches_bothdiff]++;
      } else {
	nhits_minus[hit->nmismatches_bothdiff]++;
      }
    }
  }

  debug5(
	 printf("Sizes of pieces by score level and plus/minus:\n");
	 for (score = 0; score <= cutoff_level; score++) {
	   printf("  score %d: %d plus, %d minus\n",score,nhits_plus[score],nhits_minus[score]);
	 }
	 );


  /* Reset cutoff_level */
  score = 0;
  nhits = nhits_plus[score] + nhits_minus[score];
  while (score+1 <= cutoff_level && nhits + nhits_plus[score+1] + nhits_minus[score+1] < MAX_HITS) {
    nhits += nhits_plus[score+1] + nhits_minus[score+1];
    score++;
    debug5(printf("Allowing score to go to %d, because nhits = %d\n",score,nhits));
  }
  debug5(printf("Resetting cutoff_level to be %d\n",score));
  cutoff_level = score;


  /* Store hits */
  for (score = 0; score <= cutoff_level; score++) {
    if (nhits_plus[score] == 0) {
      hits_plus[score] = (T *) NULL;
    } else {
      hits_plus[score] = (T *) MALLOC(nhits_plus[score] * sizeof(Stage3end_T));
    }
  }

  for (score = 0; score <= cutoff_level; score++) {
    if (nhits_minus[score] == 0) {
      hits_minus[score] = (T *) NULL;
    } else {
      hits_minus[score] = (T *) MALLOC(nhits_minus[score] * sizeof(Stage3end_T));
    }
  }

  for (score = 0; score <= cutoff_level; score++) {
    nhits_plus[score] = 0;
    nhits_minus[score] = 0;
  }

  for (i = 0; i < narray; i++) {
    for (q = hitarray[i]; q != NULL; q = q->rest) {
      hit = (T) q->first;
      if (hit->nmismatches_bothdiff > cutoff_level) {
	/* Skip */
      } else if (hit->plusp == true) {
	hits_plus[hit->nmismatches_bothdiff][nhits_plus[hit->nmismatches_bothdiff]++] = hit;
      } else {
	hits_minus[hit->nmismatches_bothdiff][nhits_minus[hit->nmismatches_bothdiff]++] = hit;
      }
    }
  }

  return cutoff_level;
}


/* Finds concordant pairs if nconcordant is 0 */
List_T
Stage3_pair_up_concordant (bool *abort_pairing_p, int *found_score, int *nconcordant, int *nsamechr,
			   List_T *samechr, List_T *conc_transloc,
			   List_T hitpairs, List_T *hitarray5, int narray5, List_T *hitarray3, int narray3,
			   int cutoff_level_5, int cutoff_level_3, int subopt_levels,
			   Univcoord_T *splicesites,
			   Compress_T query5_compress_fwd, Compress_T query5_compress_rev,
			   Compress_T query3_compress_fwd, Compress_T query3_compress_rev,
			   int querylength5, int querylength3, int maxpairedpaths,
			   int splicing_penalty, int genestrand) {
  T **hits5_plus, **hits5_minus, **hits3_plus, **hits3_minus;
  int *nhits5_plus, *nhits5_minus, *nhits3_plus, *nhits3_minus;
  T **terminals5_plus, **terminals5_minus, **terminals3_plus, **terminals3_minus;
  int *nterminals5_plus, *nterminals5_minus, *nterminals3_plus, *nterminals3_minus;
  int score5, score3;
  bool *sorted_hits5_p, *sorted_hits3_p, *sorted_terminals5_p = NULL, *sorted_terminals3_p = NULL;
  int cutoff_level_hits5, cutoff_level_hits3, cutoff_level_terminals5 = 0, cutoff_level_terminals3 = 0;
  int ignore_found_score;
  
  int i;
  int min_score_5 = -1, min_score_3 = -1;
  List_T q;
  T hit;

  debug5(printf("Starting Stage3_pair_up_concordant with %d concordant, narray5 %d, narray3 %d, terminals5 %p, terminals3 %p, found_score %d\n",
		*nconcordant,narray5,narray3,terminals5,terminals3,*found_score));

  /* Re-adjust cutoff_levels for best hit5 and best hit3 */
  for (i = 0; i < narray5; i++) {
    for (q = hitarray5[i]; q != NULL; q = q->rest) {
      hit = (T) q->first;
      if (min_score_5 == -1 || hit->nmismatches_bothdiff < min_score_5) {
	min_score_5 = hit->nmismatches_bothdiff;
      }
    }
  }
  debug5(printf("min_score_5 = %d\n",min_score_5));
  if (min_score_5 > cutoff_level_5) {
    cutoff_level_5 = min_score_5;
  }

  for (i = 0; i < narray3; i++) {
    for (q = hitarray3[i]; q != NULL; q = q->rest) {
      hit = (T) q->first;
      if (min_score_3 == -1 || hit->nmismatches_bothdiff < min_score_3) {
	min_score_3 = hit->nmismatches_bothdiff;
      }
    }
  }
  debug5(printf("min_score_3 = %d\n",min_score_3));
  if (min_score_3 > cutoff_level_3) {
    cutoff_level_3 = min_score_3;
  }


  /* Find sizes for allocating memory */
  debug5(printf("Sizes of 5-end pieces by score level:\n"));
  hits5_plus = (T **) MALLOCA((cutoff_level_5+1) * sizeof(T *));
  hits5_minus = (T **) MALLOCA((cutoff_level_5+1) * sizeof(T *));
  nhits5_plus = (int *) CALLOCA(cutoff_level_5+1,sizeof(int));
  nhits5_minus = (int *) CALLOCA(cutoff_level_5+1,sizeof(int));
  cutoff_level_hits5 = sort_hits_by_trimmed_score(hits5_plus,hits5_minus,nhits5_plus,nhits5_minus,
						  hitarray5,narray5,cutoff_level_5);
  sorted_hits5_p = (bool *) CALLOCA(cutoff_level_hits5+1,sizeof(bool));


  debug5(printf("Sizes of 3-end pieces by score level:\n"));
  hits3_plus = (T **) MALLOCA((cutoff_level_3+1) * sizeof(T *));
  hits3_minus = (T **) MALLOCA((cutoff_level_3+1) * sizeof(T *));
  nhits3_plus = (int *) CALLOCA(cutoff_level_3+1,sizeof(int));
  nhits3_minus = (int *) CALLOCA(cutoff_level_3+1,sizeof(int));
  cutoff_level_hits3 = sort_hits_by_trimmed_score(hits3_plus,hits3_minus,nhits3_plus,nhits3_minus,
						  hitarray3,narray3,cutoff_level_3);
  sorted_hits3_p = (bool *) CALLOCA(cutoff_level_hits3+1,sizeof(bool));

  /* Look for concordant pairs among the non-terminals */
  hitpairs = pair_up_concordant_aux(&(*abort_pairing_p),&(*found_score),&(*nconcordant),&(*nsamechr),
				    &(*samechr),&(*conc_transloc),hitpairs,
				    hits5_plus,nhits5_plus,hits5_minus,nhits5_minus,
				    hits3_plus,nhits3_plus,hits3_minus,nhits3_minus,
				    /*sorted5p*/sorted_hits5_p,/*sorted3p*/sorted_hits3_p,
				    cutoff_level_hits5,cutoff_level_hits3,subopt_levels,
				    splicesites,query5_compress_fwd,query5_compress_rev,
				    query3_compress_fwd,query3_compress_rev,querylength5,querylength3,
				    maxpairedpaths,splicing_penalty,genestrand);

#if 0
    /* Look for single terminals */
    if (terminals3 != NULL) {
      terminals3_plus = (T **) MALLOCA((cutoff_level_3+1) * sizeof(T *));
      terminals3_minus = (T **) MALLOCA((cutoff_level_3+1) * sizeof(T *));
      nterminals3_plus = (int *) CALLOCA(cutoff_level_3+1,sizeof(int));
      nterminals3_minus = (int *) CALLOCA(cutoff_level_3+1,sizeof(int));
      cutoff_level_terminals3 = sort_hits_by_trimmed_score(terminals3_plus,terminals3_minus,nterminals3_plus,nterminals3_minus,
							   &terminals3,/*narray3*/1,cutoff_level_3);
      sorted_terminals3_p = (bool *) CALLOCA(cutoff_level_terminals3+1,sizeof(bool));

      /* Do not allow terminals to alter found_score */
      ignore_found_score = *found_score;
      *with_terminal = pair_up_concordant_aux(&(*abort_pairing_p),&ignore_found_score,&(*nconcordant),&(*nsamechr),
					      &(*samechr),&(*conc_transloc),*with_terminal,
					      hits5_plus,nhits5_plus,hits5_minus,nhits5_minus,
					      terminals3_plus,nterminals3_plus,terminals3_minus,nterminals3_minus,
					      /*sorted5p*/sorted_hits5_p,/*sorted3p*/sorted_terminals3_p,
					      cutoff_level_hits5,cutoff_level_terminals3,subopt_levels,
					      splicesites,query5_compress_fwd,query5_compress_rev,
					      query3_compress_fwd,query3_compress_rev,querylength5,querylength3,
					      maxpairedpaths,splicing_penalty,genestrand);
    }

    if (terminals5 != NULL) {
      terminals5_plus = (T **) MALLOCA((cutoff_level_5+1) * sizeof(T *));
      terminals5_minus = (T **) MALLOCA((cutoff_level_5+1) * sizeof(T *));
      nterminals5_plus = (int *) CALLOCA(cutoff_level_5+1,sizeof(int));
      nterminals5_minus = (int *) CALLOCA(cutoff_level_5+1,sizeof(int));
      cutoff_level_terminals5 = sort_hits_by_trimmed_score(terminals5_plus,terminals5_minus,nterminals5_plus,nterminals5_minus,
							   &terminals5,/*narray5*/1,cutoff_level_5);
      sorted_terminals5_p = (bool *) CALLOCA(cutoff_level_terminals5+1,sizeof(bool));

      /* Do not allow terminals to alter found_score */
      ignore_found_score = *found_score;
      *with_terminal = pair_up_concordant_aux(&(*abort_pairing_p),&ignore_found_score,&(*nconcordant),&(*nsamechr),
					      &(*samechr),&(*conc_transloc),*with_terminal,
					      terminals5_plus,nterminals5_plus,terminals5_minus,nterminals5_minus,
					      hits3_plus,nhits3_plus,hits3_minus,nhits3_minus,
					      /*sorted5p*/sorted_terminals5_p,/*sorted3p*/sorted_hits3_p,
					      cutoff_level_terminals5,cutoff_level_hits3,subopt_levels,
					      splicesites,query5_compress_fwd,query5_compress_rev,
					      query3_compress_fwd,query3_compress_rev,querylength5,querylength3,
					      maxpairedpaths,splicing_penalty,genestrand);
    }

    /* Previously required *with_terminal to be NULL also */
    if (terminals3 != NULL && terminals5 != NULL) {
      /* Look for double terminals */
      /* Do not allow terminals to alter found_score */
      ignore_found_score = *found_score;
      *with_terminal = pair_up_concordant_aux(&(*abort_pairing_p),&ignore_found_score,&(*nconcordant),&(*nsamechr),
					      &(*samechr),&(*conc_transloc),*with_terminal,
					      terminals5_plus,nterminals5_plus,terminals5_minus,nterminals5_minus,
					      terminals3_plus,nterminals3_plus,terminals3_minus,nterminals3_minus,
					      /*sorted5p*/sorted_terminals5_p,/*sorted3p*/sorted_terminals3_p,
					      cutoff_level_terminals5,cutoff_level_terminals3,subopt_levels,
					      splicesites,query5_compress_fwd,query5_compress_rev,
					      query3_compress_fwd,query3_compress_rev,querylength5,querylength3,
					      maxpairedpaths,splicing_penalty,genestrand);
    }

    if (sorted_terminals3_p != NULL) {
      for (score3 = 0; score3 <= cutoff_level_terminals3; score3++) {
	FREE(terminals3_plus[score3]);
	FREE(terminals3_minus[score3]);
      }
      FREEA(sorted_terminals3_p);
      FREEA(terminals3_plus);
      FREEA(terminals3_minus);
      FREEA(nterminals3_plus);
      FREEA(nterminals3_minus);
    }

    if (sorted_terminals5_p != NULL) {
      for (score5 = 0; score5 <= cutoff_level_terminals5; score5++) {
	FREE(terminals5_plus[score5]);
	FREE(terminals5_minus[score5]);
      }
      FREEA(sorted_terminals5_p);
      FREEA(terminals5_plus);
      FREEA(terminals5_minus);
      FREEA(nterminals5_plus);
      FREEA(nterminals5_minus);
    }
  }
#endif

  for (score3 = 0; score3 <= cutoff_level_hits3; score3++) {
    FREE(hits3_plus[score3]);
    FREE(hits3_minus[score3]);
  }
  FREEA(sorted_hits3_p);
  FREEA(hits3_plus);
  FREEA(hits3_minus);
  FREEA(nhits3_plus);
  FREEA(nhits3_minus);

  for (score5 = 0; score5 <= cutoff_level_hits5; score5++) {
    FREE(hits5_plus[score5]);
    FREE(hits5_minus[score5]);
  }
  FREEA(sorted_hits5_p);
  FREEA(hits5_plus);
  FREEA(hits5_minus);
  FREEA(nhits5_plus);
  FREEA(nhits5_minus);

  debug5(printf("Finished with Stage3_pair_up_concordant: %d concordant, %d samechr, %d conc_transloc\n",
		List_length(hitpairs),List_length(*samechr),List_length(*conc_transloc)));

  return hitpairs;
}


