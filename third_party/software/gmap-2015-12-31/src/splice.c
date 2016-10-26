static char rcsid[] = "$Id: splice.c 173900 2015-09-12 00:46:34Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "splice.h"

#include <stdio.h>
#include "mem.h"
#include "assert.h"
#include "sense.h"
#include "genome128_hr.h"
#include "genome_sites.h"
#include "substring.h"
#include "maxent.h"
#include "maxent_hr.h"
#include "stage3hr.h"


#define LOWPROB_SUPPORT 20

#if 0
/* Creates issues with ambiguous substrings */
#define LOCALSPLICING_NMATCHES_SLOP 1
#else
#define LOCALSPLICING_NMATCHES_SLOP 0
#endif
#define LOCALSPLICING_PROB_SLOP 0.05


/* Splice_solve_single */
#ifdef DEBUG1
#define debug1(x) x
#else
#define debug1(x)
#endif

/* Splice_solve_double */
#ifdef DEBUG2
#define debug2(x) x
#else
#define debug2(x)
#endif

/* Group by segmentj */
#ifdef DEBUG7
#define debug7(x) x
#else
#define debug7(x)
#endif




static bool novelsplicingp = true;
static int min_shortend;

void
Splice_setup (int min_shortend_in) {
  min_shortend = min_shortend_in;
  return;
}



/* Do not compare against true or false */
/* Loosest criterion */
static int
sufficient_splice_prob_local (int support, int nmismatches, double spliceprob) {
  support -= 3*nmismatches;
  if (support < 14) {
    return (spliceprob > 0.95);
  } else if (support < 20) {
    return (spliceprob > 0.90);
  } else if (support < 26) {
    return (spliceprob > 0.85);
  } else {
    return (spliceprob > 0.70);
  }
}



/* Note: knowni holds joffset + j + 1, so 0 represents no known site
   and values greater than 0 represent a known site.  Need to subtract
   1 to obtain joffset + j. */

/* Called only by sarray-read.c, where plusp is always true */
int
Splice_resolve_sense (int *best_knowni_i, int *best_knowni_j,
		      int *best_nmismatches_i, int *best_nmismatches_j,
		      double *best_prob_i, double *best_prob_j,

		      Univcoord_T segmenti_left, Univcoord_T segmentj_left,
		      Univcoord_T segmenti_chroffset, Univcoord_T segmentj_chroffset,
		     
		      int querystart, int queryend, int querylength, Compress_T query_compress,
		      int *segmenti_donor_knownpos, int *segmentj_acceptor_knownpos,
		      int *segmentj_antidonor_knownpos, int *segmenti_antiacceptor_knownpos,
		      int *segmenti_donor_knowni, int *segmentj_acceptor_knowni,
		      int *segmentj_antidonor_knowni, int *segmenti_antiacceptor_knowni,
		      int segmenti_donor_nknown, int segmentj_acceptor_nknown,
		      int segmentj_antidonor_nknown, int segmenti_antiacceptor_nknown,
		      int splicing_penalty, int max_mismatches_allowed,
		      bool plusp, int genestrand, bool first_read_p) {
  int best_splice_pos = -1, splice_pos_start, splice_pos_end, splice_pos, i, j;

  int best_nmismatches, nmismatches;
  int best_segmenti_nmismatches, best_segmentj_nmismatches, segmenti_nmismatches, segmentj_nmismatches;
  Univcoord_T best_donor_splicecoord, best_acceptor_splicecoord;
  int best_donor_knowni, best_acceptor_knowni;
  double best_prob, best_donor_prob, best_acceptor_prob, probi, probj;
  /* bool sufficient1p, sufficient2p; */

  int donori_nsites, acceptorj_nsites, antiacceptori_nsites, antidonorj_nsites;
  int *donori_positions, *acceptorj_positions, *antiacceptori_positions, *antidonorj_positions;
  int *donori_knowni, *acceptorj_knowni, *antiacceptori_knowni, *antidonorj_knowni;

#ifdef HAVE_ALLOCA
  int *donor_positions_alloc = (int *) alloca((querylength+1)*sizeof(int));
  int *acceptor_positions_alloc = (int *) alloca((querylength+1)*sizeof(int));
  int *donor_knowni_alloc = (int *) alloca((querylength+1)*sizeof(int));
  int *acceptor_knowni_alloc = (int *) alloca((querylength+1)*sizeof(int));
#else
  int donor_positions_alloc[MAX_READLENGTH+1], acceptor_positions_alloc[MAX_READLENGTH+1];
  int donor_knowni_alloc[MAX_READLENGTH+1], acceptor_knowni_alloc[MAX_READLENGTH+1];
#endif


  debug1(printf("Splice_resolve_sense: Getting genome at lefti %u and leftj %u (diff: %d), range %d..%d\n",
		segmenti_left,segmentj_left,segmentj_left-segmenti_left,querystart,queryend));

  *best_knowni_i = *best_knowni_j = -1;
  *best_nmismatches_i = *best_nmismatches_j = 0;
  *best_prob_i = *best_prob_j = 0.0;

  splice_pos_start = querystart;
  splice_pos_end = queryend;

  if (plusp == true) {
    /* Originally from plus strand.  No complement.  */
    /* Sense (End 1 to End 2) or Antisense (End 5 to End 6) */
    if (novelsplicingp && segmenti_left + splice_pos_start >= DONOR_MODEL_LEFT_MARGIN) {
      donori_nsites = Genome_donor_positions(donor_positions_alloc,donor_knowni_alloc,
					     segmenti_donor_knownpos,segmenti_donor_knowni,
					     segmenti_left,splice_pos_start,splice_pos_end);
      donori_positions = donor_positions_alloc;
      donori_knowni = donor_knowni_alloc;
    } else {
      donori_nsites = segmenti_donor_nknown;
      donori_positions = segmenti_donor_knownpos;
      donori_knowni = segmenti_donor_knowni;
    }

#ifdef DEBUG1
    printf("Found %d donori sites:",donori_nsites);
    for (i = 0; i < donori_nsites; i++) {
      printf(" %d",donori_positions[i]);
      if (donori_knowni[i] >= 0) {
	printf(" (%d)",donori_knowni[i]);
      }
    }
    printf("\n");
#endif

    if (novelsplicingp && segmentj_left + splice_pos_start >= ACCEPTOR_MODEL_LEFT_MARGIN) {
      acceptorj_nsites = Genome_acceptor_positions(acceptor_positions_alloc,acceptor_knowni_alloc,
						   segmentj_acceptor_knownpos,segmentj_acceptor_knowni,
						   segmentj_left,splice_pos_start,splice_pos_end);
      acceptorj_positions = acceptor_positions_alloc;
      acceptorj_knowni = acceptor_knowni_alloc;
    } else {
      acceptorj_nsites = segmentj_acceptor_nknown;
      acceptorj_positions = segmentj_acceptor_knownpos;
      acceptorj_knowni = segmentj_acceptor_knowni;
    }

#ifdef DEBUG1
    printf("Found %d acceptorj sites:",acceptorj_nsites);
    for (i = 0; i < acceptorj_nsites; i++) {
      printf(" %d",acceptorj_positions[i]);
      if (acceptorj_knowni[i] >= 0) {
	printf(" (%d)",acceptorj_knowni[i]);
      }
    }
    printf("\n");
#endif

    best_nmismatches = max_mismatches_allowed;
    best_prob = 0.0;

    i = j = 0;
    while (i < donori_nsites && j < acceptorj_nsites) {
      if ((splice_pos = donori_positions[i]) < acceptorj_positions[j]) {
	i++;
      } else if (splice_pos > acceptorj_positions[j]) {
	j++;
      } else {
	segmenti_nmismatches = Genome_count_mismatches_substring(query_compress,/*left*/segmenti_left,/*pos5*/querystart,/*pos3*/splice_pos,
								 plusp,genestrand,first_read_p);
	segmentj_nmismatches = Genome_count_mismatches_substring(query_compress,/*left*/segmentj_left,/*pos5*/splice_pos,/*pos3*/queryend,
								 plusp,genestrand,first_read_p);
	if ((nmismatches = segmenti_nmismatches + segmentj_nmismatches) <= best_nmismatches) {
	  if (donori_knowni[i] >= 0) {
	    probi = 1.0; /* Needs to be 1.0 for output */
	  } else {
	    probi = Maxent_hr_donor_prob(segmenti_left + splice_pos,segmenti_chroffset);
	  }

	  if (acceptorj_knowni[j] >= 0) {
	    probj = 1.0; /* Needs to be 1.0 for output */
	  } else {
	    probj = Maxent_hr_acceptor_prob(segmentj_left + splice_pos,segmentj_chroffset);
	  }

	  debug1(
		 if (plusp == true) {
		   printf("plus sense splice_pos  %d, i.donor %f, j.acceptor %f\n",splice_pos,probi,probj);
		 } else {
		   printf("minus antisense splice_pos  %d, i.donor %f, j.acceptor %f\n",splice_pos,probi,probj);
		 });


#if 0
	  sufficient1p = sufficient_splice_prob_local(/*support*/splice_pos,segmenti_nmismatches,probi);
	  sufficient2p = sufficient_splice_prob_local(/*support*/querylength - splice_pos,segmentj_nmismatches,probj);
#endif

	  /* if (sufficient1p && sufficient2p) { */
	    if (nmismatches < best_nmismatches ||
		(nmismatches == best_nmismatches && probi + probj > best_prob)) {
	      /* Success */
	      best_nmismatches = nmismatches;
	      best_prob = probi + probj;

	      /* best_donor_splicecoord = segmenti_left + splice_pos; */
	      /* best_acceptor_splicecoord = segmentj_left + splice_pos; */
	      *best_knowni_i = donori_knowni[i];
	      *best_knowni_j = acceptorj_knowni[j];
	      *best_prob_i = probi; /* donor_prob */
	      *best_prob_j = probj; /* acceptor_prob */
	      best_splice_pos = splice_pos;
	      *best_nmismatches_i = segmenti_nmismatches;
	      *best_nmismatches_j = segmentj_nmismatches;
	    }
	    /* 	} */
	}
	i++;
	j++;
      }
    }

  } else {
    /* minus */
    /* Originally from minus strand.  Complement. */
    /* Antisense (End 7 to End 8) or Sense (End 3 to End 4) */
    if (novelsplicingp && segmenti_left + splice_pos_start >= ACCEPTOR_MODEL_RIGHT_MARGIN) {
      antiacceptori_nsites = Genome_antiacceptor_positions(acceptor_positions_alloc,acceptor_knowni_alloc,
							   segmenti_antiacceptor_knownpos,segmenti_antiacceptor_knowni,
							   segmenti_left,splice_pos_start,splice_pos_end);
      antiacceptori_positions = acceptor_positions_alloc;
      antiacceptori_knowni = acceptor_knowni_alloc;
    } else {
      antiacceptori_nsites = segmenti_antiacceptor_nknown;
      antiacceptori_positions = segmenti_antiacceptor_knownpos;
      antiacceptori_knowni = segmenti_antiacceptor_knowni;
    }

#ifdef DEBUG1
    printf("Found %d antiacceptori sites:",antiacceptori_nsites);
    for (i = 0; i < antiacceptori_nsites; i++) {
      printf(" %d",antiacceptori_positions[i]);
      if (antiacceptori_knowni[i] >= 0) {
	printf(" (%d)",antiacceptori_knowni[i]);
      }
    }
    printf("\n");
#endif

    if (novelsplicingp && segmentj_left + splice_pos_start >= DONOR_MODEL_RIGHT_MARGIN) {
      antidonorj_nsites = Genome_antidonor_positions(donor_positions_alloc,donor_knowni_alloc,
						     segmentj_antidonor_knownpos,segmentj_antidonor_knowni,
						     segmentj_left,splice_pos_start,splice_pos_end);
      antidonorj_positions = donor_positions_alloc;
      antidonorj_knowni = donor_knowni_alloc;
    } else {
      antidonorj_nsites = segmentj_antidonor_nknown;
      antidonorj_positions = segmentj_antidonor_knownpos;
      antidonorj_knowni = segmentj_antidonor_knowni;
    }

#ifdef DEBUG1
    printf("Found %d antidonorj sites:",antidonorj_nsites);
    for (i = 0; i < antidonorj_nsites; i++) {
      printf(" %d",antidonorj_positions[i]);
      if (antidonorj_knowni[i] >= 0) {
	printf(" (%d)",antidonorj_knowni[i]);
      }
    }
    printf("\n");
#endif

    best_nmismatches = max_mismatches_allowed;
    best_prob = 0.0;

    i = j = 0;
    while (i < antiacceptori_nsites && j < antidonorj_nsites) {
      if ((splice_pos = antiacceptori_positions[i]) < antidonorj_positions[j]) {
	i++;
      } else if (splice_pos > antidonorj_positions[j]) {
	j++;
      } else {
	segmenti_nmismatches = Genome_count_mismatches_substring(query_compress,/*left*/segmenti_left,/*pos5*/querystart,/*pos3*/splice_pos,
								 plusp,genestrand,first_read_p);
	segmentj_nmismatches = Genome_count_mismatches_substring(query_compress,/*left*/segmentj_left,/*pos5*/splice_pos,/*pos3*/queryend,
								 plusp,genestrand,first_read_p);
	if ((nmismatches = segmenti_nmismatches + segmentj_nmismatches) <= best_nmismatches) {
	  if (antiacceptori_knowni[i] >= 0) {
	    probi = 1.0; /* Needs to be 1.0 for output */
	  } else {
	    probi = Maxent_hr_antiacceptor_prob(segmenti_left + splice_pos,segmenti_chroffset);
	  }

	  if (antidonorj_knowni[j] >= 0) {
	    probj = 1.0; /* Needs to be 1.0 for output */
	  } else {
	    probj = Maxent_hr_antidonor_prob(segmentj_left + splice_pos,segmentj_chroffset);
	  }

	  debug1(
		 if (plusp == true) {
		   printf("plus antisense splice_pos  %d, j.donor %f, i.acceptor %f\n",splice_pos,probj,probi);
		 } else {
		   printf("minus sense splice_pos  %d, j.donor %f, i.acceptor %f\n",splice_pos,probj,probi);
		 });
	  
#if 0
	  sufficient1p = sufficient_splice_prob_local(/*support*/splice_pos,segmenti_nmismatches,probi);
	  sufficient2p = sufficient_splice_prob_local(/*support*/querylength - splice_pos,segmentj_nmismatches,probj);
#endif

	  /* if (sufficient1p && sufficient2p) { */
	    if (nmismatches < best_nmismatches ||
		(nmismatches == best_nmismatches && probi + probj > best_prob)) {
	      /* Success */
	      best_nmismatches = nmismatches;
	      best_prob = probi + probj;
	      
	      /* best_donor_splicecoord = segmentj_left + splice_pos; */
	      /* best_acceptor_splicecoord = segmenti_left + splice_pos; */
	      *best_knowni_j = antidonorj_knowni[j];
	      *best_knowni_i = antiacceptori_knowni[i];
	      *best_prob_j = probj; /* donor_prob */
	      *best_prob_i = probi;
	      best_splice_pos = splice_pos;
	      *best_nmismatches_j = segmentj_nmismatches;
	      *best_nmismatches_i = segmenti_nmismatches;
	    }
	    /* } */
	}
	i++;
	j++;
      }
    }
  }

  if (*best_prob_i > 0.95 && *best_prob_j > 0.70) {
    debug1(printf("Returning %d with probi %f and probj %f\n",best_splice_pos,*best_prob_i,*best_prob_j));
    return best_splice_pos;
  } else if (*best_prob_i > 0.70 && *best_prob_j > 0.95) {
    debug1(printf("Returning %d with probi %f and probj %f\n",best_splice_pos,*best_prob_i,*best_prob_j));
    return best_splice_pos;
  } else if (*best_prob_i > 0.80 && *best_prob_j > 0.85) {
    debug1(printf("Returning %d with probi %f and probj %f\n",best_splice_pos,*best_prob_i,*best_prob_j));
    return best_splice_pos;
  } else {
    debug1(printf("Not returning %d with probi %f and probj %f\n",best_splice_pos,*best_prob_i,*best_prob_j));
    return -1;
  }
}


/* Called only by sarray-read.c, where plusp is always true */
int
Splice_resolve_antisense (int *best_knowni_i, int *best_knowni_j,
			  int *best_nmismatches_i, int *best_nmismatches_j,
			  double *best_prob_i, double *best_prob_j,
			  
			  Univcoord_T segmenti_left, Univcoord_T segmentj_left,
			  Univcoord_T segmenti_chroffset, Univcoord_T segmentj_chroffset,
			  
			  int querystart, int queryend, int querylength, Compress_T query_compress,
			  int *segmenti_donor_knownpos, int *segmentj_acceptor_knownpos,
			  int *segmentj_antidonor_knownpos, int *segmenti_antiacceptor_knownpos,
			  int *segmenti_donor_knowni, int *segmentj_acceptor_knowni,
			  int *segmentj_antidonor_knowni, int *segmenti_antiacceptor_knowni,
			  int segmenti_donor_nknown, int segmentj_acceptor_nknown,
			  int segmentj_antidonor_nknown, int segmenti_antiacceptor_nknown,
			  int splicing_penalty, int max_mismatches_allowed,
			  bool plusp, int genestrand, bool first_read_p) {
  int best_splice_pos = -1, splice_pos_start, splice_pos_end, splice_pos, i, j;

  int best_nmismatches, nmismatches;
  int best_segmenti_nmismatches, best_segmentj_nmismatches, segmenti_nmismatches, segmentj_nmismatches;
  Univcoord_T best_donor_splicecoord, best_acceptor_splicecoord;
  int best_donor_knowni, best_acceptor_knowni;
  double best_prob, best_donor_prob, best_acceptor_prob, probi, probj;
  /* bool sufficient1p, sufficient2p; */

  int donori_nsites, acceptorj_nsites, antiacceptori_nsites, antidonorj_nsites;
  int *donori_positions, *acceptorj_positions, *antiacceptori_positions, *antidonorj_positions;
  int *donori_knowni, *acceptorj_knowni, *antiacceptori_knowni, *antidonorj_knowni;

#ifdef HAVE_ALLOCA
  int *donor_positions_alloc = (int *) alloca((querylength+1)*sizeof(int));
  int *acceptor_positions_alloc = (int *) alloca((querylength+1)*sizeof(int));
  int *donor_knowni_alloc = (int *) alloca((querylength+1)*sizeof(int));
  int *acceptor_knowni_alloc = (int *) alloca((querylength+1)*sizeof(int));
#else
  int donor_positions_alloc[MAX_READLENGTH+1], acceptor_positions_alloc[MAX_READLENGTH+1];
  int donor_knowni_alloc[MAX_READLENGTH+1], acceptor_knowni_alloc[MAX_READLENGTH+1];
#endif

  debug1(printf("Splice_resolve_antisense: Getting genome at lefti %u and leftj %u (diff: %d), range %d..%d\n",
		segmenti_left,segmentj_left,segmentj_left-segmenti_left,querystart,queryend));

  *best_knowni_i = *best_knowni_j = -1;
  *best_nmismatches_i = *best_nmismatches_j = 0;
  *best_prob_i = *best_prob_j = 0.0;

  splice_pos_start = querystart;
  splice_pos_end = queryend;

  if (plusp == false) {
    /* minus */
    /* Originally from plus strand.  No complement.  */
    /* Sense (End 1 to End 2) or Antisense (End 5 to End 6) */
    if (novelsplicingp && segmenti_left + splice_pos_start >= DONOR_MODEL_LEFT_MARGIN) {
      donori_nsites = Genome_donor_positions(donor_positions_alloc,donor_knowni_alloc,
					     segmenti_donor_knownpos,segmenti_donor_knowni,
					     segmenti_left,splice_pos_start,splice_pos_end);
      donori_positions = donor_positions_alloc;
      donori_knowni = donor_knowni_alloc;
    } else {
      donori_nsites = segmenti_donor_nknown;
      donori_positions = segmenti_donor_knownpos;
      donori_knowni = segmenti_donor_knowni;
    }

#ifdef DEBUG1
    printf("Found %d donori sites:",donori_nsites);
    for (i = 0; i < donori_nsites; i++) {
      printf(" %d",donori_positions[i]);
      if (donori_knowni[i] >= 0) {
	printf(" (%d)",donori_knowni[i]);
      }
    }
    printf("\n");
#endif

    if (novelsplicingp && segmentj_left + splice_pos_start >= ACCEPTOR_MODEL_LEFT_MARGIN) {
      acceptorj_nsites = Genome_acceptor_positions(acceptor_positions_alloc,acceptor_knowni_alloc,
						   segmentj_acceptor_knownpos,segmentj_acceptor_knowni,
						   segmentj_left,splice_pos_start,splice_pos_end);
      acceptorj_positions = acceptor_positions_alloc;
      acceptorj_knowni = acceptor_knowni_alloc;
    } else {
      acceptorj_nsites = segmentj_acceptor_nknown;
      acceptorj_positions = segmentj_acceptor_knownpos;
      acceptorj_knowni = segmentj_acceptor_knowni;
    }

#ifdef DEBUG1
    printf("Found %d acceptorj sites:",acceptorj_nsites);
    for (i = 0; i < acceptorj_nsites; i++) {
      printf(" %d",acceptorj_positions[i]);
      if (acceptorj_knowni[i] >= 0) {
	printf(" (%d)",acceptorj_knowni[i]);
      }
    }
    printf("\n");
#endif

    best_nmismatches = max_mismatches_allowed;
    best_prob = 0.0;

    i = j = 0;
    while (i < donori_nsites && j < acceptorj_nsites) {
      if ((splice_pos = donori_positions[i]) < acceptorj_positions[j]) {
	i++;
      } else if (splice_pos > acceptorj_positions[j]) {
	j++;
      } else {
	segmenti_nmismatches = Genome_count_mismatches_substring(query_compress,/*left*/segmenti_left,/*pos5*/querystart,/*pos3*/splice_pos,
								 plusp,genestrand,first_read_p);
	segmentj_nmismatches = Genome_count_mismatches_substring(query_compress,/*left*/segmentj_left,/*pos5*/splice_pos,/*pos3*/queryend,
								 plusp,genestrand,first_read_p);
	if ((nmismatches = segmenti_nmismatches + segmentj_nmismatches) <= best_nmismatches) {
	  if (donori_knowni[i] >= 0) {
	    probi = 1.0; /* Needs to be 1.0 for output */
	  } else {
	    probi = Maxent_hr_donor_prob(segmenti_left + splice_pos,segmenti_chroffset);
	  }
	  
	  if (acceptorj_knowni[j] >= 0) {
	    probj = 1.0; /* Needs to be 1.0 for output */
	  } else {
	    probj = Maxent_hr_acceptor_prob(segmentj_left + splice_pos,segmentj_chroffset);
	  }

	  debug1(
		 if (plusp == true) {
		   printf("plus sense splice_pos  %d, i.donor %f, j.acceptor %f\n",splice_pos,probi,probj);
		 } else {
		   printf("minus antisense splice_pos  %d, i.donor %f, j.acceptor %f\n",splice_pos,probi,probj);
		 });

#if 0
	  sufficient1p = sufficient_splice_prob_local(/*support*/splice_pos,segmenti_nmismatches,probi);
	  sufficient2p = sufficient_splice_prob_local(/*support*/querylength - splice_pos,segmentj_nmismatches,probj);
#endif

	  /* if (sufficient1p && sufficient2p) { */
	    if (nmismatches < best_nmismatches ||
		(nmismatches == best_nmismatches && probi + probj > best_prob)) {
	      /* Success */
	      best_nmismatches = nmismatches;
	      best_prob = probi + probj;
	      
	      /* best_donor_splicecoord = segmenti_left + splice_pos; */
	      /* best_acceptor_splicecoord = segmentj_left + splice_pos; */
	      *best_knowni_i = donori_knowni[i];
	      *best_knowni_j = acceptorj_knowni[j];
	      *best_prob_i = probi; /* donor_prob */
	      *best_prob_j = probj; /* acceptor_prob */
	      best_splice_pos = splice_pos;
	      *best_nmismatches_i = segmenti_nmismatches;
	      *best_nmismatches_j = segmentj_nmismatches;
	    }
	    /* } */
	}
	i++;
	j++;
      }
    }

  } else {
    /* plus */
    /* Originally from minus strand.  Complement. */
    /* Antisense (End 7 to End 8) or Sense (End 3 to End 4) */
    if (novelsplicingp && segmenti_left + splice_pos_start >= ACCEPTOR_MODEL_RIGHT_MARGIN) {
      antiacceptori_nsites = Genome_antiacceptor_positions(acceptor_positions_alloc,acceptor_knowni_alloc,
							   segmenti_antiacceptor_knownpos,segmenti_antiacceptor_knowni,
							   segmenti_left,splice_pos_start,splice_pos_end);
      antiacceptori_positions = acceptor_positions_alloc;
      antiacceptori_knowni = acceptor_knowni_alloc;
    } else {
      antiacceptori_nsites = segmenti_antiacceptor_nknown;
      antiacceptori_positions = segmenti_antiacceptor_knownpos;
      antiacceptori_knowni = segmenti_antiacceptor_knowni;
    }

#ifdef DEBUG1
    printf("Found %d antiacceptori sites:",antiacceptori_nsites);
    for (i = 0; i < antiacceptori_nsites; i++) {
      printf(" %d",antiacceptori_positions[i]);
      if (antiacceptori_knowni[i] >= 0) {
	printf(" (%d)",antiacceptori_knowni[i]);
      }
    }
    printf("\n");
#endif

    if (novelsplicingp && segmentj_left + splice_pos_start >= DONOR_MODEL_RIGHT_MARGIN) {
      antidonorj_nsites = Genome_antidonor_positions(donor_positions_alloc,donor_knowni_alloc,
						     segmentj_antidonor_knownpos,segmentj_antidonor_knowni,
						     segmentj_left,splice_pos_start,splice_pos_end);
      antidonorj_positions = donor_positions_alloc;
      antidonorj_knowni = donor_knowni_alloc;
    } else {
      antidonorj_nsites = segmentj_antidonor_nknown;
      antidonorj_positions = segmentj_antidonor_knownpos;
      antidonorj_knowni = segmentj_antidonor_knowni;
    }

#ifdef DEBUG1
    printf("Found %d antidonorj sites:",antidonorj_nsites);
    for (i = 0; i < antidonorj_nsites; i++) {
      printf(" %d",antidonorj_positions[i]);
      if (antidonorj_knowni[i] >= 0) {
	printf(" (%d)",antidonorj_knowni[i]);
      }
    }
    printf("\n");
#endif

    best_nmismatches = max_mismatches_allowed;
    best_prob = 0.0;

    i = j = 0;
    while (i < antiacceptori_nsites && j < antidonorj_nsites) {
      if ((splice_pos = antiacceptori_positions[i]) < antidonorj_positions[j]) {
	i++;
      } else if (splice_pos > antidonorj_positions[j]) {
	j++;
      } else {
	segmenti_nmismatches = Genome_count_mismatches_substring(query_compress,/*left*/segmenti_left,/*pos5*/querystart,/*pos3*/splice_pos,
								 plusp,genestrand,first_read_p);
	segmentj_nmismatches = Genome_count_mismatches_substring(query_compress,/*left*/segmentj_left,/*pos5*/splice_pos,/*pos3*/queryend,
								 plusp,genestrand,first_read_p);
	if ((nmismatches = segmenti_nmismatches + segmentj_nmismatches) <= best_nmismatches) {
	  if (antiacceptori_knowni[i] >= 0) {
	    probi = 1.0; /* Needs to be 1.0 for output */
	  } else {
	    probi = Maxent_hr_antiacceptor_prob(segmenti_left + splice_pos,segmenti_chroffset);
	  }

	  if (antidonorj_knowni[j] >= 0) {
	    probj = 1.0; /* Needs to be 1.0 for output */
	  } else {
	    probj = Maxent_hr_antidonor_prob(segmentj_left + splice_pos,segmentj_chroffset);
	  }

	  debug1(
		 if (plusp == true) {
		   printf("plus antisense splice_pos  %d, j.donor %f, i.acceptor %f\n",splice_pos,probj,probi);
		 } else {
		   printf("minus sense splice_pos  %d, j.donor %f, i.acceptor %f\n",splice_pos,probj,probi);
		 });
	  
#if 0
	  sufficient1p = sufficient_splice_prob_local(/*support*/splice_pos,segmenti_nmismatches,probi);
	  sufficient2p = sufficient_splice_prob_local(/*support*/querylength - splice_pos,segmentj_nmismatches,probj);
#endif

	  /* if (sufficient1p && sufficient2p) { */
	    if (nmismatches < best_nmismatches ||
		(nmismatches == best_nmismatches && probi + probj > best_prob)) {
	      /* Success */
	      best_nmismatches = nmismatches;
	      best_prob = probi + probj;
	      
	      /* best_donor_splicecoord = segmentj_left + splice_pos; */
	      /* best_acceptor_splicecoord = segmenti_left + splice_pos; */
	      *best_knowni_j = antidonorj_knowni[j];
	      *best_knowni_i = antiacceptori_knowni[i];
	      *best_prob_j = probj; /* donor_prob */
	      *best_prob_i = probi; /* acceptor_prob */
	      best_splice_pos = splice_pos;
	      *best_nmismatches_j = segmentj_nmismatches;
	      *best_nmismatches_i = segmenti_nmismatches;
	    }
	    /* } */
	}
	i++;
	j++;
      }
    }
  }

  if (*best_prob_i > 0.95 && *best_prob_j > 0.70) {
    debug1(printf("Returning %d with probi %f and probj %f\n",best_splice_pos,*best_prob_i,*best_prob_j));
    return best_splice_pos;
  } else if (*best_prob_i > 0.70 && *best_prob_j > 0.95) {
    debug1(printf("Returning %d with probi %f and probj %f\n",best_splice_pos,*best_prob_i,*best_prob_j));
    return best_splice_pos;
  } else if (*best_prob_i > 0.80 && *best_prob_j > 0.85) {
    debug1(printf("Returning %d with probi %f and probj %f\n",best_splice_pos,*best_prob_i,*best_prob_j));
    return best_splice_pos;
  } else {
    debug1(printf("Not returning %d with probi %f and probj %f\n",best_splice_pos,*best_prob_i,*best_prob_j));
    return -1;
  }
}



/* Note: knowni holds joffset + j + 1, so 0 represents no known site
   and values greater than 0 represent a known site.  Need to subtract
   1 to obtain joffset + j. */

List_T
Splice_solve_single_sense (int *found_score, int *nhits, List_T hits, List_T *lowprob,

			   bool *segmenti_usedp, bool *segmentj_usedp,
			   Univcoord_T segmenti_left, Univcoord_T segmentj_left,
			   Chrnum_T segmenti_chrnum, Univcoord_T segmenti_chroffset,
			   Univcoord_T segmenti_chrhigh, Chrpos_T segmenti_chrlength,
			   Chrnum_T segmentj_chrnum, Univcoord_T segmentj_chroffset,
			   Univcoord_T segmentj_chrhigh, Chrpos_T segmentj_chrlength,
		     
			   int querylength, Compress_T query_compress,
			   int *segmenti_donor_knownpos, int *segmentj_acceptor_knownpos,
			   int *segmentj_antidonor_knownpos, int *segmenti_antiacceptor_knownpos,
			   int *segmenti_donor_knowni, int *segmentj_acceptor_knowni,
			   int *segmentj_antidonor_knowni, int *segmenti_antiacceptor_knowni,
			   int segmenti_donor_nknown, int segmentj_acceptor_nknown,
			   int segmentj_antidonor_nknown, int segmenti_antiacceptor_nknown,
			   int splicing_penalty, int max_mismatches_allowed,
			   bool plusp, int genestrand, bool first_read_p,
			   bool subs_or_indels_p, bool sarrayp) {
  Substring_T donor, acceptor;
  int best_splice_pos, splice_pos_start, splice_pos_end, splice_pos, i, j;

  int best_nmismatches, nmismatches;
  int best_segmenti_nmismatches, best_segmentj_nmismatches, segmenti_nmismatches, segmentj_nmismatches;
  int donor_support, acceptor_support;
  Univcoord_T best_donor_splicecoord, best_acceptor_splicecoord;
  int best_donor_knowni, best_acceptor_knowni;
  double best_prob, best_donor_prob, best_acceptor_prob, probi, probj;
  bool sufficient1p, sufficient2p, orig_plusp;
  int sensedir;

  int donori_nsites, acceptorj_nsites, antiacceptori_nsites, antidonorj_nsites;
  int *donori_positions, *acceptorj_positions, *antiacceptori_positions, *antidonorj_positions;
  int *donori_knowni, *acceptorj_knowni, *antiacceptori_knowni, *antidonorj_knowni;

#ifdef HAVE_ALLOCA
  int *donor_positions_alloc = (int *) alloca((querylength+1)*sizeof(int));
  int *acceptor_positions_alloc = (int *) alloca((querylength+1)*sizeof(int));
  int *donor_knowni_alloc = (int *) alloca((querylength+1)*sizeof(int));
  int *acceptor_knowni_alloc = (int *) alloca((querylength+1)*sizeof(int));
#else
  int donor_positions_alloc[MAX_READLENGTH+1], acceptor_positions_alloc[MAX_READLENGTH+1];
  int donor_knowni_alloc[MAX_READLENGTH+1], acceptor_knowni_alloc[MAX_READLENGTH+1];
#endif


  debug1(printf("Splice_solve_single: Getting genome at lefti %u and leftj %u (diff: %d)\n",
		segmenti_left,segmentj_left,segmentj_left-segmenti_left));
  *nhits = 0;

#if 0
  int sum, lefti, righti;
  splice_pos_start = querylength;
  splice_pos_end = 0;
  for (sum = 0; sum <= max_mismatches_allowed; sum++) {
    for (lefti = 0; lefti <= sum && lefti < nmismatches_left; lefti++) {
      if ((righti = sum - lefti) < nmismatches_right &&
	  mismatch_positions_left[lefti] > mismatch_positions_right[righti]) {
	debug1(printf("At %d+%d mismatches, splice_pos using right: %d\n",lefti,righti,mismatch_positions_right[righti]+1));
	debug1(printf("At %d+%d mismatches, splice_pos using left: %d\n",lefti,righti,mismatch_positions_left[lefti]));
	if (mismatch_positions_right[righti] + 1 < splice_pos_start) {
	  splice_pos_start = mismatch_positions_right[righti] + 1;	/* This is leftmost position in righti+1 .. lefti */
	}
	if (mismatch_positions_left[lefti] > splice_pos_end) {
	  splice_pos_end = mismatch_positions_left[lefti];	/* This is rightmost position in righti+1 .. lefti */
	}
      }
    }
  }

  /* Exclude ends */
  if (splice_pos_start < min_localsplicing_end_matches) {
    splice_pos_start = min_localsplicing_end_matches;
  }
  if (splice_pos_end > querylength - min_localsplicing_end_matches) {
    splice_pos_end = querylength - min_localsplicing_end_matches;
  }
#else
  /* splice_pos_start = min_localsplicing_end_matches; */
  /* splice_pos_end = querylength - min_localsplicing_end_matches; */
  splice_pos_start = min_shortend;
  splice_pos_end = querylength - min_shortend; /* ? off by 1, so -l 3 allows only ends of up to 2 */
#endif


  if (splice_pos_start <= splice_pos_end) {
    if (plusp == true) {
      /* Originally from plus strand.  No complement.  */
      /* Sense (End 1 to End 2) or Antisense (End 5 to End 6) */
      if (novelsplicingp && segmenti_left + splice_pos_start >= DONOR_MODEL_LEFT_MARGIN) {
	donori_nsites = Genome_donor_positions(donor_positions_alloc,donor_knowni_alloc,
					       segmenti_donor_knownpos,segmenti_donor_knowni,
					       segmenti_left,splice_pos_start,splice_pos_end);
	donori_positions = donor_positions_alloc;
	donori_knowni = donor_knowni_alloc;
      } else {
	donori_nsites = segmenti_donor_nknown;
	donori_positions = segmenti_donor_knownpos;
	donori_knowni = segmenti_donor_knowni;
      }

#ifdef DEBUG1
      printf("Found %d donori sites:",donori_nsites);
      for (i = 0; i < donori_nsites; i++) {
	printf(" %d",donori_positions[i]);
	if (donori_knowni[i] >= 0) {
	  printf(" (%d)",donori_knowni[i]);
	}
      }
      printf("\n");
#endif

      if (novelsplicingp && segmentj_left + splice_pos_start >= ACCEPTOR_MODEL_LEFT_MARGIN) {
	acceptorj_nsites = Genome_acceptor_positions(acceptor_positions_alloc,acceptor_knowni_alloc,
						     segmentj_acceptor_knownpos,segmentj_acceptor_knowni,
						     segmentj_left,splice_pos_start,splice_pos_end);
	acceptorj_positions = acceptor_positions_alloc;
	acceptorj_knowni = acceptor_knowni_alloc;
      } else {
	acceptorj_nsites = segmentj_acceptor_nknown;
	acceptorj_positions = segmentj_acceptor_knownpos;
	acceptorj_knowni = segmentj_acceptor_knowni;
      }

#ifdef DEBUG1
      printf("Found %d acceptorj sites:",acceptorj_nsites);
      for (i = 0; i < acceptorj_nsites; i++) {
	printf(" %d",acceptorj_positions[i]);
	if (acceptorj_knowni[i] >= 0) {
	  printf(" (%d)",acceptorj_knowni[i]);
	}
      }
      printf("\n");
#endif

      best_nmismatches = max_mismatches_allowed;
      best_prob = 0.0;

      i = j = 0;
      while (i < donori_nsites && j < acceptorj_nsites) {
	if ((splice_pos = donori_positions[i]) < acceptorj_positions[j]) {
	  i++;
	} else if (splice_pos > acceptorj_positions[j]) {
	  j++;
	} else {
	  segmenti_nmismatches = Genome_count_mismatches_substring(query_compress,/*left*/segmenti_left,/*pos5*/0,/*pos3*/splice_pos,
								   plusp,genestrand,first_read_p);
	  segmentj_nmismatches = Genome_count_mismatches_substring(query_compress,/*left*/segmentj_left,/*pos5*/splice_pos,/*pos3*/querylength,
								   plusp,genestrand,first_read_p);
	  if ((nmismatches = segmenti_nmismatches + segmentj_nmismatches) <= best_nmismatches) {
	    if (donori_knowni[i] >= 0) {
	      probi = 1.0; /* Needs to be 1.0 for output */
	    } else {
	      probi = Maxent_hr_donor_prob(segmenti_left + splice_pos,segmenti_chroffset);
	    }

	    if (acceptorj_knowni[j] >= 0) {
	      probj = 1.0; /* Needs to be 1.0 for output */
	    } else {
	      probj = Maxent_hr_acceptor_prob(segmentj_left + splice_pos,segmentj_chroffset);
	    }

	    debug1(
		   if (plusp == true) {
		     printf("plus sense splice_pos  %d, i.donor %f, j.acceptor %f\n",splice_pos,probi,probj);
		   } else {
		     printf("minus antisense splice_pos  %d, i.donor %f, j.acceptor %f\n",splice_pos,probi,probj);
		   });

	    if (nmismatches < best_nmismatches ||
		(nmismatches == best_nmismatches && probi + probj > best_prob)) {
	      /* Success */
	      best_nmismatches = nmismatches;
	      best_prob = probi + probj;

	      best_donor_splicecoord = segmenti_left + splice_pos;
	      best_acceptor_splicecoord = segmentj_left + splice_pos;
	      best_donor_knowni = donori_knowni[i];
	      best_acceptor_knowni = acceptorj_knowni[j];
	      best_donor_prob = probi;
	      best_acceptor_prob = probj;
	      best_splice_pos = splice_pos;
	      best_segmenti_nmismatches = segmenti_nmismatches;
	      best_segmentj_nmismatches = segmentj_nmismatches;
	      orig_plusp = true;	/* for sense, require plusp to be true */
	    }
	  }
	  i++;
	  j++;
	}
      }

    } else {
      /* minus */
      /* Originally from minus strand.  Complement. */
      /* Antisense (End 7 to End 8) or Sense (End 3 to End 4) */
      if (novelsplicingp && segmenti_left + splice_pos_start >= ACCEPTOR_MODEL_RIGHT_MARGIN) {
	antiacceptori_nsites = Genome_antiacceptor_positions(acceptor_positions_alloc,acceptor_knowni_alloc,
							     segmenti_antiacceptor_knownpos,segmenti_antiacceptor_knowni,
							     segmenti_left,splice_pos_start,splice_pos_end);
	antiacceptori_positions = acceptor_positions_alloc;
	antiacceptori_knowni = acceptor_knowni_alloc;
      } else {
	antiacceptori_nsites = segmenti_antiacceptor_nknown;
	antiacceptori_positions = segmenti_antiacceptor_knownpos;
	antiacceptori_knowni = segmenti_antiacceptor_knowni;
      }

#ifdef DEBUG1
      printf("Found %d antiacceptori sites:",antiacceptori_nsites);
      for (i = 0; i < antiacceptori_nsites; i++) {
	printf(" %d",antiacceptori_positions[i]);
	if (antiacceptori_knowni[i] >= 0) {
	  printf(" (%d)",antiacceptori_knowni[i]);
	}
      }
      printf("\n");
#endif

      if (novelsplicingp && segmentj_left + splice_pos_start >= DONOR_MODEL_RIGHT_MARGIN) {
	antidonorj_nsites = Genome_antidonor_positions(donor_positions_alloc,donor_knowni_alloc,
						       segmentj_antidonor_knownpos,segmentj_antidonor_knowni,
						       segmentj_left,splice_pos_start,splice_pos_end);
	antidonorj_positions = donor_positions_alloc;
	antidonorj_knowni = donor_knowni_alloc;
      } else {
	antidonorj_nsites = segmentj_antidonor_nknown;
	antidonorj_positions = segmentj_antidonor_knownpos;
	antidonorj_knowni = segmentj_antidonor_knowni;
      }

#ifdef DEBUG1
      printf("Found %d antidonorj sites:",antidonorj_nsites);
      for (i = 0; i < antidonorj_nsites; i++) {
	printf(" %d",antidonorj_positions[i]);
	if (antidonorj_knowni[i] >= 0) {
	  printf(" (%d)",antidonorj_knowni[i]);
	}
      }
      printf("\n");
#endif

      best_nmismatches = max_mismatches_allowed;
      best_prob = 0.0;

      i = j = 0;
      while (i < antiacceptori_nsites && j < antidonorj_nsites) {
	if ((splice_pos = antiacceptori_positions[i]) < antidonorj_positions[j]) {
	  i++;
	} else if (splice_pos > antidonorj_positions[j]) {
	  j++;
	} else {
	  segmenti_nmismatches = Genome_count_mismatches_substring(query_compress,/*left*/segmenti_left,/*pos5*/0,/*pos3*/splice_pos,
								   plusp,genestrand,first_read_p);
	  segmentj_nmismatches = Genome_count_mismatches_substring(query_compress,/*left*/segmentj_left,/*pos5*/splice_pos,/*pos3*/querylength,
								   plusp,genestrand,first_read_p);
	  if ((nmismatches = segmenti_nmismatches + segmentj_nmismatches) <= best_nmismatches) {
	    if (antiacceptori_knowni[i] >= 0) {
	      probi = 1.0; /* Needs to be 1.0 for output */
	    } else {
	      probi = Maxent_hr_antiacceptor_prob(segmenti_left + splice_pos,segmenti_chroffset);
	    }

	    if (antidonorj_knowni[j] >= 0) {
	      probj = 1.0; /* Needs to be 1.0 for output */
	    } else {
	      probj = Maxent_hr_antidonor_prob(segmentj_left + splice_pos,segmentj_chroffset);
	    }

	    debug1(
		   if (plusp == true) {
		     printf("plus antisense splice_pos  %d, j.donor %f, i.acceptor %f\n",splice_pos,probj,probi);
		   } else {
		     printf("minus sense splice_pos  %d, j.donor %f, i.acceptor %f\n",splice_pos,probj,probi);
		   });
	  
	    if (nmismatches < best_nmismatches ||
		(nmismatches == best_nmismatches && probi + probj > best_prob)) {
	      /* Success */
	      best_nmismatches = nmismatches;
	      best_prob = probi + probj;

	      best_donor_splicecoord = segmentj_left + splice_pos;
	      best_acceptor_splicecoord = segmenti_left + splice_pos;
	      best_donor_knowni = antidonorj_knowni[j];
	      best_acceptor_knowni = antiacceptori_knowni[i];
	      best_donor_prob = probj;
	      best_acceptor_prob = probi;
	      best_splice_pos = splice_pos;
	      best_segmentj_nmismatches = segmentj_nmismatches;
	      best_segmenti_nmismatches = segmenti_nmismatches;
	      orig_plusp = false;	/* for sense, require plusp to be false */
	    }
	  }
	  i++;
	  j++;
	}
      }
    }

    if (best_prob > 0.0) {
      debug1(printf("best_prob = %f at splice_pos %d (%u,%u)\n",
		    best_prob,best_splice_pos,best_donor_splicecoord,best_acceptor_splicecoord));
      if (orig_plusp == true) {
	/* Originally from plus strand.  No complement. */
	sensedir = (plusp == true) ? SENSE_FORWARD : SENSE_ANTI;
	assert(sensedir == SENSE_FORWARD);

	donor = Substring_new_donor(best_donor_splicecoord,best_donor_knowni,
				    best_splice_pos,best_segmenti_nmismatches,
				    best_donor_prob,/*left*/segmenti_left,query_compress,
				    querylength,plusp,genestrand,first_read_p,sensedir,
				    segmenti_chrnum,segmenti_chroffset,segmenti_chrhigh,segmenti_chrlength);

	acceptor = Substring_new_acceptor(best_acceptor_splicecoord,best_acceptor_knowni,
					  best_splice_pos,best_segmentj_nmismatches,
					  best_acceptor_prob,/*left*/segmentj_left,query_compress,
					  querylength,plusp,genestrand,first_read_p,sensedir,
					  segmentj_chrnum,segmentj_chroffset,segmentj_chrhigh,segmentj_chrlength);

	if (donor == NULL || acceptor == NULL) {
	  if (donor != NULL) Substring_free(&donor);
	  if (acceptor != NULL) Substring_free(&acceptor);
	} else {
	  debug1(printf("Splice_solve_single_sense success\n"));
	  *segmenti_usedp = *segmentj_usedp = true;

	  donor_support = best_splice_pos;
	  acceptor_support = querylength - best_splice_pos;
	  sufficient1p = sufficient_splice_prob_local(donor_support,best_segmenti_nmismatches,best_donor_prob);
	  sufficient2p = sufficient_splice_prob_local(acceptor_support,best_segmentj_nmismatches,best_acceptor_prob);

	  if (sufficient1p && sufficient2p) {
	    *nhits += 1;
	    return List_push(hits,(void *) Stage3end_new_splice(&(*found_score),best_segmenti_nmismatches,best_segmentj_nmismatches,
								donor,acceptor,best_donor_prob,best_acceptor_prob,
								/*distance*/segmentj_left - segmenti_left,
								/*shortdistancep*/true,splicing_penalty,querylength,/*amb_length*/0,/*amb_prob*/0.0,
								/*ambcoords_donor*/NULL,/*ambcoords_acceptor*/NULL,
								/*amb_knowni_donor*/NULL,/*amb_knowni_acceptor*/NULL,
								/*amb_nmismatches_donor*/NULL,/*amb_nmismatches_acceptor*/NULL,
								/*amb_probs_donor*/NULL,/*amb_probs_acceptor*/NULL,
								/*copy_donor_p*/false,/*copy_acceptor_p*/false,first_read_p,sensedir,
								sarrayp));
	  } else if (subs_or_indels_p == true) {
	    if (donor != NULL) Substring_free(&donor);
	    if (acceptor != NULL) Substring_free(&acceptor);
	    return hits;
	  } else if (donor_support < LOWPROB_SUPPORT || acceptor_support < LOWPROB_SUPPORT) {
	    if (donor != NULL) Substring_free(&donor);
	    if (acceptor != NULL) Substring_free(&acceptor);
	    return hits;
	  } else if (sufficient1p || sufficient2p) {
	    *lowprob = List_push(*lowprob,
				 (void *) Stage3end_new_splice(&(*found_score),best_segmenti_nmismatches,best_segmentj_nmismatches,
							       donor,acceptor,best_donor_prob,best_acceptor_prob,
							       /*distance*/segmentj_left - segmenti_left,
							       /*shortdistancep*/true,splicing_penalty,querylength,/*amb_length*/0,/*amb_prob*/0.0,
							       /*ambcoords_donor*/NULL,/*ambcoords_acceptor*/NULL,
							       /*amb_knowni_donor*/NULL,/*amb_knowni_acceptor*/NULL,
							       /*amb_nmismatches_donor*/NULL,/*amb_nmismatches_acceptor*/NULL,
							       /*amb_probs_donor*/NULL,/*amb_probs_acceptor*/NULL,
							       /*copy_donor_p*/false,/*copy_acceptor_p*/false,first_read_p,sensedir,
							       sarrayp));
	    return hits;
	  } else {
	    if (donor != NULL) Substring_free(&donor);
	    if (acceptor != NULL) Substring_free(&acceptor);
	  }
	}

      } else {
	/* Originally from minus strand.  Complement. */
	sensedir = (plusp == true) ? SENSE_ANTI : SENSE_FORWARD;
	assert(sensedir == SENSE_FORWARD);

	donor = Substring_new_donor(best_donor_splicecoord,best_donor_knowni,
				    best_splice_pos,best_segmentj_nmismatches,
				    best_donor_prob,/*left*/segmentj_left,query_compress,
				    querylength,plusp,genestrand,first_read_p,sensedir,
				    segmentj_chrnum,segmentj_chroffset,segmentj_chrhigh,segmentj_chrlength);

	acceptor = Substring_new_acceptor(best_acceptor_splicecoord,best_acceptor_knowni,
					  best_splice_pos,best_segmenti_nmismatches,
					  best_acceptor_prob,/*left*/segmenti_left,query_compress,
					  querylength,plusp,genestrand,first_read_p,sensedir,
					  segmenti_chrnum,segmenti_chroffset,segmenti_chrhigh,segmenti_chrlength);

	if (donor == NULL || acceptor == NULL) {
	  if (donor != NULL) Substring_free(&donor);
	  if (acceptor != NULL) Substring_free(&acceptor);
	} else {
	  debug1(printf("Splice_solve_single_sense success\n"));
	  *segmenti_usedp = *segmentj_usedp = true;

	  acceptor_support = best_splice_pos;
	  donor_support = querylength - best_splice_pos;
	  sufficient1p = sufficient_splice_prob_local(acceptor_support,best_segmenti_nmismatches,best_acceptor_prob);
	  sufficient2p = sufficient_splice_prob_local(donor_support,best_segmentj_nmismatches,best_donor_prob);
	  if (sufficient1p && sufficient2p) {
	    *nhits += 1;
	    return List_push(hits,(void *) Stage3end_new_splice(&(*found_score),best_segmentj_nmismatches,best_segmenti_nmismatches,
								donor,acceptor,best_donor_prob,best_acceptor_prob,
								/*distance*/segmentj_left - segmenti_left,
								/*shortdistancep*/true,splicing_penalty,querylength,/*amb_length*/0,/*amb_prob*/0.0,
								/*ambcoords_donor*/NULL,/*ambcoords_acceptor*/NULL,
								/*amb_knowni_donor*/NULL,/*amb_knowni_acceptor*/NULL,
								/*amb_nmismatches_donor*/NULL,/*amb_nmismatches_acceptor*/NULL,
								/*amb_probs_donor*/NULL,/*amb_probs_acceptor*/NULL,
								/*copy_donor_p*/false,/*copy_acceptor_p*/false,first_read_p,sensedir,
								sarrayp));
	  } else if (subs_or_indels_p == true) {
	    if (donor != NULL) Substring_free(&donor);
	    if (acceptor != NULL) Substring_free(&acceptor);
	    return hits;
	  } else if (donor_support < LOWPROB_SUPPORT || acceptor_support < LOWPROB_SUPPORT) {
	    if (donor != NULL) Substring_free(&donor);
	    if (acceptor != NULL) Substring_free(&acceptor);
	    return hits;
	  } else if (sufficient1p || sufficient2p) {
	    *lowprob = List_push(*lowprob,
				 (void *) Stage3end_new_splice(&(*found_score),best_segmentj_nmismatches,best_segmenti_nmismatches,
							       donor,acceptor,best_donor_prob,best_acceptor_prob,
							       /*distance*/segmentj_left - segmenti_left,
							       /*shortdistancep*/true,splicing_penalty,querylength,/*amb_length*/0,/*amb_prob*/0.0,
							       /*ambcoords_donor*/NULL,/*ambcoords_acceptor*/NULL,
							       /*amb_knowni_donor*/NULL,/*amb_knowni_acceptor*/NULL,
							       /*amb_nmismatches_donor*/NULL,/*amb_nmismatches_acceptor*/NULL,
							       /*amb_probs_donor*/NULL,/*amb_probs_acceptor*/NULL,
							       /*copy_donor_p*/false,/*copy_acceptor_p*/false,first_read_p,sensedir,
							       sarrayp));
	    return hits;
	  } else {
	    if (donor != NULL) Substring_free(&donor);
	    if (acceptor != NULL) Substring_free(&acceptor);
	    return hits;
	  }
	}
      }
    }
  }

  debug1(printf("Splice_solve_single_sense fail\n"));
  return hits;
}


List_T
Splice_solve_single_antisense (int *found_score, int *nhits, List_T hits, List_T *lowprob,

			       bool *segmenti_usedp, bool *segmentj_usedp,
			       Univcoord_T segmenti_left, Univcoord_T segmentj_left,
			       Chrnum_T segmenti_chrnum, Univcoord_T segmenti_chroffset,
			       Univcoord_T segmenti_chrhigh, Chrpos_T segmenti_chrlength,
			       Chrnum_T segmentj_chrnum, Univcoord_T segmentj_chroffset,
			       Univcoord_T segmentj_chrhigh, Chrpos_T segmentj_chrlength,
		     
			       int querylength, Compress_T query_compress,
			       int *segmenti_donor_knownpos, int *segmentj_acceptor_knownpos,
			       int *segmentj_antidonor_knownpos, int *segmenti_antiacceptor_knownpos,
			       int *segmenti_donor_knowni, int *segmentj_acceptor_knowni,
			       int *segmentj_antidonor_knowni, int *segmenti_antiacceptor_knowni,
			       int segmenti_donor_nknown, int segmentj_acceptor_nknown,
			       int segmentj_antidonor_nknown, int segmenti_antiacceptor_nknown,
			       int splicing_penalty, int max_mismatches_allowed,
			       bool plusp, int genestrand, bool first_read_p,
			       bool subs_or_indels_p, bool sarrayp) {
  Substring_T donor, acceptor;
  int best_splice_pos, splice_pos_start, splice_pos_end, splice_pos, i, j;

  int best_nmismatches, nmismatches;
  int best_segmenti_nmismatches, best_segmentj_nmismatches, segmenti_nmismatches, segmentj_nmismatches;
  int donor_support, acceptor_support;
  Univcoord_T best_donor_splicecoord, best_acceptor_splicecoord;
  int best_donor_knowni, best_acceptor_knowni;
  double best_prob, best_donor_prob, best_acceptor_prob, probi, probj;
  bool sufficient1p, sufficient2p, orig_plusp;
  int sensedir;

  int donori_nsites, acceptorj_nsites, antiacceptori_nsites, antidonorj_nsites;
  int *donori_positions, *acceptorj_positions, *antiacceptori_positions, *antidonorj_positions;
  int *donori_knowni, *acceptorj_knowni, *antiacceptori_knowni, *antidonorj_knowni;

#ifdef HAVE_ALLOCA
  int *donor_positions_alloc = (int *) alloca((querylength+1)*sizeof(int));
  int *acceptor_positions_alloc = (int *) alloca((querylength+1)*sizeof(int));
  int *donor_knowni_alloc = (int *) alloca((querylength+1)*sizeof(int));
  int *acceptor_knowni_alloc = (int *) alloca((querylength+1)*sizeof(int));
#else
  int donor_positions_alloc[MAX_READLENGTH+1], acceptor_positions_alloc[MAX_READLENGTH+1];
  int donor_knowni_alloc[MAX_READLENGTH+1], acceptor_knowni_alloc[MAX_READLENGTH+1];
#endif

  debug1(printf("Splice_solve_single: Getting genome at lefti %u and leftj %u (diff: %d)\n",
		segmenti_left,segmentj_left,segmentj_left-segmenti_left));
  *nhits = 0;

#if 0
  int sum, lefti, righti;
  splice_pos_start = querylength;
  splice_pos_end = 0;
  for (sum = 0; sum <= max_mismatches_allowed; sum++) {
    for (lefti = 0; lefti <= sum && lefti < nmismatches_left; lefti++) {
      if ((righti = sum - lefti) < nmismatches_right &&
	  mismatch_positions_left[lefti] > mismatch_positions_right[righti]) {
	debug1(printf("At %d+%d mismatches, splice_pos using right: %d\n",lefti,righti,mismatch_positions_right[righti]+1));
	debug1(printf("At %d+%d mismatches, splice_pos using left: %d\n",lefti,righti,mismatch_positions_left[lefti]));
	if (mismatch_positions_right[righti] + 1 < splice_pos_start) {
	  splice_pos_start = mismatch_positions_right[righti] + 1;	/* This is leftmost position in righti+1 .. lefti */
	}
	if (mismatch_positions_left[lefti] > splice_pos_end) {
	  splice_pos_end = mismatch_positions_left[lefti];	/* This is rightmost position in righti+1 .. lefti */
	}
      }
    }
  }

  /* Exclude ends */
  if (splice_pos_start < min_localsplicing_end_matches) {
    splice_pos_start = min_localsplicing_end_matches;
  }
  if (splice_pos_end > querylength - min_localsplicing_end_matches) {
    splice_pos_end = querylength - min_localsplicing_end_matches;
  }
#else
  /* splice_pos_start = min_localsplicing_end_matches; */
  /* splice_pos_end = querylength - min_localsplicing_end_matches; */
  splice_pos_start = min_shortend;
  splice_pos_end = querylength - min_shortend; /* ? off by 1, so -l 3 allows only ends of up to 2 */
#endif


  if (splice_pos_start <= splice_pos_end) {
    if (plusp == false) {
      /* minus */
      /* Originally from plus strand.  No complement.  */
      /* Sense (End 1 to End 2) or Antisense (End 5 to End 6) */
      if (novelsplicingp && segmenti_left + splice_pos_start >= DONOR_MODEL_LEFT_MARGIN) {
	donori_nsites = Genome_donor_positions(donor_positions_alloc,donor_knowni_alloc,
					       segmenti_donor_knownpos,segmenti_donor_knowni,
					       segmenti_left,splice_pos_start,splice_pos_end);
	donori_positions = donor_positions_alloc;
	donori_knowni = donor_knowni_alloc;
      } else {
	donori_nsites = segmenti_donor_nknown;
	donori_positions = segmenti_donor_knownpos;
	donori_knowni = segmenti_donor_knowni;
      }

#ifdef DEBUG1
      printf("Found %d donori sites:",donori_nsites);
      for (i = 0; i < donori_nsites; i++) {
	printf(" %d",donori_positions[i]);
	if (donori_knowni[i] >= 0) {
	  printf(" (%d)",donori_knowni[i]);
	}
      }
      printf("\n");
#endif

      if (novelsplicingp && segmentj_left + splice_pos_start >= ACCEPTOR_MODEL_LEFT_MARGIN) {
	acceptorj_nsites = Genome_acceptor_positions(acceptor_positions_alloc,acceptor_knowni_alloc,
						     segmentj_acceptor_knownpos,segmentj_acceptor_knowni,
						     segmentj_left,splice_pos_start,splice_pos_end);
	acceptorj_positions = acceptor_positions_alloc;
	acceptorj_knowni = acceptor_knowni_alloc;
      } else {
	acceptorj_nsites = segmentj_acceptor_nknown;
	acceptorj_positions = segmentj_acceptor_knownpos;
	acceptorj_knowni = segmentj_acceptor_knowni;
      }

#ifdef DEBUG1
      printf("Found %d acceptorj sites:",acceptorj_nsites);
      for (i = 0; i < acceptorj_nsites; i++) {
	printf(" %d",acceptorj_positions[i]);
	if (acceptorj_knowni[i] >= 0) {
	  printf(" (%d)",acceptorj_knowni[i]);
	}
      }
      printf("\n");
#endif

      best_nmismatches = max_mismatches_allowed;
      best_prob = 0.0;

      i = j = 0;
      while (i < donori_nsites && j < acceptorj_nsites) {
	if ((splice_pos = donori_positions[i]) < acceptorj_positions[j]) {
	  i++;
	} else if (splice_pos > acceptorj_positions[j]) {
	  j++;
	} else {
	  segmenti_nmismatches = Genome_count_mismatches_substring(query_compress,/*left*/segmenti_left,/*pos5*/0,/*pos3*/splice_pos,
								   plusp,genestrand,first_read_p);
	  segmentj_nmismatches = Genome_count_mismatches_substring(query_compress,/*left*/segmentj_left,/*pos5*/splice_pos,/*pos3*/querylength,
								   plusp,genestrand,first_read_p);
	  if ((nmismatches = segmenti_nmismatches + segmentj_nmismatches) <= best_nmismatches) {
	    if (donori_knowni[i] >= 0) {
	      probi = 1.0; /* Needs to be 1.0 for output */
	    } else {
	      probi = Maxent_hr_donor_prob(segmenti_left + splice_pos,segmenti_chroffset);
	    }

	    if (acceptorj_knowni[j] >= 0) {
	      probj = 1.0; /* Needs to be 1.0 for output */
	    } else {
	      probj = Maxent_hr_acceptor_prob(segmentj_left + splice_pos,segmentj_chroffset);
	    }

	    debug1(
		   if (plusp == true) {
		     printf("plus sense splice_pos  %d, i.donor %f, j.acceptor %f\n",splice_pos,probi,probj);
		   } else {
		     printf("minus antisense splice_pos  %d, i.donor %f, j.acceptor %f\n",splice_pos,probi,probj);
		   });

	    if (nmismatches < best_nmismatches ||
		(nmismatches == best_nmismatches && probi + probj > best_prob)) {
	      /* Success */
	      best_nmismatches = nmismatches;
	      best_prob = probi + probj;

	      best_donor_splicecoord = segmenti_left + splice_pos;
	      best_acceptor_splicecoord = segmentj_left + splice_pos;
	      best_donor_knowni = donori_knowni[i];
	      best_acceptor_knowni = acceptorj_knowni[j];
	      best_donor_prob = probi;
	      best_acceptor_prob = probj;
	      best_splice_pos = splice_pos;
	      best_segmenti_nmismatches = segmenti_nmismatches;
	      best_segmentj_nmismatches = segmentj_nmismatches;
	      orig_plusp = true;	/* for antisense, require plusp to be false */
	    }
	  }
	  i++;
	  j++;
	}
      }

    } else {
      /* plus */
      /* Originally from minus strand.  Complement. */
      /* Antisense (End 7 to End 8) or Sense (End 3 to End 4) */
      if (novelsplicingp && segmenti_left + splice_pos_start >= ACCEPTOR_MODEL_RIGHT_MARGIN) {
	antiacceptori_nsites = Genome_antiacceptor_positions(acceptor_positions_alloc,acceptor_knowni_alloc,
							     segmenti_antiacceptor_knownpos,segmenti_antiacceptor_knowni,
							     segmenti_left,splice_pos_start,splice_pos_end);
	antiacceptori_positions = acceptor_positions_alloc;
	antiacceptori_knowni = acceptor_knowni_alloc;
      } else {
	antiacceptori_nsites = segmenti_antiacceptor_nknown;
	antiacceptori_positions = segmenti_antiacceptor_knownpos;
	antiacceptori_knowni = segmenti_antiacceptor_knowni;
      }

#ifdef DEBUG1
      printf("Found %d antiacceptori sites:",antiacceptori_nsites);
      for (i = 0; i < antiacceptori_nsites; i++) {
	printf(" %d",antiacceptori_positions[i]);
	if (antiacceptori_knowni[i] >= 0) {
	  printf(" (%d)",antiacceptori_knowni[i]);
	}
      }
      printf("\n");
#endif

      if (novelsplicingp && segmentj_left + splice_pos_start >= DONOR_MODEL_RIGHT_MARGIN) {
	antidonorj_nsites = Genome_antidonor_positions(donor_positions_alloc,donor_knowni_alloc,
						       segmentj_antidonor_knownpos,segmentj_antidonor_knowni,
						       segmentj_left,splice_pos_start,splice_pos_end);
	antidonorj_positions = donor_positions_alloc;
	antidonorj_knowni = donor_knowni_alloc;
      } else {
	antidonorj_nsites = segmentj_antidonor_nknown;
	antidonorj_positions = segmentj_antidonor_knownpos;
	antidonorj_knowni = segmentj_antidonor_knowni;
      }

#ifdef DEBUG1
      printf("Found %d antidonorj sites:",antidonorj_nsites);
      for (i = 0; i < antidonorj_nsites; i++) {
	printf(" %d",antidonorj_positions[i]);
	if (antidonorj_knowni[i] >= 0) {
	  printf(" (%d)",antidonorj_knowni[i]);
	}
      }
      printf("\n");
#endif

      best_nmismatches = max_mismatches_allowed;
      best_prob = 0.0;

      i = j = 0;
      while (i < antiacceptori_nsites && j < antidonorj_nsites) {
	if ((splice_pos = antiacceptori_positions[i]) < antidonorj_positions[j]) {
	  i++;
	} else if (splice_pos > antidonorj_positions[j]) {
	  j++;
	} else {
	  segmenti_nmismatches = Genome_count_mismatches_substring(query_compress,/*left*/segmenti_left,/*pos5*/0,/*pos3*/splice_pos,
								   plusp,genestrand,first_read_p);
	  segmentj_nmismatches = Genome_count_mismatches_substring(query_compress,/*left*/segmentj_left,/*pos5*/splice_pos,/*pos3*/querylength,
								   plusp,genestrand,first_read_p);
	  if ((nmismatches = segmenti_nmismatches + segmentj_nmismatches) <= best_nmismatches) {
	    if (antiacceptori_knowni[i] >= 0) {
	      probi = 1.0; /* Needs to be 1.0 for output */
	    } else {
	      probi = Maxent_hr_antiacceptor_prob(segmenti_left + splice_pos,segmenti_chroffset);
	    }

	    if (antidonorj_knowni[j] >= 0) {
	      probj = 1.0; /* Needs to be 1.0 for output */
	    } else {
	      probj = Maxent_hr_antidonor_prob(segmentj_left + splice_pos,segmentj_chroffset);
	    }

	    debug1(
		   if (plusp == true) {
		     printf("plus antisense splice_pos  %d, j.donor %f, i.acceptor %f\n",splice_pos,probj,probi);
		   } else {
		     printf("minus sense splice_pos  %d, j.donor %f, i.acceptor %f\n",splice_pos,probj,probi);
		   });
	  
	    if (nmismatches < best_nmismatches ||
		(nmismatches == best_nmismatches && probi + probj > best_prob)) {
	      /* Success */
	      best_nmismatches = nmismatches;
	      best_prob = probi + probj;

	      best_donor_splicecoord = segmentj_left + splice_pos;
	      best_acceptor_splicecoord = segmenti_left + splice_pos;
	      best_donor_knowni = antidonorj_knowni[j];
	      best_acceptor_knowni = antiacceptori_knowni[i];
	      best_donor_prob = probj;
	      best_acceptor_prob = probi;
	      best_splice_pos = splice_pos;
	      best_segmentj_nmismatches = segmentj_nmismatches;
	      best_segmenti_nmismatches = segmenti_nmismatches;
	      orig_plusp = false;	/* for antisense, require plusp to be true */
	    }
	  }
	  i++;
	  j++;
	}
      }
    }

    if (best_prob > 0.0) {
      debug1(printf("best_prob = %f at splice_pos %d (%u,%u)\n",
		    best_prob,best_splice_pos,best_donor_splicecoord,best_acceptor_splicecoord));
      if (orig_plusp == true) {
	/* Originally from plus strand.  No complement. */
	sensedir = (plusp == true) ? SENSE_FORWARD : SENSE_ANTI;
	assert(sensedir == SENSE_ANTI);

	donor = Substring_new_donor(best_donor_splicecoord,best_donor_knowni,
				    best_splice_pos,best_segmenti_nmismatches,
				    best_donor_prob,/*left*/segmenti_left,query_compress,
				    querylength,plusp,genestrand,first_read_p,sensedir,
				    segmenti_chrnum,segmenti_chroffset,segmenti_chrhigh,segmenti_chrlength);

	acceptor = Substring_new_acceptor(best_acceptor_splicecoord,best_acceptor_knowni,
					  best_splice_pos,best_segmentj_nmismatches,
					  best_acceptor_prob,/*left*/segmentj_left,query_compress,
					  querylength,plusp,genestrand,first_read_p,sensedir,
					  segmentj_chrnum,segmentj_chroffset,segmentj_chrhigh,segmentj_chrlength);

	if (donor == NULL || acceptor == NULL) {
	  if (donor != NULL) Substring_free(&donor);
	  if (acceptor != NULL) Substring_free(&acceptor);
	} else {
	  debug1(printf("Splice_solve_single_antisense success\n"));
	  *segmenti_usedp = *segmentj_usedp = true;

	  donor_support = best_splice_pos;
	  acceptor_support = querylength - best_splice_pos;
	  sufficient1p = sufficient_splice_prob_local(donor_support,best_segmenti_nmismatches,best_donor_prob);
	  sufficient2p = sufficient_splice_prob_local(acceptor_support,best_segmentj_nmismatches,best_acceptor_prob);

	  if (sufficient1p && sufficient2p) {
	    *nhits += 1;
	    return List_push(hits,(void *) Stage3end_new_splice(&(*found_score),best_segmenti_nmismatches,best_segmentj_nmismatches,
								donor,acceptor,best_donor_prob,best_acceptor_prob,
								/*distance*/segmentj_left - segmenti_left,
								/*shortdistancep*/true,splicing_penalty,querylength,/*amb_length*/0,/*amb_prob*/0.0,
								/*ambcoords_donor*/NULL,/*ambcoords_acceptor*/NULL,
								/*amb_knowni_donor*/NULL,/*amb_knowni_acceptor*/NULL,
								/*amb_nmismatches_donor*/NULL,/*amb_nmismatches_acceptor*/NULL,
								/*amb_probs_donor*/NULL,/*amb_probs_acceptor*/NULL,
								/*copy_donor_p*/false,/*copy_acceptor_p*/false,first_read_p,sensedir,
								sarrayp));
	  } else if (subs_or_indels_p == true) {
	    if (donor != NULL) Substring_free(&donor);
	    if (acceptor != NULL) Substring_free(&acceptor);
	    return hits;
	  } else if (donor_support < LOWPROB_SUPPORT || acceptor_support < LOWPROB_SUPPORT) {
	    if (donor != NULL) Substring_free(&donor);
	    if (acceptor != NULL) Substring_free(&acceptor);
	    return hits;
	  } else if (sufficient1p || sufficient2p) {
	    *lowprob = List_push(*lowprob,
				 (void *) Stage3end_new_splice(&(*found_score),best_segmenti_nmismatches,best_segmentj_nmismatches,
							       donor,acceptor,best_donor_prob,best_acceptor_prob,
							       /*distance*/segmentj_left - segmenti_left,
							       /*shortdistancep*/true,splicing_penalty,querylength,/*amb_length*/0,/*amb_prob*/0.0,
							       /*ambcoords_donor*/NULL,/*ambcoords_acceptor*/NULL,
							       /*amb_knowni_donor*/NULL,/*amb_knowni_acceptor*/NULL,
							       /*amb_nmismatches_donor*/NULL,/*amb_nmismatches_acceptor*/NULL,
							       /*amb_probs_donor*/NULL,/*amb_probs_acceptor*/NULL,
							       /*copy_donor_p*/false,/*copy_acceptor_p*/false,first_read_p,sensedir,
							       sarrayp));
	    return hits;
	  } else {
	    if (donor != NULL) Substring_free(&donor);
	    if (acceptor != NULL) Substring_free(&acceptor);
	  }
	}

      } else {
	/* Originally from minus strand.  Complement. */
	sensedir = (plusp == true) ? SENSE_ANTI : SENSE_FORWARD;
	assert(sensedir == SENSE_ANTI);

	donor = Substring_new_donor(best_donor_splicecoord,best_donor_knowni,
				    best_splice_pos,best_segmentj_nmismatches,
				    best_donor_prob,/*left*/segmentj_left,query_compress,
				    querylength,plusp,genestrand,first_read_p,sensedir,
				    segmentj_chrnum,segmentj_chroffset,segmentj_chrhigh,segmentj_chrlength);

	acceptor = Substring_new_acceptor(best_acceptor_splicecoord,best_acceptor_knowni,
					  best_splice_pos,best_segmenti_nmismatches,
					  best_acceptor_prob,/*left*/segmenti_left,query_compress,
					  querylength,plusp,genestrand,first_read_p,sensedir,
					  segmenti_chrnum,segmenti_chroffset,segmenti_chrhigh,segmenti_chrlength);

	if (donor == NULL || acceptor == NULL) {
	  if (donor != NULL) Substring_free(&donor);
	  if (acceptor != NULL) Substring_free(&acceptor);
	} else {
	  debug1(printf("Splice_solve_single_antisense success\n"));
	  *segmenti_usedp = *segmentj_usedp = true;

	  acceptor_support = best_splice_pos;
	  donor_support = querylength - best_splice_pos;
	  sufficient1p = sufficient_splice_prob_local(acceptor_support,best_segmenti_nmismatches,best_acceptor_prob);
	  sufficient2p = sufficient_splice_prob_local(donor_support,best_segmentj_nmismatches,best_donor_prob);
	  if (sufficient1p && sufficient2p) {
	    *nhits += 1;
	    return List_push(hits,(void *) Stage3end_new_splice(&(*found_score),best_segmentj_nmismatches,best_segmenti_nmismatches,
								donor,acceptor,best_donor_prob,best_acceptor_prob,
								/*distance*/segmentj_left - segmenti_left,
								/*shortdistancep*/true,splicing_penalty,querylength,/*amb_length*/0,/*amb_prob*/0.0,
								/*ambcoords_donor*/NULL,/*ambcoords_acceptor*/NULL,
								/*amb_knowni_donor*/NULL,/*amb_knowni_acceptor*/NULL,
								/*amb_nmismatches_donor*/NULL,/*amb_nmismatches_acceptor*/NULL,
								/*amb_probs_donor*/NULL,/*amb_probs_acceptor*/NULL,
								/*copy_donor_p*/false,/*copy_acceptor_p*/false,first_read_p,sensedir,
								sarrayp));
	  } else if (subs_or_indels_p == true) {
	    if (donor != NULL) Substring_free(&donor);
	    if (acceptor != NULL) Substring_free(&acceptor);
	    return hits;
	  } else if (donor_support < LOWPROB_SUPPORT || acceptor_support < LOWPROB_SUPPORT) {
	    if (donor != NULL) Substring_free(&donor);
	    if (acceptor != NULL) Substring_free(&acceptor);
	    return hits;
	  } else if (sufficient1p || sufficient2p) {
	    *lowprob = List_push(*lowprob,
				 (void *) Stage3end_new_splice(&(*found_score),best_segmentj_nmismatches,best_segmenti_nmismatches,
							       donor,acceptor,best_donor_prob,best_acceptor_prob,
							       /*distance*/segmentj_left - segmenti_left,
							       /*shortdistancep*/true,splicing_penalty,querylength,/*amb_length*/0,/*amb_prob*/0.0,
							       /*ambcoords_donor*/NULL,/*ambcoords_acceptor*/NULL,
							       /*amb_knowni_donor*/NULL,/*amb_knowni_acceptor*/NULL,
							       /*amb_nmismatches_donor*/NULL,/*amb_nmismatches_acceptor*/NULL,
							       /*amb_probs_donor*/NULL,/*amb_probs_acceptor*/NULL,
							       /*copy_donor_p*/false,/*copy_acceptor_p*/false,first_read_p,sensedir,
							       sarrayp));
	    return hits;
	  } else {
	    if (donor != NULL) Substring_free(&donor);
	    if (acceptor != NULL) Substring_free(&acceptor);
	    return hits;
	  }
	}
      }
    }
  }

  debug1(printf("Splice_solve_single_antisense fail\n"));
  return hits;
}


List_T
Splice_solve_double (int *found_score, int *nhits, List_T hits, List_T *lowprob,

		     bool *segmenti_usedp, bool *segmentm_usedp, bool *segmentj_usedp,
		     Univcoord_T segmenti_left, Univcoord_T segmentm_left, Univcoord_T segmentj_left,
		     Chrnum_T segmenti_chrnum, Univcoord_T segmenti_chroffset,
		     Univcoord_T segmenti_chrhigh, Chrpos_T segmenti_chrlength,
		     Chrnum_T segmentm_chrnum, Univcoord_T segmentm_chroffset,
		     Univcoord_T segmentm_chrhigh, Chrpos_T segmentm_chrlength,
		     Chrnum_T segmentj_chrnum, Univcoord_T segmentj_chroffset,
		     Univcoord_T segmentj_chrhigh, Chrpos_T segmentj_chrlength,

		     int querylength, Compress_T query_compress,
		     int *segmenti_donor_knownpos, int *segmentm_acceptor_knownpos, int *segmentm_donor_knownpos, int *segmentj_acceptor_knownpos,
		     int *segmentj_antidonor_knownpos, int *segmentm_antiacceptor_knownpos, int *segmentm_antidonor_knownpos, int *segmenti_antiacceptor_knownpos,
		     int *segmenti_donor_knowni, int *segmentm_acceptor_knowni, int *segmentm_donor_knowni, int *segmentj_acceptor_knowni,
		     int *segmentj_antidonor_knowni, int *segmentm_antiacceptor_knowni, int *segmentm_antidonor_knowni, int *segmenti_antiacceptor_knowni,
		     int segmenti_donor_nknown, int segmentm_acceptor_nknown, int segmentm_donor_nknown, int segmentj_acceptor_nknown,
		     int segmentj_antidonor_nknown, int segmentm_antiacceptor_nknown, int segmentm_antidonor_nknown, int segmenti_antiacceptor_nknown,
		     int splicing_penalty, int max_mismatches_allowed, bool plusp, int genestrand, bool first_read_p,
		     bool subs_or_indels_p, bool sarrayp) {
  Substring_T donor, shortexon, acceptor;
  int best_splice_pos_1, best_splice_pos_2, splice_pos_start, splice_pos_end, splice_pos_1, splice_pos_2;
  int i, a, b, j;

  int best_nmismatches, nmismatches;
  int best_segmenti_nmismatches, best_segmentm_nmismatches, best_segmentj_nmismatches,
    segmenti_nmismatches, segmentm_nmismatches, segmentj_nmismatches;
  int donor_support, acceptor_support, middle_support;
  Univcoord_T best_donor1_splicecoord, best_acceptor1_splicecoord, best_donor2_splicecoord, best_acceptor2_splicecoord;
  int best_donor1_knowni, best_acceptor1_knowni, best_donor2_knowni, best_acceptor2_knowni;
  double best_prob, best_donor1_prob, best_acceptor1_prob, best_donor2_prob, best_acceptor2_prob,
    probi, proba, probb, probj;
  bool sufficient1p, sufficient2p, sufficient3p, sufficient4p, orig_plusp, matchp;
  int sensedir;

  int donori_nsites, acceptora_nsites, donorb_nsites, acceptorj_nsites,
    antiacceptori_nsites, antidonora_nsites, antiacceptorb_nsites, antidonorj_nsites;
  int *donori_positions, *acceptora_positions, *donorb_positions, *acceptorj_positions,
    *antiacceptori_positions, *antidonora_positions, *antiacceptorb_positions, *antidonorj_positions;
  int *donori_knowni, *acceptora_knowni, *donorb_knowni, *acceptorj_knowni,
    *antiacceptori_knowni, *antidonora_knowni, *antiacceptorb_knowni, *antidonorj_knowni;

#ifdef HAVE_ALLOCA
  int *donor1_positions_alloc = (int *) alloca((querylength+1)*sizeof(int));
  int *acceptor1_positions_alloc = (int *) alloca((querylength+1)*sizeof(int));
  int *donor2_positions_alloc = (int *) alloca((querylength+1)*sizeof(int));
  int *acceptor2_positions_alloc = (int *) alloca((querylength+1)*sizeof(int));
  int *donor1_knowni_alloc = (int *) alloca((querylength+1)*sizeof(int));
  int *acceptor1_knowni_alloc = (int *) alloca((querylength+1)*sizeof(int));
  int *donor2_knowni_alloc = (int *) alloca((querylength+1)*sizeof(int));
  int *acceptor2_knowni_alloc = (int *) alloca((querylength+1)*sizeof(int));
#else  
  int donor1_positions_alloc[MAX_READLENGTH+1], acceptor1_positions_alloc[MAX_READLENGTH+1],
    donor2_positions_alloc[MAX_READLENGTH+1], acceptor2_positions_alloc[MAX_READLENGTH+1];
  int donor1_knowni_alloc[MAX_READLENGTH+1], acceptor1_knowni_alloc[MAX_READLENGTH+1],
    donor2_knowni_alloc[MAX_READLENGTH+1], acceptor2_knowni_alloc[MAX_READLENGTH+1];
#endif


  debug2(printf("Splice_solve_double: Getting genome at lefti %u, leftm %u, and leftj %u\n",
		segmenti_left,segmentm_left,segmentj_left));

  *nhits = 0;
  splice_pos_start = min_shortend;
  splice_pos_end = querylength - min_shortend; /* ? off by 1, so -l 3 allows only ends of up to 2 */

  if (splice_pos_start <= splice_pos_end) {
    /* Originally from plus strand.  No complement. */
    /* Sense (End 1 to End 2) or Antisense (End 5 to End 6) */

    /* Segment i */
    if (novelsplicingp && segmenti_left + splice_pos_start >= DONOR_MODEL_LEFT_MARGIN) {
      donori_nsites = Genome_donor_positions(donor1_positions_alloc,donor1_knowni_alloc,
					     segmenti_donor_knownpos,segmenti_donor_knowni,
					     segmenti_left,splice_pos_start,splice_pos_end);
      donori_positions = donor1_positions_alloc;
      donori_knowni = donor1_knowni_alloc;
    } else {
      donori_nsites = segmenti_donor_nknown;
      donori_positions = segmenti_donor_knownpos;
      donori_knowni = segmenti_donor_knowni;
    }

#ifdef DEBUG2
    printf("Found %d donori sites:",donori_nsites);
    for (i = 0; i < donori_nsites; i++) {
      printf(" %d",donori_positions[i]);
      if (donori_knowni[i] >= 0) {
	printf(" (%d)",donori_knowni[i]);
      }
    }
    printf("\n");
#endif

    /* Segment m1 */
    if (novelsplicingp && segmentm_left + splice_pos_start >= ACCEPTOR_MODEL_LEFT_MARGIN) {
      acceptora_nsites = Genome_acceptor_positions(acceptor1_positions_alloc,acceptor1_knowni_alloc,
						   segmentm_acceptor_knownpos,segmentm_acceptor_knowni,
						   segmentm_left,splice_pos_start,splice_pos_end);
      acceptora_positions = acceptor1_positions_alloc;
      acceptora_knowni = acceptor1_knowni_alloc;
    } else {
      acceptora_nsites = segmentm_acceptor_nknown;
      acceptora_positions = segmentm_acceptor_knownpos;
      acceptora_knowni = segmentm_acceptor_knowni;
    }

#ifdef DEBUG2
    printf("Found %d acceptora sites:",acceptora_nsites);
    for (i = 0; i < acceptora_nsites; i++) {
      printf(" %d",acceptora_positions[i]);
      if (acceptora_knowni[i] >= 0) {
	printf(" (%d)",acceptora_knowni[i]);
      }
    }
    printf("\n");
#endif

    /* Segment m2 */
    if (novelsplicingp && segmentm_left + splice_pos_start >= DONOR_MODEL_LEFT_MARGIN) {
      donorb_nsites = Genome_donor_positions(donor2_positions_alloc,donor2_knowni_alloc,
					     segmentm_donor_knownpos,segmentm_donor_knowni,
					     segmentm_left,splice_pos_start,splice_pos_end);
      donorb_positions = donor2_positions_alloc;
      donorb_knowni = donor2_knowni_alloc;
    } else {
      donorb_nsites = segmentm_donor_nknown;
      donorb_positions = segmentm_donor_knownpos;
      donorb_knowni = segmentm_donor_knowni;
    }

#ifdef DEBUG2
    printf("Found %d donorb sites:",donorb_nsites);
    for (i = 0; i < donorb_nsites; i++) {
      printf(" %d",donorb_positions[i]);
      if (donorb_knowni[i] >= 0) {
	printf(" (%d)",donorb_knowni[i]);
      }
    }
    printf("\n");
#endif

    /* Segment j */
    if (novelsplicingp && segmentj_left + splice_pos_start >= ACCEPTOR_MODEL_LEFT_MARGIN) {
      acceptorj_nsites = Genome_acceptor_positions(acceptor2_positions_alloc,acceptor2_knowni_alloc,
						   segmentj_acceptor_knownpos,segmentj_acceptor_knowni,
						   segmentj_left,splice_pos_start,splice_pos_end);
      acceptorj_positions = acceptor2_positions_alloc;
      acceptorj_knowni = acceptor2_knowni_alloc;
    } else {
      acceptorj_nsites = segmentj_acceptor_nknown;
      acceptorj_positions = segmentj_acceptor_knownpos;
      acceptorj_knowni = segmentj_acceptor_knowni;
    }

#ifdef DEBUG2
    printf("Found %d acceptorj sites:",acceptorj_nsites);
    for (i = 0; i < acceptorj_nsites; i++) {
      printf(" %d",acceptorj_positions[i]);
      if (acceptorj_knowni[i] >= 0) {
	printf(" (%d)",acceptorj_knowni[i]);
      }
    }
    printf("\n");
#endif

    best_nmismatches = max_mismatches_allowed;
    best_prob = 0.0;
    orig_plusp = true;

    i = a = b = j = 0;
    while (i < donori_nsites && a < acceptora_nsites) {
      if ((splice_pos_1 = donori_positions[i]) < acceptora_positions[a]) {
	i++;
      } else if (splice_pos_1 > acceptora_positions[a]) {
	a++;
      } else {
	while (b < donorb_nsites && donorb_positions[b] <= splice_pos_1) {
	  b++;
	}
	while (j < acceptorj_nsites && acceptorj_positions[j] <= splice_pos_1) {
	  j++;
	}
	matchp = false;
	while (b < donorb_nsites && j < acceptorj_nsites && matchp == false) {
	  if ((splice_pos_2 = donorb_positions[b]) < acceptorj_positions[j]) {
	    b++;
	  } else if (splice_pos_2 > acceptorj_positions[j]) {
	    j++;
	  } else {
	    segmenti_nmismatches = Genome_count_mismatches_substring(query_compress,/*left*/segmenti_left,/*pos5*/0,/*pos3*/splice_pos_1,
								     plusp,genestrand,first_read_p);
	    segmentm_nmismatches = Genome_count_mismatches_substring(query_compress,/*left*/segmentm_left,/*pos5*/splice_pos_1,/*pos3*/splice_pos_2,
								     plusp,genestrand,first_read_p);
	    segmentj_nmismatches = Genome_count_mismatches_substring(query_compress,/*left*/segmentj_left,/*pos5*/splice_pos_2,/*pos3*/querylength,
								     plusp,genestrand,first_read_p);
	    if ((nmismatches = segmenti_nmismatches + segmentm_nmismatches + segmentj_nmismatches) <= best_nmismatches) {
	      if (donori_knowni[i] >= 0) {
		probi = 1.0; /* Needs to be 1.0 for output */
	      } else {
		probi = Maxent_hr_donor_prob(segmenti_left + splice_pos_1,segmenti_chroffset);
	      }

	      if (acceptora_knowni[a] >= 0) {
		proba = 1.0; /* Needs to be 1.0 for output */
	      } else {
		proba = Maxent_hr_acceptor_prob(segmentm_left + splice_pos_1,segmentm_chroffset);
	      }

	      if (donorb_knowni[b] >= 0) {
		probb = 1.0; /* Needs to be 1.0 for output */
	      } else {
		probb = Maxent_hr_donor_prob(segmentm_left + splice_pos_2,segmentm_chroffset);
	      }
	      
	      if (acceptorj_knowni[j] >= 0) {
		probj = 1.0; /* Needs to be 1.0 for output */
	      } else {
		probj = Maxent_hr_acceptor_prob(segmentj_left + splice_pos_2,segmentj_chroffset);
	      }

	      debug2(
		     if (plusp == true) {
		       printf("plus sense splice_pos  %d, %d, i.donor %f, m.acceptor %f, m.donor %f, j.acceptor %f\n",
			      splice_pos_1,splice_pos_2,probi,proba,probb,probj);
		     } else {
		       printf("minus antisense splice_pos  %d %d, i.donor %f, m.acceptor %f, m.donor %f, j.acceptor %f\n",
			      splice_pos_1,splice_pos_2,probi,proba,probb,probj);
		     });

	      if (nmismatches < best_nmismatches ||
		  (nmismatches == best_nmismatches && probi + proba + probb + probj > best_prob)) {
		/* Success */
		best_nmismatches = nmismatches;
		best_prob = probi + proba + probb + probj;

		best_donor1_splicecoord = segmenti_left + splice_pos_1;
		best_acceptor1_splicecoord = segmentm_left + splice_pos_1;
		best_donor2_splicecoord = segmentm_left + splice_pos_2;
		best_acceptor2_splicecoord = segmentj_left + splice_pos_2;
		best_donor1_knowni = donori_knowni[i];
		best_acceptor1_knowni = acceptora_knowni[a];
		best_donor2_knowni = donorb_knowni[b];
		best_acceptor2_knowni = acceptorj_knowni[j];
		best_donor1_prob = probi;
		best_acceptor1_prob = proba;
		best_donor2_prob = probb;
		best_acceptor2_prob = probj;
		best_splice_pos_1 = splice_pos_1;
		best_splice_pos_2 = splice_pos_2;
		best_segmenti_nmismatches = segmenti_nmismatches;
		best_segmentm_nmismatches = segmentm_nmismatches;
		best_segmentj_nmismatches = segmentj_nmismatches;
	      }
	    }
	    /* b++; j++; Don't advance b or j, so next i/a can match */
	    matchp = true;
	  }
	}
	i++;
	a++;
      }
    }


    /* Originally from minus strand.  Complement. */
    /* Antisense (End 7 to End 8) or Sense (End 3 to End 4) */

    /* Segment i */
    if (novelsplicingp && segmenti_left + splice_pos_start >= ACCEPTOR_MODEL_RIGHT_MARGIN) {
      antiacceptori_nsites = Genome_antiacceptor_positions(acceptor1_positions_alloc,acceptor1_knowni_alloc,
							   segmenti_antiacceptor_knownpos,segmenti_antiacceptor_knowni,
							   segmenti_left,splice_pos_start,splice_pos_end);
      antiacceptori_positions = acceptor1_positions_alloc;
      antiacceptori_knowni = acceptor1_knowni_alloc;
    } else {
      antiacceptori_nsites = segmenti_antiacceptor_nknown;
      antiacceptori_positions = segmenti_antiacceptor_knownpos;
      antiacceptori_knowni = segmenti_antiacceptor_knowni;
    }

#ifdef DEBUG2
    printf("Found %d antiacceptori sites:",antiacceptori_nsites);
    for (i = 0; i < antiacceptori_nsites; i++) {
      printf(" %d",antiacceptori_positions[i]);
      if (antiacceptori_knowni[i] >= 0) {
	printf(" (%d)",antiacceptori_knowni[i]);
      }
    }
    printf("\n");
#endif

    /* Segment m1 */
    if (novelsplicingp && segmentm_left + splice_pos_start >= DONOR_MODEL_RIGHT_MARGIN) {
      antidonora_nsites = Genome_antidonor_positions(donor1_positions_alloc,donor1_knowni_alloc,
						     segmentm_antidonor_knownpos,segmentm_antidonor_knowni,
						     segmentm_left,splice_pos_start,splice_pos_end);
      antidonora_positions = donor1_positions_alloc;
      antidonora_knowni = donor1_knowni_alloc;
    } else {
      antidonora_nsites = segmentm_antidonor_nknown;
      antidonora_positions = segmentm_antidonor_knownpos;
      antidonora_knowni = segmentm_antidonor_knowni;
    }

#ifdef DEBUG2
    printf("Found %d antidonora sites:",antidonora_nsites);
    for (i = 0; i < antidonora_nsites; i++) {
      printf(" %d",antidonora_positions[i]);
      if (antidonora_knowni[i] >= 0) {
	printf(" (%d)",antidonora_knowni[i]);
      }
    }
    printf("\n");
#endif

    /* Segment m2 */
    if (novelsplicingp && segmentm_left + splice_pos_start >= ACCEPTOR_MODEL_RIGHT_MARGIN) {
      antiacceptorb_nsites = Genome_antiacceptor_positions(acceptor2_positions_alloc,acceptor2_knowni_alloc,
							   segmentm_antiacceptor_knownpos,segmentm_antiacceptor_knowni,
							   segmentm_left,splice_pos_start,splice_pos_end);
      antiacceptorb_positions = acceptor2_positions_alloc;
      antiacceptorb_knowni = acceptor2_knowni_alloc;
    } else {
      antiacceptorb_nsites = segmentm_antiacceptor_nknown;
      antiacceptorb_positions = segmentm_antiacceptor_knownpos;
      antiacceptorb_knowni = segmentm_antiacceptor_knowni;
    }

#ifdef DEBUG2
    printf("Found %d antiacceptorb sites:",antiacceptorb_nsites);
    for (i = 0; i < antiacceptorb_nsites; i++) {
      printf(" %d",antiacceptorb_positions[i]);
      if (antiacceptorb_knowni[i] >= 0) {
	printf(" (%d)",antiacceptorb_knowni[i]);
      }
    }
    printf("\n");
#endif

    /* Segment j */
    if (novelsplicingp && segmentj_left + splice_pos_start >= DONOR_MODEL_RIGHT_MARGIN) {
      antidonorj_nsites = Genome_antidonor_positions(donor2_positions_alloc,donor2_knowni_alloc,
						     segmentj_antidonor_knownpos,segmentj_antidonor_knowni,
						     segmentj_left,splice_pos_start,splice_pos_end);
      antidonorj_positions = donor2_positions_alloc;
      antidonorj_knowni = donor2_knowni_alloc;
    } else {
      antidonorj_nsites = segmentj_antidonor_nknown;
      antidonorj_positions = segmentj_antidonor_knownpos;
      antidonorj_knowni = segmentj_antidonor_knowni;
    }

#ifdef DEBUG2
    printf("Found %d antidonorj sites:",antidonorj_nsites);
    for (i = 0; i < antidonorj_nsites; i++) {
      printf(" %d",antidonorj_positions[i]);
      if (antidonorj_knowni[i] >= 0) {
	printf(" (%d)",antidonorj_knowni[i]);
      }
    }
    printf("\n");
#endif


    i = a = b = j = 0;
    while (i < antiacceptori_nsites && a < antidonora_nsites) {
      if ((splice_pos_1 = antiacceptori_positions[i]) < antidonora_positions[a]) {
	i++;
      } else if (splice_pos_1 > antidonora_positions[a]) {
	a++;
      } else {
	while (b < antiacceptorb_nsites && antiacceptorb_positions[b] <= splice_pos_1) {
	  b++;
	}
	while (j < antidonorj_nsites && antidonorj_positions[j] <= splice_pos_1) {
	  j++;
	}
	matchp = false;
	while (b < antiacceptorb_nsites && j < antidonorj_nsites && matchp == false) {
	  if ((splice_pos_2 = antiacceptorb_positions[b]) < antidonorj_positions[j]) {
	    b++;
	  } else if (splice_pos_2 > antidonorj_positions[j]) {
	    j++;
	  } else {
	    segmenti_nmismatches = Genome_count_mismatches_substring(query_compress,/*left*/segmenti_left,/*pos5*/0,/*pos3*/splice_pos_1,
								     plusp,genestrand,first_read_p);
	    segmentm_nmismatches = Genome_count_mismatches_substring(query_compress,/*left*/segmentm_left,/*pos5*/splice_pos_1,/*pos3*/splice_pos_2,
								     plusp,genestrand,first_read_p);
	    segmentj_nmismatches = Genome_count_mismatches_substring(query_compress,/*left*/segmentj_left,/*pos5*/splice_pos_2,/*pos3*/querylength,
								     plusp,genestrand,first_read_p);
	    
	    if ((nmismatches = segmenti_nmismatches + segmentm_nmismatches + segmentj_nmismatches) <= best_nmismatches) {
	      if (antiacceptori_knowni[i] >= 0) {
		probi = 1.0; /* Needs to be 1.0 for output */
	      } else {
		probi = Maxent_hr_antiacceptor_prob(segmenti_left + splice_pos_1,segmenti_chroffset);
	      }
	    
	      if (antidonora_knowni[a] >= 0) {
		proba = 1.0; /* Needs to be 1.0 for output */
	      } else {
		proba = Maxent_hr_antidonor_prob(segmentm_left + splice_pos_1,segmentm_chroffset);
	      }

	      if (antiacceptorb_knowni[b] >= 0) {
		probb = 1.0; /* Needs to be 1.0 for output */
	      } else {
		probb = Maxent_hr_antiacceptor_prob(segmentm_left + splice_pos_2,segmentm_chroffset);
	      }

	      if (antidonorj_knowni[j] >= 0) {
		probj = 1.0; /* Needs to be 1.0 for output */
	      } else {
		probj = Maxent_hr_antidonor_prob(segmentj_left + splice_pos_2,segmentj_chroffset);
	      }

	      debug2(
		     if (plusp == true) {
		       printf("plus antisense splice_pos  %d, %d, i.antiacceptor %f, m.antidonor %f, m.antiacceptor %f, j.antidonor %f\n",
			      splice_pos_1,splice_pos_2,probi,proba,probb,probj);
		     } else {
		       printf("minus sense splice_pos  %d, %d, i.antiacceptor %f, m.antidonor %f, m.antiacceptor %f, j.antidonor %f\n",
			      splice_pos_1,splice_pos_2,probi,proba,probb,probj);
		     });

	      if (nmismatches < best_nmismatches ||
		  (nmismatches == best_nmismatches && probi + proba + probb + probj > best_prob)) {
		/* Success */
		best_nmismatches = nmismatches;
		best_prob = probi + proba + probb + probj;

		best_acceptor1_splicecoord = segmenti_left + splice_pos_1;
		best_donor1_splicecoord = segmentm_left + splice_pos_1;
		best_acceptor2_splicecoord = segmentm_left + splice_pos_2;
		best_donor2_splicecoord = segmentj_left + splice_pos_2;
		best_acceptor1_knowni = antiacceptori_knowni[i];
		best_donor1_knowni = antidonora_knowni[a];
		best_acceptor2_knowni = antiacceptorb_knowni[b];
		best_donor2_knowni = antidonorj_knowni[j];
		best_acceptor1_prob = probi;
		best_donor1_prob = proba;
		best_acceptor2_prob = probb;
		best_donor2_prob = probj;
		best_splice_pos_1 = splice_pos_1;
		best_splice_pos_2 = splice_pos_2;
		best_segmenti_nmismatches = segmenti_nmismatches;
		best_segmentm_nmismatches = segmentm_nmismatches;
		best_segmentj_nmismatches = segmentj_nmismatches;
		orig_plusp = false;
	      }
	    }
	    /* b++; j++; Don't advance b or j, so next i/a can match */
	    matchp = true;
	  }
	}
	i++;
	a++;
      }
    }


    if (best_prob > 0.0) {
      debug2(printf("best_prob = %f at splice_pos %d and %d\n",best_prob,best_splice_pos_1,best_splice_pos_2));
      if (orig_plusp == true) {
	/* Originally from plus strand.  No complement. */
	sensedir = (plusp == true) ? SENSE_FORWARD : SENSE_ANTI;

	donor = Substring_new_donor(best_donor1_splicecoord,best_donor1_knowni,
				    best_splice_pos_1,best_segmenti_nmismatches,
				    best_donor1_prob,/*left*/segmenti_left,query_compress,
				    querylength,plusp,genestrand,first_read_p,sensedir,
				    segmenti_chrnum,segmenti_chroffset,segmenti_chrhigh,segmenti_chrlength);

	shortexon = Substring_new_shortexon(best_acceptor1_splicecoord,best_acceptor1_knowni,
					    best_donor2_splicecoord,best_donor2_knowni,
					    /*acceptor_pos*/best_splice_pos_1,/*donor_pos*/best_splice_pos_2,best_segmentm_nmismatches,
					    /*acceptor_prob*/best_acceptor1_prob,/*donor_prob*/best_donor2_prob,
					    /*left*/segmentm_left,query_compress,
					    querylength,plusp,genestrand,first_read_p,
					    sensedir,/*acceptor_ambp*/false,/*donor_ambp*/false,
					    segmentm_chrnum,segmentm_chroffset,segmentm_chrhigh,segmentm_chrlength);

	acceptor = Substring_new_acceptor(best_acceptor2_splicecoord,best_acceptor2_knowni,
					  best_splice_pos_2,best_segmentj_nmismatches,
					  best_acceptor2_prob,/*left*/segmentj_left,query_compress,
					  querylength,plusp,genestrand,first_read_p,sensedir,
					  segmentj_chrnum,segmentj_chroffset,segmentj_chrhigh,segmentj_chrlength);

	if (donor == NULL || shortexon == NULL || acceptor == NULL) {
	  if (donor != NULL) Substring_free(&donor);
	  if (shortexon != NULL) Substring_free(&shortexon);
	  if (acceptor != NULL) Substring_free(&acceptor);
	} else {
	  *segmenti_usedp = *segmentm_usedp = *segmentj_usedp = true;

	  donor_support = best_splice_pos_1;
	  middle_support = best_splice_pos_2 - best_splice_pos_1;
	  acceptor_support = querylength - best_splice_pos_2;
	  sufficient1p = sufficient_splice_prob_local(donor_support,best_segmenti_nmismatches,best_donor1_prob);
	  sufficient2p = sufficient_splice_prob_local(middle_support,best_segmentm_nmismatches,best_acceptor1_prob);
	  sufficient3p = sufficient_splice_prob_local(middle_support,best_segmentm_nmismatches,best_donor2_prob);
	  sufficient4p = sufficient_splice_prob_local(acceptor_support,best_segmentj_nmismatches,best_acceptor2_prob);
	  if (sufficient1p && sufficient2p && sufficient3p && sufficient4p) {
	    *nhits += 1;
	    hits = List_push(hits,(void *) Stage3end_new_shortexon(&(*found_score),donor,acceptor,shortexon,
								   best_donor1_prob,/*shortexonA_prob*/best_acceptor1_prob,
								   /*shortexonD_prob*/best_donor2_prob,best_acceptor2_prob,
								   /*amb_length_donor*/0,/*amb_length_acceptor*/0,
								   /*amb_prob_donor*/0.0,/*amb_prob_acceptor*/0.0,
								   /*ambcoords_donor*/NULL,/*ambcoords_acceptor*/NULL,
								   /*amb_knowni_donor*/NULL,/*amb_knowni_acceptor*/NULL,
								   /*amb_nmismatches_donor*/NULL,/*amb_nmismatches_acceptor*/NULL,
								   /*amb_probs_donor*/NULL,/*amb_probs_acceptor*/NULL,
								   /*copy_donor_p*/false,/*copy_acceptor_p*/false,/*copy_shortexon_p*/false,
								   splicing_penalty,querylength,first_read_p,sensedir,sarrayp));
	  } else if (subs_or_indels_p == true) {
	    /* Don't alter hits */
	    if (donor != NULL) Substring_free(&donor);
	    if (shortexon != NULL) Substring_free(&shortexon);
	    if (acceptor != NULL) Substring_free(&acceptor);
	  } else if (donor_support < LOWPROB_SUPPORT || acceptor_support < LOWPROB_SUPPORT) {
	    if (donor != NULL) Substring_free(&donor);
	    if (shortexon != NULL) Substring_free(&shortexon);
	    if (acceptor != NULL) Substring_free(&acceptor);
	  } else if ((sufficient1p || sufficient2p) && (sufficient3p || sufficient4p)) {
	    *lowprob = List_push(*lowprob,
				 (void *) Stage3end_new_shortexon(&(*found_score),donor,acceptor,shortexon,
								   best_donor1_prob,/*shortexonA_prob*/best_acceptor1_prob,
								   /*shortexonD_prob*/best_donor2_prob,best_acceptor2_prob,
								  /*amb_length_donor*/0,/*amb_length_acceptor*/0,
								  /*amb_prob_donor*/0.0,/*amb_prob_acceptor*/0.0,
								  /*ambcoords_donor*/NULL,/*ambcoords_acceptor*/NULL,
								  /*amb_knowni_donor*/NULL,/*amb_knowni_acceptor*/NULL,
								  /*amb_nmismatches_donor*/NULL,/*amb_nmismatches_acceptor*/NULL,
								  /*amb_probs_donor*/NULL,/*amb_probs_acceptor*/NULL,
								  /*copy_donor_p*/false,/*copy_acceptor_p*/false,/*copy_shortexon_p*/false,
								  splicing_penalty,querylength,first_read_p,sensedir,sarrayp));
	  } else {
	    if (donor != NULL) Substring_free(&donor);
	    if (shortexon != NULL) Substring_free(&shortexon);
	    if (acceptor != NULL) Substring_free(&acceptor);
	  }
	}

      } else {
	/* Originally from minus strand.  Complement. */
	sensedir = (plusp == true) ? SENSE_ANTI : SENSE_FORWARD;

	donor = Substring_new_donor(best_donor2_splicecoord,best_donor2_knowni,
				    best_splice_pos_2,best_segmentj_nmismatches,
				    best_donor2_prob,/*left*/segmentj_left,query_compress,
				    querylength,plusp,genestrand,first_read_p,sensedir,
				    segmentj_chrnum,segmentj_chroffset,segmentj_chrhigh,segmentj_chrlength);

	shortexon = Substring_new_shortexon(best_acceptor2_splicecoord,best_acceptor2_knowni,
					    best_donor1_splicecoord,best_donor1_knowni,
					    /*acceptor_pos*/best_splice_pos_2,/*donor_pos*/best_splice_pos_1,best_segmentm_nmismatches,
					    /*acceptor_prob*/best_acceptor2_prob,/*donor_prob*/best_donor1_prob,
					    /*left*/segmentm_left,query_compress,querylength,
					    plusp,genestrand,first_read_p,sensedir,/*acceptor_ambp*/false,/*donor_ambp*/false,
					    segmentm_chrnum,segmentm_chroffset,segmentm_chrhigh,segmentm_chrlength);

	acceptor = Substring_new_acceptor(best_acceptor1_splicecoord,best_acceptor1_knowni,
					  best_splice_pos_1,best_segmenti_nmismatches,
					  best_acceptor1_prob,/*left*/segmenti_left,query_compress,
					  querylength,plusp,genestrand,first_read_p,sensedir,
					  segmenti_chrnum,segmenti_chroffset,segmenti_chrhigh,segmenti_chrlength);

	if (donor == NULL || shortexon == NULL || acceptor == NULL) {
	  if (donor != NULL) Substring_free(&donor);
	  if (shortexon != NULL) Substring_free(&shortexon);
	  if (acceptor != NULL) Substring_free(&acceptor);
	} else {
	  *segmenti_usedp = *segmentm_usedp = *segmentj_usedp = true;

	  acceptor_support = best_splice_pos_1;
	  middle_support = best_splice_pos_2 - best_splice_pos_1;
	  donor_support = querylength - best_splice_pos_2;
	  sufficient1p = sufficient_splice_prob_local(acceptor_support,best_segmenti_nmismatches,best_acceptor1_prob);
	  sufficient2p = sufficient_splice_prob_local(middle_support,best_segmentm_nmismatches,best_donor1_prob);
	  sufficient3p = sufficient_splice_prob_local(middle_support,best_segmentm_nmismatches,best_acceptor2_prob);
	  sufficient4p = sufficient_splice_prob_local(donor_support,best_segmentj_nmismatches,best_donor2_prob);
	  if (sufficient1p && sufficient2p && sufficient3p && sufficient4p) {
	    *nhits += 1;
	    hits = List_push(hits,(void *) Stage3end_new_shortexon(&(*found_score),donor,acceptor,shortexon,
								   best_donor2_prob,/*shortexonA_prob*/best_acceptor2_prob,
								   /*shortexonD_prob*/best_donor1_prob,best_acceptor1_prob,
								   /*amb_length_donor*/0,/*amb_length_acceptor*/0,
								   /*amb_prob_donor*/0.0,/*amb_prob_acceptor*/0.0,
								   /*ambcoords_donor*/NULL,/*ambcoords_acceptor*/NULL,
								   /*amb_knowni_donor*/NULL,/*amb_knowni_acceptor*/NULL,
								   /*amb_nmismatches_donor*/NULL,/*amb_nmismatches_acceptor*/NULL,
								   /*amb_probs_donor*/NULL,/*amb_probs_acceptor*/NULL,
								   /*copy_donor_p*/false,/*copy_acceptor_p*/false,/*copy_shortexon_p*/false,
								   splicing_penalty,querylength,first_read_p,sensedir,sarrayp));
	  } else if (subs_or_indels_p == true) {
	    /* Don't alter hits */
	    if (donor != NULL) Substring_free(&donor);
	    if (shortexon != NULL) Substring_free(&shortexon);
	    if (acceptor != NULL) Substring_free(&acceptor);
	  } else if (donor_support < LOWPROB_SUPPORT || acceptor_support < LOWPROB_SUPPORT) {
	    if (donor != NULL) Substring_free(&donor);
	    if (shortexon != NULL) Substring_free(&shortexon);
	    if (acceptor != NULL) Substring_free(&acceptor);
	  } else if ((sufficient1p || sufficient2p) && (sufficient3p || sufficient4p)) {
	    *lowprob = List_push(*lowprob,
				 (void *) Stage3end_new_shortexon(&(*found_score),donor,acceptor,shortexon,
								   best_donor2_prob,/*shortexonA_prob*/best_acceptor2_prob,
								   /*shortexonD_prob*/best_donor1_prob,best_acceptor1_prob,
								  /*amb_length_donor*/0,/*amb_length_acceptor*/0,
								  /*amb_prob_donor*/0.0,/*amb_prob_acceptor*/0.0,
								  /*ambcoords_donor*/NULL,/*ambcoords_acceptor*/NULL,
								  /*amb_knowni_donor*/NULL,/*amb_knowni_acceptor*/NULL,
								  /*amb_nmismatches_donor*/NULL,/*amb_nmismatches_acceptor*/NULL,
								  /*amb_probs_donor*/NULL,/*amb_probs_acceptor*/NULL,
								  /*copy_donor_p*/false,/*copy_acceptor_p*/false,/*copy_shortexon_p*/false,
								  splicing_penalty,querylength,first_read_p,sensedir,sarrayp));
	  } else {
	    if (donor != NULL) Substring_free(&donor);
	    if (shortexon != NULL) Substring_free(&shortexon);
	    if (acceptor != NULL) Substring_free(&acceptor);
	  }
	}
      }
    }
  }

  return hits;
}


static int
donor_match_length_cmp (const void *a, const void *b) {
  Stage3end_T x = * (Stage3end_T *) a;
  Stage3end_T y = * (Stage3end_T *) b;
  
  int x_length = Substring_match_length_orig(Stage3end_substring_donor(x));
  int y_length = Substring_match_length_orig(Stage3end_substring_donor(y));

  if (x_length < y_length) {
    return -1;
  } else if (y_length < x_length) {
    return +1;
  } else {
    return 0;
  }
}

static int
acceptor_match_length_cmp (const void *a, const void *b) {
  Stage3end_T x = * (Stage3end_T *) a;
  Stage3end_T y = * (Stage3end_T *) b;
  
  int x_length = Substring_match_length_orig(Stage3end_substring_acceptor(x));
  int y_length = Substring_match_length_orig(Stage3end_substring_acceptor(y));

  if (x_length < y_length) {
    return -1;
  } else if (y_length < x_length) {
    return +1;
  } else {
    return 0;
  }
}


static List_T
group_by_segmenti_aux (int *found_score, List_T winners, List_T *ambiguous,
		       Stage3end_T *hitarray, int n, int querylength, bool first_read_p, bool sarrayp) {
  Stage3end_T hit, *subarray;
  int i, j, k, ii, jj, kk, nn;
  int n_good_spliceends;
  Univcoord_T segmenti_left;
  Substring_T donor, acceptor;
  int best_nmismatches, nmismatches, nmismatches_donor, nmismatches_acceptor;
  double best_prob, prob, donor_prob, acceptor_prob;
  List_T accepted_hits, rejected_hits, donor_hits, acceptor_hits, p;

  int sensedir;
#ifdef LARGE_GENOMES
  Uint8list_T ambcoords;
#else
  Uintlist_T ambcoords;
#endif
  Intlist_T amb_knowni, amb_nmismatches;
  Doublelist_T amb_probs;
  int donor_length, acceptor_length;
  bool plusp;

  i = 0;
  while (i < n) {
    hit = hitarray[i];
    segmenti_left = Stage3end_chimera_segmenti_left(hit);
    j = i + 1;
    while (j < n && Stage3end_chimera_segmenti_left(hitarray[j]) == segmenti_left) {
      j++;
    }
    if (j == i + 1) {
      /* Singleton */
      debug7(printf("Saving hit %d\n",i));
      winners = List_push(winners,(void *) hit);

    } else {
      plusp = Stage3end_plusp(hit);
      best_nmismatches = querylength;
      best_prob = 0.0;
      for (k = i; k < j; k++) {
	hit = hitarray[k];
	debug7(printf("analyzing distance %d, donor length %d (%llu..%llu) and acceptor length %d (%llu..%llu), nmismatches %d, probabilities %f and %f\n",
		      Stage3end_distance(hit),Substring_match_length_orig(Stage3end_substring_donor(hit)),
		      Substring_genomicstart(Stage3end_substring_donor(hit)),Substring_genomicend(Stage3end_substring_donor(hit)),
		      Substring_match_length_orig(Stage3end_substring_acceptor(hit)),
		      Substring_genomicstart(Stage3end_substring_acceptor(hit)),Substring_genomicend(Stage3end_substring_acceptor(hit)),
		      Stage3end_nmismatches_whole(hit),Substring_chimera_prob(Stage3end_substring_donor(hit)),
		      Substring_chimera_prob(Stage3end_substring_acceptor(hit))));
	if ((nmismatches = Stage3end_nmismatches_whole(hit)) < best_nmismatches) {
	  best_nmismatches = nmismatches;
	}
	if ((prob = Stage3end_chimera_prob(hit)) > best_prob) {
	  best_prob = prob;
	}
      }

      n_good_spliceends = 0;
      accepted_hits = rejected_hits = (List_T) NULL;
      for (k = i; k < j; k++) {
	hit = hitarray[k];
	if (Stage3end_nmismatches_whole(hit) <= best_nmismatches + LOCALSPLICING_NMATCHES_SLOP &&
	    Stage3end_chimera_prob(hit) >= best_prob - LOCALSPLICING_PROB_SLOP) {
	  debug7(printf("accepting distance %d, probabilities %f and %f\n",
			Stage3end_distance(hit),Substring_chimera_prob(Stage3end_substring_donor(hit)),
			Substring_chimera_prob(Stage3end_substring_acceptor(hit))));
	  n_good_spliceends += 1;
	  accepted_hits = List_push(accepted_hits,(void *) hit);
	} else {
	  rejected_hits = List_push(rejected_hits,(void *) hit);
	}
      }

      if (n_good_spliceends == 0) {
	/* Conjunction is too strict.  Allow for disjunction instead. */
	List_free(&rejected_hits);
	for (k = i; k < j; k++) {
	  hit = hitarray[k];
	  if (Stage3end_nmismatches_whole(hit) <= best_nmismatches + LOCALSPLICING_NMATCHES_SLOP ||
	      Stage3end_chimera_prob(hit) >= best_prob - LOCALSPLICING_PROB_SLOP) {
	    debug7(printf("accepting distance %d, probabilities %f and %f\n",
			  Stage3end_distance(hit),Substring_chimera_prob(Stage3end_substring_donor(hit)),
			  Substring_chimera_prob(Stage3end_substring_acceptor(hit))));
	    n_good_spliceends += 1;
	    accepted_hits = List_push(accepted_hits,(void *) hit);
	  } else {
	    rejected_hits = List_push(rejected_hits,(void *) hit);
	  }
	}
      }
	
      for (p = rejected_hits; p != NULL; p = List_next(p)) {
	hit = (Stage3end_T) List_head(p);
	Stage3end_free(&hit);
      }
      List_free(&rejected_hits);

      if (n_good_spliceends == 1) {
	winners = List_push(winners,List_head(accepted_hits));
	List_free(&accepted_hits);

      } else {
	/* Multiple hits */
	donor_hits = acceptor_hits = (List_T) NULL;
	for (p = accepted_hits; p != NULL; p = List_next(p)) {
	  hit = (Stage3end_T) List_head(p);
	  donor = Stage3end_substring_donor(hit);
	  acceptor = Stage3end_substring_acceptor(hit);
	  if (Stage3end_plusp(hit) == true) {
	    if (Substring_genomicstart(donor) == segmenti_left) {
	      donor_hits = List_push(donor_hits,(void *) hit);
	    } else if (Substring_genomicstart(acceptor) == segmenti_left) {
	      acceptor_hits = List_push(acceptor_hits,(void *) hit);
	    } else {
	      Stage3end_free(&hit);
	    }
	  } else {
	    if (Substring_genomicend(donor) == segmenti_left) {
	      donor_hits = List_push(donor_hits,(void *) hit);
	    } else if (Substring_genomicend(acceptor) == segmenti_left) {
	      acceptor_hits = List_push(acceptor_hits,(void *) hit);
	    } else {
	      Stage3end_free(&hit);
	    }
	  }
	}

	if (donor_hits != NULL) {
	  subarray = (Stage3end_T *) List_to_array_n(&nn,donor_hits);
	  qsort(subarray,nn,sizeof(Stage3end_T),donor_match_length_cmp);
	  ii = 0;
	  while (ii < nn) {
	    hit = subarray[ii];
	    donor = Stage3end_substring_donor(hit);
	    donor_length = Substring_match_length_orig(donor);
	    jj = ii + 1;
	    while (jj < nn && Substring_match_length_orig(Stage3end_substring_donor(subarray[jj])) == donor_length) {
	      jj++;
	    }
	    if (jj == ii + 1) {
	      winners = List_push(winners,(void *) hit);
	    } else {
	      sensedir = Stage3end_sensedir(hit);

	      ambcoords = NULL;
	      amb_knowni = (Intlist_T) NULL;
	      amb_nmismatches = (Intlist_T) NULL;
	      amb_probs = (Doublelist_T) NULL;

	      for (kk = ii; kk < jj; kk++) {
		acceptor = Stage3end_substring_acceptor(subarray[kk]);
#ifdef LARGE_GENOMES
		ambcoords = Uint8list_push(ambcoords,Substring_splicecoord(acceptor));
#else
		ambcoords = Uintlist_push(ambcoords,Substring_splicecoord(acceptor));
#endif
		amb_knowni = Intlist_push(amb_knowni,-1);
		amb_nmismatches = Intlist_push(amb_nmismatches,Substring_nmismatches_whole(acceptor));
		amb_probs = Doublelist_push(amb_probs,Substring_chimera_prob(acceptor));
	      }

	      nmismatches_acceptor = best_nmismatches - Substring_nmismatches_whole(donor);
	      donor_prob = Junction_donor_prob(Stage3end_junctionA(hit));
	      prob = best_prob - donor_prob;
	      *ambiguous = List_push(*ambiguous,
				     (void *) Stage3end_new_splice(&(*found_score),
								   /*nmismatches_donor*/Substring_nmismatches_whole(donor),nmismatches_acceptor,
								   donor,/*acceptor*/NULL,donor_prob,/*acceptor_prob*/prob,/*distance*/0U,
								   /*shortdistancep*/false,/*penalty*/0,querylength,
								   /*amb_length*/Substring_match_length_orig(acceptor),/*amb_prob*/prob,
								   /*ambcoords_donor*/NULL,ambcoords,
								   /*amb_knowni_donor*/NULL,amb_knowni,
								   /*amb_nmismatches_donor*/NULL,amb_nmismatches,
								   /*amb_probs_donor*/NULL,amb_probs,
								   /*copy_donor_p*/true,/*copy_acceptor_p*/false,first_read_p,
								   sensedir,sarrayp));
	      Doublelist_free(&amb_probs);
	      Intlist_free(&amb_knowni);
	      Intlist_free(&amb_nmismatches);
#ifdef LARGE_GENOMES
	      Uint8list_free(&ambcoords);
#else
	      Uintlist_free(&ambcoords);
#endif
	      for (kk = ii; kk < jj; kk++) {
		hit = subarray[kk];
		Stage3end_free(&hit);
	      }
	    }

	    ii = jj;
	  }
	  FREE(subarray);
	  List_free(&donor_hits);
	}

	if (acceptor_hits != NULL) {
	  subarray = (Stage3end_T *) List_to_array_n(&nn,acceptor_hits);
	  qsort(subarray,nn,sizeof(Stage3end_T),acceptor_match_length_cmp);
	  ii = 0;
	  while (ii < nn) {
	    hit = subarray[ii];
	    acceptor = Stage3end_substring_acceptor(hit);
	    acceptor_length = Substring_match_length_orig(acceptor);
	    jj = ii + 1;
	    while (jj < nn && Substring_match_length_orig(Stage3end_substring_acceptor(subarray[jj])) == acceptor_length) {
	      jj++;
	    }
	    if (jj == ii + 1) {
	      winners = List_push(winners,(void *) hit);
	    } else {
	      sensedir = Stage3end_sensedir(hit);

	      ambcoords = NULL;
	      amb_knowni = (Intlist_T) NULL;
	      amb_nmismatches = (Intlist_T) NULL;
	      amb_probs = (Doublelist_T) NULL;

	      for (kk = ii; kk < jj; kk++) {
		donor = Stage3end_substring_donor(subarray[kk]);
#ifdef LARGE_GENOMES
		ambcoords = Uint8list_push(ambcoords,Substring_splicecoord(donor));
#else
		ambcoords = Uintlist_push(ambcoords,Substring_splicecoord(donor));
#endif
		amb_knowni = Intlist_push(amb_knowni,-1);
		amb_nmismatches = Intlist_push(amb_nmismatches,Substring_nmismatches_whole(donor));
		amb_probs = Doublelist_push(amb_probs,Substring_chimera_prob(donor));
	      }

	      nmismatches_donor = best_nmismatches - Substring_nmismatches_whole(acceptor);
	      acceptor_prob = Junction_acceptor_prob(Stage3end_junctionD(hit));
	      prob = best_prob - acceptor_prob;
	      *ambiguous = List_push(*ambiguous,
				     (void *) Stage3end_new_splice(&(*found_score),
								   nmismatches_donor,/*nmismatches_acceptor*/Substring_nmismatches_whole(acceptor),
								   /*donor*/NULL,acceptor,/*donor_prob*/prob,acceptor_prob,/*distance*/0U,
								   /*shortdistancep*/false,/*penalty*/0,querylength,
								   /*amb_length*/Substring_match_length_orig(donor),/*amb_prob*/prob,
								   ambcoords,/*ambcoords_acceptor*/NULL,
								   amb_knowni,/*amb_knowni_acceptor*/NULL,
								   amb_nmismatches,/*amb_nmismatches_acceptor*/NULL,
								   amb_probs,/*amb_probs_acceptor*/NULL,
								   /*copy_donor_p*/false,/*copy_acceptor_p*/true,first_read_p,
								   sensedir,sarrayp));
	      Doublelist_free(&amb_probs);
	      Intlist_free(&amb_knowni);
	      Intlist_free(&amb_nmismatches);
#ifdef LARGE_GENOMES
	      Uint8list_free(&ambcoords);
#else
	      Uintlist_free(&ambcoords);
#endif
	      for (kk = ii; kk < jj; kk++) {
		hit = subarray[kk];
		Stage3end_free(&hit);
	      }
	    }

	    ii = jj;
	  }
	  FREE(subarray);
	  List_free(&acceptor_hits);
	}

	List_free(&accepted_hits);
      }
    }

    i = j;
  }

  return winners;
}

List_T
Splice_group_by_segmenti (int *found_score, List_T localsplicing, List_T *ambiguous,
			  int querylength, bool first_read_p, bool sarrayp) {
  List_T winners = NULL, p;
  Stage3end_T *array_forward, *array_anti, hit;
  int n_sense_forward = 0, n_sense_anti = 0, k_forward, k_anti;

  for (p = localsplicing; p != NULL; p = List_next(p)) {
    hit = (Stage3end_T) List_head(p);
    if (Stage3end_sensedir(hit) == SENSE_FORWARD) {
      n_sense_forward++;
    } else {
      assert(Stage3end_sensedir(hit) == SENSE_ANTI);
      n_sense_anti++;
    }
  }

  if (n_sense_forward > 0) {
    array_forward = (Stage3end_T *) MALLOCA(n_sense_forward * sizeof(Stage3end_T));
    k_forward = 0;
  }
  if (n_sense_anti > 0) {
    array_anti = (Stage3end_T *) MALLOCA(n_sense_anti * sizeof(Stage3end_T));
    k_anti = 0;
  }

  for (p = localsplicing; p != NULL; p = List_next(p)) {
    hit = (Stage3end_T) List_head(p);
    if (Stage3end_sensedir(hit) == SENSE_FORWARD) {
      array_forward[k_forward++] = (Stage3end_T) List_head(p);
    } else {
      array_anti[k_anti++] = (Stage3end_T) List_head(p);
    }
  }

  if (n_sense_forward > 0) {
    qsort(array_forward,n_sense_forward,sizeof(Stage3end_T),Stage3end_chimera_segmenti_cmp);
    winners = group_by_segmenti_aux(&(*found_score),winners,&(*ambiguous),array_forward,n_sense_forward,
				    querylength,first_read_p,sarrayp);
    FREEA(array_forward);
  }

  if (n_sense_anti > 0) {
    qsort(array_anti,n_sense_anti,sizeof(Stage3end_T),Stage3end_chimera_segmenti_cmp);
    winners = group_by_segmenti_aux(&(*found_score),winners,&(*ambiguous),array_anti,n_sense_anti,
				    querylength,first_read_p,sarrayp);
    FREEA(array_anti);
  }

  List_free(&localsplicing);

  return winners;
}




static List_T
group_by_segmentj_aux (int *found_score, List_T winners, List_T *ambiguous, 
		       Stage3end_T *hitarray, int n, int querylength, bool first_read_p, bool sarrayp) {
  Stage3end_T hit, *subarray;
  int i, j, k, ii, jj, kk, nn;
  int n_good_spliceends;
  Univcoord_T segmentj_left;
  Substring_T donor, acceptor;
  int best_nmismatches, nmismatches, nmismatches_donor, nmismatches_acceptor;
  double best_prob, prob, donor_prob, acceptor_prob;
  List_T accepted_hits, rejected_hits, donor_hits, acceptor_hits, p;
  int donor_length, acceptor_length;
  bool plusp;

  int sensedir;
#ifdef LARGE_GENOMES
  Uint8list_T ambcoords;
#else
  Uintlist_T ambcoords;
#endif
  Intlist_T amb_knowni, amb_nmismatches;
  Doublelist_T amb_probs;

  i = 0;
  while (i < n) {
    hit = hitarray[i];
    segmentj_left = Stage3end_chimera_segmentj_left(hit);
    j = i + 1;
    while (j < n && Stage3end_chimera_segmentj_left(hitarray[j]) == segmentj_left) {
      j++;
    }
    if (j == i + 1) {
      /* Singleton */
      debug7(printf("Saving hit %d\n",i));
      winners = List_push(winners,(void *) hit);

    } else {
      plusp = Stage3end_plusp(hit);
      best_nmismatches = querylength;
      best_prob = 0.0;
      for (k = i; k < j; k++) {
	hit = hitarray[k];
	debug7(printf("analyzing distance %d, donor length %d (%llu..%llu) and acceptor length %d (%llu..%llu), nmismatches %d, probabilities %f and %f\n",
		      Stage3end_distance(hit),Substring_match_length_orig(Stage3end_substring_donor(hit)),
		      Substring_genomicstart(Stage3end_substring_donor(hit)),Substring_genomicend(Stage3end_substring_donor(hit)),
		      Substring_match_length_orig(Stage3end_substring_acceptor(hit)),
		      Substring_genomicstart(Stage3end_substring_acceptor(hit)),Substring_genomicend(Stage3end_substring_acceptor(hit)),
		      Stage3end_nmismatches_whole(hit),Substring_chimera_prob(Stage3end_substring_donor(hit)),
		      Substring_chimera_prob(Stage3end_substring_acceptor(hit))));
	if ((nmismatches = Stage3end_nmismatches_whole(hit)) < best_nmismatches) {
	  best_nmismatches = nmismatches;
	}
	if ((prob = Stage3end_chimera_prob(hit)) > best_prob) {
	  best_prob = prob;
	}
      }

      n_good_spliceends = 0;
      accepted_hits = rejected_hits = (List_T) NULL;
      for (k = i; k < j; k++) {
	hit = hitarray[k];
	if (Stage3end_nmismatches_whole(hit) <= best_nmismatches + LOCALSPLICING_NMATCHES_SLOP &&
	    Stage3end_chimera_prob(hit) >= best_prob - LOCALSPLICING_PROB_SLOP) {
	  debug7(printf("accepting distance %d, probabilities %f and %f\n",
			Stage3end_distance(hit),Substring_chimera_prob(Stage3end_substring_donor(hit)),
			Substring_chimera_prob(Stage3end_substring_acceptor(hit))));
	  n_good_spliceends += 1;
	  accepted_hits = List_push(accepted_hits,(void *) hit);
	} else {
	  rejected_hits = List_push(rejected_hits,(void *) hit);
	}
      }

      if (n_good_spliceends == 0) {
	/* Conjunction is too strict.  Allow for disjunction instead. */
	List_free(&rejected_hits);
	for (k = i; k < j; k++) {
	  hit = hitarray[k];
	  if (Stage3end_nmismatches_whole(hit) <= best_nmismatches + LOCALSPLICING_NMATCHES_SLOP ||
	      Stage3end_chimera_prob(hit) >= best_prob - LOCALSPLICING_PROB_SLOP) {
	    debug7(printf("accepting distance %d, probabilities %f and %f\n",
			  Stage3end_distance(hit),Substring_chimera_prob(Stage3end_substring_donor(hit)),
			  Substring_chimera_prob(Stage3end_substring_acceptor(hit))));
	    n_good_spliceends += 1;
	    accepted_hits = List_push(accepted_hits,(void *) hit);
	  } else {
	    rejected_hits = List_push(rejected_hits,(void *) hit);
	  }
	}
      }
	
      for (p = rejected_hits; p != NULL; p = List_next(p)) {
	hit = (Stage3end_T) List_head(p);
	Stage3end_free(&hit);
      }
      List_free(&rejected_hits);

      if (n_good_spliceends == 1) {
	assert(List_length(accepted_hits) == 1);
	winners = List_push(winners,List_head(accepted_hits));
	List_free(&accepted_hits);

      } else {
	/* Multiple hits */
	donor_hits = acceptor_hits = (List_T) NULL;
	for (p = accepted_hits; p != NULL; p = List_next(p)) {
	  hit = (Stage3end_T) List_head(p);
	  donor = Stage3end_substring_donor(hit);
	  acceptor = Stage3end_substring_acceptor(hit);
	  if (Stage3end_plusp(hit) == true) {
	    if (Substring_genomicstart(donor) == segmentj_left) {
	      donor_hits = List_push(donor_hits,(void *) hit);
	    } else if (Substring_genomicstart(acceptor) == segmentj_left) {
	      acceptor_hits = List_push(acceptor_hits,(void *) hit);
	    } else {
	      abort();
	      Stage3end_free(&hit);
	    }
	  } else {
	    if (Substring_genomicend(donor) == segmentj_left) {
	      donor_hits = List_push(donor_hits,(void *) hit);
	    } else if (Substring_genomicend(acceptor) == segmentj_left) {
	      acceptor_hits = List_push(acceptor_hits,(void *) hit);
	    } else {
	      abort();
	      Stage3end_free(&hit);
	    }
	  }
	}

	if (donor_hits != NULL) {
	  subarray = (Stage3end_T *) List_to_array_n(&nn,donor_hits);
	  qsort(subarray,nn,sizeof(Stage3end_T),donor_match_length_cmp);
	  ii = 0;
	  while (ii < nn) {
	    hit = subarray[ii];
	    donor = Stage3end_substring_donor(hit);
	    donor_length = Substring_match_length_orig(donor);
	    jj = ii + 1;
	    while (jj < nn && Substring_match_length_orig(Stage3end_substring_donor(subarray[jj])) == donor_length) {
	      jj++;
	    }
	    if (jj == ii + 1) {
	      winners = List_push(winners,(void *) hit);
	    } else {
	      sensedir = Stage3end_sensedir(hit);

	      ambcoords = NULL;
	      amb_knowni = (Intlist_T) NULL;
	      amb_nmismatches = (Intlist_T) NULL;
	      amb_probs = (Doublelist_T) NULL;

	      for (kk = ii; kk < jj; kk++) {
		acceptor = Stage3end_substring_acceptor(subarray[kk]);
#ifdef LARGE_GENOMES
		ambcoords = Uint8list_push(ambcoords,Substring_splicecoord(acceptor));
#else
		ambcoords = Uintlist_push(ambcoords,Substring_splicecoord(acceptor));
#endif
		amb_knowni = Intlist_push(amb_knowni,-1);
		amb_nmismatches = Intlist_push(amb_nmismatches,Substring_nmismatches_whole(acceptor));
		amb_probs = Doublelist_push(amb_probs,Substring_chimera_prob(acceptor));
	      }

	      nmismatches_acceptor = best_nmismatches - Substring_nmismatches_whole(donor);
	      donor_prob = Junction_donor_prob(Stage3end_junctionA(hit));
	      prob = best_prob - donor_prob;
	      *ambiguous = List_push(*ambiguous,
				     (void *) Stage3end_new_splice(&(*found_score),
								   /*nmismatches_donor*/Substring_nmismatches_whole(donor),nmismatches_acceptor,
								   donor,/*acceptor*/NULL,donor_prob,/*acceptor_prob*/prob,/*distance*/0U,
								   /*shortdistancep*/false,/*penalty*/0,querylength,
								   /*amb_length*/Substring_match_length_orig(acceptor),/*amb_prob*/prob,
								   /*ambcoords_donor*/NULL,ambcoords,
								   /*amb_knowni_donor*/NULL,amb_knowni,
								   /*amb_nmismatches_donor*/NULL,amb_nmismatches,
								   /*amb_probs_donor*/NULL,amb_probs,
								   /*copy_donor_p*/true,/*copy_acceptor_p*/false,first_read_p,
								   sensedir,sarrayp));
	      Doublelist_free(&amb_probs);
	      Intlist_free(&amb_knowni);
	      Intlist_free(&amb_nmismatches);
#ifdef LARGE_GENOMES
	      Uint8list_free(&ambcoords);
#else
	      Uintlist_free(&ambcoords);
#endif
	      for (kk = ii; kk < jj; kk++) {
		hit = subarray[kk];
		Stage3end_free(&hit);
	      }
	    }

	    ii = jj;
	  }
	  FREE(subarray);
	  List_free(&donor_hits);
	}

	if (acceptor_hits != NULL) {
	  subarray = (Stage3end_T *) List_to_array_n(&nn,acceptor_hits);
	  qsort(subarray,nn,sizeof(Stage3end_T),acceptor_match_length_cmp);
	  ii = 0;
	  while (ii < nn) {
	    hit = subarray[ii];
	    acceptor = Stage3end_substring_acceptor(hit);
	    acceptor_length = Substring_match_length_orig(acceptor);
	    jj = ii + 1;
	    while (jj < nn && Substring_match_length_orig(Stage3end_substring_acceptor(subarray[jj])) == acceptor_length) {
	      jj++;
	    }
	    if (jj == ii + 1) {
	      winners = List_push(winners,(void *) hit);
	    } else {
	      sensedir = Stage3end_sensedir(hit);

	      ambcoords = NULL;
	      amb_knowni = (Intlist_T) NULL;
	      amb_nmismatches = (Intlist_T) NULL;
	      amb_probs = (Doublelist_T) NULL;

	      for (kk = ii; kk < jj; kk++) {
		donor = Stage3end_substring_donor(subarray[kk]);
#ifdef LARGE_GENOMES
		ambcoords = Uint8list_push(ambcoords,Substring_splicecoord(donor));
#else
		ambcoords = Uintlist_push(ambcoords,Substring_splicecoord(donor));
#endif
		amb_knowni = Intlist_push(amb_knowni,-1);
		amb_nmismatches = Intlist_push(amb_nmismatches,Substring_nmismatches_whole(donor));
		amb_probs = Doublelist_push(amb_probs,Substring_chimera_prob(donor));
	      }

	      nmismatches_donor = best_nmismatches - Substring_nmismatches_whole(acceptor);
	      acceptor_prob = Junction_acceptor_prob(Stage3end_junctionD(hit));
	      prob = best_prob - acceptor_prob;
	      *ambiguous = List_push(*ambiguous,
				     (void *) Stage3end_new_splice(&(*found_score),
								   nmismatches_donor,/*nmismatches_acceptor*/Substring_nmismatches_whole(acceptor),
								   /*donor*/NULL,acceptor,/*donor_prob*/prob,acceptor_prob,/*distance*/0U,
								   /*shortdistancep*/false,/*penalty*/0,querylength,
								   /*amb_length*/Substring_match_length_orig(donor),/*amb_prob*/prob,
								   ambcoords,/*ambcoords_acceptor*/NULL,
								   amb_knowni,/*amb_knowni_acceptor*/NULL,
								   amb_nmismatches,/*amb_nmismatches_acceptor*/NULL,
								   amb_probs,/*amb_probs_acceptor*/NULL,
								   /*copy_donor_p*/false,/*copy_acceptor_p*/true,first_read_p,
								   sensedir,sarrayp));
	      Doublelist_free(&amb_probs);
	      Intlist_free(&amb_knowni);
	      Intlist_free(&amb_nmismatches);
#ifdef LARGE_GENOMES
	      Uint8list_free(&ambcoords);
#else
	      Uintlist_free(&ambcoords);
#endif
	      for (kk = ii; kk < jj; kk++) {
		hit = subarray[kk];
		Stage3end_free(&hit);
	      }
	    }

	    ii = jj;
	  }
	  FREE(subarray);
	  List_free(&acceptor_hits);
	}

	List_free(&accepted_hits);
      }
    }

    i = j;
  }

  return winners;
}

List_T
Splice_group_by_segmentj (int *found_score, List_T localsplicing, List_T *ambiguous,
			  int querylength, bool first_read_p, bool sarrayp) {
  List_T winners = NULL, p;
  Stage3end_T *array_forward, *array_anti, hit;
  int n_sense_forward = 0, n_sense_anti = 0, k_forward, k_anti;

  for (p = localsplicing; p != NULL; p = List_next(p)) {
    hit = (Stage3end_T) List_head(p);
    if (Stage3end_sensedir(hit) == SENSE_FORWARD) {
      n_sense_forward++;
    } else {
      assert(Stage3end_sensedir(hit) == SENSE_ANTI);
      n_sense_anti++;
    }
  }

  if (n_sense_forward > 0) {
    array_forward = (Stage3end_T *) MALLOCA(n_sense_forward * sizeof(Stage3end_T));
    k_forward = 0;
  }
  if (n_sense_anti > 0) {
    array_anti = (Stage3end_T *) MALLOCA(n_sense_anti * sizeof(Stage3end_T));
    k_anti = 0;
  }

  for (p = localsplicing; p != NULL; p = List_next(p)) {
    hit = (Stage3end_T) List_head(p);
    if (Stage3end_sensedir(hit) == SENSE_FORWARD) {
      array_forward[k_forward++] = (Stage3end_T) List_head(p);
    } else {
      array_anti[k_anti++] = (Stage3end_T) List_head(p);
    }
  }

  if (n_sense_forward > 0) {
    qsort(array_forward,n_sense_forward,sizeof(Stage3end_T),Stage3end_chimera_segmentj_cmp);
    winners = group_by_segmentj_aux(&(*found_score),winners,&(*ambiguous),array_forward,n_sense_forward,
				    querylength,first_read_p,sarrayp);
    FREEA(array_forward);
  }

  if (n_sense_anti > 0) {
    qsort(array_anti,n_sense_anti,sizeof(Stage3end_T),Stage3end_chimera_segmentj_cmp);
    winners = group_by_segmentj_aux(&(*found_score),winners,&(*ambiguous),array_anti,n_sense_anti,
				    querylength,first_read_p,sarrayp);
    FREEA(array_anti);
  }

  List_free(&localsplicing);

  return winners;
}

