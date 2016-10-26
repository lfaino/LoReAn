static char rcsid[] = "$Id: smooth.c 145990 2014-08-25 21:47:32Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "smooth.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>		/* For pow(), log() */
#include "bool.h"
#include "mem.h"
#include "except.h"
#include "comp.h"
#include "pair.h"
#include "pairdef.h"
#include "listdef.h"
#include "intlist.h"
#include "doublelist.h"
#include "pbinom.h"


/* Will examine internal exons smaller than this prior to single gaps */
#define ZERONETGAP 9
#define SHORTEXONLEN_NETGAP 15

/* Will delete/mark internal exons smaller than this for dual genome gap */
#define DELETE_THRESHOLD 0.1
#define MARK_THRESHOLD 1e-7

/* Will delete internal exons smaller than this at ends */
#define SHORTEXONLEN_END 10


#define LOOSE_EXON_PVALUE 0.05
#define STRICT_EXON_PVALUE 1e-4


#ifdef GSNAP
/* Allow more intron predictions in ends of short reads */
#define SHORTEXONPROB_END 0.10
#else
#define SHORTEXONPROB_END 0.05
#endif

#define BIGGAP 9


#ifdef DEBUG
#define debug(x) x
#else 
#define debug(x)
#endif


typedef enum {KEEP, DELETE, MARK} Exonstatus_T;


static bool
big_gap_p (Pair_T pair, bool bysizep) {

#if 0
  /* No longer following the logic below */
  /* When bysizep is true, single gaps should have been solved already
     (in pass 2), so any gap remaining must be a poor one (must have
     been removed because of a negative score), and we should consider
     all gaps, regardless of size.  When bysizep is false, we must be
     doing initial smoothing, before any single gaps have been solved,
     so we want to consider only big gaps.  This change was motivated
     by alignment of mouse HER2 against human genome, where a 58/58
     section was deleted in stage 3, but subsequent smoothing failed
     to recognize the short exon between this section and the next
     intron. */
#endif

  if (bysizep == true) {
    return true;
  } else if (pair->genomejump > pair->queryjump) {
    return (pair->genomejump - pair->queryjump) > BIGGAP ? true : false;
  } else {
    return false;
  }
}


static int *
get_exonlengths (int **exonmatches, int **exonprotects, int *nexons, List_T pairs, bool bysizep) {
  int *exonlengths;
  Intlist_T list = NULL, matchlist = NULL, protectlist = NULL;
  Pair_T firstpair, pair;
  int querypos, nmatches = 0, nprotected = 0;
#ifdef DEBUG
  int i;
#endif

  firstpair = List_head(pairs);
  querypos = firstpair->querypos;

  while (pairs != NULL) {
    pair = List_head(pairs);
    if (pair->gapp == true && big_gap_p(pair,bysizep) == true) {
      list = Intlist_push(list,querypos-firstpair->querypos+1);
      matchlist = Intlist_push(matchlist,nmatches);
      protectlist = Intlist_push(protectlist,nprotected);

      pairs = List_next(pairs);
      if (pairs != NULL) {
	firstpair = List_head(pairs);
	querypos = firstpair->querypos;
	nmatches = nprotected = 0;
      }
    } else {
      querypos = pair->querypos;
      if (pair->comp == MATCH_COMP || pair->comp == DYNPROG_MATCH_COMP) {
	nmatches++;
      }
      if (pair->protectedp == true) {
	nprotected++;
      }
      pairs = List_next(pairs);
    }
  }

  list = Intlist_push(list,querypos-firstpair->querypos+1);
  list = Intlist_reverse(list);
  matchlist = Intlist_push(matchlist,nmatches);
  matchlist = Intlist_reverse(matchlist);
  protectlist = Intlist_push(protectlist,nprotected);
  protectlist = Intlist_reverse(protectlist);

  exonlengths = Intlist_to_array(&(*nexons),list);
  *exonmatches = Intlist_to_array(&(*nexons),matchlist);
  *exonprotects = Intlist_to_array(&(*nexons),protectlist);
  debug(
	printf("%d exons: ",*nexons);
	for (i = 0; i < *nexons; i++) {
	  printf("%d (%d matches),",exonlengths[i],(*exonmatches)[i]);
	}
	printf("\n");
	);
  Intlist_free(&list);
  Intlist_free(&matchlist);
  Intlist_free(&protectlist);

  return exonlengths;
}


static void
get_exonmatches (int *total_matches, int *total_denominator, int **exonmatches, int **exon_denominator,
		 int *nexons, List_T pairs) {
  Intlist_T matchlist = NULL, denominatorlist = NULL;
  Pair_T firstpair, pair;
  int querypos, nmatches = 0, nmismatches = 0;
#ifdef DEBUG
  int i;
#endif

  *total_matches = *total_denominator = 0;

  firstpair = List_head(pairs);
  querypos = firstpair->querypos;

  while (pairs != NULL) {
    pair = List_head(pairs);
    if (pair->gapp == true && big_gap_p(pair,/*bysizep*/false) == true) {
      matchlist = Intlist_push(matchlist,nmatches);
      denominatorlist = Intlist_push(denominatorlist,nmatches + nmismatches);
      *total_matches += nmatches;
      *total_denominator += (nmatches + nmismatches);

      pairs = List_next(pairs);
      if (pairs != NULL) {
	firstpair = List_head(pairs);
	querypos = firstpair->querypos;
	nmatches = nmismatches = 0;
      }
    } else {
      querypos = pair->querypos;
      if (pair->comp == MATCH_COMP || pair->comp == DYNPROG_MATCH_COMP) {
	nmatches++;
      } else {
	nmismatches++;
      }
      pairs = List_next(pairs);
    }
  }

  matchlist = Intlist_push(matchlist,nmatches);
  matchlist = Intlist_reverse(matchlist);
  denominatorlist = Intlist_push(denominatorlist,nmatches + nmismatches);
  denominatorlist = Intlist_reverse(denominatorlist);
  *total_matches += nmatches;
  *total_denominator += (nmatches + nmismatches);

  *exonmatches = Intlist_to_array(&(*nexons),matchlist);
  *exon_denominator = Intlist_to_array(&(*nexons),denominatorlist);
  debug(
	printf("%d exons: ",*nexons);
	for (i = 0; i < *nexons; i++) {
	  printf("%d matches, %d denominator,",(*exonmatches)[i],(*exon_denominator)[i]);
	}
	printf("\n");
	);
  Intlist_free(&denominatorlist);
  Intlist_free(&matchlist);

  return;
}



static List_T
get_intron_neighborhoods (int **intron_matches_left, int **intron_denominator_left,
			  int **intron_matches_right, int **intron_denominator_right,
			  int nintrons, List_T pairs) {
  List_T path, p;
  Pair_T pair;
  int introni;
  int i;

  *intron_matches_left = (int *) CALLOC(nintrons,sizeof(int));
  *intron_denominator_left = (int *) CALLOC(nintrons,sizeof(int));
  *intron_matches_right = (int *) CALLOC(nintrons,sizeof(int));
  *intron_denominator_right = (int *) CALLOC(nintrons,sizeof(int));

  introni = -1;
  for (p = pairs; p != NULL; p = List_next(p)) {
    pair = List_head(p);
    if (pair->gapp == true && big_gap_p(pair,/*bysizep*/false) == true) {
      introni += 1;
      i = 0;
      
    } else if (introni >= 0 && i < 25) {
      if (pair->comp == MATCH_COMP || pair->comp == DYNPROG_MATCH_COMP) {
	(*intron_matches_left)[introni] += 1;
      }
      (*intron_denominator_left)[introni] += 1;
      i++;

    } else {
      i++;
    }
  }

  path = List_reverse(pairs);
  introni = nintrons;
  for (p = path; p != NULL; p = List_next(p)) {
    pair = List_head(p);
    if (pair->gapp == true && big_gap_p(pair,/*bysizep*/false) == true) {
      introni -= 1;
      i = 0;
      
    } else if (introni < nintrons && i < 25) {
      if (pair->comp == MATCH_COMP || pair->comp == DYNPROG_MATCH_COMP) {
	(*intron_matches_right)[introni] += 1;
      }
      (*intron_denominator_right)[introni] += 1;
      i++;

    } else {
      i++;
    }
  }

  pairs = List_reverse(path);
  return pairs;
}


static int *
get_intronlengths (int *nintrons, List_T pairs, bool bysizep) {
  int *intronlengths, length;
  Intlist_T list = NULL;
  Pair_T pair;
#ifdef DEBUG
  int i;
#endif

  while (pairs != NULL) {
    pair = List_head(pairs);
    if (pair->gapp == true && big_gap_p(pair,bysizep) == true) {
#if 0
      if (pair->genomejump > pair->queryjump) {
	list = Intlist_push(list,pair->genomejump);
      } else {
	list = Intlist_push(list,pair->queryjump);
      }
#else
      if ((length = pair->genomejump - pair->queryjump) < 0) {
	/* cDNA insertion */
	list = Intlist_push(list,-length);
      } else {
	list = Intlist_push(list,length);
      }
#endif
    }
    pairs = List_next(pairs);
  }

  list = Intlist_reverse(list);
  intronlengths = Intlist_to_array(&(*nintrons),list);
  debug(
	printf("%d introns: ",*nintrons);
	for (i = 0; i < *nintrons; i++) {
	  printf("%d,",intronlengths[i]);
	}
	printf("\n");
	);
  Intlist_free(&list);

  return intronlengths;
}


static void
get_intronprobs (double **donor_probs, double **acceptor_probs, int *nintrons, List_T pairs) {
  Doublelist_T donor_prob_list = NULL, acceptor_prob_list = NULL;
  Pair_T pair;
#ifdef DEBUG
  int i;
#endif

  while (pairs != NULL) {
    pair = List_head(pairs);
    if (pair->gapp == true && big_gap_p(pair,/*bysizep*/false) == true) {
      donor_prob_list = Doublelist_push(donor_prob_list,pair->donor_prob);
      acceptor_prob_list = Doublelist_push(acceptor_prob_list,pair->acceptor_prob);
    }
    pairs = List_next(pairs);
  }

  donor_prob_list = Doublelist_reverse(donor_prob_list);
  acceptor_prob_list = Doublelist_reverse(acceptor_prob_list);
  *donor_probs = Doublelist_to_array(&(*nintrons),donor_prob_list);
  *acceptor_probs = Doublelist_to_array(&(*nintrons),acceptor_prob_list);
  debug(
	printf("%d introns: ",*nintrons);
	for (i = 0; i < *nintrons; i++) {
	  printf("%f+%f,",(*donor_probs)[i],(*acceptor_probs)[i]);
	}
	printf("\n");
	);
  Doublelist_free(&acceptor_prob_list);
  Doublelist_free(&donor_prob_list);

  return;
}


static const Except_T length_error = {"Negative exon or intron length"};

static double
compute_prob (int exonlen, int intronlen, int indexsize) {
  double prob;

  if (exonlen < indexsize) {
    prob = 1.0;
  } else {
    prob = 1 - pow(1.0-pow(4.0,(double) -exonlen),(double) intronlen);
  }
  debug(printf("Probability of exon of length %d (indexsize %d) with intron of length %d is %g\n",
	       exonlen,indexsize,intronlen,prob));
  return prob;
}

static bool
short_exon_byprob (int exonlen, int intronlen, int indexsize, double prob_threshold) {
  double prob;

  prob = compute_prob(exonlen,intronlen,indexsize);
  if (prob < prob_threshold) {
    return false;
  } else {
    return true;
  }
}

static bool
short_exon_bylength (int exonlen, int length_threshold) {

  if (exonlen >= length_threshold) {
    return false;
  } else {
    return true;
  }
}

static void
zero_net_gap (int *starti, int *startj, int i, int j, int *intronlengths) {
  int netgap, bestnetgap = 1000000;
  int k, l, adji;

  if (i == 0) {
    adji = 0;
  } else {
    adji = i-1;
  }

  for (k = adji; k < j; k++) {
    netgap = intronlengths[k];
    for (l = k+1; l < j; l++) {
      netgap += intronlengths[l];
      debug(printf("zero_net_gap: netgap from %d to %d is %d\n",k,l,netgap));
      if (abs(netgap) < bestnetgap) {
	bestnetgap = abs(netgap);
	*starti = k+1;
	*startj = l;
      }
    }
  }

  debug(printf("zero_net_gap: best result is %d from %d to %d\n",bestnetgap,*starti,*startj));

  if (0 && bestnetgap > ZERONETGAP) {
    debug(printf("zero_net_gap: not recommending any deletions because bestnetgap %d > zeronetgap %d\n",
		 bestnetgap,ZERONETGAP));
    *starti = *startj = -1;
  }

  return;
}


static void
find_internal_shorts_by_netgap (bool *deletep, int *exonstatus, int *exonmatches, int nexons,
				int *intronlengths) {
  int starti, startj, i, j;
  int exonlen;

  *deletep = false;
  for (i = 0; i < nexons; i++) {
    exonstatus[i] = KEEP;
  }
  
  /* Mark short middle exons */
  for (i = 1; i < nexons - 1; i++) {
    exonlen = exonmatches[i];
    if (exonlen < SHORTEXONLEN_NETGAP) {
      exonstatus[i] = MARK;
    }
  }

  /* Find internal shorts */
  i = 0;
  while (i < nexons) {
    if (exonstatus[i] == MARK) {
      j = i;
      while (j < nexons && exonstatus[j] == MARK) {
	j++;
      }
      debug(printf("Calling zero_net_gap with %d exons\n",j-i));
      zero_net_gap(&starti,&startj,i,j,intronlengths);
      if (starti >= 0) {
	for (j = starti; j <= startj; j++) {
	  *deletep = true;
	  exonstatus[j] = DELETE;
	}
      } else if (j - i == 1) {
	exonstatus[i] = MARK;
      }
      i = j;
    } else {
      i++;
    }
  }

  return;
}


static void
find_internal_shorts_by_size (bool *shortp, bool *deletep, int *exonstatus,
			      int *exonmatches, int *exonprotects, int nexons,
			      int *intronlengths, int stage2_indexsize) {
  int i;
  int exonlen, intronlen;
  double prob;

  *shortp = *deletep = false;
  for (i = 0; i < nexons; i++) {
    exonstatus[i] = KEEP;
  }
  
  /* Mark short middle exons */
  for (i = 1; i < nexons - 1; i++) {
    exonlen = exonmatches[i];
    intronlen = intronlengths[i-1]+intronlengths[i];
    if ((double) exonprotects[i] > 0.50 * (double) exonlen) {
      /* Keep protected exon */
    } else {
      prob = compute_prob(exonlen,intronlen,stage2_indexsize); /* Hack: add 4 for possible canonical dinucleotides */
      if (prob > DELETE_THRESHOLD) {
	*deletep = true;
	exonstatus[i] = DELETE;
      } else if (prob > MARK_THRESHOLD) {
	*shortp = true;
	exonstatus[i] = MARK;
      }
    }
  }

  return;
}


static bool
bad_intron_p (double donor_prob, double acceptor_prob, int intron_matches_left, int intron_denominator_left,
	      int intron_matches_right, int intron_denominator_right, int total_matches, int total_denominator) {
  double theta_left, theta_right;

  if (donor_prob < 0.9) {
    debug(printf("Intron with donor_prob %f is bad\n",donor_prob));
    return true;
  } else if (acceptor_prob < 0.9) {
    debug(printf("Intron with acceptor_prob %f is bad\n",acceptor_prob));
    return true;
  } else {
    theta_left = (double) (total_matches - intron_matches_left)/(double) (total_denominator - intron_denominator_left + 1);
    if (theta_left > 1.0) {
      theta_left = 1.0;
    }
    if (Pbinom(intron_matches_left,intron_denominator_left,theta_left) < 1e-3) {
      debug(printf("Intron with matches %d/%d is bad with pvalue %f\n",
		   intron_matches_left,intron_denominator_left,Pbinom(intron_matches_left,intron_denominator_left,theta_left)));
      return true;
    } else {
      theta_right = (double) (total_matches - intron_matches_right)/(double) (total_denominator - intron_denominator_right + 1);
      if (theta_right > 1.0) {
	theta_right = 1.0;
      }
      if (Pbinom(intron_matches_right,intron_denominator_right,theta_right) < 1e-3) {
	debug(printf("Intron with matches %d/%d is bad with pvalue %f\n",
		     intron_matches_right,intron_denominator_right,Pbinom(intron_matches_right,intron_denominator_right,theta_right)));
	return true;
      } else {
	debug(printf("Intron with matches %d/%d and %d/%d is okay with pvalues %f and %f\n",
		     intron_matches_left,intron_denominator_left,intron_matches_right,intron_denominator_right,
		     Pbinom(intron_matches_left,intron_denominator_left,theta_left),
		     Pbinom(intron_matches_right,intron_denominator_right,theta_right)));
	return false;
      }
    }
  }
}


struct Smoothcell_T {
  int exoni;
  double pvalue;
  int exonstatus;
};

static int
Smoothcell_cmp (const void *x, const void *y) {
  struct Smoothcell_T a = * (struct Smoothcell_T *) x;
  struct Smoothcell_T b = * (struct Smoothcell_T *) y;

  if (a.pvalue < b.pvalue) {
    return -1;
  } else if (b.pvalue < a.pvalue) {
    return +1;
  } else {
    return 0;
  }
}


static void
find_internal_bads_by_prob (bool *deletep, int *exonstatus, int *exonmatches, int *exon_denominator, int nexons,
			    double *donor_probs, double *acceptor_probs,
			    int *intron_matches_left, int *intron_denominator_left,
			    int *intron_matches_right, int *intron_denominator_right,
			    int total_matches, int total_denominator) {
  struct Smoothcell_T *cells;
  int starti, startj, i, j;
  int exonlen;
  double theta_1, theta0, theta1;
  int numerator_1, denominator_1, numerator0, denominator0, numerator1, denominator1;
  bool intron1_bad_p, intron2_bad_p;
  
  double worst_exon_pvalue = 1.0, pvalue;
  int worst_exoni = -1;

  *deletep = false;
  cells = (struct Smoothcell_T *) MALLOCA(nexons * sizeof(struct Smoothcell_T));
  for (i = 0; i < nexons; i++) {
    exonstatus[i] = KEEP;
    cells[i].exoni = i;
    cells[i].pvalue = 1.0;
    cells[i].exonstatus = KEEP;
  }
  
  /* Mark short middle exons */
  intron1_bad_p = bad_intron_p(donor_probs[0],acceptor_probs[0],intron_matches_left[0],intron_denominator_left[0],
			       intron_matches_right[0],intron_denominator_right[0],total_matches,total_denominator);
  for (i = 1; i < nexons - 1; i++) {
    /* exonlen = exonmatches[i]; */
    intron2_bad_p = bad_intron_p(donor_probs[i],acceptor_probs[i],intron_matches_left[i],intron_denominator_left[i],
				 intron_matches_right[i],intron_denominator_right[i],total_matches,total_denominator);
    debug(printf("For exon %d, left intron bad %d, right intron bad %d\n",i,intron1_bad_p,intron2_bad_p));

    if (intron1_bad_p == true && intron2_bad_p == true) {
      numerator0 = exonmatches[i];
      denominator0 = exon_denominator[i];
      theta0 = (double) (total_matches - numerator0 + 1)/(double) (total_denominator - denominator0 + 1);
      debug(printf("  Checking %d matches, %d denominator, theta %f => pvalue %g\n",
		   numerator0,denominator0,theta0,Pbinom(numerator0,denominator0,theta0)));
      if ((pvalue = Pbinom(numerator0,denominator0,theta0)) < STRICT_EXON_PVALUE) {
	cells[i].pvalue = pvalue;
	cells[i].exonstatus = DELETE;
      }

    } else {
      /* Do nothing */
    }

    intron1_bad_p = intron2_bad_p;
  }

  qsort(cells,nexons,sizeof(struct Smoothcell_T),Smoothcell_cmp);
  i = 0;
  while (i < nexons && cells[i].pvalue < STRICT_EXON_PVALUE) {
    if (cells[i].exonstatus == DELETE) {
      *deletep = true;
      exonstatus[cells[i].exoni] = DELETE;
      exonstatus[cells[i].exoni - 1] = KEEP; /* Prevent consecutive deletes */
      exonstatus[cells[i].exoni + 1] = KEEP; /* Prevent consecutive deletes */
    }
    i++;
  }
  FREEA(cells);

  return;
}



/* For ends, we turn off the indexsize parameter */
static void
find_end_shorts (bool *deletep, int *exonstatus, int *exonmatches, int nexons, int *intronlengths) {
  int i;
  bool shortp;

  *deletep = false;
  for (i = 0; i < nexons; i++) {
    exonstatus[i] = KEEP;
  }

  shortp = true;
  i = 0;
  while (i < nexons - 1 && shortp == true) {
    if (short_exon_bylength(exonmatches[i],SHORTEXONLEN_END) == true &&
	short_exon_byprob(exonmatches[i],intronlengths[i],/*indexsize*/0,SHORTEXONPROB_END) == true) {
      *deletep = true;
      exonstatus[i] = DELETE;
    } else {
      shortp = false;
    }
    i++;
  }
    
  shortp = true;
  i = nexons - 1;
  while (i > 0 && shortp == true) {
    if (short_exon_bylength(exonmatches[i],SHORTEXONLEN_END) == true &&
	short_exon_byprob(exonmatches[i],intronlengths[i-1],/*indexsize*/0,SHORTEXONPROB_END) == true) {
      *deletep = true;
      exonstatus[i] = DELETE;
    } else {
      shortp = false;
    }
    --i;
  }
    
  return;
}


#ifdef DEBUG
static void
print_exon_status (int *exonstatus, int *exonmatches, int *exonprotects, int nexons) {
  int i;

  if (exonprotects == NULL) {
    for (i = 0; i < nexons; i++) {
      if (exonstatus[i] == KEEP) {
	printf("Long exon %d of %d matches => keep\n",i,exonmatches[i]);
      } else if (exonstatus[i] == MARK) {
	printf("Short exon %d of %d matches => mark\n",i,exonmatches[i]);
      } else if (exonstatus[i] == DELETE) {
	printf("Exon %d of %d matches => delete\n",i,exonmatches[i]);
      } else {
	abort();
      }
    }
  } else {
    for (i = 0; i < nexons; i++) {
      if (exonstatus[i] == KEEP) {
	printf("Long exon %d of %d matches (%d protected) => keep\n",i,exonmatches[i],exonprotects[i]);
      } else if (exonstatus[i] == MARK) {
	printf("Short exon %d of %d matches (%d protected) => mark\n",i,exonmatches[i],exonprotects[i]);
      } else if (exonstatus[i] == DELETE) {
	printf("Exon %d of %d matches (%d protected) => delete\n",i,exonmatches[i],exonprotects[i]);
      } else {
	abort();
      }
    }
  }

  return;
}
#endif


static List_T
delete_and_mark_exons (List_T pairs, Pairpool_T pairpool,
		       int *exonstatus, bool markp, bool bysizep) {
  List_T newpairs = NULL;
  Pair_T pair;
  int currstatus, prevstatus = KEEP;
  int i;

  i = 0;
  currstatus = exonstatus[i];
  while (pairs != NULL) {
    /* pairptr = pairs; */
    /* pairs = Pairpool_pop(pairs,&pair); */
    pair = (Pair_T) pairs->first;
    if (pair->gapp == true && big_gap_p(pair,bysizep) == true) {
      prevstatus = currstatus;
      currstatus = exonstatus[++i];
      debug(printf("Gap observed\n"));
      if (prevstatus != DELETE && currstatus != DELETE) {
#ifdef WASTE
	newpairs = Pairpool_push_existing(newpairs,pairpool,pair);
#else
	newpairs = List_transfer_one(newpairs,&pairs);
#endif
      } else {
	pairs = Pairpool_pop(pairs,&pair);
      }

    } else if (currstatus == KEEP) {
      /* debug(printf("Marking position %d as keep\n",pair->querypos)); */
      if (prevstatus == DELETE) {
	if (newpairs != NULL) {
	  newpairs = Pairpool_push_gapholder(newpairs,pairpool,/*queryjump*/0,/*genomejump*/0,
					     /*leftpair*/newpairs->first,/*rightpair*/pair,/*knownp*/false);
	}
	prevstatus = /*currstatus*/ KEEP;
      }

#ifdef WASTE
      newpairs = Pairpool_push_existing(newpairs,pairpool,pair);
#else
      newpairs = List_transfer_one(newpairs,&pairs);
#endif
    } else if (currstatus == MARK) {
      debug(printf("Marking position %d as short in pair %p\n",pair->querypos,pair));
      if (prevstatus == DELETE) {
	if (newpairs != NULL) {
	  newpairs = Pairpool_push_gapholder(newpairs,pairpool,/*queryjump*/0,/*genomejump*/0,
					     /*leftpair*/newpairs->first,/*rightpair*/pair,/*knownp*/false);
	}
	prevstatus = /*currstatus*/ MARK;
      }

      if (markp == true) {
	pair->shortexonp = true;
      }
#ifdef WASTE
      newpairs = Pairpool_push_existing(newpairs,pairpool,pair);
#else
      newpairs = List_transfer_one(newpairs,&pairs);
#endif
    } else {
      debug(printf("Marking position %d for deletion\n",pair->querypos));
      pairs = Pairpool_pop(pairs,&pair);
    }
  }

  /* Remove gaps at end */
  if (newpairs != NULL) {
    pair = List_head(newpairs);
  }
  while (newpairs != NULL && pair->gapp == true) {
    debug(printf("Popping gap at end\n"));
    newpairs = Pairpool_pop(newpairs,&pair);
  }

  /* Remove gaps at beginning */
  newpairs = List_reverse(newpairs);
  if (newpairs != NULL) {
    pair = List_head(newpairs);
  }
  while (newpairs != NULL && pair->gapp == true) {
    debug(printf("Popping gap at beginning\n"));
    newpairs = Pairpool_pop(newpairs,&pair);
  }

  debug(printf("Result of delete_and_mark_exons:\n"));
  debug(Pair_dump_list(newpairs,/*zerobasedp*/true));
  debug(printf("\n"));

  return newpairs;
}


#if 0
static List_T
mark_exons (List_T pairs,
#ifdef WASTE
	    Pairpool_T pairpool,
#endif
	    int *exonstatus) {
  List_T newpairs = NULL, pairptr;
  Pair_T pair;
  int currstatus, prevstatus;
  int i;

  debug(
	for (i = 0; i < nexons; i++) {
	  if (exonstatus[i] == KEEP) {
	    printf("Long exon %d of %d matches => keep\n",i,exonmatches[i]);
	  } else if (exonstatus[i] == MARK) {
	    printf("Short exon %d of %d matches => mark\n",i,exonmatches[i]);
	  } else if (exonstatus[i] == DELETE) {
	    printf("Exon %d of %d matches => delete\n",i,exonmatches[i]);
	  } else {
	    abort();
	  }
	}
	);

  i = 0;
  currstatus = exonstatus[i];
  while (pairs != NULL) {
    pairptr = pairs;
    pairs = Pairpool_pop(pairs,&pair);
    if (pair->gapp == true && big_gap_p(pair,/*bysizep*/true) == true) {
      prevstatus = currstatus;
      currstatus = exonstatus[++i];
      debug(printf("Gap observed\n"));
      if (prevstatus != DELETE && currstatus != DELETE) {
#ifdef WASTE
	newpairs = Pairpool_push_existing(newpairs,pairpool,pair);
#else
	newpairs = List_push_existing(newpairs,pairptr);
#endif
      }

    } else if (currstatus == KEEP) {
      /* debug(printf("Marking position %d as keep\n",pair->querypos)); */
#ifdef WASTE
      newpairs = Pairpool_push_existing(newpairs,pairpool,pair);
#else
      newpairs = List_push_existing(newpairs,pairptr);
#endif
    } else if (currstatus == MARK) {
      debug(printf("Marking position %d as short in pair %p\n",pair->querypos,pair));
      pair->shortexonp = true;

#ifdef WASTE
      newpairs = Pairpool_push_existing(newpairs,pairpool,pair);
#else
      newpairs = List_push_existing(newpairs,pairptr);
#endif
    } else {
      /* Normally would delete */
#ifdef WASTE
      newpairs = Pairpool_push_existing(newpairs,pairpool,pair);
#else
      newpairs = List_push_existing(newpairs,pairptr);
#endif
    }
  }

  debug(printf("Result of mark_exons:\n"));
  debug(Pair_dump_list(newpairs,/*zerobasedp*/true));
  debug(printf("\n"));

  return List_reverse(newpairs);
}
#endif



void
Smooth_reset (List_T pairs) {
  List_T p;
  Pair_T pair;

  for (p = pairs; p != NULL; p = List_next(p)) {
    pair = List_head(p);
    pair->shortexonp = false;
  }
  return;
}



/* Needed for low-identity alignments */
/* Assumes pairs are from 1..querylength.  Reverses the pairs to be querylength..1 */
List_T
Smooth_pairs_by_netgap (bool *deletep, List_T pairs, Pairpool_T pairpool) {
  int *exonstatus;
  int *exonlengths, *exonmatches, *exonprotects, *intronlengths, nexons, nintrons;
#ifdef DEBUG
  int i;
#endif

  *deletep = false;
  /* smooth_reset(pairs); */
  if (pairs != NULL) {
    /* Remove internal shorts */
    exonlengths = get_exonlengths(&exonmatches,&exonprotects,&nexons,pairs,/*bysizep*/false);
    intronlengths = get_intronlengths(&nintrons,pairs,/*bysizep*/false);

    debug(
	  printf("Beginning of smoothing.  Initial structure:\n");
	  for (i = 0; i < nexons-1; i++) {
	    printf("Exon %d of length %d (%d matches, %d protects)\n",i,exonlengths[i],exonmatches[i],exonprotects[i]);
	    printf("Intron %d of length %d\n",i,intronlengths[i]);
	  }
	  printf("Exon %d of length %d (%d matches)\n",nexons-1,exonlengths[nexons-1],exonmatches[nexons-1]);
	  );

    debug(printf("\nFind internal shorts\n"));
    exonstatus = (int *) MALLOCA(nexons * sizeof(int));
    find_internal_shorts_by_netgap(&(*deletep),exonstatus,exonmatches,nexons,intronlengths);
    debug(printf("\nRemove internal shorts\n"));
    if (*deletep == true) {
      debug(print_exon_status(exonstatus,exonmatches,exonprotects,nexons));
      pairs = delete_and_mark_exons(pairs,pairpool,exonstatus,/*markp*/false,/*bysizep*/false);
    }
    
#if 0
    /* This is not correct */
    debug(
	  printf("After removing internal shorts:\n");
	  for (i = 0; i < nexons-1; i++) {
	    printf("Exon %d of length %d (%d matches)\n",i,exonlengths[i],exonmatches[i]);
	    printf("Intron %d of length %d\n",i,intronlengths[i]);
	  }
	  printf("Exon %d of length %d\n",nexons-1,exonlengths[nexons-1]);
	  );
#endif

    FREEA(exonstatus);
    FREE(intronlengths);
    FREE(exonprotects);
    FREE(exonmatches);
    FREE(exonlengths);

    debug(printf("Ending of smoothing\n\n"));
  }

  return pairs;
}


/* Assumes pairs are from 1..querylength.  Reverses the pairs to be querylength..1 */
List_T
Smooth_pairs_by_size (bool *shortp, bool *deletep, List_T pairs, Pairpool_T pairpool, int stage2_indexsize) {
  int *exonstatus;
  int *exonlengths, *exonmatches, *exonprotects, *intronlengths, nexons, nintrons;
  bool delete1p, delete2p;
#ifdef DEBUG
  int i;
#endif

  *shortp = *deletep = false;

#if 0
  /* Should be handled instead by trim_noncanonical_end5 and trim_noncanonical_end3 */
  /* smooth_reset(pairs); */
  if (pairs != NULL) {
    /* Trim ends */
    exonlengths = get_exonlengths(&exonmatches,&exonprotects,&nexons,pairs,/*bysizep*/true);
    intronlengths = get_intronlengths(&nintrons,pairs,/*bysizep*/true);

    debug(printf("\nFind end shorts\n"));
    exonstatus = (int *) MALLOCA(nexons * sizeof(int));
    find_end_shorts(&delete1p,exonstatus,exonmatches,nexons,intronlengths);
    if (delete1p == true) {
      *deletep = true;
      debug(print_exon_status(exonstatus,exonmatches,exonprotects,nexons));
      pairs = delete_and_mark_exons(pairs,pairpool,exonstatus,/*markp*/false,/*bysizep*/true);
    }

    FREEA(exonstatus);

    FREE(intronlengths);
    FREE(exonprotects);
    FREE(exonmatches);
    FREE(exonlengths);
  }
#endif

  if (pairs != NULL) {
    /* Remove internal shorts */
    exonlengths = get_exonlengths(&exonmatches,&exonprotects,&nexons,pairs,/*bysizep*/true);
    intronlengths = get_intronlengths(&nintrons,pairs,/*bysizep*/true);

    debug(
	  printf("Beginning of smoothing.  Initial structure:\n");
	  for (i = 0; i < nexons-1; i++) {
	    printf("Exon %d of length %d (%d matches, %d protects)\n",i,exonlengths[i],exonmatches[i],exonprotects[i]);
	    printf("Intron %d of length %d\n",i,intronlengths[i]);
	  }
	  printf("Exon %d of length %d (%d matches)\n",nexons-1,exonlengths[nexons-1],exonmatches[nexons-1]);
	  );

    debug(printf("\nFind internal shorts\n"));
    exonstatus = (int *) MALLOCA(nexons * sizeof(int));
    find_internal_shorts_by_size(&(*shortp),&delete2p,exonstatus,exonmatches,exonprotects,nexons,intronlengths,stage2_indexsize);
    debug(printf("\nRemove internal shorts\n"));
    if (delete2p == true) {
      *deletep = true;
    }
    if (delete2p == true || *shortp == true) {
      debug(print_exon_status(exonstatus,exonmatches,exonprotects,nexons));
      pairs = delete_and_mark_exons(pairs,pairpool,exonstatus,/*markp*/true,/*bysizep*/true);
    }

    debug(
	  printf("After removing internal shorts:\n");
	  for (i = 0; i < nexons-1; i++) {
	    printf("Exon %d of length %d (%d matches)\n",i,exonlengths[i],exonmatches[i]);
	    printf("Intron %d of length %d\n",i,intronlengths[i]);
	  }
	  printf("Exon %d of length %d\n",nexons-1,exonlengths[nexons-1]);
	  );

    FREEA(exonstatus);
    FREE(intronlengths);
    FREE(exonprotects);
    FREE(exonmatches);
    FREE(exonlengths);

  }

  debug(printf("Ending of smoothing\n\n"));

  return pairs;
}


#if 0
List_T
Smooth_mark_short_exons (List_T pairs, Pairpool_T pairpool, int stage2_indexsize) {
  int *exonstatus;
  int *exonlengths, *exonmatches, *intronlengths, nexons, nintrons;
  bool shortp;
  bool delete1p, delete2p;
#ifdef DEBUG
  int i;
#endif

  shortp = false;
  /* smooth_reset(pairs); */
  if (pairs != NULL) {
    /* Trim ends */
    exonlengths = get_exonlengths(&exonmatches,&nexons,pairs,/*bysizep*/true);
    intronlengths = get_intronlengths(&nintrons,pairs,/*bysizep*/true);

    debug(printf("\nFind end shorts\n"));
    exonstatus = (int *) MALLOCA(nexons * sizeof(int));
    find_end_shorts(&delete1p,exonstatus,exonmatches,nexons,intronlengths);
    pairs = mark_exons(pairs,
#ifdef WASTE
		       pairpool,
#endif
		       exonstatus);

    FREEA(exonstatus);

    FREE(intronlengths);
    FREE(exonmatches);
    FREE(exonlengths);
  }

  if (pairs != NULL) {
    /* Find internal shorts */
    exonlengths = get_exonlengths(&exonmatches,&exonprotects,&nexons,pairs,/*bysizep*/true);
    intronlengths = get_intronlengths(&nintrons,pairs,/*bysizep*/true);

    debug(
	  printf("Beginning of smoothing.  Initial structure:\n");
	  for (i = 0; i < nexons-1; i++) {
	    printf("Exon %d of length %d (%d matches, %d protects)\n",i,exonlengths[i],exonmatches[i],exonprotects[i]);
	    printf("Intron %d of length %d\n",i,intronlengths[i]);
	  }
	  printf("Exon %d of length %d (%d matches)\n",nexons-1,exonlengths[nexons-1],exonmatches[nexons-1]);
	  );

    debug(printf("\nFind internal shorts\n"));
    exonstatus = (int *) MALLOCA(nexons * sizeof(int));
    find_internal_shorts_by_size(&shortp,&delete2p,exonstatus,exonmatches,nexons,intronlengths,stage2_indexsize);
    debug(printf("\nMark internal shorts\n"));
    pairs = mark_exons(pairs,
#ifdef WASTE
		       pairpool,
#endif
		       exonstatus);

    debug(
	  printf("After marking internal shorts:\n");
	  for (i = 0; i < nexons-1; i++) {
	    printf("Exon %d of length %d (%d matches)\n",i,exonlengths[i],exonmatches[i]);
	    printf("Intron %d of length %d\n",i,intronlengths[i]);
	  }
	  printf("Exon %d of length %d\n",nexons-1,exonlengths[nexons-1]);
	  );

    FREEA(exonstatus);
    FREE(intronlengths);
    FREE(exonmatches);
    FREE(exonlengths);

  }

  debug(printf("End of Smooth_mark_short_exons\n\n"));

  return pairs;
}
#endif



List_T
Smooth_pairs_by_intronprobs (bool *badp, List_T pairs, Pairpool_T pairpool) {
  int *exonstatus;
  int *exonmatches, *exon_denominator, nexons, nintrons;
  int *intron_matches_left, *intron_denominator_left, *intron_matches_right, *intron_denominator_right;
  double *donor_probs, *acceptor_probs;
  int total_matches, total_denominator;
#ifdef DEBUG
  int *intronlengths;
  int i;
#endif

  debug(Pair_dump_list(pairs,true));

  *badp = false;
  /* smooth_reset(pairs); */
  if (pairs != NULL) {
    /* Remove internal shorts */
    get_intronprobs(&donor_probs,&acceptor_probs,&nintrons,pairs);
    if (nintrons > 0) {
      get_exonmatches(&total_matches,&total_denominator,&exonmatches,&exon_denominator,&nexons,pairs);
      debug(intronlengths = get_intronlengths(&nintrons,pairs,/*bysizep*/false));

      pairs = get_intron_neighborhoods(&intron_matches_left,&intron_denominator_left,
				       &intron_matches_right,&intron_denominator_right,nintrons,pairs);

      debug(
	    printf("Beginning of smoothing.  Initial structure:\n");
	    for (i = 0; i < nexons-1; i++) {
	      printf("Exon %d with %d matches and %d denominator\n",i,exonmatches[i],exon_denominator[i]);
	      printf("Intron %d, length %d, with probs %f+%f and matches %d/%d and %d/%d\n",
		     i,intronlengths[i],donor_probs[i],acceptor_probs[i],intron_matches_left[i],intron_denominator_left[i],
		     intron_matches_right[i],intron_denominator_right[i]);
	    }
	    printf("Exon %d with %d matches and %d denominator\n",nexons-1,exonmatches[nexons-1],exon_denominator[nexons-1]);
	    );

      debug(printf("\nFind internal bads\n"));
      exonstatus = (int *) MALLOCA(nexons * sizeof(int));
      find_internal_bads_by_prob(&(*badp),exonstatus,exonmatches,exon_denominator,nexons,
				 donor_probs,acceptor_probs,
				 intron_matches_left,intron_denominator_left,
				 intron_matches_right,intron_denominator_right,
				 total_matches,total_denominator);
      debug(printf("\nRemove internal bads\n"));
      if (*badp == false) {
	debug(printf("No internal bads\n"));
      } else {
	debug(print_exon_status(exonstatus,exonmatches,/*exonprotects*/NULL,nexons));
	pairs = delete_and_mark_exons(pairs,pairpool,exonstatus,/*markp*/true,/*bysizep*/false);
      }
    
#if 0
      /* This is not correct */
      debug(
	    printf("After removing internal shorts:\n");
	    for (i = 0; i < nexons-1; i++) {
	      printf("Exon %d of length %d (%d matches)\n",i,exonlengths[i],exonmatches[i]);
	      printf("Intron %d of length %d\n",i,intronlengths[i]);
	    }
	    printf("Exon %d of length %d\n",nexons-1,exonlengths[nexons-1]);
	    );
#endif

      debug(FREE(intronlengths));
      FREEA(exonstatus);
      FREE(intron_denominator_right);
      FREE(intron_matches_right);
      FREE(intron_denominator_left);
      FREE(intron_matches_left);
      FREE(exon_denominator);
      FREE(exonmatches);

      debug(printf("Ending of smoothing\n\n"));
    }

    FREE(acceptor_probs);
    FREE(donor_probs);
  }

  return pairs;
}
