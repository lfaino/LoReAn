static char rcsid[] = "$Id: indel.c 167164 2015-06-09 20:54:17Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "indel.h"

#include "assert.h"
#include "mem.h"
#include "genome128_hr.h"
#include "stage3hr.h"


/* Indels */ 
#ifdef DEBUG2
#define debug2(x) x
#else
#define debug2(x)
#endif


static int min_indel_end_matches;
static int indel_penalty_middle;


void
Indel_setup (int min_indel_end_matches_in, int indel_penalty_middle_in) {
  min_indel_end_matches = min_indel_end_matches_in;
  indel_penalty_middle = indel_penalty_middle_in;
  return;
}


/* Called only by sarray-read.c, where plusp is always true */
/* indels is positive here */
int
Indel_resolve_middle_insertion (int *best_nmismatches_i, int *best_nmismatches_j,
				Univcoord_T left, int indels, Compress_T query_compress,
				int querystart, int queryend, int querylength,
				int max_mismatches_allowed,
				bool plusp, int genestrand, bool first_read_p) {
  int best_indel_pos = -1, indel_pos;
#ifdef DEBUG2
  int i;
  char *gbuffer;
#endif
  int nmismatches_left, nmismatches_right;
  int best_sum, sum, nmismatches_lefti, nmismatches_righti, lefti, righti;
  int nmismatches1, nmismatches2;

#ifdef HAVE_ALLOCA
  int *mismatch_positions_left = (int *) ALLOCA(querylength * sizeof(int));
  int *mismatch_positions_right = (int *) ALLOCA(querylength * sizeof(int));
#else
  int mismatch_positions_left[MAX_READLENGTH], mismatch_positions_right[MAX_READLENGTH];

  if (max_mismatches_allowed > MAX_READLENGTH) {
    max_mismatches_allowed = MAX_READLENGTH;
  }
#endif


  /* query has insertion.  Get |indels| less from genome; trim from left. */
  /* left = ptr->diagonal - querylength; */

  assert(indels > 0);
  debug2(gbuffer = (char *) CALLOC(querylength-indels+1,sizeof(char)));
  debug2(Genome_fill_buffer_blocks(left+indels,querylength-indels,gbuffer));
  debug2(printf("solve_middle_indel, plus, insertion: Getting genome at diagonal - querylength %d + indels %d = %llu\n",
		querylength,indels,(unsigned long long) left+indels));
  debug2(printf("g1: %s\n",gbuffer));
  debug2(printf("g2: %s\n",&(gbuffer[indels])));

  /* No need to check chromosome bounds */
  debug2(printf("max_mismatches_allowed is %d\n",max_mismatches_allowed));
  nmismatches_left = Genome_mismatches_left(mismatch_positions_left,max_mismatches_allowed,
					    query_compress,left,/*pos5*/querystart,/*pos3*/queryend,
					    plusp,genestrand,first_read_p);

  debug2(
	 printf("%d mismatches on left at:",nmismatches_left);
	 for (i = 0; i <= nmismatches_left; i++) {
	   printf(" %d",mismatch_positions_left[i]);
	 }
	 printf("\n");
	 );


  /* No need to check chromosome bounds */
  debug2(printf("max_mismatches_allowed is %d\n",max_mismatches_allowed));
  nmismatches_right = Genome_mismatches_right(mismatch_positions_right,max_mismatches_allowed,
					      query_compress,left-indels,/*pos5*/querystart,/*pos3*/queryend,
					      plusp,genestrand,first_read_p);

  debug2(
	 printf("%d mismatches on right at:",nmismatches_right);
	 for (i = 0; i <= nmismatches_right; i++) {
	   printf(" %d",mismatch_positions_right[i]);
	 }
	 printf("\n");
	 );

  best_sum = querylength + querylength;

  /* Modeled after end D to get lowest possible coordinate */
  righti = 0;
  lefti = nmismatches_left - 1;
  nmismatches_righti = /*righti*/ 0;
  nmismatches_lefti = /*lefti+1*/ nmismatches_left;

  while (righti < nmismatches_right) {
    while (lefti >= 0 && mismatch_positions_left[lefti] > mismatch_positions_right[righti] - indels) {
      lefti--;
    }
    sum = righti + lefti + 1;
    debug2(printf("  (Case D) sum %d=%d+%d at indel_pos %d.",
		  sum,righti,lefti+1,mismatch_positions_right[righti]-indels+1));
    if (sum <= best_sum) {
      indel_pos = mismatch_positions_right[righti] - indels + 1;
      if (indel_pos >= min_indel_end_matches && indel_pos + indels <= querylength - min_indel_end_matches) {
	best_indel_pos = indel_pos;
	nmismatches_righti = righti;
	nmismatches_lefti = lefti + 1;
	debug2(printf("**"));
	best_sum = sum;
      }
    }
    righti++;
  }
  debug2(printf("\n"));


  /* Try from other side to see if we missed anything */
  lefti = 0;
  righti = nmismatches_right - 1;

  while (lefti < nmismatches_left) {
    while (righti >= 0 && mismatch_positions_right[righti] < mismatch_positions_left[lefti] + indels) {
      righti--;
    }
    sum = lefti + righti + 1;
    debug2(printf("  (Case D2) sum %d=%d+%d at indel_pos %d.",
		  sum,lefti,righti+1,mismatch_positions_left[lefti]));
    if (sum < best_sum) {
      indel_pos = mismatch_positions_left[lefti];
      if (indel_pos >= min_indel_end_matches && indel_pos + indels <= querylength - min_indel_end_matches) {
	best_indel_pos = indel_pos;
	nmismatches_righti = righti + 1;
	nmismatches_lefti = lefti;
	debug2(printf("**"));
	best_sum = sum;
      }
    } else if (sum == best_sum) {
      indel_pos = mismatch_positions_left[lefti];
      if (indel_pos < best_indel_pos) {
	if (indel_pos >= min_indel_end_matches && indel_pos + indels <= querylength - min_indel_end_matches) {
	  best_indel_pos = indel_pos;
	  nmismatches_righti = righti + 1;
	  nmismatches_lefti = lefti;
	  debug2(printf("**"));
	  /* best_sum = sum; */
	}
      }
    }
    lefti++;
  }
  debug2(printf("\n"));

  *best_nmismatches_i = nmismatches_lefti;
  *best_nmismatches_j = nmismatches_righti;

  if (best_sum > max_mismatches_allowed) {
    debug2(printf("Returning -1\n"));
    return -1;
#if 0
  } else if (plusp == true) {
    return best_indel_pos;
  } else {
    return querylength - best_indel_pos - indels;
#else
  } else {
    debug2(printf("Returning %d\n",best_indel_pos));
    return best_indel_pos;
#endif
  }
}


/* Called only by sarray-read.c, where plusp is always true */
/* indels is negative here */
int
Indel_resolve_middle_deletion (int *best_nmismatches_i, int *best_nmismatches_j,
			       Univcoord_T left, int indels, Compress_T query_compress,
			       int querystart, int queryend, int querylength,
			       int max_mismatches_allowed,
			       bool plusp, int genestrand, bool first_read_p) {
  int best_indel_pos = -1, indel_pos;
#ifdef DEBUG2
  int i;
  char *gbuffer;
#endif
  int nmismatches_left, nmismatches_right, nmismatches_lefti, nmismatches_righti;
  int best_sum, sum, lefti, righti;

#ifdef HAVE_ALLOCA
  int *mismatch_positions_left = (int *) ALLOCA(querylength * sizeof(int));
  int *mismatch_positions_right = (int *) ALLOCA(querylength * sizeof(int));
#else
  int mismatch_positions_left[MAX_READLENGTH], mismatch_positions_right[MAX_READLENGTH];

  if (max_mismatches_allowed > MAX_READLENGTH) {
    max_mismatches_allowed = MAX_READLENGTH;
  }
#endif


  /* query has deletion.  Get |indels| more from genome; add to right. */
  /* left = ptr->diagonal - querylength; */

  assert(indels < 0);
  debug2(gbuffer = (char *) CALLOC(querylength-indels+1,sizeof(char)));
  debug2(Genome_fill_buffer_blocks(left,querylength-indels,gbuffer));
  debug2(printf("solve_middle_indel, plus, deletion (indels %d), max_mismatches_allowed %d: Getting genome at diagonal - querylength %d = %llu\n",
		indels,max_mismatches_allowed,querylength,(unsigned long long) left));
  debug2(printf("g1: %s\n",gbuffer));
  debug2(printf("g2: %s\n",&(gbuffer[-indels])));
  debug2(FREE(gbuffer));

  /* No need to check chromosome bounds */
  nmismatches_left = Genome_mismatches_left(mismatch_positions_left,max_mismatches_allowed,
					    query_compress,left,/*pos5*/querystart,/*pos3*/queryend,
					    plusp,genestrand,first_read_p);

  debug2(
	 printf("%d mismatches on left at:",nmismatches_left);
	 for (i = 0; i <= nmismatches_left; i++) {
	   printf(" %d",mismatch_positions_left[i]);
	 }
	 printf("\n");
	 );

  /* No need to check chromosome bounds */
  nmismatches_right = Genome_mismatches_right(mismatch_positions_right,max_mismatches_allowed,
					      query_compress,left-indels,/*pos5*/querystart,/*pos3*/queryend,
					      plusp,genestrand,first_read_p);

  debug2(
	 printf("%d mismatches on right at:",nmismatches_right);
	 for (i = 0; i <= nmismatches_right; i++) {
	   printf(" %d",mismatch_positions_right[i]);
	 }
	 printf("\n");
	 );

  best_sum = querylength + querylength;

  /* Modeled after end C to get lowest possible coordinate */
  righti = 0;
  lefti = nmismatches_left - 1;
  nmismatches_righti = /*righti*/ 0;
  nmismatches_lefti = /*lefti+1*/ nmismatches_left;

  while (righti < nmismatches_right) {
    while (lefti >= 0 && mismatch_positions_left[lefti] > mismatch_positions_right[righti]) {
      lefti--;
    }
    sum = righti + lefti + 1;
    debug2(printf("  (Case C1) sum %d=%d+%d at indel_pos %d.",
		  sum,righti,lefti+1,mismatch_positions_right[righti]+1));
    if (sum <= best_sum) {
      indel_pos = mismatch_positions_right[righti] + 1;
      if (indel_pos >= min_indel_end_matches && indel_pos <= querylength - min_indel_end_matches) {
	best_indel_pos = indel_pos;
	nmismatches_righti = righti;
	nmismatches_lefti = lefti + 1;
	debug2(printf("**"));
	best_sum = sum;
      }
    }
    righti++;
  }
  debug2(printf("\n"));

  /* Try from other side to see if we missed anything */
  lefti = 0;
  righti = nmismatches_right - 1;

  while (lefti < nmismatches_left) {
    while (righti >= 0 && mismatch_positions_right[righti] < mismatch_positions_left[lefti]) {
      righti--;
    }
    sum = lefti + righti + 1;
    debug2(printf("  (Case C2) sum %d=%d+%d at indel_pos %d.",
		  sum,lefti,righti+1,mismatch_positions_left[lefti]));
    if (sum < best_sum) {
      indel_pos = mismatch_positions_left[lefti];
      if (indel_pos >= min_indel_end_matches && indel_pos <= querylength - min_indel_end_matches) {
	best_indel_pos = indel_pos;
	nmismatches_lefti = lefti;
	nmismatches_righti = righti + 1;
	debug2(printf("**"));
	best_sum = sum;
      }
    } else if (sum == best_sum) {
      indel_pos = mismatch_positions_left[lefti];
      if (indel_pos < best_indel_pos) {
	if (indel_pos >= min_indel_end_matches && indel_pos <= querylength - min_indel_end_matches) {
	  best_indel_pos = indel_pos;
	  nmismatches_lefti = lefti;
	  nmismatches_righti = righti + 1;
	  debug2(printf("**"));
	  /* best_sum = sum; */
	}
      }
    }
    lefti++;
  }
  debug2(printf("\n"));

  *best_nmismatches_i = nmismatches_lefti;
  *best_nmismatches_j = nmismatches_righti;

  if (best_sum > max_mismatches_allowed) {
    debug2(printf("Returning -1\n"));
    return -1;
#if 0
  } else if (plusp == true) {
    return best_indel_pos;
  } else {
    return querylength - best_indel_pos;
#else
  } else {
    debug2(printf("Returning %d\n",best_indel_pos));
    return best_indel_pos;
#endif
  }
}


/* indels is positive here */
List_T
Indel_solve_middle_insertion (bool *foundp, int *found_score, int *nhits, List_T hits,
			      Univcoord_T left, Chrnum_T chrnum, Univcoord_T chroffset,
			      Univcoord_T chrhigh, Chrpos_T chrlength,
			      int indels, Compress_T query_compress,
			      int querylength, int max_mismatches_allowed,
			      bool plusp, int genestrand, bool first_read_p, bool sarrayp) {
#ifdef DEBUG2
  int i;
  char *gbuffer;
#endif
  Stage3end_T hit;
  int best_indel_pos, query_indel_pos, indel_pos;
  int nmismatches_left, nmismatches_right;
  int best_sum, sum, nmismatches_lefti, nmismatches_righti, lefti, righti;
  int nmismatches1, nmismatches2;

#ifdef HAVE_ALLOCA
  int *mismatch_positions_left = (int *) ALLOCA(querylength * sizeof(int));
  int *mismatch_positions_right = (int *) ALLOCA(querylength * sizeof(int));
#else
  int mismatch_positions_left[MAX_READLENGTH], mismatch_positions_right[MAX_READLENGTH];
#endif


  *foundp = false;

  /* query has insertion.  Get |indels| less from genome; trim from left. */
  /* left = ptr->diagonal - querylength; */

  assert(indels > 0);
  debug2(gbuffer = (char *) CALLOC(querylength-indels+1,sizeof(char)));
  debug2(Genome_fill_buffer_blocks(left+indels,querylength-indels,gbuffer));
  debug2(printf("solve_middle_indel, plus, insertion: Getting genome at diagonal - querylength %d + indels %d = %llu\n",
		querylength,indels,(unsigned long long) left+indels));
  debug2(printf("g1: %s\n",gbuffer));
  debug2(printf("g2: %s\n",&(gbuffer[indels])));

  /* No need to check chromosome bounds */
  nmismatches_left = Genome_mismatches_left(mismatch_positions_left,max_mismatches_allowed,
					    query_compress,left+indels,/*pos5*/0,/*pos3*/querylength,
					    plusp,genestrand,first_read_p);

  debug2(
	 printf("%d mismatches on left at:",nmismatches_left);
	 for (i = 0; i <= nmismatches_left; i++) {
	   printf(" %d",mismatch_positions_left[i]);
	 }
	 printf("\n");
	 );


  /* No need to check chromosome bounds */
  nmismatches_right = Genome_mismatches_right(mismatch_positions_right,max_mismatches_allowed,
					      query_compress,left,/*pos5*/0,/*pos3*/querylength,
					      plusp,genestrand,first_read_p);

  debug2(
	 printf("%d mismatches on right at:",nmismatches_right);
	 for (i = 0; i <= nmismatches_right; i++) {
	   printf(" %d",mismatch_positions_right[i]);
	 }
	 printf("\n");
	 );

  best_sum = querylength + querylength;

  /* Modeled after end D to get lowest possible coordinate */
  righti = 0;
  lefti = nmismatches_left - 1;
  nmismatches_righti = /*righti*/ 0;
  nmismatches_lefti = /*lefti+1*/ nmismatches_left;

  while (righti < nmismatches_right) {
    while (lefti >= 0 && mismatch_positions_left[lefti] > mismatch_positions_right[righti] - indels) {
      lefti--;
    }
    sum = righti + lefti + 1;
    debug2(printf("  (Case D) sum %d=%d+%d at indel_pos %d.",
		  sum,righti,lefti+1,mismatch_positions_right[righti]-indels+1));
    if (sum <= best_sum) {
      indel_pos = mismatch_positions_right[righti] - indels + 1;
      if (indel_pos >= min_indel_end_matches && indel_pos + indels <= querylength - min_indel_end_matches) {
	best_indel_pos = indel_pos;
	nmismatches_righti = righti;
	nmismatches_lefti = lefti + 1;
	debug2(printf("**"));
	best_sum = sum;
      }
    }
    righti++;
  }
  debug2(printf("\n"));


  /* Try from other side to see if we missed anything */
  lefti = 0;
  righti = nmismatches_right - 1;

  while (lefti < nmismatches_left) {
    while (righti >= 0 && mismatch_positions_right[righti] < mismatch_positions_left[lefti] + indels) {
      righti--;
    }
    sum = lefti + righti + 1;
    debug2(printf("  (Case D2) sum %d=%d+%d at indel_pos %d.",
		  sum,lefti,righti+1,mismatch_positions_left[lefti]));
    if (sum < best_sum) {
      indel_pos = mismatch_positions_left[lefti];
      if (indel_pos >= min_indel_end_matches && indel_pos + indels <= querylength - min_indel_end_matches) {
	best_indel_pos = indel_pos;
	nmismatches_righti = righti + 1;
	nmismatches_lefti = lefti;
	debug2(printf("**"));
	best_sum = sum;
      }
    } else if (sum == best_sum) {
      indel_pos = mismatch_positions_left[lefti];
      if (indel_pos < best_indel_pos) {
	if (indel_pos >= min_indel_end_matches && indel_pos + indels <= querylength - min_indel_end_matches) {
	  best_indel_pos = indel_pos;
	  nmismatches_righti = righti + 1;
	  nmismatches_lefti = lefti;
	  debug2(printf("**"));
	  /* best_sum = sum; */
	}
      }
    }
    lefti++;
  }
  debug2(printf("\n"));

  if (best_sum <= max_mismatches_allowed) {
    if (plusp == true) {
      query_indel_pos = best_indel_pos;
      nmismatches1 = nmismatches_lefti;
      nmismatches2 = nmismatches_righti;
    } else {
      query_indel_pos = querylength - best_indel_pos - indels;
      nmismatches1 = nmismatches_righti;
      nmismatches2 = nmismatches_lefti;
    }

    if ((hit = Stage3end_new_insertion(&(*found_score),indels,query_indel_pos,
				       nmismatches1,nmismatches2,
				       /*left*/left+indels,/*genomiclength*/querylength-indels,
				       query_compress,querylength,plusp,genestrand,first_read_p,
				       chrnum,chroffset,chrhigh,chrlength,
				       indel_penalty_middle,sarrayp)) != NULL) {
      debug2(printf("successful insertion with %d=%d+%d mismatches and indel_pos at %d\n",
		    sum,nmismatches_lefti,nmismatches_righti,best_indel_pos));
      /* ptr->usedp = ptr2->usedp = true; */
      *foundp = true;
      *nhits += 1;
      hits = List_push(hits,(void *) hit);
    }
  }

  return hits;
}



/* indels is negative here */
List_T
Indel_solve_middle_deletion (bool *foundp, int *found_score, int *nhits, List_T hits,
			     Univcoord_T left, Chrnum_T chrnum, Univcoord_T chroffset,
			     Univcoord_T chrhigh, Chrpos_T chrlength,
			     int indels, Compress_T query_compress, int querylength,
			     int max_mismatches_allowed,
			     bool plusp, int genestrand, bool first_read_p, bool sarrayp) {
#ifdef DEBUG2
  int i;
  char *gbuffer;
#endif
  Stage3end_T hit;
  int best_indel_pos, query_indel_pos, indel_pos;
  int nmismatches_left, nmismatches_right;
  int best_sum, sum, nmismatches_lefti, nmismatches_righti, lefti, righti;
  int nmismatches1, nmismatches2;

#ifdef HAVE_ALLOCA
  int *mismatch_positions_left = (int *) ALLOCA(querylength * sizeof(int));
  int *mismatch_positions_right = (int *) ALLOCA(querylength * sizeof(int));
#else
  int mismatch_positions_left[MAX_READLENGTH], mismatch_positions_right[MAX_READLENGTH];
#endif


  *foundp = false;

  /* query has deletion.  Get |indels| more from genome; add to right. */
  /* left = ptr->diagonal - querylength; */

  assert(indels < 0);
  debug2(gbuffer = (char *) CALLOC(querylength-indels+1,sizeof(char)));
  debug2(Genome_fill_buffer_blocks(left,querylength-indels,gbuffer));
  debug2(printf("solve_middle_indel, plus, deletion (indels %d): Getting genome at diagonal - querylength %d = %llu\n",
		indels,querylength,(unsigned long long) left));
  debug2(printf("g1: %s\n",gbuffer));
  debug2(printf("g2: %s\n",&(gbuffer[-indels])));
  debug2(FREE(gbuffer));

  /* No need to check chromosome bounds */
  nmismatches_left = Genome_mismatches_left(mismatch_positions_left,max_mismatches_allowed,
					    query_compress,left,/*pos5*/0,/*pos3*/querylength,
					    plusp,genestrand,first_read_p);

  debug2(
	 printf("%d mismatches on left at:",nmismatches_left);
	 for (i = 0; i <= nmismatches_left; i++) {
	   printf(" %d",mismatch_positions_left[i]);
	 }
	 printf("\n");
	 );

  /* No need to check chromosome bounds */
  nmismatches_right = Genome_mismatches_right(mismatch_positions_right,max_mismatches_allowed,
					      query_compress,left-indels,/*pos5*/0,/*pos3*/querylength,
					      plusp,genestrand,first_read_p);

  debug2(
	 printf("%d mismatches on right at:",nmismatches_right);
	 for (i = 0; i <= nmismatches_right; i++) {
	   printf(" %d",mismatch_positions_right[i]);
	 }
	 printf("\n");
	 );

  best_sum = querylength + querylength;

  /* Modeled after end C to get lowest possible coordinate */
  righti = 0;
  lefti = nmismatches_left - 1;
  nmismatches_righti = /*righti*/ 0;
  nmismatches_lefti = /*lefti+1*/ nmismatches_left;

  while (righti < nmismatches_right) {
    while (lefti >= 0 && mismatch_positions_left[lefti] > mismatch_positions_right[righti]) {
      lefti--;
    }
    sum = righti + lefti + 1;
    debug2(printf("  (Case C1) sum %d=%d+%d at indel_pos %d.",
		  sum,righti,lefti+1,mismatch_positions_right[righti]+1));
    if (sum <= best_sum) {
      indel_pos = mismatch_positions_right[righti] + 1;
      if (indel_pos >= min_indel_end_matches && indel_pos <= querylength - min_indel_end_matches) {
	best_indel_pos = indel_pos;
	nmismatches_righti = righti;
	nmismatches_lefti = lefti + 1;
	debug2(printf("**"));
	best_sum = sum;
      }
    }
    righti++;
  }
  debug2(printf("\n"));

  /* Try from other side to see if we missed anything */
  lefti = 0;
  righti = nmismatches_right - 1;

  while (lefti < nmismatches_left) {
    while (righti >= 0 && mismatch_positions_right[righti] < mismatch_positions_left[lefti]) {
      righti--;
    }
    sum = lefti + righti + 1;
    debug2(printf("  (Case C2) sum %d=%d+%d at indel_pos %d.",
		  sum,lefti,righti+1,mismatch_positions_left[lefti]));
    if (sum < best_sum) {
      indel_pos = mismatch_positions_left[lefti];
      if (indel_pos >= min_indel_end_matches && indel_pos <= querylength - min_indel_end_matches) {
	best_indel_pos = indel_pos;
	nmismatches_lefti = lefti;
	nmismatches_righti = righti + 1;
	debug2(printf("**"));
	best_sum = sum;
      }
    } else if (sum == best_sum) {
      indel_pos = mismatch_positions_left[lefti];
      if (indel_pos < best_indel_pos) {
	if (indel_pos >= min_indel_end_matches && indel_pos <= querylength - min_indel_end_matches) {
	  best_indel_pos = indel_pos;
	  nmismatches_lefti = lefti;
	  nmismatches_righti = righti + 1;
	  debug2(printf("**"));
	  /* best_sum = sum; */
	}
      }
    }
    lefti++;
  }
  debug2(printf("\n"));


  if (best_sum <= max_mismatches_allowed) {
    if (plusp == true) {
      query_indel_pos = best_indel_pos;
      nmismatches1 = nmismatches_lefti;
      nmismatches2 = nmismatches_righti;
    } else {
      query_indel_pos = querylength - best_indel_pos;
      nmismatches1 = nmismatches_righti;
      nmismatches2 = nmismatches_lefti;
    }

    if ((hit = Stage3end_new_deletion(&(*found_score),-indels,query_indel_pos,
				      nmismatches1,nmismatches2,
				      left,/*genomiclength*/querylength-indels,
				      query_compress,querylength,plusp,genestrand,first_read_p,
				      chrnum,chroffset,chrhigh,chrlength,
				      indel_penalty_middle,sarrayp)) != NULL) {
      debug2(printf("successful middle deletion with %d=%d+%d mismatches and indel_pos at %d and nindels %d\n",
		    best_sum,nmismatches_lefti,nmismatches_righti,best_indel_pos,-indels));
      *foundp = true;
      /* ptr->usedp = ptr2->usedp = true; */
      *nhits += 1;
      hits = List_push(hits,(void *) hit);
    }
  }

  return hits;
}


