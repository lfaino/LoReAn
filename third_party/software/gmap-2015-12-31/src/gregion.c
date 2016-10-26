static char rcsid[] = "$Id: gregion.c 151507 2014-10-22 19:37:48Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "gregion.h"
#include <stdlib.h>
#include "assert.h"
#include "mem.h"
#include "indexdb.h"		/* For SUFFICIENT_SUPPORT */
#include "univinterval.h"
#include "chrnum.h"
#include "listdef.h"
#include "uinttable.h"


#define MAX_GENOMICLENGTH 2000000
#define MAX_NCHRS_FOR_ALLOCA 100
#define MAX_GREGIONS_FOR_ALLOCA 100


#define EXTRA_SHORTEND  30000
#define EXTRA_LONGEND   100000

#define USE_CLEAN 1

#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif

/* filter_clean */
#ifdef DEBUG2
#define debug2(x) x
#else
#define debug2(x)
#endif

/* Sufficient support */
#ifdef DEBUG5
#define debug5(x) x
#else
#define debug5(x)
#endif


#define T Gregion_T
struct T {
  int nexons;

  Chrpos_T genomiclength;
  bool plusp;
  int genestrand;

  Chrpos_T extension5;
  Chrpos_T extension3;
  bool extendedp;

  Chrnum_T chrnum;

  Chrpos_T chrstart;
  Chrpos_T chrend;

  Chrpos_T extentstart;
  Chrpos_T extentend;

  Univcoord_T chroffset;	/* This is for chr, not the segment */
  Univcoord_T chrhigh;
  Chrpos_T chrlength;

  int querystart;		/* Used only for maponly mode */
  int queryend;			/* Used only for maponly mode */
  int matchsize;		/* Used only for maponly mode */

  int trimstart;
  int trimend;
  bool sufficient_support_p;

  double weight;		/* Product of match weights */
  int support;

  int ncovered;
  int source;			/* Oligoindex source in stage 2 */

  /* For cleaning */
  

  int total_up;
  int total_down;
  Chrpos_T bestprev_ceiling;
  Chrpos_T bestprev_floor;
  int score_ceiling;
  int score_floor;

#ifdef USE_CLEAN
  bool bounded_low_p;
  bool bounded_high_p;
#endif
};



void
Gregion_print (T this) {

#if 1
  /* Off for debugging */
  printf(" %d..%d ",this->querystart,this->queryend);
#endif
  printf("%u %u %u %u #%d:%u..%u  length:%u  weight:%.2f  support:%d",
	 this->extentstart,this->extentend,this->chrstart,this->chrend,
	 this->chrnum,this->chrstart,this->chrend,this->genomiclength,
	 this->weight,this->support);
#ifdef USE_CLEAN
  printf("  bounded_low:%d, bounded_high:%d",this->bounded_low_p,this->bounded_high_p);
#endif

  printf("\n");

  return;
}


void
Gregion_free (T *old) {
  FREE(*old);
  return;
}

Univcoord_T
Gregion_genomicstart (T this) {
  return this->chroffset + this->chrstart;
}

Univcoord_T
Gregion_genomicend (T this) {
  return this->chroffset + this->chrend;
}

Chrpos_T
Gregion_chrstart (T this) {
  return this->chrstart;
}

Chrpos_T
Gregion_chrend (T this) {
  return this->chrend;
}

Chrpos_T
Gregion_genomiclength (T this) {
  return this->genomiclength;
}

bool
Gregion_plusp (T this) {
  return this->plusp;
}

bool
Gregion_revcompp (T this) {
  return !(this->plusp);
}

int
Gregion_genestrand (T this) {
  return this->genestrand;
}

Chrnum_T
Gregion_chrnum (T this) {
  return this->chrnum;
}

/* Used only for debugging.  String is allocated and should be freed. */
char *
Gregion_chr (T this, Univ_IIT_T chromosome_iit) {
  return Chrnum_to_string(this->chrnum,chromosome_iit);
}

Univcoord_T
Gregion_chroffset (T this) {
  return this->chroffset;
}

Univcoord_T
Gregion_chrhigh (T this) {
  return this->chrhigh;
}

Chrpos_T
Gregion_chrlength (T this) {
  return this->chrlength;
}

int
Gregion_querystart (T this) {
  return this->querystart;
}

int
Gregion_queryend (T this) {
  return this->queryend;
}

int
Gregion_matchsize (T this) {
  return this->matchsize;
}

double
Gregion_weight (T this) {
  return this->weight;
}

int
Gregion_support (T this) {
  return this->support;
}

bool 
Gregion_extendedp (T this) {
  return this->extendedp;
}

void
Gregion_set_ncovered (T this, int ncovered, int source) {
  this->ncovered = ncovered;
  this->source = source;
  return;
}

int
Gregion_ncovered (T this) {
  return this->ncovered;
}



T
Gregion_new (int nexons, Univcoord_T genomicstart, Univcoord_T genomicend,
	     bool plusp, int genestrand, Univ_IIT_T chromosome_iit, int querystart, int queryend, int querylength,
	     int matchsize, int trimstart, int trimend, int circular_typeint) {
  T new = (T) MALLOC(sizeof(*new));

  debug(printf("Creating gregion with genomicstart %u, genomicend %u, trimstart %d, trimend %d\n",
	       genomicstart,genomicend,trimstart,trimend));

  new->nexons = nexons;
  if (chromosome_iit == NULL) {
    /* Should not reach here, since user_genomicseg does not call stage 1 */
    new->chrnum = 0;
    new->chroffset = 0U;
    new->chrlength = 0U;
  } else {
    new->chrnum = Univ_IIT_get_one(chromosome_iit,genomicstart,genomicstart);
    /* new->chroffset = Univinterval_low(Univ_IIT_interval(chromosome_iit,new->chrnum)); */
    /* new->chrhigh = Univinterval_high(Univ_IIT_interval(chromosome_iit,new->chrnum)); */
    Univ_IIT_interval_bounds(&new->chroffset,&new->chrhigh,&new->chrlength,
			     chromosome_iit,new->chrnum,circular_typeint);
  }
  
  assert(genomicstart < genomicend);
  new->genomiclength = genomicend - genomicstart;

  new->chrstart = genomicstart - new->chroffset;
  new->chrend = new->chrstart + new->genomiclength;

  new->plusp = plusp;
  new->genestrand = genestrand;

  if (plusp == true) {
    if (new->chrstart < (Chrpos_T) querystart) {
      /* Extends past beginning of genome */
      new->extentstart = 0;
    } else {
      new->extentstart = new->chrstart - querystart;
    }
    new->extentend = new->chrend + (querylength - queryend);
  } else {
    if (new->chrstart < (Chrpos_T) (querylength - queryend)) {
      /* Extends past beginning of genome */
      new->extentstart = 0;
    } else {
      new->extentstart = new->chrstart - (querylength - queryend);
    }
    new->extentend = new->chrend + querystart;
  }

  new->extension5 = 0U;
  new->extension3 = 0U;
  new->extendedp = false;

  new->querystart = querystart;
  new->queryend = queryend;
  new->matchsize = matchsize;

  new->trimstart = trimstart;
  new->trimend = trimend;

#ifdef PMAP
  debug5(printf("  Testing bound5+extension5 = %d - %u < %d + %d, bound3+extension3 = %d + %u > %d - %d\n",
		new->querystart,new->extension5,trimstart,SUFFICIENT_SUPPORT,
		new->queryend,new->extension3,trimend,SUFFICIENT_SUPPORT));
  if (new->querystart - new->extension5 < trimstart + SUFFICIENT_SUPPORT && 
      new->queryend + new->extension3 > trimend - SUFFICIENT_SUPPORT) {
    new->sufficient_support_p = true;
  } else {
    new->sufficient_support_p = false;
  }
#else
  debug5(printf("  Testing bound5+extension5 = %d - %u < %d + %d, bound3 = %d + %u > %d - %d\n",
		new->querystart,new->extension5,trimstart,SUFFICIENT_SUPPORT,
		new->queryend,new->extension3,trimend,SUFFICIENT_SUPPORT));
  if (new->querystart - new->extension5 < trimstart + SUFFICIENT_SUPPORT && 
      new->queryend + new->extension3 > trimend - SUFFICIENT_SUPPORT) {
    new->sufficient_support_p = true;
  } else {
    new->sufficient_support_p = false;
  }
#endif

  new->weight = 0.0;
  new->support = queryend - querystart + matchsize;

  new->total_up = 0;
  new->total_down = 0;
  new->bestprev_ceiling = 0;
  new->bestprev_floor = (Chrpos_T) -1;
  new->score_ceiling = 0;
  new->score_floor = 0;

#ifdef USE_CLEAN
  new->bounded_low_p = false;
  new->bounded_high_p = false;
#endif

  return new;
}


T
Gregion_new_from_matches (Match_T match5, Match_T match3, int genestrand, Univ_IIT_T chromosome_iit, 
			  int querylength, int matchsize, int trimstart, int trimend, int circular_typeint) {
  T gregion;
  Univcoord_T genomicstart, genomicend;
  int querystart, queryend;

  querystart = Match_querypos(match5);
  queryend = Match_querypos(match3);

  if (Match_forwardp(match5)) {
    genomicstart = Match_position(match5);
    genomicend = Match_position(match3) + 1U;
    /* chrnum = Match_chrnum(match5); */

  } else {
    genomicstart = Match_position(match3);
    genomicend = Match_position(match5) + 1U;
    /* chrnum = Match_chrnum(match3); */
  }

#if 0
  if (chromosome_iit == NULL) {
    chroffset = 0U;
  } else {
    chroffset = Univinterval_low(Univ_IIT_interval(chromosome_iit,chrnum));
  }
#endif

  debug(printf("Coordinates are %u .. %u\n",genomicstart,genomicend));

  gregion = Gregion_new(/*nexons*/0,genomicstart,genomicend,Match_forwardp(match5),genestrand,
			chromosome_iit,querystart,queryend,querylength,matchsize,trimstart,trimend,
			circular_typeint);

  gregion->weight = Match_weight(match5) * Match_weight(match3);
  Match_incr_npairings(match5);
  Match_incr_npairings(match3);

  return gregion;
}



static int
weight_cmp (const void *x, const void *y) {
  T a = * (T *) x;
  T b = * (T *) y;

  if (a->weight > b->weight) {
    return -1;
  } else if (a->weight < b->weight) {
    return +1;
  } else if (a->support > b->support) {
    return -1;
  } else if (a->support < b->support) {
    return +1;
  } else if (a->genomiclength > b->genomiclength) {
    return -1;
  } else if (a->genomiclength < b->genomiclength) {
    return +1;
  } else {
    return 0;
  }
}


/* Not intended for qsort.  Returns 0 when not comparable. */
static bool
gregion_overlap_p (T x, T y) {
  Univcoord_T x_genomicstart, x_genomicend, y_genomicstart, y_genomicend;
  Univcoord_T overlap;
  double fraction;

  if (x->plusp != y->plusp) {
    return false;			/* Different strands */

  } else {
    x_genomicstart = x->chroffset + x->chrstart;
    x_genomicend = x->chroffset + x->chrend;
    y_genomicstart = y->chroffset + y->chrstart;
    y_genomicend = y->chroffset + y->chrend;

    if (y_genomicstart > x_genomicend || x_genomicstart > y_genomicend) {
      /*
      /-- x --/ /-- y --/ or /-- y --/ /-- x --/
      */
      return false;		/* No overlap */

    } else if (y_genomicstart < x_genomicstart) {
      debug(printf("x %u..%u, y %u..%u",x_genomicstart,x_genomicend,y_genomicstart,y_genomicend));
      if (y_genomicend < x_genomicend) {
	/*
	    /-- x --/
	/-- y --/
	*/
	overlap = y_genomicend - x_genomicstart;
	if (y_genomicend - y_genomicstart < x_genomicend - x_genomicstart) {
	  fraction = (double) overlap/(double) (y_genomicend - y_genomicstart);
	} else {
	  fraction = (double) overlap/(double) (x_genomicend - x_genomicstart);
	}
	debug(printf(" => fraction %f",fraction));
	if (fraction > 0.5) {
	  debug(printf(" => overlap\n",fraction));
	  return true;
	} else {
	  debug(printf(" => no overlap\n",fraction));
	  return false;
	}

      } else {
	/*
	    /-- x --/
	/----- y -----/
	*/
	debug(printf(" => subsumption\n"));
	return true;
      }

    } else {
      debug(printf("x %u..%u, y %u..%u\n",x_genomicstart,x_genomicend,y_genomicstart,y_genomicend));
      if (y_genomicend < x_genomicend) {
	/*
	/----- x -----/
	  /-- y --/
	*/
	debug(printf(" => subsumption\n"));
	return true;

      } else {
	/*
	/-- x --/
	  /-- y --/
	*/
	overlap = x_genomicend - y_genomicstart;
	if (y_genomicend - y_genomicstart < x_genomicend - x_genomicstart) {
	  fraction = (double) overlap/(double) (y_genomicend - y_genomicstart);
	} else {
	  fraction = (double) overlap/(double) (x_genomicend - x_genomicstart);
	}
	debug(printf(" => fraction %f",fraction));
	if (fraction > 0.5) {
	  debug(printf(" => overlap\n",fraction));
	  return true;
	} else {
	  debug(printf(" => no overlap\n",fraction));
	  return false;
	}

      }
    }
  }
}


List_T
Gregion_filter_unique (List_T gregionlist) {
  List_T unique = NULL;
  T x, y, gregion, *array;
  int n, i, j;
  bool *eliminate;
#ifdef DEBUG
  List_T p, q;
#endif

  n = List_length(gregionlist);
  if (n == 0) {
    return NULL;
  }

  debug(
	for (p = gregionlist, i = 0; p != NULL; p = p->rest, i++) {
	  gregion = (T) p->first;
	  printf("  Initial %d: %d..%d #%d:%u-%u (plusp = %d)\n",
		 i,gregion->querystart,gregion->queryend,gregion->chrnum,gregion->chrstart,gregion->chrend,gregion->plusp);
	}
	);

  if (n < MAX_GREGIONS_FOR_ALLOCA) {
    eliminate = (bool *) CALLOCA(n,sizeof(bool));
    array = (T *) MALLOCA(n * sizeof(T));
    List_fill_array_and_free((void **) array,&gregionlist);
  } else {
    eliminate = (bool *) CALLOC(n,sizeof(bool));
    array = (T *) List_to_array(gregionlist,NULL);
    List_free(&gregionlist);
  }

  qsort(array,n,sizeof(T),weight_cmp);

  for (i = 0; i < n; i++) {
    x = array[i];
    for (j = i+1; j < n; j++) {
      y = array[j];
      if (gregion_overlap_p(x,y) == true) {
#ifdef DEBUG
	printf("Found overlap between these regions:\n");
	printf("   ");
	Gregion_print(x);
	printf("   ");
	Gregion_print(y);
	printf("\n");
#endif
	eliminate[j] = true;
      }
    }
  }

  for (i = n-1; i >= 0; i--) {
    gregion = array[i];
    if (eliminate[i] == false) {
      debug(printf("  Keeping #%d:%u-%u (plusp = %d)\n",
		   gregion->chrnum,gregion->chrstart,gregion->chrend,gregion->plusp));
      unique = List_push(unique,(void *) gregion);
    } else {
      debug(printf("  Eliminating #%d:%u-%u (plusp = %d)\n",
		   gregion->chrnum,gregion->chrstart,gregion->chrend,gregion->plusp));
      /*
      if (gregion->match5 != NULL) {
	Match_decr_npairings(gregion->match5);
	Match_decr_npairings(gregion->match3);
      }
      */
      Gregion_free(&gregion);
    }
  }

  if (n < MAX_GREGIONS_FOR_ALLOCA) {
    FREEA(eliminate);
    FREEA(array);
  } else {
    FREE(eliminate);
    FREE(array);
  }

  debug(
	for (p = unique, i = 0; p != NULL; p = p->rest, i++) {
	  gregion = (T) p->first;
	  printf("  Final %d: #%d:%u-%u (plusp = %d)\n",
		 i,gregion->chrnum,gregion->chrstart,gregion->chrend,gregion->plusp);
	}
	);

  return unique;
}


List_T
Gregion_filter_support (List_T gregionlist, int boundary, double pct_max, int diff_max) {
  List_T good = NULL, p;
  int threshold, maxsupport = 0;
  T gregion;

  if (gregionlist == NULL) {
    return NULL;
  }

  for (p = gregionlist; p != NULL; p = List_next(p)) {
    gregion = (T) List_head(p);
    if (gregion->support > maxsupport) {
      maxsupport = gregion->support;
    }
  }

  if (maxsupport > boundary) {
    /* Use diff_max */
    for (p = gregionlist; p != NULL; p = List_next(p)) {
      gregion = (T) List_head(p);
      if (maxsupport - gregion->support < diff_max) {
	good = List_push(good,gregion);
      } else {
	Gregion_free(&gregion);
      }
    }
  } else {
    threshold = (int) ((double) maxsupport * pct_max);
    /* Use threshold only */
    for (p = gregionlist; p != NULL; p = List_next(p)) {
      gregion = (T) List_head(p);
      if (gregion->support >= threshold) {
	good = List_push(good,gregion);
      } else {
	Gregion_free(&gregion);
      }
    }
  }

  List_free(&gregionlist);
  return List_reverse(good);
}


double
Gregion_best_weight (List_T gregionlist) {
  double best_weight = 0.0;
  T gregion;
  List_T p;

  for (p = gregionlist; p != NULL; p = List_next(p)) {
    gregion = (T) List_head(p);
    if (gregion->weight > best_weight) {
      best_weight = gregion->weight;
    }
  }

  return best_weight;
}


bool
Gregion_sufficient_support (T this) {
  return this->sufficient_support_p;
}


void
Gregion_extend (T this, Chrpos_T extension5, Chrpos_T extension3, int querylength,
		int min_extra_end) {
  Chrpos_T left, right;
  int extra_end;

  debug(printf("Entering Gregion_extend with extension5 %u and extension3 %u\n",extension5,extension3));
  debug(printf("  #%d:%u..%u\n",this->chrnum,this->chrstart,this->chrlength));
  this->extension5 = extension5;
  this->extension3 = extension3;
  this->extendedp = true;


  if (this->nexons == 1 || Gregion_sufficient_support(this) == true || this->support < 100) {
    extra_end = EXTRA_SHORTEND;
#if 0
    /* Should no longer be necessary for known splicesites */
    if (extra_end < min_extra_end) {
      extra_end = min_extra_end;
    }
#endif
    if (this->plusp == true) {
      left = extension5 + querylength + extra_end;
      right = extension3 + querylength + extra_end;
    } else {
      left = extension3 + querylength + extra_end;
      right = extension5 + querylength + extra_end;
    }
  } else {
    extra_end = EXTRA_LONGEND;
#if 0
    /* Should no longer be necessary for known splicesites */
    if (extra_end < min_extra_end) {
      extra_end = min_extra_end;
    }
#endif
    if (this->plusp == true) {
      left = extension5 + extra_end;
      right = extension3 + extra_end;
    } else {
      left = extension3 + extra_end;
      right = extension5 + extra_end;
    }
  }

  /* printf("chrstart %u vs left %u\n",this->chrstart,left); */
  if (this->chrstart < left) {
    /* At beginning of chromosome */
    this->chrstart = 0U;
  } else {
    this->chrstart -= left;
  }

  /* printf("chroffset %u + chrend %u + right %u vs chrhigh %u\n",this->chroffset,this->chrend,right,this->chrhigh); */
  if (this->chroffset + this->chrend + right >= this->chrhigh) {
    /* At end of chromosome */
    this->chrend = this->chrlength - 1;
  } else {
    this->chrend += right;
  }

  /* Prevent very large genomic segments */
  if (this->chrend > this->chrstart + MAX_GENOMICLENGTH) {
    this->chrend = this->chrstart + MAX_GENOMICLENGTH;
  }

  this->genomiclength = this->chrend - this->chrstart + 1U;
  /* printf("chrstart %u, chrend %u, genomiclength %u\n",this->chrstart,this->chrend,this->genomiclength); */

#ifdef PMAP
  debug5(printf("  Testing bound5+extension5 = %d - %u < %d + %d, bound3+extension3 = %d + %u > %d - %d\n",
		this->querystart,this->extension5,this->trimstart,SUFFICIENT_SUPPORT,
		this->queryend,this->extension3,this->trimend,SUFFICIENT_SUPPORT));
  if (this->querystart - this->extension5 < this->trimstart + SUFFICIENT_SUPPORT && 
      this->queryend + this->extension3 > this->trimend - SUFFICIENT_SUPPORT) {
    this->sufficient_support_p = true;
  } else {
    this->sufficient_support_p = false;
  }
#else
  debug5(printf("  Testing bound5+extension5 = %d - %u < %d + %d, bound3 = %d + %u > %d - %d\n",
		this->querystart,this->extension5,this->trimstart,SUFFICIENT_SUPPORT,
		this->queryend,this->extension3,this->trimend,SUFFICIENT_SUPPORT));
  if (this->querystart - this->extension5 < this->trimstart + SUFFICIENT_SUPPORT && 
      this->queryend + this->extension3 > this->trimend - SUFFICIENT_SUPPORT) {
    this->sufficient_support_p = true;
  } else {
    this->sufficient_support_p = false;
  }
#endif

  debug(printf("  #%d:%u..%u\n",this->chrnum,this->chrstart,this->chrend));

  return;
}


int
Gregion_cmp (const void *a, const void *b) {
  T x = * (T *) a;
  T y = * (T *) b;

  if (x->ncovered > y->ncovered) {
    return -1;
  } else if (y->ncovered > x->ncovered) {
    return +1;
  } else {
    return 0;
  }
}


/************************************************************************
 *  Filtering, based on spliceclean
 ************************************************************************/

#ifdef USE_CLEAN

#define MAXLOOKBACK 60

static void
compute_total_up (List_T gregions) {
  long int total_up;
  List_T p;
  Gregion_T gregion;

  total_up = 0;
  for (p = gregions; p != NULL; p = List_next(p)) {
    gregion = (Gregion_T) List_head(p);
    total_up += 1;
    gregion->total_up = total_up;
  }

  return;
}

static void
compute_total_down (List_T gregions) {
  long int total_down;
  List_T p;
  Gregion_T gregion;

  total_down = 0;
  for (p = gregions; p != NULL; p = List_next(p)) {
    gregion = (Gregion_T) List_head(p);
    total_down += 1;
    gregion->total_down = total_down;
  }

  return;
}


static int
Gregion_low_descending_cmp (const void *a, const void *b) {
  Gregion_T x = * (Gregion_T *) a;
  Gregion_T y = * (Gregion_T *) b;

  if (x->extentstart > y->extentstart) {
    return -1;
  } else if (y->extentstart > x->extentstart) {
    return +1;
  } else {
    return 0;
  }
}

static int
Gregion_high_ascending_cmp (const void *a, const void *b) {
  Gregion_T x = * (Gregion_T *) a;
  Gregion_T y = * (Gregion_T *) b;

  if (x->extentend < y->extentend) {
    return -1;
  } else if (y->extentend < x->extentend) {
    return +1;
  } else {
    return 0;
  }
}

static List_T
Gregion_sort_low_descending (List_T gregions) {
  List_T sorted = NULL;
  Gregion_T *array;
  int n, i;

  if ((n = List_length(gregions)) == 0) {
    return (List_T) NULL;

  } else if (n < MAX_GREGIONS_FOR_ALLOCA) {
    array = (Gregion_T *) MALLOCA(n * sizeof(Gregion_T));
    List_fill_array((void **) array,gregions);

    qsort(array,n,sizeof(Gregion_T),Gregion_low_descending_cmp);
    for (i = n-1; i >= 0; i--) {
      sorted = List_push(sorted,array[i]);
    }

    FREEA(array);
    return sorted;

  } else {
    array = (Gregion_T *) List_to_array(gregions,NULL);

    qsort(array,n,sizeof(Gregion_T),Gregion_low_descending_cmp);
    for (i = n-1; i >= 0; i--) {
      sorted = List_push(sorted,array[i]);
    }

    FREE(array);
    return sorted;
  }
}



static List_T
Gregion_sort_high_ascending (List_T gregions) {
  List_T sorted = NULL;
  Gregion_T *array;
  int n, i;

  if ((n = List_length(gregions)) == 0) {
    return (List_T) NULL;

  } else if (n < MAX_GREGIONS_FOR_ALLOCA) {
    array = (Gregion_T *) MALLOCA(n * sizeof(Gregion_T));
    List_fill_array((void **) array,gregions);

    qsort(array,n,sizeof(Gregion_T),Gregion_high_ascending_cmp);
    for (i = n-1; i >= 0; i--) {
      sorted = List_push(sorted,array[i]);
    }

    FREEA(array);
    return sorted;

  } else {
    array = (Gregion_T *) List_to_array(gregions,NULL);

    qsort(array,n,sizeof(Gregion_T),Gregion_high_ascending_cmp);
    for (i = n-1; i >= 0; i--) {
      sorted = List_push(sorted,array[i]);
    }

    FREE(array);
    return sorted;
  }
}



typedef struct Base_T *Base_T;
struct Base_T {
  Chrpos_T prevpos;
  Chrpos_T minextent;
  Chrpos_T maxextent;
  List_T gregions;

  bool usedp;
};

#if 0
static void
Base_free (Base_T *old) {
  if ((*old)->gregions != NULL) {
    List_free(&(*old)->gregions);
  }
  FREE(*old);
  return;
}
#endif


static Base_T
Base_new () {
  Base_T new = (Base_T) MALLOC(sizeof(*new));

  new->minextent = -1U;
  new->maxextent = -1U;
  new->gregions = (List_T) NULL;
  new->usedp = false;

  return new;
}



/* Assumes gregions are arranged low to high */
static Gregion_T
apply_ceiling (List_T gregions, Chrpos_T ceiling) {
  List_T p;
  Gregion_T prevgregion = NULL, gregion;

  for (p = gregions; p != NULL; p = List_next(p)) {
    gregion = (Gregion_T) List_head(p);
    if (gregion->extentend > ceiling) {
      return prevgregion;
    }
    prevgregion = gregion;
  }

  return prevgregion;
}

/* Assumes gregions are arranged high to low */
static Gregion_T
apply_floor (List_T gregions, Chrpos_T floor) {
  List_T p;
  Gregion_T prevgregion = NULL, gregion;

  for (p = gregions; p != NULL; p = List_next(p)) {
    gregion = (Gregion_T) List_head(p);
    if (gregion->extentstart < floor) {
      return prevgregion;
    }
    prevgregion = gregion;
  }

  return prevgregion;
}


static Chrpos_T
compute_ceilings (Uinttable_T low_basetable) {
  Chrpos_T ceiling, bestprevpos, prevpos;
  long int bestscore, score;
  int nlookback;
  Chrpos_T *keys;
  Base_T base, prevbase;
  Gregion_T gregion, prevgregion;
#ifdef DEBUG2
  Gregion_T bestprevgregion;
#endif
  List_T p;
  int n, i;
  
  n = Uinttable_length(low_basetable);
  if (n < MAX_NCHRS_FOR_ALLOCA) {
    keys = (Chrpos_T *) MALLOCA(n * sizeof(Chrpos_T));
    Uinttable_fill_keys(keys,low_basetable,/*sortp*/true);
  } else {
    keys = Uinttable_keys(low_basetable,/*sortp*/true);
  }
  debug2(printf("low_basetable has %d entries\n",n));
  
  prevpos = 0U;
  for (i = 0; i < n; i++) {
    base = (Base_T) Uinttable_get(low_basetable,keys[i]);
    base->gregions = Gregion_sort_high_ascending(base->gregions);
    compute_total_up(base->gregions);

    base->prevpos = prevpos;
    prevpos = keys[i];
  }

  /* Initialize minbaselow */
  base = (Base_T) Uinttable_get(low_basetable,keys[0]);
  for (p = base->gregions; p != NULL; p = List_next(p)) {
    gregion = (Gregion_T) List_head(p);
    gregion->score_ceiling = gregion->total_up;
    gregion->bestprev_ceiling = 0U;
  }

  for (i = 1; i < n; i++) {
    base = (Base_T) Uinttable_get(low_basetable,keys[i]);
    
    debug2(printf("At base %u, have %d gregions\n",keys[i],List_length(base->gregions)));

    for (p = base->gregions; p != NULL; p = List_next(p)) {
      gregion = (Gregion_T) List_head(p);

      bestscore = 0;
      bestprevpos = keys[0];
      debug2(bestprevgregion = NULL);

      prevpos = base->prevpos;
      nlookback = 0;

      ceiling = gregion->extentend;
      debug2(printf("  Gregion %u",gregion->extentend));
      while (prevpos >= keys[0] && nlookback < MAXLOOKBACK) {
	prevbase = (Base_T) Uinttable_get(low_basetable,prevpos);
	if ((prevgregion = apply_ceiling(prevbase->gregions,ceiling)) != NULL) {
	  debug2(printf(" ... prev:%u score:%d+%d = %d",
			prevpos,prevgregion->score_ceiling,gregion->total_up,
			prevgregion->score_ceiling+gregion->total_up));

	  if ((score = prevgregion->score_ceiling + gregion->total_up) > bestscore) {
	    debug2(printf("*"));
	    bestscore = score;
	    bestprevpos = prevpos;
#ifdef DEBUG2
	    bestprevgregion = prevgregion;
#endif
	  }
	}

	prevpos = prevbase->prevpos;
	nlookback++;
      }
      debug2(printf("\n"));

      gregion->score_ceiling = bestscore;
      gregion->bestprev_ceiling = bestprevpos;

#ifdef DEBUG2
      if (bestprevgregion == NULL) {
	printf(" no prevgregion, bestprev is %u, score is %d\n",
	       bestprevpos,gregion->score_ceiling);
      } else {
	printf(" bestprev is pos %u, gregion %u, score is %d\n",
	       bestprevpos,bestprevgregion->extentend,gregion->score_ceiling);
      }
#endif
    }
  }

  /* Get best overall gregion */
  bestscore = 0;
  bestprevpos = 0U;

  for (i = 0; i < n; i++) {
    base = (Base_T) Uinttable_get(low_basetable,keys[i]);
    for (p = base->gregions; p != NULL; p = List_next(p)) {
      gregion = (Gregion_T) List_head(p);
      
      if (gregion->score_ceiling > bestscore) {
	bestscore = gregion->score_ceiling;
	bestprevpos = keys[i];
      }
    }
  }

  if (n < MAX_NCHRS_FOR_ALLOCA) {
    FREEA(keys);
  } else {
    FREE(keys);
  }

  return bestprevpos;
}



static Chrpos_T
compute_floors (Uinttable_T high_basetable) {
  Chrpos_T floor, bestprevpos, prevpos;
  long int bestscore, score;
  int nlookback;
  Chrpos_T *keys;
  Base_T base, prevbase;
  Gregion_T gregion, prevgregion;
#ifdef DEBUG2
  Gregion_T bestprevgregion;
#endif
  List_T p;
  int n, i;

  n = Uinttable_length(high_basetable);
  if (n < MAX_NCHRS_FOR_ALLOCA) {
    keys = (Chrpos_T *) MALLOCA(n * sizeof(Chrpos_T));
    Uinttable_fill_keys(keys,high_basetable,/*sortp*/true);
  } else {
    keys = Uinttable_keys(high_basetable,/*sortp*/true);
  }
  debug2(printf("high_basetable has %d entries\n",n));

  prevpos = (Chrpos_T) -1U;
  for (i = n-1; i >= 0; --i) {
    base = (Base_T) Uinttable_get(high_basetable,keys[i]);
    base->gregions = Gregion_sort_low_descending(base->gregions);
    compute_total_down(base->gregions);

    base->prevpos = prevpos;
    prevpos = keys[i];
  }

  /* Initialize maxbasehigh */
  base = (Base_T) Uinttable_get(high_basetable,keys[n-1]);
  for (p = base->gregions; p != NULL; p = List_next(p)) {
    gregion = (Gregion_T) List_head(p);
    gregion->score_floor = gregion->total_down;
    gregion->bestprev_floor = (Chrpos_T) -1U;
  }

  for (i = n-2; i >= 0; --i) {
    base = (Base_T) Uinttable_get(high_basetable,keys[i]);

    debug2(printf("At base %u, have %d gregions\n",keys[i],List_length(base->gregions)));
    for (p = base->gregions; p != NULL; p = List_next(p)) {
      gregion = (Gregion_T) List_head(p);
      
      bestscore = 0;
      bestprevpos = keys[n-1];
      debug2(bestprevgregion = NULL);

      prevpos = base->prevpos;
      nlookback = 0;

      floor = gregion->extentstart /*+1*/;
      debug2(printf("  Gregion %u",gregion->extentstart));
      while (prevpos <= keys[n-1] && nlookback < MAXLOOKBACK) {
	prevbase = (Base_T) Uinttable_get(high_basetable,prevpos);
	if ((prevgregion = apply_floor(prevbase->gregions,floor)) != NULL) {
	  debug2(printf(" ... prev:%u, score:%d+%d = %d",
			prevpos,prevgregion->score_floor,gregion->total_down,
			prevgregion->score_floor+gregion->total_down));

	  if ((score = prevgregion->score_floor + gregion->total_down) > bestscore) {
	    debug2(printf("*"));
	    bestscore = score;
	    bestprevpos = prevpos;
#ifdef DEBUG2
	    bestprevgregion = prevgregion;
#endif
	  }
	}

	prevpos = prevbase->prevpos;
	nlookback++;
      }
      debug2(printf("\n"));

      gregion->score_floor = bestscore;
      gregion->bestprev_floor = bestprevpos;

#ifdef DEBUG2
      if (bestprevgregion == NULL) {
	printf(" no prevgregion, bestprev is %u, score is %d\n",
	       bestprevpos,gregion->score_floor);
      } else {
	printf(" bestprev is pos %u, gregion %u, score is %d\n",
	       bestprevpos,bestprevgregion->extentstart,gregion->score_floor);
      }
#endif
    }
  }

  /* Get best overall gregion */
  bestscore = 0;
  bestprevpos = (Chrpos_T) -1U;

  for (i = n-1; i >= 0; --i) {
    base = (Base_T) Uinttable_get(high_basetable,keys[i]);
    for (p = base->gregions; p != NULL; p = List_next(p)) {
      gregion = (Gregion_T) List_head(p);
      
      if (gregion->score_floor > bestscore) {
	bestscore = gregion->score_floor;
	bestprevpos = keys[i];
      }
    }
  }

  if (n < MAX_NCHRS_FOR_ALLOCA) {
    FREEA(keys);
  } else {
    FREE(keys);
  }

  return bestprevpos;
}


static void
traceback_ceilings (Uinttable_T low_basetable, Chrpos_T prevpos) {
  Chrpos_T ceiling;
  Gregion_T end_gregion;
  Base_T base, prevbase;
  Chrpos_T *keys;
  int n, i;
  
  n = Uinttable_length(low_basetable);
  if (n < MAX_NCHRS_FOR_ALLOCA) {
    keys = (Chrpos_T *) MALLOCA(n * sizeof(Chrpos_T));
    Uinttable_fill_keys(keys,low_basetable,/*sortp*/true);
  } else {
    keys = Uinttable_keys(low_basetable,/*sortp*/true);
  }

  ceiling = (Chrpos_T) -1U;

  i = n-1;
  while (prevpos > keys[0]) {
    debug2(printf("traceback from endpos %u, back to %u\n",keys[i],prevpos));
    while (/*startpos*/keys[i] > prevpos) {
      base = (Base_T) Uinttable_get(low_basetable,/*startpos*/keys[i]);
      base->maxextent = ceiling;
      debug2(printf("At low %u, maxextent is %u\n",/*startpos*/keys[i],ceiling));
      i--;
    }

    prevbase = (Base_T) Uinttable_get(low_basetable,prevpos);
    if ((end_gregion = apply_ceiling(prevbase->gregions,ceiling)) == NULL) {
      prevpos = keys[0];	/* Ends loop */
    } else {
      ceiling = end_gregion->extentend /*-1*/;
      prevpos = end_gregion->bestprev_ceiling;
    }
  }

  debug2(printf("End of loop\n"));
  while (i >= 0) {
    base = (Base_T) Uinttable_get(low_basetable,keys[i]);
    base->maxextent = ceiling;
    debug2(printf("At low %u, maxextent is %u\n",keys[i],ceiling));
    i--;
  }

  if (n < MAX_NCHRS_FOR_ALLOCA) {
    FREEA(keys);
  } else {
    FREE(keys);
  }

  return;
}


static void
traceback_floors (Uinttable_T high_basetable, Chrpos_T prevpos) {
  Chrpos_T floor;
  Gregion_T start_gregion;
  Base_T base, prevbase;
  Chrpos_T *keys;
  int n, i;

  n = Uinttable_length(high_basetable);
  if (n < MAX_NCHRS_FOR_ALLOCA) {
    keys = (Chrpos_T *) MALLOCA(n * sizeof(Chrpos_T));
    Uinttable_fill_keys(keys,high_basetable,/*sortp*/true);
  } else {
    keys = Uinttable_keys(high_basetable,/*sortp*/true);
  }

  floor = 0U;

  i = 0;
  while (prevpos < keys[n-1]) {
    debug2(printf("traceback from startpos %u, forward to %u\n",keys[i],prevpos));
    while (/*endpos*/keys[i] < prevpos) {
      base = (Base_T) Uinttable_get(high_basetable,/*endpos*/keys[i]);
      base->minextent = floor;
      debug2(printf("At high %u, minextent is %u\n",/*endpos*/keys[i],floor));
      i++;
    }

    prevbase = (Base_T) Uinttable_get(high_basetable,prevpos);
    if ((start_gregion = apply_floor(prevbase->gregions,floor)) == NULL) {
      prevpos = keys[n-1];	/* Ends loop */
    } else {
      floor = start_gregion->extentstart /*+1*/;
      prevpos = start_gregion->bestprev_floor;
    }
  }

  debug2(printf("End of loop\n"));
  while (i < n) {
    base = (Base_T) Uinttable_get(high_basetable,keys[i]);
    base->minextent = floor;
    debug2(printf("At high %u, minextent is %u\n",/*endpos*/keys[i],floor));
    i++;
  }

  if (n < MAX_NCHRS_FOR_ALLOCA) {
    FREEA(keys);
  } else {
    FREE(keys);
  }

  return;
}


static void
bound_gregions (Uinttable_T low_basetable, Uinttable_T high_basetable) {
  Chrpos_T minextent, maxextent;
  Base_T base_low, base_high;
  Gregion_T gregion;
  List_T p;
  Chrpos_T *keys;
  int n, i;

  n = Uinttable_length(low_basetable);
  if (n < MAX_NCHRS_FOR_ALLOCA) {
    keys = (Chrpos_T *) MALLOCA(n * sizeof(Chrpos_T));
    Uinttable_fill_keys(keys,low_basetable,/*sortp*/true);
  } else {
    keys = Uinttable_keys(low_basetable,/*sortp*/true);
  }

  for (i = 0; i < n; i++) {
    base_low = (Base_T) Uinttable_get(low_basetable,keys[i]);
    maxextent = base_low->maxextent;
    for (p = base_low->gregions; p != NULL; p = List_next(p)) {
      gregion = (Gregion_T) List_head(p);
      base_high = (Base_T) Uinttable_get(high_basetable,gregion->extentend);
      minextent = base_high->minextent;
      if (gregion->extentstart < minextent || gregion->extentend > maxextent) {
	debug2(printf("Not bounded low: #%d:%u..%u\n",gregion->chrnum,gregion->extentstart,gregion->extentend));
	gregion->bounded_low_p = false;
      } else {
	debug2(printf("Bounded low: #%d:%u..%u\n",gregion->chrnum,gregion->extentstart,gregion->extentend));
	gregion->bounded_low_p = true;
      }
    }
  }

  if (n < MAX_NCHRS_FOR_ALLOCA) {
    FREEA(keys);
  } else {
    FREE(keys);
  }


  n = Uinttable_length(high_basetable);
  if (n < MAX_NCHRS_FOR_ALLOCA) {
    keys = (Chrpos_T *) MALLOCA(n * sizeof(Chrpos_T));
    Uinttable_fill_keys(keys,high_basetable,/*sortp*/true);
  } else {
    keys = Uinttable_keys(high_basetable,/*sortp*/true);
  }

  for (i = n-1; i >= 0; i--) {
    base_high = (Base_T) Uinttable_get(high_basetable,keys[i]);
    minextent = base_high->minextent;
    for (p = base_high->gregions; p != NULL; p = List_next(p)) {
      gregion = (Gregion_T) List_head(p);
      base_low = (Base_T) Uinttable_get(low_basetable,gregion->extentstart);
      maxextent = base_low->maxextent;
      if (gregion->extentstart < minextent || gregion->extentend > maxextent) {
	debug2(printf("Not bounded high: #%d:%u..%u\n",gregion->chrnum,gregion->extentstart,gregion->extentend));
	gregion->bounded_high_p = false;
      } else {
	debug2(printf("Bounded high: #%d:%u..%u\n",gregion->chrnum,gregion->extentstart,gregion->extentend));
	gregion->bounded_high_p = true;
      }
    }
  }

  if (n < MAX_NCHRS_FOR_ALLOCA) {
    FREEA(keys);
  } else {
    FREE(keys);
  }

  return;
}




/* Called only if USE_CLEAN is defined in stage1.c */
void
Gregion_filter_clean (List_T gregionlist, int nchrs) {
  Uinttable_T *low_basetables, *high_basetables, basetable;
  Base_T base;
  Chrpos_T prevpos;
  Chrnum_T chrnum;

  List_T p;
  T gregion;
  int n;
#if 0
  List_T unique = NULL;
  T x, y, *array;
  int i, j;
  bool *eliminate;
#endif
#ifdef DEBUG
  List_T q;
  int i;
#endif

  if ((n = List_length(gregionlist)) == 0) {
    return;
  }

  debug(
	for (p = gregionlist, i = 0; p != NULL; p = p->rest, i++) {
	  gregion = (T) p->first;
	  printf("  Initial %d: %d..%d %u-%u (plusp = %d)\n",
		 i,gregion->querystart,gregion->queryend,gregion->extentstart,gregion->extentend,gregion->plusp);
	}
	);

  if (nchrs < MAX_NCHRS_FOR_ALLOCA) {
    low_basetables = (Uinttable_T *) CALLOCA(nchrs+1,sizeof(Uinttable_T));
    high_basetables = (Uinttable_T *) CALLOCA(nchrs+1,sizeof(Uinttable_T));
  } else {
    low_basetables = (Uinttable_T *) CALLOC(nchrs+1,sizeof(Uinttable_T));
    high_basetables = (Uinttable_T *) CALLOC(nchrs+1,sizeof(Uinttable_T));
  }

  for (p = gregionlist; p != NULL; p = List_next(p)) {
    gregion = (T) List_head(p);

    if (low_basetables[gregion->chrnum] == NULL) {
      low_basetables[gregion->chrnum] = Uinttable_new(n);
      high_basetables[gregion->chrnum] = Uinttable_new(n);
    }
    
    basetable = low_basetables[gregion->chrnum];
    if ((base = (Base_T) Uinttable_get(basetable,gregion->extentstart)) == NULL) {
      base = Base_new();
      Uinttable_put(basetable,gregion->extentstart,(void *) base);
    }
    base->gregions = List_push(base->gregions,(void *) gregion);

    basetable = high_basetables[gregion->chrnum];
    if ((base = (Base_T) Uinttable_get(basetable,gregion->extentend)) == NULL) {
      base = Base_new();
      Uinttable_put(basetable,gregion->extentend,(void *) base);
    }
    base->gregions = List_push(base->gregions,(void *) gregion);
  }

  for (chrnum = 0; chrnum <= nchrs; chrnum++) {
    if ((basetable = low_basetables[chrnum]) != NULL) {
      debug2(printf("Processing gregions for chrnum %d\n",chrnum));
      prevpos = compute_ceilings(basetable);
      traceback_ceilings(basetable,prevpos);
  
      basetable = high_basetables[chrnum];
      prevpos = compute_floors(basetable);
      traceback_floors(basetable,prevpos);

      bound_gregions(low_basetables[chrnum],high_basetables[chrnum]);
    }
  }
  
  /* Todo: Free each table */

  if (nchrs < MAX_NCHRS_FOR_ALLOCA) {
    FREEA(high_basetables);
    FREEA(low_basetables);
  } else {
    FREE(high_basetables);
    FREE(low_basetables);
  }

#if 0
  eliminate = (bool *) CALLOC(n,sizeof(bool));

  /* Not necessary if false is zero */
  /*
  for (i = 0; i < n; i++) {
    eliminate[i] = false;
  }
  */

  array = (T *) List_to_array(gregionlist,NULL);
  List_free(&gregionlist);
  qsort(array,n,sizeof(T),weight_cmp);

  for (i = 0; i < n; i++) {
    x = array[i];
    for (j = i+1; j < n; j++) {
      y = array[j];
      if (gregion_overlap_p(x,y) == true) {
	eliminate[j] = true;
      }
    }
  }

  for (i = n-1; i >= 0; i--) {
    gregion = array[i];
    if (eliminate[i] == false) {
      debug(printf("  Keeping %u-%u (plusp = %d)\n",
		   gregion->extentstart,gregion->extentend,gregion->plusp));
      unique = List_push(unique,(void *) gregion);
    } else {
      debug(printf("  Eliminating %u-%u (plusp = %d)\n",
		   gregion->extentstart,gregion->extentend,gregion->plusp));
      /*
      if (gregion->match5 != NULL) {
	Match_decr_npairings(gregion->match5);
	Match_decr_npairings(gregion->match3);
      }
      */
      Gregion_free(&gregion);
    }
  }

  FREE(eliminate);
  FREE(array);

  debug(
	for (p = unique, i = 0; p != NULL; p = p->rest, i++) {
	  gregion = (T) p->first;
	  printf("  Final %d: %u-%u (plusp = %d)\n",
		 i,gregion->extentstart,gregion->extentend,gregion->plusp);
	}
	);

  return unique;

#endif

  return;
}

#endif

