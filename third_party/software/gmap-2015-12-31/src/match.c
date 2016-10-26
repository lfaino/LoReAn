static char rcsid[] = "$Id: match.c 99737 2013-06-27 19:33:03Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#ifndef HAVE_MEMCPY
# define memcpy(d,s,n) bcopy((s),(d),(n))
#endif
#ifndef HAVE_MEMMOVE
# define memmove(d,s,n) bcopy((s),(d),(n))
#endif

#include "match.h"
#include "matchdef.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>		/* For memcpy */
#include "mem.h"
#include "segmentpos.h"
#include "separator.h"
#include "sequence.h"
#include "indexdb.h"
#include "univinterval.h"

#define MIN_STAGE1_FSUPPORT 0.20
#define MAX_STAGE1_STRETCH 2000.0


#define T Match_T

#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif

/* Filtering of poor matchpairs */
#ifdef DEBUG3
#define debug3(x) x
#else
#define debug3(x)
#endif


int
Match_querypos (T this) {
  return this->querypos;
}

bool
Match_forwardp (T this) {
  return this->forwardp;
}

bool
Match_fivep (T this) {
  return this->fivep;
}

Univcoord_T
Match_position (T this) {
  return this->position;
}

Chrnum_T
Match_chrnum (T this) {
  return this->chrnum;
}

/* Used only for debugging.  String is allocated and should be freed. */
char *
Match_chr (T this, Univ_IIT_T chromosome_iit) {
  return Chrnum_to_string(this->chrnum,chromosome_iit);
}

Chrpos_T
Match_chrpos (T this) {
  return this->chrpos;
}

int
Match_incr_npairings (T this) {
  this->npairings += 1;
  return this->npairings;
}

int
Match_decr_npairings (T this) {
  this->npairings -= 1;
  return this->npairings;
}

int
Match_npairings (T this) {
  return this->npairings;
}

void
Match_set_weight (T this, double weight) {
  this->weight = weight;
  this->has_weight_p = true;
  return;
}

double
Match_weight (T this) {
  return this->weight;
}

bool
Match_has_weight_p (T this) {
  return this->has_weight_p;
}


int
Match_cmp (const void *a, const void *b) {
  T x = * (T *) a;
  T y = * (T *) b;

  if (x->position < y->position) {
    return -1;
  } else if (x->position > y->position) {
    return 1;
  } else {
    return 0;
  }
}

#ifndef USE_MATCHPOOL
/* Matches now made in matchpool.c */
T
Match_new (int querypos, bool forwardp, bool fivep,
	   Univcoord_T position, Univ_IIT_T chromosome_iit) {
  T new = (T) MALLOC(sizeof(*new));
  int index;

  new->querypos = querypos;
  new->weight = 0.0;		/* Will be entered later */
  new->has_weight_p = false;
  new->position = position;
  new->forwardp = forwardp;
  new->fivep = fivep;
  new->npairings = 0;

  if (chromosome_iit == NULL) {
    new->chrnum = 0;
    new->chrpos = position;
  } else {
    index = Univ_IIT_get_one(chromosome_iit,position,position);
    new->chrpos = position - Univinterval_low(Univ_IIT_interval(chromosome_iit,index));
    new->chrnum = index;
  }

  return new;
}

void
Match_free (T *old) {
  if (*old) {
    FREE(*old);
  }
  return;
}
#endif


/* Static gbuffer1, gbuffer2 are allowable only for debugging */
#ifdef PMAP
#define MAXSTAGE1SIZE 36
#else
#define MAXSTAGE1SIZE 24
#endif

void
Match_print_mer (T this, char *queryseq_ptr, Genome_T genome, Univ_IIT_T chromosome_iit, int stage1size) {
  char *genomicseg_ptr;
  Sequence_T genomicseg;
  int querypos;
  Univcoord_T position;

  querypos = this->querypos;
  position = this->position;

#ifdef PMAP
  if (this->forwardp == true) {
    genomicseg = Genome_get_segment(genome,position,3*stage1size,chromosome_iit,/*revcomp*/false);
  } else {
    genomicseg = Genome_get_segment(genome,position-(3*stage1size-1U),3*stage1size,chromosome_iit,/*revcomp*/true);
  }
#else
  if (this->forwardp == true) {
    genomicseg = Genome_get_segment(genome,position,stage1size,chromosome_iit,/*revcomp*/false);
  } else {
    genomicseg = Genome_get_segment(genome,position-(stage1size-1U),stage1size,chromosome_iit,/*revcomp*/true);
  }
#endif
  genomicseg_ptr = Sequence_fullpointer(genomicseg);

  printf("query:%.*s ",stage1size,&(queryseq_ptr[querypos]));
#ifdef PMAP
  printf("genomic:%.*s",3*stage1size,genomicseg_ptr);
#else
  printf("genomic:%.*s",stage1size,genomicseg_ptr);
#endif

  Sequence_free(&genomicseg);

  return;
}


void
Match_print (T this, Univ_IIT_T chromosome_iit) {
  char *chr;

  chr = Match_chr(this,chromosome_iit);
  printf("  Match at %d: #%d(chr%s):%u (forwardp:%d, npairings:%d), weight %.3f ",
	 this->querypos,this->chrnum,chr,this->chrpos,this->forwardp,
	 this->npairings,this->weight);
  FREE(chr);

  return;
}



static double
compute_fsupport (T bound5, T bound3, int trimlength, int stage1size) {
  return (double) (bound3->querypos - bound5->querypos + stage1size)/(double) trimlength;
}

static double
compute_stretch (T bound5, T bound3) {
  Univcoord_T position5, position3;
  int querypos5, querypos3;

  querypos5 = bound5->querypos;
  querypos3 = bound3->querypos;
  
  if (querypos5 == querypos3) {
    return 1.0;
  } else {
    position5 = bound5->position;
    position3 = bound3->position;
    if (position3 > position5) {
#ifdef PMAP
      return (double) (position3 - position5)/(double) (querypos3 - querypos5)/3.0;
#else
      return (double) (position3 - position5)/(double) (querypos3 - querypos5);
#endif
    } else {
#ifdef PMAP
      return (double) (position5 - position3)/(double) (querypos3 - querypos5)/3.0;
#else
      return (double) (position5 - position3)/(double) (querypos3 - querypos5);
#endif
    }
  }
}

bool
Match_acceptable_pair (T match5, T match3, int trimlength, int stage1size) {
  debug3(printf("support = %d/trimlength = %d\n",
		match3->querypos - match5->querypos + stage1size,trimlength));
  debug3(printf("stretch = %f\n",compute_stretch(match5,match3)));

  if (compute_fsupport(match5,match3,trimlength,stage1size) < MIN_STAGE1_FSUPPORT) {
    debug3(printf("Insufficient coverage of the query sequence\n"));
    return false;
  } else if (compute_stretch(match5,match3) > MAX_STAGE1_STRETCH) {
    debug3(printf("Genomic region is too large relative to matching cDNA region\n"));
    return false;
  } else {
    return true;
  }
}


#if 0

bool
Match_sufficient_support (T match5, T match3, int trimstart, int trimend) {
#ifdef PMAP
  debug(printf("  Testing bound5 = %d < %d + %d, bound3 = %d (+%d) > %d - %d\n",
	       Match_querypos(match5),trimstart,SUFFICIENT_SUPPORT,
	       Match_querypos(match3),INDEX1PART_AA,trimend,SUFFICIENT_SUPPORT));
  if (match5->querypos < trimstart + SUFFICIENT_SUPPORT && 
      match3->querypos + INDEX1PART_AA > trimend - SUFFICIENT_SUPPORT) {
    return true;
  } else {
    return false;
  }
#else
  debug(printf("  Testing bound5 = %d < %d + %d, bound3 = %d (+%d) > %d - %d\n",
	       Match_querypos(match5),trimstart,SUFFICIENT_SUPPORT,
	       Match_querypos(match3),INDEX1PART,trimend,SUFFICIENT_SUPPORT));
  if (match5->querypos < trimstart + SUFFICIENT_SUPPORT && 
      match3->querypos + INDEX1PART > trimend - SUFFICIENT_SUPPORT) {
    return true;
  } else {
    return false;
  }
#endif
}

#endif
