static char rcsid[] = "$Id: splicetrie_build.c 135656 2014-05-09 01:47:27Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "splicetrie_build.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>		/* For qsort */

#include "assert.h"
#include "mem.h"
#include "iitdef.h"
#include "interval.h"


#define MININTRONLEN 9
#define SPLICEDIST_EXTRA 10 /* Reported intron length does not include dinucleotides */

#define EMPTY_POINTER 0U
#define INTRON_HIGH_TO_LOW 1U /* Because splicesites based on Interval_low (x..x+1 and y..y+1), but introns go from low to high (x..y+1) */


/* Generating trie */
#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif

/* Generating trie, details */
#ifdef DEBUG0
#define debug0(x) x
#else
#define debug0(x)
#endif

/* Full dump */
#ifdef DEBUG2
#define debug2(x) x
#else
#define debug2(x)
#endif


char *
Splicetype_string (Splicetype_T splicetype) {
  switch (splicetype) {
  case DONOR: return "donor";
  case ACCEPTOR: return "acceptor";
  case ANTIDONOR: return "antidonor";
  case ANTIACCEPTOR: return "antiacceptor";
  default: abort();
  }
}


#define RIGHT_A 0x00
#define RIGHT_C 0x01
#define RIGHT_G 0x02
#define RIGHT_T 0x03


/*                      87654321 */
#define LOW_TWO_BITS  0x00000003

static void
splicefrag_nt_leftward (char *nt, Genomecomp_T splicefrag) {
  int i, j;
  Genomecomp_T highbits, lowbits;

  highbits = splicefrag >> 16;
  lowbits = splicefrag /* & 0x0000FFFF */;

  j = 15;
  for (i = 0; i < 16; i++) {
    switch (((highbits & 0x01) << 1) | (lowbits & 0x01)) {
    case RIGHT_A: nt[j] = 'A'; break;
    case RIGHT_C: nt[j] = 'C'; break;
    case RIGHT_G: nt[j] = 'G'; break;
    case RIGHT_T: nt[j] = 'T'; break;
    }
    highbits >>= 1;
    lowbits >>= 1;
    j--;
  }

  return;
}

static void
splicefrag_nt_rightward (char *nt, Genomecomp_T splicefrag) {
  int i, j;
  Genomecomp_T highbits, lowbits;

  highbits = splicefrag >> 16;
  lowbits = splicefrag /* & 0x0000FFFF */;

  j = 0;
  for (i = 0; i < 16; i++) {
    switch (((highbits & 0x01) << 1) | (lowbits & 0x01)) {
    case RIGHT_A: nt[j] = 'A'; break;
    case RIGHT_C: nt[j] = 'C'; break;
    case RIGHT_G: nt[j] = 'G'; break;
    case RIGHT_T: nt[j] = 'T'; break;
    }
    highbits >>= 1;
    lowbits >>= 1;
    j++;
  }

  return;
}



/* Splicetype is for the anchor splice */
void
Splicetrie_dump (Trieoffset_T *triestart, Univcoord_T *splicesites,
		 Splicetype_T splicetype, Genomecomp_T *splicefrags_ref) {
  Triecontent_T leaf;
  int nleaves, i;
  Univcoord_T position;
  Triecontent_T offseta, offsetc, offsetg, offsett;
  char gbuffer[17];

  gbuffer[16] = '\0';

  if (single_leaf_p(leaf = triestart[0])) {
    position = splicesites[leaf];
    printf("%d %u",(int) leaf,position);
    if (splicetype == DONOR || splicetype == ANTIACCEPTOR) {
      splicefrag_nt_rightward(gbuffer,splicefrags_ref[leaf]);
      printf(" %s (rightward)\n",gbuffer);
    } else {
      splicefrag_nt_leftward(gbuffer,splicefrags_ref[leaf]);
      printf(" %s (leftward)\n",gbuffer);
    }

  } else if (multiple_leaf_p(leaf)) {
    nleaves = (int) (-leaf);
    for (i = 1; i <= nleaves; i++) {
      leaf = triestart[i];
      position = splicesites[leaf];
      printf("%d %u",(int) leaf,position);
      if (splicetype == DONOR || splicetype == ANTIACCEPTOR) {
	splicefrag_nt_rightward(gbuffer,splicefrags_ref[leaf]);
	printf(" %s (rightward)\n",gbuffer);
      } else {
	splicefrag_nt_leftward(gbuffer,splicefrags_ref[leaf]);
	printf(" %s (leftward)\n",gbuffer);
      }
    }

  } else {
#ifdef USE_2BYTE_RELOFFSETS
    get_offsets(&offseta,&offsetc,&offsetg,&offsett,triestart[1],triestart[2]);
#else
    offseta = (int) triestart[1];
    offsetc = (int) triestart[2];
    offsetg = (int) triestart[3];
    offsett = (int) triestart[4];
#endif

    if (offseta > 0) {
      Splicetrie_dump(&(triestart[-offseta]),splicesites,splicetype,splicefrags_ref);
    }
    if (offsetc > 0) {
      Splicetrie_dump(&(triestart[-offsetc]),splicesites,splicetype,splicefrags_ref);
    }
    if (offsetg > 0) {
      Splicetrie_dump(&(triestart[-offsetg]),splicesites,splicetype,splicefrags_ref);
    }
    if (offsett > 0) {
      Splicetrie_dump(&(triestart[-offsett]),splicesites,splicetype,splicefrags_ref);
    }
  }

  return;
}


/*               87654321 */
#define LEFT_A 0x00000000
#define LEFT_C 0x40000000
#define LEFT_G 0x80000000
#define LEFT_T 0xC0000000
#define HIGH2  0xC0000000


/* Puts leftmost character into lowest bits */
/* For right splicestrings, we want the leftmost character in the highest bits */

static Genomecomp_T
compress16 (bool *saw_n_p, char *buffer) {
  Genomecomp_T low = 0U;
  int c;
  int i;

  /* *saw_n_p = false; -- Want to check both ref and alt, so rely on caller to set */
  for (i = 0; i < 16; i++) {
    c = buffer[i];
    low >>= 2;
    switch (c) {
    case 'A': break;
    case 'C': low |= LEFT_C; break;
    case 'G': low |= LEFT_G; break;
    case 'T': low |= LEFT_T; break;
    default: *saw_n_p = true; break;
    }
  }

  return low;
}


static Genomecomp_T
uint4_reverse (Genomecomp_T forward) {
  Genomecomp_T reverse = 0U;
  int c;
  int i;

  for (i = 0; i < 16; i++) {
    c = forward & 0x03;
    reverse <<= 2;
    reverse |= c;
    forward >>= 2;
  }

  return reverse;
}



/*                   87654321 */
#define LEFT_SET   0x80000000
#define LEFT_CLEAR 0x00000000

/* Puts the odd bits in the top 16 bits, and the even bits in the
   bottom 16 bits.  Needed to match genomebits format, although that
   has 32-bit blocks instead of 16 bits. */
static Genomecomp_T
unshuffle16 (bool *saw_n_p, char *buffer) {
  Genomecomp_T bits = 0U;
  int c;
  int i;

  /* *saw_n_p = false; -- Want to check both ref and alt, so rely on caller to set */
  /* Low bits */
  for (i = 0; i < 16; i++) {
    c = buffer[i];
    bits >>= 1;
    switch (c) {
    case 'A': /* bits |= LEFT_CLEAR; */ break;
    case 'C':    bits |= LEFT_SET; break;
    case 'G': /* bits |= LEFT_CLEAR; */ break;
    case 'T':    bits |= LEFT_SET; break;
    default: *saw_n_p = true; break;
    }
  }

  /* High bits */
  for (i = 0; i < 16; i++) {
    c = buffer[i];
    bits >>= 1;
    switch (c) {
    case 'A': /* bits |= LEFT_CLEAR; */ break;
    case 'C': /* bits |= LEFT_CLEAR; */ break;
    case 'G':    bits |= LEFT_SET; break;
    case 'T':    bits |= LEFT_SET; break;
    default: *saw_n_p = true; break;
    }
  }

  return bits;
}


#if 0
/* Replaced by procedures in splicestringpool.c */
static void
Splicestring_free (Splicestring_T *old) {
  FREE(*old);
  return;
}

void
Splicestring_gc (List_T *splicestrings, int nsplicesites) {
  int i;
  List_T list, p;
  Splicestring_T splicestring;

  for (i = 0; i < nsplicesites; i++) {
    list = splicestrings[i];
    for (p = list; p != NULL; p = List_next(p)) {
      splicestring = (Splicestring_T) List_head(p);
      Splicestring_free(&splicestring);
    }
    List_free(&list);
  }
  FREE(splicestrings);
  return;
}

static Splicestring_T
Splicestring_new (Genomecomp_T string, Genomecomp_T splicesite, int splicesite_i) {
  Splicestring_T new = (Splicestring_T) MALLOC(sizeof(*new));

  new->string = string;
  new->splicesite = splicesite;
  new->splicesite_i = (Genomecomp_T) splicesite_i;
  return new;
}
#endif

static int
Splicestring_cmp (const void *a, const void *b) {
  Splicestring_T x = * (Splicestring_T *) a;
  Splicestring_T y = * (Splicestring_T *) b;

#ifdef USE_STRINGS
  return strcmp(x->string,y->string);
#else
  if (x->string < y->string) {
    return -1;
  } else if (x->string > y->string) {
    return +1;
  } else {
    return 0;
  }
#endif
}


static List_T
allelic_combinations (Genomecomp_T refstring, Genomecomp_T altstring, Genomecomp_T splicesite, int splicesite_i,
		      Splicestringpool_T splicestringpool) {
  List_T splicestrings = NULL;
  Uintlist_T combinations, newcombinations, temp, p;
  int refc, altc;
  int i;

  combinations = Uintlist_push(NULL,0U);
  for (i = 0; i < 16; i++) {
    newcombinations = NULL;

    refc = (refstring & HIGH2) >> 30;
    altc = (altstring & HIGH2) >> 30;
    refstring <<= 2;
    altstring <<= 2;
    if (refc != altc) {
      for (p = combinations; p != NULL; p = Uintlist_next(p)) {
	newcombinations = Uintlist_push(newcombinations,(Uintlist_head(p) << 2) | refc);
	newcombinations = Uintlist_push(newcombinations,(Uintlist_head(p) << 2) | altc);
      }
    } else {
      for (p = combinations; p != NULL; p = Uintlist_next(p)) {
	newcombinations = Uintlist_push(newcombinations,(Uintlist_head(p) << 2) | refc);
      }
    }

    temp = combinations;
    combinations = Uintlist_reverse(newcombinations);
    Uintlist_free(&temp);
  }

  for (p = combinations; p != NULL; p = Uintlist_next(p)) {
    splicestrings = Splicestringpool_push(splicestrings,splicestringpool,
					  Uintlist_head(p),splicesite,splicesite_i);
  }

  Uintlist_free(&combinations);
  return splicestrings;
}






/* If distances are provided, splicedists are the observed distances */
Univcoord_T *
Splicetrie_retrieve_via_splicesites (bool *distances_observed_p,
#ifdef GSNAP
				     Genomecomp_T **splicecomp,
#endif
				     Splicetype_T **splicetypes, Chrpos_T **splicedists,
				     List_T **splicestrings, Genomecomp_T **splicefrags_ref, Genomecomp_T **splicefrags_alt,
				     int *nsplicesites, IIT_T splicing_iit, int *splicing_divint_crosstable,
				     int donor_typeint, int acceptor_typeint, Univ_IIT_T chromosome_iit,
				     Genome_T genome, Genome_T genomealt, Chrpos_T shortsplicedist,
				     Splicestringpool_T splicestringpool) {
  Univcoord_T *splicesites, chroffset, chrhigh, position;
  Chrpos_T chrlength, chrpos;
  Univcoord_T last_donor, last_antidonor, last_acceptor, last_antiacceptor;
  int last_donor_k, last_antidonor_k, last_acceptor_k, last_antiacceptor_k;
  int distance;
  Genomecomp_T refstring, altstring;
  int *splicesites1;
  int divno, nsplicesites1, i, k;
  Chrnum_T chrnum;
  Interval_T *intervals, interval;
  struct Interval_T *interval_structs;
  char gbuffer_ref[17], gbuffer_alt[17], *chr;
  char *restofheader, *annot;
  bool firstp = true, saw_n_p, allocp, alloc_header_p;
  int ntoolong = 0;

#ifdef DEBUG2
  List_T p;
  Splicestring_T splicestring;
#endif

#ifdef GSNAP
  int nblocks;
#endif

  k = 0;
  for (chrnum = 1; chrnum <= Univ_IIT_total_nintervals(chromosome_iit); chrnum++) {
    if ((divno = splicing_divint_crosstable[chrnum]) > 0) {
      Univ_IIT_interval_bounds(&chroffset,&chrhigh,&chrlength,chromosome_iit,chrnum,/*circular_typeint*/-1);
      /* chrlength = Univ_IIT_length(chromosome_iit,chrnum); */
      splicesites1 = IIT_get_with_divno(&nsplicesites1,splicing_iit,divno,
					0U,chrlength-1U,/*sortp*/false);
      if (nsplicesites1 > 0) {
	if (firstp == true) {
	  annot = IIT_annotation(&restofheader,splicing_iit,splicesites1[0],&alloc_header_p);
	  if (restofheader[0] == '\0') {
	    fprintf(stderr,"splice distances absent...");
	    *distances_observed_p = false;
	  } else if (sscanf(restofheader,"%d",&distance) < 1) {
	    fprintf(stderr,"splice distances absent...");
	    *distances_observed_p = false;
	  } else {
	    fprintf(stderr,"splice distances present...");
	    *distances_observed_p = true;
	  }
	  if (alloc_header_p == true) {
	    FREE(restofheader);
	  }
	  firstp = false;
	}

	intervals = (Interval_T *) CALLOC(nsplicesites1,sizeof(Interval_T));
	for (i = 0; i < nsplicesites1; i++) {
	  intervals[i] = &(splicing_iit->intervals[divno][i]);
	}
	qsort(intervals,nsplicesites1,sizeof(Interval_T),Interval_cmp_low);

	last_donor = last_antidonor = last_acceptor = last_antiacceptor = 0U;
	for (i = 0; i < nsplicesites1; i++) {
	  interval = intervals[i];
	  chrpos = Interval_low(interval);
	  position = chrpos + chroffset;

	  if (position >= chrhigh) {
	    chr = Univ_IIT_label(chromosome_iit,chrnum,&allocp);
	    fprintf(stderr,"\nSplice site %s:%u extends beyond chromosome length %u.  Discarding...",
		    chr,chrpos,chrlength);
	    if (allocp) FREE(chr);

	  } else if (Interval_type(interval) == donor_typeint) {
	    if (Interval_sign(interval) > 0) {
	      if (position != last_donor) {
		last_donor = position;
		k++;
	      }
	    } else {
	      if (position != last_antidonor) {
		last_antidonor = position;
		k++;
	      }
	    }
	  } else if (Interval_type(interval) == acceptor_typeint) {
	    if (Interval_sign(interval) > 0) {
	      if (position != last_acceptor) {
		last_acceptor = position;
		k++;
	      }
	    } else {
	      if (position != last_antiacceptor) {
		last_antiacceptor = position;
		k++;
	      }
	    }
	  }
	}
	FREE(intervals);
	FREE(splicesites1);
      }
    }
  }

  *nsplicesites = k;
  debug(printf("total unique splicesites: %d\n",*nsplicesites));
  fprintf(stderr,"%d unique splicesites...",*nsplicesites);

  /* The above procedure determines number of unique splicesites */

  if (*nsplicesites == 0) {
#ifdef GSNAP
    *splicecomp = (Genomecomp_T *) NULL;
#endif
    splicesites = (Univcoord_T *) NULL;
    *splicetypes = (Splicetype_T *) NULL;
    *splicedists = (Chrpos_T *) NULL;
    *splicestrings = (List_T *) NULL;
    *splicefrags_ref = (Genomecomp_T *) NULL;
    *splicefrags_alt = (Genomecomp_T *) NULL;
    return splicesites;
  }

#ifdef GSNAP
  nblocks = (Genome_totallength(genome)+31)/32U;
  *splicecomp = (Genomecomp_T *) CALLOC(nblocks,sizeof(Genomecomp_T));
#endif
  splicesites = (Univcoord_T *) CALLOC((*nsplicesites) + 1,sizeof(Univcoord_T));
  *splicetypes = (Splicetype_T *) CALLOC(*nsplicesites,sizeof(Splicetype_T));
  *splicedists = (Chrpos_T *) CALLOC(*nsplicesites,sizeof(Chrpos_T));
  *splicestrings = (List_T *) CALLOC(*nsplicesites,sizeof(List_T));
  *splicefrags_ref = (Genomecomp_T *) CALLOC(*nsplicesites,sizeof(Genomecomp_T));
  if (genomealt == NULL) {
    *splicefrags_alt = *splicefrags_ref;
  } else {
    *splicefrags_alt = (Genomecomp_T *) CALLOC(*nsplicesites,sizeof(Genomecomp_T));
  }

  /* Use interval_structs instead of intervals, because we want to
     copy information and avoid creating many Interval_T objects */

  k = 0;
  for (chrnum = 1; chrnum <= Univ_IIT_total_nintervals(chromosome_iit); chrnum++) {
    if ((divno = splicing_divint_crosstable[chrnum]) > 0) {
      Univ_IIT_interval_bounds(&chroffset,&chrhigh,&chrlength,chromosome_iit,chrnum,/*circular_typeint*/-1);
      /* chrlength = Univ_IIT_length(chromosome_iit,chrnum); */
      splicesites1 = IIT_get_with_divno(&nsplicesites1,splicing_iit,divno,
					0U,chrlength-1U,/*sortp*/false);
      if (nsplicesites1 > 0) {
	interval_structs = (struct Interval_T *) CALLOC(nsplicesites1,sizeof(struct Interval_T));
	for (i = 0; i < nsplicesites1; i++) {
	  /* intervals[i] = &(splicing_iit->intervals[divno][i]); */
	  /* Copy so we can store distance information in Interval_high */
	  Interval_copy_existing(&(interval_structs[i]),&(splicing_iit->intervals[divno][i]));
	  if (*distances_observed_p == false) {
	    /* No, want to have essentially zero distance */
	    /* Interval_store_length(intervals[i],shortsplicedist); */
	  } else {
	    annot = IIT_annotation(&restofheader,splicing_iit,splicesites1[i],&alloc_header_p);
	    if (sscanf(restofheader,"%d",&distance) != 1) {
	      fprintf(stderr,"splicesites file missing distance in entry %s...exiting\n",
		      IIT_label(splicing_iit,splicesites1[i],&allocp));
	      exit(9);
	    } else if (distance < 0) {
	      fprintf(stderr,"splicesites file has a negative distance %d in entry %s...exiting\n",
		      distance,IIT_label(splicing_iit,splicesites1[i],&allocp));
	      exit(9);
	    } else if (distance > (int) shortsplicedist) {
	      ntoolong++;
	      Interval_store_length(&(interval_structs[i]),distance + SPLICEDIST_EXTRA);  /* Previously stored shortsplicedist */
	    } else {
	      Interval_store_length(&(interval_structs[i]),distance + SPLICEDIST_EXTRA);
	    }
	    if (alloc_header_p == true) {
	      FREE(restofheader);
	    }
	  }
	}

	qsort(interval_structs,nsplicesites1,sizeof(struct Interval_T),Interval_cmp_low_struct);

	last_donor = last_antidonor = last_acceptor = last_antiacceptor = 0U;
	for (i = 0; i < nsplicesites1; i++) {
	  interval = &(interval_structs[i]);
	  chrpos = Interval_low(interval);
	  position = chrpos + chroffset;

	  if (position >= chrhigh) {
#if 0
	    /* Warning given previously */
	    chr = Univ_IIT_label(chromosome_iit,chrnum,&allocp);
	    fprintf(stderr,"\nSplice site %s:%u extends beyond chromosome length %u.  Discarding...",
		    chr,chrpos,chrlength);
	    if (allocp) FREE(chr);
#endif

	  } else if (Interval_type(interval) == donor_typeint) {
	    if (Interval_sign(interval) > 0) {
	      if (position == last_donor) {
		if (Interval_length(interval) > (*splicedists)[last_donor_k]) {
		  (*splicedists)[last_donor_k] = Interval_length(interval);
		}

	      } else {
		last_donor_k = k;
		last_donor = splicesites[k] = position;
#ifdef GSNAP
		assert(position/32U < (unsigned int) nblocks);
		(*splicecomp)[position/32U] |= (1 << (position % 32));
#endif
		(*splicetypes)[k] = DONOR;
		(*splicedists)[k] = Interval_length(interval);

		saw_n_p = false;
		Genome_fill_buffer_simple(genome,position-16,16,gbuffer_ref);
		refstring = compress16(&saw_n_p,gbuffer_ref);
		(*splicefrags_ref)[k] = unshuffle16(&saw_n_p,gbuffer_ref);
		if (genomealt) {
		  Genome_fill_buffer_simple_alt(genome,genomealt,position-16,16,gbuffer_alt);
		  altstring = compress16(&saw_n_p,gbuffer_alt);
		  (*splicefrags_alt)[k] = unshuffle16(&saw_n_p,gbuffer_alt);
		}

		if (saw_n_p == true) {
		  chr = Univ_IIT_label(chromosome_iit,chrnum,&allocp);
		  fprintf(stderr,"\nNon-standard nucleotide N near splice site %s:%u.  Discarding...",
			  chr,chrpos);
		  if (allocp) FREE(chr);

		} else if (genomealt) {
		  (*splicestrings)[k] = allelic_combinations(refstring,altstring,position,k,splicestringpool);
		  k++;
		} else {
		  (*splicestrings)[k] = Splicestringpool_push(NULL,splicestringpool,refstring,position,k);
		  k++;
		}
	      }

	    } else {
	      if (position == last_antidonor) {
		if (Interval_length(interval) > (*splicedists)[last_antidonor_k]) {
		  (*splicedists)[last_antidonor_k] = Interval_length(interval);
		}

	      } else {
		last_antidonor_k = k;
		last_antidonor = splicesites[k] = position;
#ifdef GSNAP
		assert(position/32U < (unsigned int) nblocks);
		(*splicecomp)[position/32U] |= (1 << (position % 32));
#endif
		(*splicetypes)[k] = ANTIDONOR;
		(*splicedists)[k] = Interval_length(interval);

		saw_n_p = false;
		Genome_fill_buffer_simple(genome,position,16,gbuffer_ref);
		refstring = uint4_reverse(compress16(&saw_n_p,gbuffer_ref));
		(*splicefrags_ref)[k] = unshuffle16(&saw_n_p,gbuffer_ref);
		if (genomealt) {
		  Genome_fill_buffer_simple_alt(genome,genomealt,position,16,gbuffer_alt);
		  altstring = uint4_reverse(compress16(&saw_n_p,gbuffer_alt));
		  (*splicefrags_alt)[k] = unshuffle16(&saw_n_p,gbuffer_alt);
		}

		if (saw_n_p == true) {
		  chr = Univ_IIT_label(chromosome_iit,chrnum,&allocp);
		  fprintf(stderr,"\nNon-standard nucleotide N near splice site %s:%u.  Discarding...",
			  chr,chrpos);
		  if (allocp) FREE(chr);
		} else if (genomealt) {
		  (*splicestrings)[k] = allelic_combinations(refstring,altstring,position,k,splicestringpool);
		  k++;
		} else {
		  (*splicestrings)[k] = Splicestringpool_push(NULL,splicestringpool,refstring,position,k);
		  k++;
		}
	      }
	    }

	  } else if (Interval_type(interval) == acceptor_typeint) {
	    if (Interval_sign(interval) > 0) {
	      if (position == last_acceptor) {
		if (Interval_length(interval) > (*splicedists)[last_acceptor_k]) {
		  (*splicedists)[last_acceptor_k] = Interval_length(interval);
		}

	      } else {
		last_acceptor_k = k;
		last_acceptor = splicesites[k] = position;
#ifdef GSNAP
		assert(position/32U < (unsigned int) nblocks);
		(*splicecomp)[position/32U] |= (1 << (position % 32));
#endif
		(*splicetypes)[k] = ACCEPTOR;
		(*splicedists)[k] = Interval_length(interval);

		saw_n_p = false;
		Genome_fill_buffer_simple(genome,position,16,gbuffer_ref);
		refstring = uint4_reverse(compress16(&saw_n_p,gbuffer_ref));
		(*splicefrags_ref)[k] = unshuffle16(&saw_n_p,gbuffer_ref);
		if (genomealt) {
		  Genome_fill_buffer_simple_alt(genome,genomealt,position,16,gbuffer_alt);
		  altstring = uint4_reverse(compress16(&saw_n_p,gbuffer_alt));
		  (*splicefrags_alt)[k] = unshuffle16(&saw_n_p,gbuffer_alt);
		}

		if (saw_n_p == true) {
		  chr = Univ_IIT_label(chromosome_iit,chrnum,&allocp);
		  fprintf(stderr,"\nNon-standard nucleotide N near splice site %s:%u.  Discarding...",
			  chr,chrpos);
		  if (allocp) FREE(chr);
		} else if (genomealt) {
		  (*splicestrings)[k] = allelic_combinations(refstring,altstring,position,k,splicestringpool);
		  k++;
		} else {
		  (*splicestrings)[k] = Splicestringpool_push(NULL,splicestringpool,refstring,position,k);
		  k++;
		}
	      }

	    } else {
	      if (position == last_antiacceptor) {
		if (Interval_length(interval) > (*splicedists)[last_antiacceptor_k]) {
		  (*splicedists)[last_antiacceptor_k] = Interval_length(interval);
		}

	      } else {
		last_antiacceptor_k = k;
		last_antiacceptor = splicesites[k] = position;
#ifdef GSNAP
		assert(position/32U < (unsigned int) nblocks);
		(*splicecomp)[position/32U] |= (1 << (position % 32));
#endif
		(*splicetypes)[k] = ANTIACCEPTOR;
		(*splicedists)[k] = Interval_length(interval);

		saw_n_p = false;
		Genome_fill_buffer_simple(genome,position-16,16,gbuffer_ref);
		refstring = compress16(&saw_n_p,gbuffer_ref);
		(*splicefrags_ref)[k] = unshuffle16(&saw_n_p,gbuffer_ref);
		if (genomealt) {
		  Genome_fill_buffer_simple_alt(genome,genomealt,position-16,16,gbuffer_alt);
		  altstring = compress16(&saw_n_p,gbuffer_alt);
		  (*splicefrags_alt)[k] = unshuffle16(&saw_n_p,gbuffer_alt);
		}

		if (saw_n_p == true) {
		  chr = Univ_IIT_label(chromosome_iit,chrnum,&allocp);
		  fprintf(stderr,"\nNon-standard nucleotide N near splice site %s:%u.  Discarding...",
			  chr,chrpos);
		  if (allocp) FREE(chr);
		} else if (genomealt) {
		  (*splicestrings)[k] = allelic_combinations(refstring,altstring,position,k,splicestringpool);
		  k++;
		} else {
		  (*splicestrings)[k] = Splicestringpool_push(NULL,splicestringpool,refstring,position,k);
		  k++;
		}
	      }

	    }
	  }
	}

	FREE(interval_structs);
	FREE(splicesites1);
      }
    }
  }

  *nsplicesites = k;
  splicesites[*nsplicesites] = (Univcoord_T) -1; /* Marker for comparison in identify_all_segments */
  fprintf(stderr,"\n%d splicesites are valid...",*nsplicesites);

#ifdef DEBUG2
  for (k = 0; k < *nsplicesites; k++) {
    printf("%d: %u %s %08X\n",k,splicesites[k],Splicetype_string((*splicetypes)[k]),(*splicefrags_ref)[k]);
    for (p = (*splicestrings)[k]; p != NULL; p = List_next(p)) {
      splicestring = (Splicestring_T) List_head(p);
      printf("  %u %u %u\n",splicestring->string,splicestring->splicesite,splicestring->splicesite_i);
    }
  }
#endif

  if (ntoolong > 0) {
    fprintf(stderr,"%d entries with distance > %d specified for local splice distance...",ntoolong,shortsplicedist);
  }

  return splicesites;
}


struct Cell_T {
  int k;			/* original k */
  Univcoord_T pos;		/* splicesite pos */
};

static int
Cell_position_cmp (const void *a, const void *b) {
  struct Cell_T x = * (struct Cell_T *) a;
  struct Cell_T y = * (struct Cell_T *) b;

  if (x.pos < y.pos) {
    return -1;
  } else if (y.pos < x.pos) {
    return +1;
  } else {
    return 0;
  }
}


Univcoord_T *
Splicetrie_retrieve_via_introns (
#ifdef GSNAP
				 Genomecomp_T **splicecomp,
#endif
				 Splicetype_T **splicetypes, Chrpos_T **splicedists,
				 List_T **splicestrings, Genomecomp_T **splicefrags_ref, Genomecomp_T **splicefrags_alt,
				 int *nsplicesites, IIT_T splicing_iit, int *splicing_divint_crosstable,
				 Univ_IIT_T chromosome_iit, Genome_T genome, Genome_T genomealt,
				 Splicestringpool_T splicestringpool) {
  Univcoord_T *splicesites, chroffset, chrhigh, position;
  Chrpos_T chrlength, chrpos;
  Univcoord_T last_donor, last_antidonor, last_acceptor, last_antiacceptor;
  int last_donor_k, last_antidonor_k, last_acceptor_k, last_antiacceptor_k;
  Genomecomp_T refstring, altstring;
  int *introns1;
  int divno, nintrons1, i, k, j;
  Chrnum_T chrnum;
  Interval_T *intervals, interval;
  char gbuffer_ref[17], gbuffer_alt[17], *chr;
  bool saw_n_p, allocp;

  struct Cell_T *cells;
  Univcoord_T *temp_splicesites;
  Splicetype_T *temp_splicetypes;
  Chrpos_T *temp_splicedists;
  List_T *temp_splicestrings, p;
  UINT4 *temp_splicefrags_ref, *temp_splicefrags_alt;
  Splicestring_T splicestring;

#ifdef GSNAP
  int nblocks;
#endif

  k = 0;
  for (chrnum = 1; chrnum <= Univ_IIT_total_nintervals(chromosome_iit); chrnum++) {
    if ((divno = splicing_divint_crosstable[chrnum]) > 0) {
      Univ_IIT_interval_bounds(&chroffset,&chrhigh,&chrlength,chromosome_iit,chrnum,/*circular_typeint*/-1);
      /* chrlength = Univ_IIT_length(chromosome_iit,chrnum); */
      introns1 = IIT_get_with_divno(&nintrons1,splicing_iit,divno,0U,chrlength-1U,/*sortp*/false);
      if (nintrons1 > 0) {
	intervals = (Interval_T *) CALLOC(nintrons1,sizeof(Interval_T));
	for (i = 0; i < nintrons1; i++) {
	  intervals[i] = &(splicing_iit->intervals[divno][i]);
	}

	qsort(intervals,nintrons1,sizeof(Interval_T),Interval_cmp_low);
	last_donor = last_antiacceptor = 0U;
	for (i = 0; i < nintrons1; i++) {
	  interval = intervals[i];
	  chrpos = Interval_low(interval);
	  position = chrpos + chroffset;

	  if (position >= chrhigh) {
	    chr = Univ_IIT_label(chromosome_iit,chrnum,&allocp);
	    fprintf(stderr,"\nSplice site %s:%u extends beyond chromosome length %u.  Discarding...",
		    chr,chrpos,chrlength);
	    if (allocp) FREE(chr);

	  } else if (Interval_sign(interval) > 0) {
	    if (position != last_donor) {
	      last_donor = position;
	      k++;
	    }
	  } else {
	    if (position != last_antiacceptor) {
	      last_antiacceptor = position;
	      k++;
	    }
	  }
	}

	qsort(intervals,nintrons1,sizeof(Interval_T),Interval_cmp_high);
	last_acceptor = last_antidonor = 0U;
	for (i = 0; i < nintrons1; i++) {
	  interval = intervals[i];
	  chrpos = Interval_high(interval) - INTRON_HIGH_TO_LOW;
	  position = chrpos + chroffset;

	  if (position >= chrhigh) {
	    chr = Univ_IIT_label(chromosome_iit,chrnum,&allocp);
	    fprintf(stderr,"\nSplice site %s:%u extends beyond chromosome length %u.  Discarding...",
		    chr,chrpos,chrlength);
	    if (allocp) FREE(chr);

	  } else if (Interval_sign(interval) > 0) {
	    if (position != last_acceptor) {
	      last_acceptor = position;
	      k++;
	    }
	  } else {
	    if (position != last_antidonor) {
	      last_antidonor = position;
	      k++;
	    }
	  }
	}

	FREE(intervals);
	FREE(introns1);
      }
    }
  }

  *nsplicesites = k;
  debug(printf("total unique splicesites: %d\n",*nsplicesites));
  fprintf(stderr,"%d unique splicesites...",*nsplicesites);

  /* The above procedure determines number of unique splicesites */

  if (*nsplicesites == 0) {
#ifdef GSNAP
    *splicecomp = (Genomecomp_T *) NULL;
#endif
    splicesites = (Univcoord_T *) NULL;
    *splicetypes = (Splicetype_T *) NULL;
    *splicedists = (Chrpos_T *) NULL;
    *splicestrings = (List_T *) NULL;
    *splicefrags_ref = (Genomecomp_T *) NULL;
    *splicefrags_alt = (Genomecomp_T *) NULL;

    return splicesites;
  }

#ifdef GSNAP
  nblocks = (Genome_totallength(genome)+31)/32U;
  *splicecomp = (Genomecomp_T *) CALLOC(nblocks,sizeof(Genomecomp_T));
#endif
  splicesites = (Univcoord_T *) CALLOC((*nsplicesites) + 1,sizeof(Univcoord_T));
  *splicetypes = (Splicetype_T *) CALLOC(*nsplicesites,sizeof(Splicetype_T));
  *splicedists = (Chrpos_T *) CALLOC(*nsplicesites,sizeof(Chrpos_T));
  *splicestrings = (List_T *) CALLOC(*nsplicesites,sizeof(List_T));
  *splicefrags_ref = (Genomecomp_T *) CALLOC(*nsplicesites,sizeof(Genomecomp_T));
  if (genomealt == NULL) {
    *splicefrags_alt = *splicefrags_ref;
  } else {
    *splicefrags_alt = (Genomecomp_T *) CALLOC(*nsplicesites,sizeof(Genomecomp_T));
  }


  k = 0;
  for (chrnum = 1; chrnum <= Univ_IIT_total_nintervals(chromosome_iit); chrnum++) {
    if ((divno = splicing_divint_crosstable[chrnum]) > 0) {
      Univ_IIT_interval_bounds(&chroffset,&chrhigh,&chrlength,chromosome_iit,chrnum,/*circular_typeint*/-1);
      /* chrlength = Univ_IIT_length(chromosome_iit,chrnum); */
      introns1 = IIT_get_with_divno(&nintrons1,splicing_iit,divno,0U,chrlength-1U,/*sortp*/false);
      if (nintrons1 > 0) {
	intervals = (Interval_T *) CALLOC(nintrons1,sizeof(Interval_T));
	for (i = 0; i < nintrons1; i++) {
	  intervals[i] = &(splicing_iit->intervals[divno][i]);
	}

	qsort(intervals,nintrons1,sizeof(Interval_T),Interval_cmp_low);
	last_donor = last_antiacceptor = 0U;
	for (i = 0; i < nintrons1; i++) {
	  interval = intervals[i];
	  chrpos = Interval_low(interval);
	  position = chrpos + chroffset;

	  if (position >= chrhigh) {
#if 0
	    /* Warning given previously */
	    chr = Univ_IIT_label(chromosome_iit,chrnum,&allocp);
	    fprintf(stderr,"\nSplice site %s:%u extends beyond chromosome length %u.  Discarding...",
		    chr,chrpos,chrlength);
	    if (allocp) FREE(chr);
#endif

	  } else if (Interval_sign(interval) > 0) {
	    if (position == last_donor) {
	      if (Interval_length(interval) > (*splicedists)[last_donor_k]) {
		(*splicedists)[last_donor_k] = Interval_length(interval);
	      }

	    } else {
	      last_donor_k = k;
	      last_donor = splicesites[k] = position;
#ifdef GSNAP
	      assert(position/32U < (unsigned int) nblocks);
	      (*splicecomp)[position/32U] |= (1 << (position % 32));
#endif
	      (*splicetypes)[k] = DONOR;
	      (*splicedists)[k] = Interval_length(interval);

	      saw_n_p = false;
	      Genome_fill_buffer_simple(genome,position-16,16,gbuffer_ref);
	      refstring = compress16(&saw_n_p,gbuffer_ref);
	      (*splicefrags_ref)[k] = unshuffle16(&saw_n_p,gbuffer_ref);
	      if (genomealt) {
		Genome_fill_buffer_simple_alt(genome,genomealt,position-16,16,gbuffer_alt);
		altstring = compress16(&saw_n_p,gbuffer_alt);
		(*splicefrags_alt)[k] = unshuffle16(&saw_n_p,gbuffer_alt);
	      }

	      if (saw_n_p == true) {
		chr = Univ_IIT_label(chromosome_iit,chrnum,&allocp);
		fprintf(stderr,"\nNon-standard nucleotide N near splice site %s:%u.  Discarding...",
			chr,chrpos);
		if (allocp) FREE(chr);
	      } else if (genomealt) {
		(*splicestrings)[k] = allelic_combinations(refstring,altstring,position,k,splicestringpool);
		k++;
	      } else {
		(*splicestrings)[k] = Splicestringpool_push(NULL,splicestringpool,refstring,position,k);
		k++;
	      }
	    }

	  } else {
	    if (position == last_antiacceptor) {
	      if (Interval_length(interval) > (*splicedists)[last_antiacceptor_k]) {
		(*splicedists)[last_antiacceptor_k] = Interval_length(interval);
	      }

	    } else {
	      last_antiacceptor_k = k;
	      last_antiacceptor = splicesites[k] = position;
#ifdef GSNAP
	      assert(position/32U < (unsigned int) nblocks);
	      (*splicecomp)[position/32U] |= (1 << (position % 32));
#endif
	      (*splicetypes)[k] = ANTIACCEPTOR;
	      (*splicedists)[k] = Interval_length(interval);

	      saw_n_p = false;
	      Genome_fill_buffer_simple(genome,position-16,16,gbuffer_ref);
	      refstring = compress16(&saw_n_p,gbuffer_ref);
	      (*splicefrags_ref)[k] = unshuffle16(&saw_n_p,gbuffer_ref);
	      if (genomealt) {
		Genome_fill_buffer_simple_alt(genome,genomealt,position-16,16,gbuffer_alt);
		altstring = compress16(&saw_n_p,gbuffer_alt);
		(*splicefrags_alt)[k] = unshuffle16(&saw_n_p,gbuffer_alt);
	      }

	      if (saw_n_p == true) {
		chr = Univ_IIT_label(chromosome_iit,chrnum,&allocp);
		fprintf(stderr,"\nNon-standard nucleotide N near splice site %s:%u.  Discarding...",
			chr,chrpos);
		if (allocp) FREE(chr);
	      } else if (genomealt) {
		(*splicestrings)[k] = allelic_combinations(refstring,altstring,position,k,splicestringpool);
		k++;
	      } else {
		(*splicestrings)[k] = Splicestringpool_push(NULL,splicestringpool,refstring,position,k);
		k++;
	      }
	    }

	  }
	}

	qsort(intervals,nintrons1,sizeof(Interval_T),Interval_cmp_high);
	last_acceptor = last_antidonor = 0U;
	for (i = 0; i < nintrons1; i++) {
	  interval = intervals[i];
	  chrpos = Interval_high(interval) - INTRON_HIGH_TO_LOW;
	  position = chrpos + chroffset;

	  if (position >= chrhigh) {
#if 0
	    /* Warning given previously */
	    chr = Univ_IIT_label(chromosome_iit,chrnum,&allocp);
	    fprintf(stderr,"\nSplice site %s:%u extends beyond chromosome length %u.  Discarding...",
		    chr,chrpos,chrlength);
	    if (allocp) FREE(chr);
#endif

	  } else if (Interval_sign(interval) > 0) {
	    if (position == last_acceptor) {
	      if (Interval_length(interval) > (*splicedists)[last_acceptor_k]) {
		(*splicedists)[last_acceptor_k] = Interval_length(interval);
	      }

	    } else {
	      last_acceptor_k = k;
	      last_acceptor = splicesites[k] = position;
#ifdef GSNAP
	      assert(position/32U < (unsigned int) nblocks);
	      (*splicecomp)[position/32U] |= (1 << (position % 32));
#endif
	      (*splicetypes)[k] = ACCEPTOR;
	      (*splicedists)[k] = Interval_length(interval);

	      saw_n_p = false;
	      Genome_fill_buffer_simple(genome,position,16,gbuffer_ref);
	      refstring = uint4_reverse(compress16(&saw_n_p,gbuffer_ref));
	      (*splicefrags_ref)[k] = unshuffle16(&saw_n_p,gbuffer_ref);
	      if (genomealt) {
		Genome_fill_buffer_simple_alt(genome,genomealt,position,16,gbuffer_alt);
		altstring = uint4_reverse(compress16(&saw_n_p,gbuffer_alt));
		(*splicefrags_alt)[k] = unshuffle16(&saw_n_p,gbuffer_alt);
	      }

	      if (saw_n_p == true) {
		chr = Univ_IIT_label(chromosome_iit,chrnum,&allocp);
		fprintf(stderr,"\nNon-standard nucleotide N near splice site %s:%u.  Discarding...",
			chr,chrpos);
		if (allocp) FREE(chr);
	      } else if (genomealt) {
		(*splicestrings)[k] = allelic_combinations(refstring,altstring,position,k,splicestringpool);
		k++;
	      } else {
		(*splicestrings)[k] = Splicestringpool_push(NULL,splicestringpool,refstring,position,k);
		k++;
	      }
	    }

	  } else {
	    if (position == last_antidonor) {
	      if (Interval_length(interval) > (*splicedists)[last_antidonor_k]) {
		(*splicedists)[last_antidonor_k] = Interval_length(interval);
	      }

	    } else {
	      last_antidonor_k = k;
	      last_antidonor = splicesites[k] = position;
#ifdef GSNAP
	      assert(position/32U < (unsigned int) nblocks);
	      (*splicecomp)[position/32U] |= (1 << (position % 32));
#endif
	      (*splicetypes)[k] = ANTIDONOR;
	      (*splicedists)[k] = Interval_length(interval);

	      saw_n_p = false;
	      Genome_fill_buffer_simple(genome,position,16,gbuffer_ref);
	      refstring = uint4_reverse(compress16(&saw_n_p,gbuffer_ref));
	      (*splicefrags_ref)[k] = unshuffle16(&saw_n_p,gbuffer_ref);
	      if (genomealt) {
		Genome_fill_buffer_simple_alt(genome,genomealt,position,16,gbuffer_alt);
		altstring = uint4_reverse(compress16(&saw_n_p,gbuffer_alt));
		(*splicefrags_alt)[k] = unshuffle16(&saw_n_p,gbuffer_alt);
	      }

	      if (saw_n_p == true) {
		chr = Univ_IIT_label(chromosome_iit,chrnum,&allocp);
		fprintf(stderr,"\nNon-standard nucleotide N near splice site %s:%u.  Discarding...",
			chr,chrpos);
		if (allocp) FREE(chr);
	      } else if (genomealt) {
		(*splicestrings)[k] = allelic_combinations(refstring,altstring,position,k,splicestringpool);
		k++;
	      } else {
		(*splicestrings)[k] = Splicestringpool_push(NULL,splicestringpool,refstring,position,k);
		k++;
	      }
	    }
	  }
	}

	FREE(intervals);
	FREE(introns1);
      }
    }
  }

  *nsplicesites = k;
  fprintf(stderr,"\n%d splicesites are valid...",*nsplicesites);


  /* Need to sort by individual splicesites */
  cells = (struct Cell_T *) CALLOC(*nsplicesites,sizeof(struct Cell_T));
  for (k = 0; k < *nsplicesites; k++) {
    cells[k].k = k;
    cells[k].pos = splicesites[k];
  }
  qsort(cells,*nsplicesites,sizeof(struct Cell_T),Cell_position_cmp);


  /* Save unordered information */
  temp_splicesites = splicesites;
  temp_splicetypes = *splicetypes;
  temp_splicedists = *splicedists;
  temp_splicestrings = *splicestrings;
  temp_splicefrags_ref = *splicefrags_ref;
  if (genomealt == NULL) {
    temp_splicefrags_alt = temp_splicefrags_ref;
  } else {
    temp_splicefrags_alt = *splicefrags_alt;
  }

  /* Allocate ordered information */
  splicesites = (Univcoord_T *) CALLOC((*nsplicesites) + 1,sizeof(Univcoord_T));
  *splicetypes = (Splicetype_T *) CALLOC(*nsplicesites,sizeof(Splicetype_T));
  *splicedists = (Chrpos_T *) CALLOC(*nsplicesites,sizeof(Chrpos_T));
  *splicestrings = (List_T *) CALLOC(*nsplicesites,sizeof(List_T));
  *splicefrags_ref = (UINT4 *) CALLOC(*nsplicesites,sizeof(UINT4));
  if (genomealt == NULL) {
    *splicefrags_alt = *splicefrags_ref;
  } else {
    *splicefrags_alt = (UINT4 *) CALLOC(*nsplicesites,sizeof(UINT4));
  }


  /* Order all information */
  for (j = 0; j < *nsplicesites; j++) {
    k = cells[j].k;
    splicesites[j] = temp_splicesites[k];
    (*splicetypes)[j] = temp_splicetypes[k];
    (*splicedists)[j] = temp_splicedists[k];
    (*splicestrings)[j] = temp_splicestrings[k];
    for (p = (*splicestrings)[j]; p != NULL; p = List_next(p)) {
      /* Fix reference back to splicesite */
      splicestring = (Splicestring_T) List_head(p);
      splicestring->splicesite_i = j;
    }
    (*splicefrags_ref)[j] = temp_splicefrags_ref[k];
    (*splicefrags_alt)[j] = temp_splicefrags_alt[k];
  }

  splicesites[*nsplicesites] = (Univcoord_T) -1; /* Marker for comparison in identify_all_segments */

  FREE(cells);
  FREE(temp_splicesites);
  FREE(temp_splicetypes);
  FREE(temp_splicedists);
  FREE(temp_splicestrings);
  FREE(temp_splicefrags_ref);
  if (genomealt != NULL) {
    FREE(temp_splicefrags_alt);
  }


#ifdef DEBUG2
  for (k = 0; k < *nsplicesites; k++) {
    printf("%d: %u %s %08X\n",k,splicesites[k],Splicetype_string((*splicetypes)[k]),(*splicefrags_ref)[k]);
    for (p = (*splicestrings)[k]; p != NULL; p = List_next(p)) {
      splicestring = (Splicestring_T) List_head(p);
      printf("  %u %u %u\n",splicestring->string,splicestring->splicesite,splicestring->splicesite_i);
    }
  }
#endif

  return splicesites;
}


static bool
get_exons (Uintlist_T *donors, Uintlist_T *acceptors, char *annot) {
  Chrpos_T exonstart, exonend;
  char *p;

  /* Skip header */
  p = annot;
  while (*p != '\0' && *p != '\n') {
    p++;
  }
  if (*p == '\n') p++;

  *donors = *acceptors = (Uintlist_T) NULL;
  while (*p != '\0') {
    if (sscanf(p,"%u %u",&exonstart,&exonend) != 2) {
      return false;
    } else {
      *donors = Uintlist_push(*donors,exonend);
      *acceptors = Uintlist_push(*acceptors,exonstart);
    }
      
    /* Advance to next exon to see if this is the last one */
    while (*p != '\0' && *p != '\n') p++;
    if (*p == '\n') p++;
  }

  if (*donors == NULL || *acceptors == NULL) {
    return false;
  } else {
    *donors = Uintlist_pop(*donors,&exonend);
    *donors = Uintlist_reverse(*donors);

    *acceptors = Uintlist_reverse(*acceptors);
    *acceptors = Uintlist_pop(*acceptors,&exonstart);

    return true;
  }
}




/************************************************************************/



typedef struct Trie_T *Trie_T;
struct Trie_T {
  Splicestring_T leaf;

  int nsites;
  int na;
  int nc;
  int ng;
  int nt;

  Trie_T triea;
  Trie_T triec;
  Trie_T trieg;
  Trie_T triet;
};


#if 0
static void
Trie_free (Trie_T *old) {
  if ((*old) == NULL) {
    return;
  } else if ((*old)->leaf != NULL) {
    FREE(*old);
    return;
  } else {
    Trie_free(&(*old)->triea);
    Trie_free(&(*old)->triec);
    Trie_free(&(*old)->trieg);
    Trie_free(&(*old)->triet);
    FREE(*old);
    return;
  }
}
#endif

/************************************************************************
 *   Building splicetries
 ************************************************************************/

#if 0
static Trie_T
Trie_new (Splicestring_T *sites, int nsites, int charpos) {
  Trie_T trie;
  int i, ptr;

  if (nsites == 0) {
    return (Trie_T) NULL;

  } else {

    if (nsites == 1) {
      trie = (Trie_T) MALLOC(sizeof(*trie));
      trie->nsites = 1;
      trie->na = trie->nc = trie->ng = trie->nt = 0;
      trie->triea = trie->triec = trie->trieg = trie->triet = (Trie_T) NULL;
      trie->leaf = sites[0];
      /* trie->remainder = &(sites[0]->string[charpos]); */
      return trie;

    } else if (charpos >= 16) {
      return (Trie_T) NULL;

    } else {
      trie = (Trie_T) MALLOC(sizeof(*trie));
      trie->nsites = nsites;
      trie->na = trie->nc = trie->ng = trie->nt = 0;
      trie->leaf = (Splicestring_T) NULL;
      /* trie->remainder = (char *) NULL; */

      for (i = 0; i < nsites; i++) {
	switch ((sites[i]->string >> (30 - 2*charpos)) & 0x03) {
	case RIGHT_A: trie->na++; break;
	case RIGHT_C: trie->nc++; break;
	case RIGHT_G: trie->ng++; break;
	case RIGHT_T: trie->nt++; break;
	default: abort();
	}
      }

      ptr = 0;
      trie->triea = Trie_new(&(sites[ptr]),trie->na,charpos+1);
      ptr += trie->na;

      trie->triec = Trie_new(&(sites[ptr]),trie->nc,charpos+1);
      ptr += trie->nc;

      trie->trieg = Trie_new(&(sites[ptr]),trie->ng,charpos+1);
      ptr += trie->ng;

      trie->triet = Trie_new(&(sites[ptr]),trie->nt,charpos+1);
      /* ptr += trie->nt; */

      return trie;
    }
  }
}
#endif


/* Based on original Trie_output_new, which used repeated calls to Uintlist_push */
static int
Trie_output_size (int size, Splicestring_T *sites, int nsites, int charpos) {
  int i, k;
  unsigned int posa, posc, posg, post;
  int na, nc, ng, nt;

  if (nsites == 0) {
    debug0(printf("nsites == 0, so not adding anything\n"));

  } else if (nsites == 1) {
    debug0(printf("nsites == 1, so adding 1\n"));
    size += 1;
    
  } else if (charpos >= 16) {
    if (nsites > MAX_DUPLICATES) {
      fprintf(stderr,"Warning: Splicetrie exceeded max duplicates value of %d\n",MAX_DUPLICATES);
    } else {
      size += 1;

      for (i = 0; i < nsites; i++) {
	debug0(printf(" %d",sites[i]->splicesite_i));
	size += 1;
      }
      debug0(printf("\n"));
    }
      
  } else {
    na = nc = ng = nt = 0;
    for (i = 0; i < nsites; i++) {
      switch ((sites[i]->string >> (30 - 2*charpos)) & 0x03) {
      case RIGHT_A: na++; break;
      case RIGHT_C: nc++; break;
      case RIGHT_G: ng++; break;
      case RIGHT_T: nt++; break;
      default: abort();
      }
    }
    debug0(printf("%d A, %d C, %d G, %d T\n",na,nc,ng,nt));

    k = 0;
    size = Trie_output_size(size,&(sites[k]),na,charpos+1);

    k += na;
    size = Trie_output_size(size,&(sites[k]),nc,charpos+1);

    k += nc;
    size = Trie_output_size(size,&(sites[k]),ng,charpos+1);

    k += ng;
    size = Trie_output_size(size,&(sites[k]),nt,charpos+1);
    /* k += nt; */

    size += 5;
  }

  return size;
}



#ifdef USE_LIST
/* Combination of Trie_new and Trie_output */
/* Note that when *ptr = NULL_POINTER, we do not need to advance
   nprinted, because no entry is made into triecontents_list */
static Uintlist_T
Trie_output_list (unsigned int *ptr, int *nprinted, Uintlist_T triecontents_list,
		  Splicestring_T *sites, int nsites, int charpos) {
  int i, k;
  unsigned int posa, posc, posg, post;
  int na, nc, ng, nt;

  if (nsites == 0) {
    debug0(printf("nsites == 0, so NULL\n"));
    *ptr = NULL_POINTER;

  } else if (nsites == 1) {
    debug0(printf("nsites == 1, so pushing %d\n",sites[0]->splicesite_i));
    *ptr = *nprinted;
    triecontents_list = Uintlist_push(triecontents_list,sites[0]->splicesite_i);
    *nprinted += 1;
    
  } else if (charpos >= 16) {
    if (nsites > MAX_DUPLICATES) {
      fprintf(stderr,"Warning: Splicetrie exceeded max duplicates value of %d\n",MAX_DUPLICATES);
      *ptr = NULL_POINTER;
    } else {
      *ptr = *nprinted;

      triecontents_list = Uintlist_push(triecontents_list,(Triecontent_T) (-nsites));
      *nprinted += 1;
      for (i = 0; i < nsites; i++) {
	debug0(printf(" %d",sites[i]->splicesite_i));
	triecontents_list = Uintlist_push(triecontents_list,sites[i]->splicesite_i);
	*nprinted += 1;
      }
      debug0(printf("\n"));
    }
      
  } else {
    na = nc = ng = nt = 0;
    for (i = 0; i < nsites; i++) {
      switch ((sites[i]->string >> (30 - 2*charpos)) & 0x03) {
      case RIGHT_A: na++; break;
      case RIGHT_C: nc++; break;
      case RIGHT_G: ng++; break;
      case RIGHT_T: nt++; break;
      default: abort();
      }
    }
    debug0(printf("%d A, %d C, %d G, %d T\n",na,nc,ng,nt));

    k = 0;
    triecontents_list = Trie_output_list(&posa,&(*nprinted),triecontents_list,
					&(sites[k]),na,charpos+1);
    k += na;
    triecontents_list = Trie_output_list(&posc,&(*nprinted),triecontents_list,
					&(sites[k]),nc,charpos+1);
    k += nc;
    triecontents_list = Trie_output_list(&posg,&(*nprinted),triecontents_list,
					&(sites[k]),ng,charpos+1);
    k += ng;
    triecontents_list = Trie_output_list(&post,&(*nprinted),triecontents_list,
					&(sites[k]),nt,charpos+1);
    /* k += nt; */

    *ptr = *nprinted;
    triecontents_list = Uintlist_push(triecontents_list,INTERNAL_NODE);

    if (posa == NULL_POINTER) {
      triecontents_list = Uintlist_push(triecontents_list,EMPTY_POINTER);
    } else {
      triecontents_list = Uintlist_push(triecontents_list,*ptr - posa);
    }
    if (posc == NULL_POINTER) {
      triecontents_list = Uintlist_push(triecontents_list,EMPTY_POINTER);
    } else {
      triecontents_list = Uintlist_push(triecontents_list,*ptr - posc);
    }
    if (posg == NULL_POINTER) {
      triecontents_list = Uintlist_push(triecontents_list,EMPTY_POINTER);
    } else {
      triecontents_list = Uintlist_push(triecontents_list,*ptr - posg);
    }
    if (post == NULL_POINTER) {
      triecontents_list = Uintlist_push(triecontents_list,EMPTY_POINTER);
    } else {
      triecontents_list = Uintlist_push(triecontents_list,*ptr - post);
    }
    *nprinted += 5;
  }

  return triecontents_list;
}

#else

/* ptr is for offsets */
static unsigned int *
Trie_output_array (unsigned int *ptr, int *nprinted, unsigned int *contents,
		   Splicestring_T *sites, int nsites, int charpos) {
  int i, k;
  unsigned int posa, posc, posg, post;
  int na, nc, ng, nt;

  if (nsites == 0) {
    debug0(printf("nsites == 0, so NULL\n"));
    *ptr = NULL_POINTER;

  } else if (nsites == 1) {
    debug0(printf("nsites == 1, so pushing %d\n",sites[0]->splicesite_i));
    *ptr = *nprinted;
    *contents++ = sites[0]->splicesite_i;
    *nprinted += 1;
    
  } else if (charpos >= 16) {
    if (nsites > MAX_DUPLICATES) {
      fprintf(stderr,"Warning: Splicetrie exceeded max duplicates value of %d\n",MAX_DUPLICATES);
      *ptr = NULL_POINTER;
    } else {
      *ptr = *nprinted;

      *contents++ = (Triecontent_T) (-nsites);
      *nprinted += 1;
      for (i = 0; i < nsites; i++) {
	debug0(printf(" %d",sites[i]->splicesite_i));
	*contents++ = sites[i]->splicesite_i;
	*nprinted += 1;
      }
      debug0(printf("\n"));
    }
      
  } else {
    na = nc = ng = nt = 0;
    for (i = 0; i < nsites; i++) {
      switch ((sites[i]->string >> (30 - 2*charpos)) & 0x03) {
      case RIGHT_A: na++; break;
      case RIGHT_C: nc++; break;
      case RIGHT_G: ng++; break;
      case RIGHT_T: nt++; break;
      default: abort();
      }
    }
    debug0(printf("%d A, %d C, %d G, %d T\n",na,nc,ng,nt));

    k = 0;
    contents = Trie_output_array(&posa,&(*nprinted),contents,&(sites[k]),na,charpos+1);

    k += na;
    contents = Trie_output_array(&posc,&(*nprinted),contents,&(sites[k]),nc,charpos+1);

    k += nc;
    contents = Trie_output_array(&posg,&(*nprinted),contents,&(sites[k]),ng,charpos+1);

    k += ng;
    contents = Trie_output_array(&post,&(*nprinted),contents,&(sites[k]),nt,charpos+1);
    /* k += nt; */

    *ptr = *nprinted;
    *contents++ = INTERNAL_NODE;

    if (posa == NULL_POINTER) {
      *contents++ = EMPTY_POINTER;
    } else {
      *contents++ = *ptr - posa;
    }
    if (posc == NULL_POINTER) {
      *contents++ = EMPTY_POINTER;
    } else {
      *contents++ = *ptr - posc;
    }
    if (posg == NULL_POINTER) {
      *contents++ = EMPTY_POINTER;
    } else {
      *contents++ = *ptr - posg;
    }
    if (post == NULL_POINTER) {
      *contents++ = EMPTY_POINTER;
    } else {
      *contents++ = *ptr - post;
    }
    *nprinted += 5;
  }

  return contents;
}
#endif


#if 0
static void
Trie_print (Trie_T trie, int level) {
  int i;

  if (trie == NULL) {
    printf("\n");

  } else {
    if (trie->leaf != NULL) {
      printf(" %d@%u\n",(int) trie->leaf->splicesite_i,trie->leaf->splicesite);

    } else {
      printf("\n");
      for (i = 0; i < level; i++) {
	printf("\t");
      }
      printf("A (%d)",trie->na);
      Trie_print(trie->triea,level+1);

      for (i = 0; i < level; i++) {
	printf("\t");
      }
      printf("C (%d)",trie->nc);
      Trie_print(trie->triec,level+1);

      for (i = 0; i < level; i++) {
	printf("\t");
      }
      printf("G (%d)",trie->ng);
      Trie_print(trie->trieg,level+1);

      for (i = 0; i < level; i++) {
	printf("\t");
      }
      printf("T (%d)",trie->nt);
      Trie_print(trie->triet,level+1);
    }
  }

  return;
}
#endif

#if 0
static int
Trie_print_compact (int *ptr, FILE *fp, Trie_T trie, int nprinted) {
  int posa, posc, posg, post;

  if (trie == NULL) {
    *ptr = NULL_POINTER;
    return nprinted;

  } else {
    if (trie->leaf != NULL) {
      printf("%d. %d@%u\n",nprinted,(int) trie->leaf->splicesite_i,trie->leaf->splicesite);
      *ptr = nprinted;
      return nprinted + 1;

    } else {
      nprinted = Trie_print_compact(&posa,fp,trie->triea,nprinted);
      nprinted = Trie_print_compact(&posc,fp,trie->triec,nprinted);
      nprinted = Trie_print_compact(&posg,fp,trie->trieg,nprinted);
      nprinted = Trie_print_compact(&post,fp,trie->triet,nprinted);

      *ptr = nprinted;

      printf("%d. %u\n",nprinted,INTERNAL_NODE);

      printf("relative: %d %d %d %d\n",
	     *ptr - posa,*ptr - posc,*ptr - posg,*ptr - post);
      printf("absolute: %d %d %d %d\n",
	     posa,posc,posg,post);

      return nprinted + 3;
    }
  }
}
#endif


#ifdef USE_2BYTE_RELOFFSETS
#define TRIE_EMPTY_SIZE 3
#else
#define TRIE_EMPTY_SIZE 5
#endif


#ifdef USE_LIST
static Uintlist_T
Trie_output_empty_list (unsigned int *ptr, int *nprinted, Uintlist_T triecontents_list) {
  *ptr = (unsigned int) *nprinted;
  triecontents_list = Uintlist_push(triecontents_list,INTERNAL_NODE);
#ifdef USE_2BYTE_RELOFFSETS
  triecontents_list = Uintlist_push(triecontents_list,EMPTY_POINTER);
  triecontents_list = Uintlist_push(triecontents_list,EMPTY_POINTER);
  *nprinted += 3;
#else
  triecontents_list = Uintlist_push(triecontents_list,EMPTY_POINTER);
  triecontents_list = Uintlist_push(triecontents_list,EMPTY_POINTER);
  triecontents_list = Uintlist_push(triecontents_list,EMPTY_POINTER);
  triecontents_list = Uintlist_push(triecontents_list,EMPTY_POINTER);
  *nprinted += 5;
#endif
  return triecontents_list;
}

#else
static unsigned int *
Trie_output_empty_array (unsigned int *ptr, int *nprinted, unsigned int *contents) {
  *ptr = (unsigned int) *nprinted;
  *contents++ = INTERNAL_NODE;
#ifdef USE_2BYTE_RELOFFSETS
  *contents++ = EMPTY_POINTER;
  *contents++ = EMPTY_POINTER;
  *nprinted += 3;
#else
  *contents++ = EMPTY_POINTER;
  *contents++ = EMPTY_POINTER;
  *contents++ = EMPTY_POINTER;
  *contents++ = EMPTY_POINTER;
  *nprinted += 5;
#endif
  return contents;
}

#endif


#if 0
static Uintlist_T
Trie_output (int *ptr, int *nprinted, Uintlist_T triecontents_list, Trie_T trie) {
  int posa, posc, posg, post;
#ifdef USE_2BYTE_RELOFFSETS
  unsigned int pointersAC, pointersGT;
  int reloffset;
#endif

  if (trie == NULL) {
    *ptr = NULL_POINTER;

  } else if (trie->leaf != NULL) {
    *ptr = *nprinted;
    /* Entering splicesite_i, rather than position */
    triecontents_list = Uintlist_push(triecontents_list,trie->leaf->splicesite_i);
    *nprinted += 1;

  } else {
    triecontents_list = Trie_output(&posa,&(*nprinted),triecontents_list,trie->triea);
    triecontents_list = Trie_output(&posc,&(*nprinted),triecontents_list,trie->triec);
    triecontents_list = Trie_output(&posg,&(*nprinted),triecontents_list,trie->trieg);
    triecontents_list = Trie_output(&post,&(*nprinted),triecontents_list,trie->triet);

    *ptr = *nprinted;
    triecontents_list = Uintlist_push(triecontents_list,INTERNAL_NODE);


#ifdef USE_2BYTE_RELOFFSETS
    if (posa == NULL_POINTER) {
      pointersAC = 0;
    } else if ((reloffset = *ptr - posa) >= 65536) {
      fprintf(stderr,"Reloffset A %d = %d - %d is too big for 2 bytes\n",reloffset,*ptr,posa);
      abort();
    } else {
      pointersAC = reloffset;
    }

    if (posc == NULL_POINTER) {
      pointersAC = (pointersAC << 16);
    } else if ((reloffset = *ptr - posc) >= 65536) {
      fprintf(stderr,"Reloffset C %d = %d - %d is too big for 2 bytes\n",reloffset,*ptr,posc);
      abort();
    } else {
      pointersAC = (pointersAC << 16) + reloffset;
    }

    if (posg == NULL_POINTER) {
      pointersGT = 0;
    } else if ((reloffset = *ptr - posg) >= 65536) {
      fprintf(stderr,"Reloffset G %d = %d - %d is too big for 2 bytes\n",reloffset,*ptr,posg);
      abort();
    } else {
      pointersGT = reloffset;
    }

    if (post == NULL_POINTER) {
      pointersGT = (pointersGT << 16);
    } else if ((reloffset = *ptr - post) >= 65536) {
      fprintf(stderr,"Reloffset T %d = %d - %d is too big for 2 bytes\n",reloffset,*ptr,post);
      abort();
    } else {
      pointersGT = (pointersGT << 16) + reloffset;
    }

    triecontents_list = Uintlist_push(triecontents_list,pointersAC);
    triecontents_list = Uintlist_push(triecontents_list,pointersGT);

    *nprinted += 3;
#else
    if (posa == NULL_POINTER) {
      triecontents_list = Uintlist_push(triecontents_list,EMPTY_POINTER);
    } else {
      triecontents_list = Uintlist_push(triecontents_list,*ptr - posa);
    }
    if (posc == NULL_POINTER) {
      triecontents_list = Uintlist_push(triecontents_list,EMPTY_POINTER);
    } else {
      triecontents_list = Uintlist_push(triecontents_list,*ptr - posc);
    }
    if (posg == NULL_POINTER) {
      triecontents_list = Uintlist_push(triecontents_list,EMPTY_POINTER);
    } else {
      triecontents_list = Uintlist_push(triecontents_list,*ptr - posg);
    }
    if (post == NULL_POINTER) {
      triecontents_list = Uintlist_push(triecontents_list,EMPTY_POINTER);
    } else {
      triecontents_list = Uintlist_push(triecontents_list,*ptr - post);
    }
    *nprinted += 5;
#endif

  }
  return triecontents_list;
}
#endif


void
Splicetrie_npartners (int **nsplicepartners_skip, int **nsplicepartners_obs, int **nsplicepartners_max,
		      Univcoord_T *splicesites, Splicetype_T *splicetypes,
		      Chrpos_T *splicedists, List_T *splicestrings, int nsplicesites,
		      Univ_IIT_T chromosome_iit, Chrpos_T max_distance, bool distances_observed_p) {
  int nsites_skip, nsites_obs, nsites_max, j, j1;
  Univcoord_T chroffset, chrhigh, leftbound_obs, leftbound_max, rightbound_obs, rightbound_max;
  Univcoord_T leftbound_min, rightbound_min;
  Chrpos_T chrlength;
  Chrnum_T chrnum;


  (*nsplicepartners_skip) = (int *) CALLOC(nsplicesites,sizeof(int));
  if (distances_observed_p == true) {
    (*nsplicepartners_obs) = (int *) CALLOC(nsplicesites,sizeof(int));
  } else {
    (*nsplicepartners_obs) = (int *) NULL;
  }
  (*nsplicepartners_max) = (int *) CALLOC(nsplicesites,sizeof(int));

  chrhigh = 0U;
  for (j = 0; j < nsplicesites; j++) {
    if (splicesites[j] > chrhigh) {
      chrnum = Univ_IIT_get_one(chromosome_iit,splicesites[j],splicesites[j]);
      Univ_IIT_interval_bounds(&chroffset,&chrhigh,&chrlength,chromosome_iit,chrnum,/*circular_typeint*/-1);
      /* chrhigh += 1U; */
    }

    switch (splicetypes[j]) {
    case DONOR:
      nsites_skip = nsites_obs = nsites_max = 0;
      j1 = j + 1;

      if ((rightbound_min = splicesites[j] + MININTRONLEN) > chrhigh) {
	rightbound_min = chrhigh;
      }
      while (j1 < nsplicesites && splicesites[j1] < rightbound_min) {
	nsites_skip++;
	j1++;
      }

      if (distances_observed_p) {
	if ((rightbound_obs = splicesites[j] + splicedists[j]) > chrhigh) {
	  rightbound_obs = chrhigh;
	}
	while (j1 < nsplicesites && splicesites[j1] < rightbound_obs) {
	  if (splicetypes[j1] == ACCEPTOR) {
	    nsites_obs += List_length(splicestrings[j1]);
	  }
	  j1++;
	}
      }

      if ((rightbound_max = splicesites[j] + max_distance) > chrhigh) {
	rightbound_max = chrhigh;
      }
      while (j1 < nsplicesites && splicesites[j1] < rightbound_max) {
	if (splicetypes[j1] == ACCEPTOR) {
	  nsites_max += List_length(splicestrings[j1]);
	}
	j1++;
      }

      break;

    case ACCEPTOR:
      nsites_skip = nsites_obs = nsites_max = 0;
      j1 = j - 1;

      if (splicesites[j] < chroffset + MININTRONLEN) {
	leftbound_min = chroffset;
      } else {
	leftbound_min = splicesites[j] - MININTRONLEN;
      }
      while (j1 >= 0 && splicesites[j1] > leftbound_min) {
	nsites_skip++;
	j1--;
      }

      if (distances_observed_p) {
	if (splicesites[j] < chroffset + splicedists[j]) {
	  leftbound_obs = chroffset;
	} else {
	  leftbound_obs = splicesites[j] - splicedists[j];
	}
	while (j1 >= 0 && splicesites[j1] > leftbound_obs) {
	  if (splicetypes[j1] == DONOR) {
	    nsites_obs += List_length(splicestrings[j1]);
	  }
	  j1--;
	}
      }

      if (splicesites[j] < chroffset + max_distance) {
	leftbound_max = chroffset;
      } else {
	leftbound_max = splicesites[j] - max_distance;
      }
      while (j1 >= 0 && splicesites[j1] > leftbound_max) {
	if (splicetypes[j1] == DONOR) {
	  nsites_max += List_length(splicestrings[j1]);
	}
	j1--;
      }

      break;

    case ANTIDONOR:
      nsites_skip = nsites_obs = nsites_max = 0;
      j1 = j - 1;

      if (splicesites[j] < chroffset + MININTRONLEN) {
	leftbound_min = chroffset;
      } else {
	leftbound_min = splicesites[j] - MININTRONLEN;
      }
      while (j1 >= 0 && splicesites[j1] > leftbound_min) {
	nsites_skip++;
	j1--;
      }

      if (distances_observed_p) {
	if (splicesites[j] < chroffset + splicedists[j]) {
	  leftbound_obs = chroffset;
	} else {
	  leftbound_obs = splicesites[j] - splicedists[j];
	}
	while (j1 >= 0 && splicesites[j1] > leftbound_obs) {
	  if (splicetypes[j1] == ANTIACCEPTOR) {
	    nsites_obs += List_length(splicestrings[j1]);
	  }
	  j1--;
	}
      }

      if (splicesites[j] < chroffset + max_distance) {
	leftbound_max = chroffset;
      } else {
	leftbound_max = splicesites[j] - max_distance;
      }
      while (j1 >= 0 && splicesites[j1] > leftbound_max) {
	if (splicetypes[j1] == ANTIACCEPTOR) {
	  nsites_max += List_length(splicestrings[j1]);
	}
	j1--;
      }

      break;

    case ANTIACCEPTOR:
      nsites_skip = nsites_obs = nsites_max = 0;
      j1 = j + 1;

      if ((rightbound_min = splicesites[j] + MININTRONLEN) > chrhigh) {
	rightbound_min = chrhigh;
      }
      while (j1 < nsplicesites && splicesites[j1] < rightbound_min) {
	nsites_skip++;
	j1++;
      }

      if (distances_observed_p) {
	if ((rightbound_obs = splicesites[j] + splicedists[j]) > chrhigh) {
	  rightbound_obs = chrhigh;
	}
	while (j1 < nsplicesites && splicesites[j1] < rightbound_obs) {
	  if (splicetypes[j1] == ANTIDONOR) {
	    nsites_obs += List_length(splicestrings[j1]);
	  }
	  j1++;
	}
      }

      if ((rightbound_max = splicesites[j] + max_distance) > chrhigh) {
	rightbound_max = chrhigh;
      }
      while (j1 < nsplicesites && splicesites[j1] < rightbound_max) {
	if (splicetypes[j1] == ANTIDONOR) {
	  nsites_max += List_length(splicestrings[j1]);
	}
	j1++;
      }

      break;

    default: 
      fprintf(stderr,"Unexpected splicetype %d\n",splicetypes[j]);
      abort();
    }

    (*nsplicepartners_skip)[j] = nsites_skip;
    if (distances_observed_p) {
      (*nsplicepartners_obs)[j] = nsites_obs;
    }
    (*nsplicepartners_max)[j] = nsites_max;
  }

  return;
}


#define MAX_SITES_ALLOCATED 1000

static void
build_via_splicesites_size (int *trie_obs_size, int *trie_max_size,
			    int *nsplicepartners_skip, int *nsplicepartners_obs, int *nsplicepartners_max,
			    Splicetype_T *splicetypes, List_T *splicestrings, int nsplicesites,
			    bool distances_observed_p) {
  List_T p;
  int nsites, j, j1;
  Splicestring_T *sites, sites_allocated[MAX_SITES_ALLOCATED];
  int npartners;

  for (j = 0; j < nsplicesites; j++) {
    switch (splicetypes[j]) {
    case DONOR:
      j1 = j + 1 + nsplicepartners_skip[j];

      if (distances_observed_p) {
	if ((npartners = nsplicepartners_obs[j]) == 0) {
#if 0
	  triecontents_obs_list = Trie_output_empty_list(&((*trieoffsets_obs)[j]),&nprinted_obs,triecontents_obs_list);
#else
	  (*trie_obs_size) += TRIE_EMPTY_SIZE;
#endif

	} else {
	  if (npartners > MAX_SITES_ALLOCATED) {
	    sites = (Splicestring_T *) CALLOC(npartners,sizeof(Splicestring_T));
	  } else {
	    sites = &(sites_allocated[0]);
	  }
	  nsites = 0;
	  while (nsites < nsplicepartners_obs[j]) {
	    if (splicetypes[j1] == ACCEPTOR) {
	      for (p = splicestrings[j1]; p != NULL; p = List_next(p)) {
		sites[nsites++] = (Splicestring_T) List_head(p);
	      }
	    }
	    j1++;
	  }
	  assert(nsites == nsplicepartners_obs[j]);
	  qsort(sites,nsites,sizeof(Splicestring_T),Splicestring_cmp);
#if 0
	  triecontents_obs_list = Trie_output_list(&((*trieoffsets_obs)[j]),&nprinted_obs,triecontents_obs_list,
						  sites,nsites,/*charpos*/0);
#else
	  *trie_obs_size = Trie_output_size(*trie_obs_size,sites,nsites,/*charpos*/0);
#endif
	  if (npartners > MAX_SITES_ALLOCATED) {
	    FREE(sites);
	  }
	}
      }

      if ((npartners = nsplicepartners_max[j]) == 0) {
#if 0
	triecontents_max_list = Trie_output_empty_list(&((*trieoffsets_max)[j]),&nprinted_max,triecontents_max_list);
#else
	(*trie_max_size) += TRIE_EMPTY_SIZE;
#endif

      } else {
	if (npartners > MAX_SITES_ALLOCATED) {
	  sites = (Splicestring_T *) CALLOC(npartners,sizeof(Splicestring_T));
	} else {
	  sites = &(sites_allocated[0]);
	}
	nsites = 0;
	while (nsites < nsplicepartners_max[j]) {
	  if (splicetypes[j1] == ACCEPTOR) {
	    for (p = splicestrings[j1]; p != NULL; p = List_next(p)) {
	      sites[nsites++] = (Splicestring_T) List_head(p);
	    }
	  }
	  j1++;
	}
	assert(nsites == nsplicepartners_max[j]);
	qsort(sites,nsites,sizeof(Splicestring_T),Splicestring_cmp);
#if 0
	triecontents_max_list = Trie_output_list(&((*trieoffsets_max)[j]),&nprinted_max,triecontents_max_list,
						sites,nsites,/*charpos*/0);
#else
	*trie_max_size = Trie_output_size(*trie_max_size,sites,nsites,/*charpos*/0);
#endif
	if (npartners > MAX_SITES_ALLOCATED) {
	  FREE(sites);
	}
      }

      break;

    case ACCEPTOR:

      j1 = j - 1 - nsplicepartners_skip[j];

      if (distances_observed_p) {
	if ((npartners = nsplicepartners_obs[j]) == 0) {
#if 0
	  triecontents_obs_list = Trie_output_empty_list(&((*trieoffsets_obs)[j]),&nprinted_obs,triecontents_obs_list);
#else
	  (*trie_obs_size) += TRIE_EMPTY_SIZE;
#endif

	} else {
	  if (npartners > MAX_SITES_ALLOCATED) {
	    sites = (Splicestring_T *) CALLOC(npartners,sizeof(Splicestring_T));
	  } else {
	    sites = &(sites_allocated[0]);
	  }
	  nsites = 0;
	  while (nsites < nsplicepartners_obs[j]) {
	    if (splicetypes[j1] == DONOR) {
	      for (p = splicestrings[j1]; p != NULL; p = List_next(p)) {
		sites[nsites++] = (Splicestring_T) List_head(p);
	      }
	    }
	    j1--;
	  }
	  assert(nsites == nsplicepartners_obs[j]);
	  qsort(sites,nsites,sizeof(Splicestring_T),Splicestring_cmp);
#if 0
	  triecontents_obs_list = Trie_output_list(&((*trieoffsets_obs)[j]),&nprinted_obs,triecontents_obs_list,
						  sites,nsites,/*charpos*/0);
#else
	  *trie_obs_size = Trie_output_size(*trie_obs_size,sites,nsites,/*charpos*/0);
#endif
	  if (npartners > MAX_SITES_ALLOCATED) {
	    FREE(sites);
	  }
	}
      }

      if ((npartners = nsplicepartners_max[j]) == 0) {
#if 0
	triecontents_max_list = Trie_output_empty_list(&((*trieoffsets_max)[j]),&nprinted_max,triecontents_max_list);
#else
	(*trie_max_size) += TRIE_EMPTY_SIZE;
#endif

      } else {
	if (npartners > MAX_SITES_ALLOCATED) {
	  sites = (Splicestring_T *) CALLOC(npartners,sizeof(Splicestring_T));
	} else {
	  sites = &(sites_allocated[0]);
	}
	nsites = 0;
	while (nsites < nsplicepartners_max[j]) {
	  if (splicetypes[j1] == DONOR) {
	    for (p = splicestrings[j1]; p != NULL; p = List_next(p)) {
	      sites[nsites++] = (Splicestring_T) List_head(p);
	    }
	  }
	  j1--;
	}
	assert(nsites == nsplicepartners_max[j]);
	qsort(sites,nsites,sizeof(Splicestring_T),Splicestring_cmp);
#if 0
	triecontents_max_list = Trie_output_list(&((*trieoffsets_max)[j]),&nprinted_max,triecontents_max_list,
						sites,nsites,/*charpos*/0);
#else
	*trie_max_size = Trie_output_size(*trie_max_size,sites,nsites,/*charpos*/0);
#endif
	if (npartners > MAX_SITES_ALLOCATED) {
	  FREE(sites);
	}
      }

      break;

    case ANTIDONOR:

      j1 = j - 1 - nsplicepartners_skip[j];

      if (distances_observed_p) {
	if ((npartners = nsplicepartners_obs[j]) == 0) {
#if 0
	  triecontents_obs_list = Trie_output_empty_list(&((*trieoffsets_obs)[j]),&nprinted_obs,triecontents_obs_list);
#else
	  (*trie_obs_size) += TRIE_EMPTY_SIZE;
#endif

	} else {
	  if (npartners > MAX_SITES_ALLOCATED) {
	    sites = (Splicestring_T *) CALLOC(npartners,sizeof(Splicestring_T));
	  } else {
	    sites = &(sites_allocated)[0];
	  }
	  nsites = 0;
	  while (nsites < nsplicepartners_obs[j]) {
	    if (splicetypes[j1] == ANTIACCEPTOR) {
	      for (p = splicestrings[j1]; p != NULL; p = List_next(p)) {
		sites[nsites++] = (Splicestring_T) List_head(p);
	      }
	    }
	    j1--;
	  }
	  assert(nsites == nsplicepartners_obs[j]);
	  qsort(sites,nsites,sizeof(Splicestring_T),Splicestring_cmp);
#if 0
	  triecontents_obs_list = Trie_output_list(&((*trieoffsets_obs)[j]),&nprinted_obs,triecontents_obs_list,
						  sites,nsites,/*charpos*/0);
#else
	  *trie_obs_size = Trie_output_size(*trie_obs_size,sites,nsites,/*charpos*/0);
#endif
	  if (npartners > MAX_SITES_ALLOCATED) {
	    FREE(sites);
	  }
	}
      }

      if ((npartners = nsplicepartners_max[j]) == 0) {
#if 0
	triecontents_max_list = Trie_output_empty_list(&((*trieoffsets_max)[j]),&nprinted_max,triecontents_max_list);
#else
	(*trie_max_size) += TRIE_EMPTY_SIZE;
#endif

      } else {
	if (npartners > MAX_SITES_ALLOCATED) {
	  sites = (Splicestring_T *) CALLOC(npartners,sizeof(Splicestring_T));
	} else {
	  sites = &(sites_allocated[0]);
	}
	nsites = 0;
	while (nsites < nsplicepartners_max[j]) {
	  if (splicetypes[j1] == ANTIACCEPTOR) {
	    for (p = splicestrings[j1]; p != NULL; p = List_next(p)) {
	      sites[nsites++] = (Splicestring_T) List_head(p);
	    }
	  }
	  j1--;
	}
	assert(nsites == nsplicepartners_max[j]);
	qsort(sites,nsites,sizeof(Splicestring_T),Splicestring_cmp);
#if 0
	triecontents_max_list = Trie_output_list(&((*trieoffsets_max)[j]),&nprinted_max,triecontents_max_list,
						sites,nsites,/*charpos*/0);
#else
	*trie_max_size = Trie_output_size(*trie_max_size,sites,nsites,/*charpos*/0);
#endif
	if (npartners > MAX_SITES_ALLOCATED) {
	  FREE(sites);
	}
      }

      break;

    case ANTIACCEPTOR:

      j1 = j + 1 + nsplicepartners_skip[j];

      if (distances_observed_p) {
	if ((npartners = nsplicepartners_obs[j]) == 0) {
#if 0
	  triecontents_obs_list = Trie_output_empty_list(&((*trieoffsets_obs)[j]),&nprinted_obs,triecontents_obs_list);
#else
	  (*trie_obs_size) += TRIE_EMPTY_SIZE;
#endif

	} else {
	  if (npartners > MAX_SITES_ALLOCATED) {
	    sites = (Splicestring_T *) CALLOC(npartners,sizeof(Splicestring_T));
	  } else {
	    sites = &(sites_allocated[0]);
	  }
	  nsites = 0;
	  while (nsites < nsplicepartners_obs[j]) {
	    if (splicetypes[j1] == ANTIDONOR) {
	      for (p = splicestrings[j1]; p != NULL; p = List_next(p)) {
		sites[nsites++] = (Splicestring_T) List_head(p);
	      }
	    }
	    j1++;
	  }
	  assert(nsites == nsplicepartners_obs[j]);
	  qsort(sites,nsites,sizeof(Splicestring_T),Splicestring_cmp);
#if 0
	  triecontents_obs_list = Trie_output_list(&((*trieoffsets_obs)[j]),&nprinted_obs,triecontents_obs_list,
						  sites,nsites,/*charpos*/0);
#else
	  *trie_obs_size = Trie_output_size(*trie_obs_size,sites,nsites,/*charpos*/0);
#endif
	  if (npartners > MAX_SITES_ALLOCATED) {
	    FREE(sites);
	  }
	}
      }

      if ((npartners = nsplicepartners_max[j]) == 0) {
#if 0
	triecontents_max_list = Trie_output_empty_list(&((*trieoffsets_max)[j]),&nprinted_max,triecontents_max_list);
#else
	(*trie_max_size) += TRIE_EMPTY_SIZE;
#endif

      } else {
	if (npartners > MAX_SITES_ALLOCATED) {
	  sites = (Splicestring_T *) CALLOC(npartners,sizeof(Splicestring_T));
	} else {
	  sites = &(sites_allocated[0]);
	}
	nsites = 0;
	while (nsites < nsplicepartners_max[j]) {
	  if (splicetypes[j1] == ANTIDONOR) {
	    for (p = splicestrings[j1]; p != NULL; p = List_next(p)) {
	      sites[nsites++] = (Splicestring_T) List_head(p);
	    }
	  }
	  j1++;
	}
	assert(nsites == nsplicepartners_max[j]);
	qsort(sites,nsites,sizeof(Splicestring_T),Splicestring_cmp);
#if 0
	triecontents_max_list = Trie_output_list(&((*trieoffsets_max)[j]),&nprinted_max,triecontents_max_list,
						sites,nsites,/*charpos*/0);
#else
	*trie_max_size = Trie_output_size(*trie_max_size,sites,nsites,/*charpos*/0);
#endif
	if (npartners > MAX_SITES_ALLOCATED) {
	  FREE(sites);
	}
      }

      break;

    default:
      fprintf(stderr,"Unexpected splicetype %d\n",splicetypes[j]);
      abort();
    }
  }

  return;
}



void
Splicetrie_build_via_splicesites (Triecontent_T **triecontents_obs, Trieoffset_T **trieoffsets_obs,
				  Triecontent_T **triecontents_max, Trieoffset_T **trieoffsets_max,
				  int *nsplicepartners_skip, int *nsplicepartners_obs, int *nsplicepartners_max,
				  Splicetype_T *splicetypes, List_T *splicestrings, int nsplicesites) {
#ifdef USE_LIST
  Uintlist_T triecontents_obs_list = NULL, triecontents_max_list = NULL;
#else
  unsigned int *triecontents_obs_ptr, *triecontents_max_ptr;
#endif
  List_T p;
  int nsites, j, j1;
  Splicestring_T *sites, sites_allocated[MAX_SITES_ALLOCATED];
  int npartners;
  int nprinted_obs = 0, nprinted_max = 0;
  int trie_obs_size = 0, trie_max_size = 0;
  bool distances_observed_p;

  if (nsplicepartners_obs == NULL) {
    distances_observed_p = false;
  } else {
    distances_observed_p = true;
  }

  if (distances_observed_p == true) {
    *trieoffsets_obs = (Trieoffset_T *) CALLOC(nsplicesites,sizeof(Trieoffset_T));
  } else {
    *trieoffsets_obs = (Trieoffset_T *) NULL;
  }
  *trieoffsets_max = (Trieoffset_T *) CALLOC(nsplicesites,sizeof(Trieoffset_T));

  build_via_splicesites_size(&trie_obs_size,&trie_max_size,
			     nsplicepartners_skip,nsplicepartners_obs,nsplicepartners_max,
			     splicetypes,splicestrings,nsplicesites,distances_observed_p);
#ifndef USE_LIST
  if (distances_observed_p) {
    triecontents_obs_ptr = *triecontents_obs = (Triecontent_T *) MALLOC(trie_obs_size * sizeof(Triecontent_T));
  } else {
    *triecontents_obs = (Triecontent_T *) NULL;
  }
  triecontents_max_ptr = *triecontents_max = (Triecontent_T *) MALLOC(trie_max_size * sizeof(Triecontent_T));
#endif


  for (j = 0; j < nsplicesites; j++) {
    switch (splicetypes[j]) {
    case DONOR:
      debug(
	    if (distances_observed_p == true) {
	      printf("donor #%d (%d partners obs, %d partners max):",
		     j,nsplicepartners_obs[j],nsplicepartners_max[j]);
	    } else {
	      printf("donor #%d (%d partners max):",
		     j,nsplicepartners_max[j]);
	    });

      j1 = j + 1 + nsplicepartners_skip[j];

      if (distances_observed_p) {
	if ((npartners = nsplicepartners_obs[j]) == 0) {
#ifdef USE_LIST
	  triecontents_obs_list = Trie_output_empty_list(&((*trieoffsets_obs)[j]),&nprinted_obs,triecontents_obs_list);
#else
	  triecontents_obs_ptr = Trie_output_empty_array(&((*trieoffsets_obs)[j]),&nprinted_obs,triecontents_obs_ptr);
#endif

	} else {
	  if (npartners > MAX_SITES_ALLOCATED) {
	    sites = (Splicestring_T *) CALLOC(npartners,sizeof(Splicestring_T));
	  } else {
	    sites = &(sites_allocated[0]);
	  }
	  nsites = 0;
	  while (nsites < nsplicepartners_obs[j]) {
	    if (splicetypes[j1] == ACCEPTOR) {
	      for (p = splicestrings[j1]; p != NULL; p = List_next(p)) {
		debug(printf(" %d",j1));
		sites[nsites++] = (Splicestring_T) List_head(p);
	      }
	    }
	    j1++;
	  }
	  assert(nsites == nsplicepartners_obs[j]);
	  qsort(sites,nsites,sizeof(Splicestring_T),Splicestring_cmp);
#ifdef USE_LIST
	  triecontents_obs_list = Trie_output_list(&((*trieoffsets_obs)[j]),&nprinted_obs,triecontents_obs_list,
						  sites,nsites,/*charpos*/0);
#else
	  triecontents_obs_ptr = Trie_output_array(&((*trieoffsets_obs)[j]),&nprinted_obs,triecontents_obs_ptr,
						  sites,nsites,/*charpos*/0);
#endif
	  if (npartners > MAX_SITES_ALLOCATED) {
	    FREE(sites);
	  }
	}
      }

      if ((npartners = nsplicepartners_max[j]) == 0) {
#ifdef USE_LIST
	triecontents_max_list = Trie_output_empty_list(&((*trieoffsets_max)[j]),&nprinted_max,triecontents_max_list);
#else
	triecontents_max_ptr = Trie_output_empty_array(&((*trieoffsets_max)[j]),&nprinted_max,triecontents_max_ptr);
#endif

      } else {
	if (npartners > MAX_SITES_ALLOCATED) {
	  sites = (Splicestring_T *) CALLOC(npartners,sizeof(Splicestring_T));
	} else {
	  sites = &(sites_allocated[0]);
	}
	nsites = 0;
	while (nsites < nsplicepartners_max[j]) {
	  if (splicetypes[j1] == ACCEPTOR) {
	    for (p = splicestrings[j1]; p != NULL; p = List_next(p)) {
	      debug(printf(" %d",j1));
	      sites[nsites++] = (Splicestring_T) List_head(p);
	    }
	  }
	  j1++;
	}
	assert(nsites == nsplicepartners_max[j]);
	qsort(sites,nsites,sizeof(Splicestring_T),Splicestring_cmp);
#ifdef USE_LIST
	triecontents_max_list = Trie_output_list(&((*trieoffsets_max)[j]),&nprinted_max,triecontents_max_list,
						sites,nsites,/*charpos*/0);
#else
	triecontents_max_ptr = Trie_output_array(&((*trieoffsets_max)[j]),&nprinted_max,triecontents_max_ptr,
						sites,nsites,/*charpos*/0);
#endif
	if (npartners > MAX_SITES_ALLOCATED) {
	  FREE(sites);
	}
      }

      debug(printf("\n"));
      break;

    case ACCEPTOR:
      debug(
	    if (distances_observed_p == true) {
	      printf("acceptor #%d (%d partners obs, %d partners max):",
		     j,nsplicepartners_obs[j],nsplicepartners_max[j]);
	    } else {
	      printf("acceptor #%d (%d partners max):",
		     j,nsplicepartners_max[j]);
	    });

      j1 = j - 1 - nsplicepartners_skip[j];

      if (distances_observed_p) {
	if ((npartners = nsplicepartners_obs[j]) == 0) {
#ifdef USE_LIST
	  triecontents_obs_list = Trie_output_empty_list(&((*trieoffsets_obs)[j]),&nprinted_obs,triecontents_obs_list);
#else
	  triecontents_obs_ptr = Trie_output_empty_array(&((*trieoffsets_obs)[j]),&nprinted_obs,triecontents_obs_ptr);
#endif

	} else {
	  if (npartners > MAX_SITES_ALLOCATED) {
	    sites = (Splicestring_T *) CALLOC(npartners,sizeof(Splicestring_T));
	  } else {
	    sites = &(sites_allocated[0]);
	  }
	  nsites = 0;
	  while (nsites < nsplicepartners_obs[j]) {
	    if (splicetypes[j1] == DONOR) {
	      for (p = splicestrings[j1]; p != NULL; p = List_next(p)) {
		debug(printf(" %d",j1));
		sites[nsites++] = (Splicestring_T) List_head(p);
	      }
	    }
	    j1--;
	  }
	  assert(nsites == nsplicepartners_obs[j]);
	  qsort(sites,nsites,sizeof(Splicestring_T),Splicestring_cmp);
#ifdef USE_LIST
	  triecontents_obs_list = Trie_output_list(&((*trieoffsets_obs)[j]),&nprinted_obs,triecontents_obs_list,
   	                                           sites,nsites,/*charpos*/0);
#else
	  triecontents_obs_ptr = Trie_output_array(&((*trieoffsets_obs)[j]),&nprinted_obs,triecontents_obs_ptr,
	                                           sites,nsites,/*charpos*/0);
#endif
	  if (npartners > MAX_SITES_ALLOCATED) {
	    FREE(sites);
	  }
	}
      }

      if ((npartners = nsplicepartners_max[j]) == 0) {
#ifdef USE_LIST
	triecontents_max_list = Trie_output_empty_list(&((*trieoffsets_max)[j]),&nprinted_max,triecontents_max_list);
#else
	triecontents_max_ptr = Trie_output_empty_array(&((*trieoffsets_max)[j]),&nprinted_max,triecontents_max_ptr);
#endif

      } else {
	if (npartners > MAX_SITES_ALLOCATED) {
	  sites = (Splicestring_T *) CALLOC(npartners,sizeof(Splicestring_T));
	} else {
	  sites = &(sites_allocated[0]);
	}
	nsites = 0;
	while (nsites < nsplicepartners_max[j]) {
	  if (splicetypes[j1] == DONOR) {
	    for (p = splicestrings[j1]; p != NULL; p = List_next(p)) {
	      debug(printf(" %d",j1));
	      sites[nsites++] = (Splicestring_T) List_head(p);
	    }
	  }
	  j1--;
	}
	assert(nsites == nsplicepartners_max[j]);
	qsort(sites,nsites,sizeof(Splicestring_T),Splicestring_cmp);
#ifdef USE_LIST
	triecontents_max_list = Trie_output_list(&((*trieoffsets_max)[j]),&nprinted_max,triecontents_max_list,
						sites,nsites,/*charpos*/0);
#else
	triecontents_max_ptr = Trie_output_array(&((*trieoffsets_max)[j]),&nprinted_max,triecontents_max_ptr,
						sites,nsites,/*charpos*/0);
#endif
	if (npartners > MAX_SITES_ALLOCATED) {
	  FREE(sites);
	}
      }

      debug(printf("\n"));
      break;

    case ANTIDONOR:
      debug(
	    if (distances_observed_p == true) {
	      printf("antidonor #%d (%d partners obs, %d partners max):",
		     j,nsplicepartners_obs[j],nsplicepartners_max[j]);
	    } else {
	      printf("antidonor #%d (%d partners max):",
		     j,nsplicepartners_max[j]);
	    });

      j1 = j - 1 - nsplicepartners_skip[j];

      if (distances_observed_p) {
	if ((npartners = nsplicepartners_obs[j]) == 0) {
#ifdef USE_LIST
	  triecontents_obs_list = Trie_output_empty_list(&((*trieoffsets_obs)[j]),&nprinted_obs,triecontents_obs_list);
#else
	  triecontents_obs_ptr = Trie_output_empty_array(&((*trieoffsets_obs)[j]),&nprinted_obs,triecontents_obs_ptr);
#endif

	} else {
	  if (npartners > MAX_SITES_ALLOCATED) {
	    sites = (Splicestring_T *) CALLOC(npartners,sizeof(Splicestring_T));
	  } else {
	    sites = &(sites_allocated[0]);
	  }
	  nsites = 0;
	  while (nsites < nsplicepartners_obs[j]) {
	    if (splicetypes[j1] == ANTIACCEPTOR) {
	      for (p = splicestrings[j1]; p != NULL; p = List_next(p)) {
		debug(printf(" %d",j1));
		sites[nsites++] = (Splicestring_T) List_head(p);
	      }
	    }
	    j1--;
	  }
	  assert(nsites == nsplicepartners_obs[j]);
	  qsort(sites,nsites,sizeof(Splicestring_T),Splicestring_cmp);
#ifdef USE_LIST
	  triecontents_obs_list = Trie_output_list(&((*trieoffsets_obs)[j]),&nprinted_obs,triecontents_obs_list,
						  sites,nsites,/*charpos*/0);
#else
	  triecontents_obs_ptr = Trie_output_array(&((*trieoffsets_obs)[j]),&nprinted_obs,triecontents_obs_ptr,
						  sites,nsites,/*charpos*/0);
#endif
	  if (npartners > MAX_SITES_ALLOCATED) {
	    FREE(sites);
	  }
	}
      }

      if ((npartners = nsplicepartners_max[j]) == 0) {
#ifdef USE_LIST
	triecontents_max_list = Trie_output_empty_list(&((*trieoffsets_max)[j]),&nprinted_max,triecontents_max_list);
#else
	triecontents_max_ptr = Trie_output_empty_array(&((*trieoffsets_max)[j]),&nprinted_max,triecontents_max_ptr);
#endif

      } else {
	if (npartners > MAX_SITES_ALLOCATED) {
	  sites = (Splicestring_T *) CALLOC(npartners,sizeof(Splicestring_T));
	} else {
	  sites = &(sites_allocated[0]);
	}
	nsites = 0;
	while (nsites < nsplicepartners_max[j]) {
	  if (splicetypes[j1] == ANTIACCEPTOR) {
	    for (p = splicestrings[j1]; p != NULL; p = List_next(p)) {
	      debug(printf(" %d",j1));
	      sites[nsites++] = (Splicestring_T) List_head(p);
	    }
	  }
	  j1--;
	}
	assert(nsites == nsplicepartners_max[j]);
	qsort(sites,nsites,sizeof(Splicestring_T),Splicestring_cmp);
#ifdef USE_LIST
	triecontents_max_list = Trie_output_list(&((*trieoffsets_max)[j]),&nprinted_max,triecontents_max_list,
						sites,nsites,/*charpos*/0);
#else
	triecontents_max_ptr = Trie_output_array(&((*trieoffsets_max)[j]),&nprinted_max,triecontents_max_ptr,
						sites,nsites,/*charpos*/0);
#endif
	if (npartners > MAX_SITES_ALLOCATED) {
	  FREE(sites);
	}
      }

      debug(printf("\n"));
      break;

    case ANTIACCEPTOR:
      debug(
	    if (distances_observed_p == true) {
	      printf("antiacceptor #%d (%d partners obs, %d partners max):",
		     j,nsplicepartners_obs[j],nsplicepartners_max[j]);
	    } else {
	      printf("antiacceptor #%d (%d partners max):",
		     j,nsplicepartners_max[j]);
	    });

      j1 = j + 1 + nsplicepartners_skip[j];

      if (distances_observed_p) {
	if ((npartners = nsplicepartners_obs[j]) == 0) {
#ifdef USE_LIST
	  triecontents_obs_list = Trie_output_empty_list(&((*trieoffsets_obs)[j]),&nprinted_obs,triecontents_obs_list);
#else
	  triecontents_obs_ptr = Trie_output_empty_array(&((*trieoffsets_obs)[j]),&nprinted_obs,triecontents_obs_ptr);
#endif

	} else {
	  if (npartners > MAX_SITES_ALLOCATED) {
	    sites = (Splicestring_T *) CALLOC(npartners,sizeof(Splicestring_T));
	  } else {
	    sites = &(sites_allocated[0]);
	  }
	  nsites = 0;
	  while (nsites < nsplicepartners_obs[j]) {
	    if (splicetypes[j1] == ANTIDONOR) {
	      for (p = splicestrings[j1]; p != NULL; p = List_next(p)) {
		debug(printf(" %d",j1));
		sites[nsites++] = (Splicestring_T) List_head(p);
	      }
	    }
	    j1++;
	  }
	  assert(nsites == nsplicepartners_obs[j]);
	  qsort(sites,nsites,sizeof(Splicestring_T),Splicestring_cmp);
#ifdef USE_LIST
	  triecontents_obs_list = Trie_output_list(&((*trieoffsets_obs)[j]),&nprinted_obs,triecontents_obs_list,
						  sites,nsites,/*charpos*/0);
#else
	  triecontents_obs_ptr = Trie_output_array(&((*trieoffsets_obs)[j]),&nprinted_obs,triecontents_obs_ptr,
						  sites,nsites,/*charpos*/0);
#endif
	  if (npartners > MAX_SITES_ALLOCATED) {
	    FREE(sites);
	  }
	}
      }

      if ((npartners = nsplicepartners_max[j]) == 0) {
#ifdef USE_LIST
	triecontents_max_list = Trie_output_empty_list(&((*trieoffsets_max)[j]),&nprinted_max,triecontents_max_list);
#else
	triecontents_max_ptr = Trie_output_empty_array(&((*trieoffsets_max)[j]),&nprinted_max,triecontents_max_ptr);
#endif

      } else {
	if (npartners > MAX_SITES_ALLOCATED) {
	  sites = (Splicestring_T *) CALLOC(npartners,sizeof(Splicestring_T));
	} else {
	  sites = &(sites_allocated[0]);
	}
	nsites = 0;
	while (nsites < nsplicepartners_max[j]) {
	  if (splicetypes[j1] == ANTIDONOR) {
	    for (p = splicestrings[j1]; p != NULL; p = List_next(p)) {
	      debug(printf(" %d",j1));
	      sites[nsites++] = (Splicestring_T) List_head(p);
	    }
	  }
	  j1++;
	}
	assert(nsites == nsplicepartners_max[j]);
	qsort(sites,nsites,sizeof(Splicestring_T),Splicestring_cmp);
#ifdef USE_LIST
	triecontents_max_list = Trie_output_list(&((*trieoffsets_max)[j]),&nprinted_max,triecontents_max_list,
						sites,nsites,/*charpos*/0);
#else
	triecontents_max_ptr = Trie_output_array(&((*trieoffsets_max)[j]),&nprinted_max,triecontents_max_ptr,
						sites,nsites,/*charpos*/0);
#endif
	if (npartners > MAX_SITES_ALLOCATED) {
	  FREE(sites);
	}
      }

      debug(printf("\n"));
      break;

    default:
      fprintf(stderr,"Unexpected splicetype %d\n",splicetypes[j]);
      abort();
    }
  }


  assert(nprinted_obs == trie_obs_size);
  assert(nprinted_max == trie_max_size);


  if (distances_observed_p) {
    fprintf(stderr,"splicetrie_obs has %d entries...",nprinted_obs);
#ifdef USE_LIST
    triecontents_obs_list = Uintlist_reverse(triecontents_obs_list);
    *triecontents_obs = Uintlist_to_array(&nprinted_obs,triecontents_obs_list);
    Uintlist_free(&triecontents_obs_list);
#endif
  } else {
    *triecontents_obs = (Triecontent_T *) NULL;
  }


  fprintf(stderr,"splicetrie_max has %d entries...",nprinted_max);
#ifdef USE_LIST
  triecontents_max_list = Uintlist_reverse(triecontents_max_list);
  *triecontents_max = Uintlist_to_array(&nprinted_max,triecontents_max_list);
  Uintlist_free(&triecontents_max_list);
#endif

  return;
}


#ifdef DEBUG
/* Splicetype is for the anchor splice */
static void
dump_sites (Splicestring_T *sites, int nsites, Splicetype_T splicetype) {
  int i;
  char gbuffer[17];

  gbuffer[16] = '\0';

  for (i = 0; i < nsites; i++) {
    printf("%d %u",sites[i]->splicesite_i,sites[i]->splicesite);
    if (splicetype == DONOR || splicetype == ANTIACCEPTOR) {
      splicefrag_nt_rightward(gbuffer,sites[i]->string);
      printf(" %s (rightward)\n",gbuffer);
    } else {
      splicefrag_nt_leftward(gbuffer,sites[i]->string);
      printf(" %s (leftward)\n",gbuffer);
    }
  }
  return;
}
#endif



static void
build_via_introns_size (int *trie_obs_size, Univcoord_T *splicesites, Splicetype_T *splicetypes,
			List_T *splicestrings, int nsplicesites,
			Univ_IIT_T chromosome_iit, IIT_T splicing_iit, int *splicing_divint_crosstable) {
  List_T p;
  int nsites, j, j1;
  Splicestring_T *sites, sites_allocated[MAX_SITES_ALLOCATED];
  
  Univcoord_T chroffset, chrhigh;
  Chrpos_T chrlength;
  Chrpos_T *coords;
  int ncoords, k;
  Chrnum_T chrnum;
  int divno;


  chrhigh = 0U;
  for (j = 0; j < nsplicesites; j++) {
    if (splicesites[j] > chrhigh) {
      chrnum = Univ_IIT_get_one(chromosome_iit,splicesites[j],splicesites[j]);
      Univ_IIT_interval_bounds(&chroffset,&chrhigh,&chrlength,chromosome_iit,chrnum,/*circular_typeint*/-1);
      /* chrhigh += 1U; */

      divno = splicing_divint_crosstable[chrnum];
      assert(divno > 0);
    }

    switch (splicetypes[j]) {
    case DONOR:
      nsites = 0;
      coords = IIT_get_highs_for_low(&ncoords,splicing_iit,divno,splicesites[j]-chroffset);

      j1 = j + 1;
      k = 0;
      while (j1 < nsplicesites && k < ncoords) {
	if (splicetypes[j1] == ACCEPTOR && splicesites[j1] == coords[k] - INTRON_HIGH_TO_LOW + chroffset) {
	  nsites += List_length(splicestrings[j1]);
	  k++;
	}
	j1++;
      }
      /* assertion may not hold if last few introns are invalidated */
      /* assert(k == ncoords); */

      if (nsites == 0) {
#if 0
	(*trieoffsets_obs)[j] = NULL_POINTER;
#endif

      } else {
	if (nsites > MAX_SITES_ALLOCATED) {
	  sites = (Splicestring_T *) CALLOC(nsites,sizeof(Splicestring_T));
	} else {
	  sites = &(sites_allocated[0]);
	}
	nsites = 0;
	j1 = j + 1;
	k = 0;
	while (j1 < nsplicesites && k < ncoords) {
	  if (splicetypes[j1] == ACCEPTOR && splicesites[j1] == coords[k] - INTRON_HIGH_TO_LOW + chroffset) {
	    for (p = splicestrings[j1]; p != NULL; p = List_next(p)) {
	      sites[nsites++] = (Splicestring_T) List_head(p);
	    }
	    k++;
	  }
	  j1++;
	}

	qsort(sites,nsites,sizeof(Splicestring_T),Splicestring_cmp);
#if 0
	triecontents_obs_list = Trie_output_list(&((*trieoffsets_obs)[j]),&nprinted_obs,triecontents_obs_list,
						sites,nsites,/*charpos*/0);
#else
	*trie_obs_size = Trie_output_size(*trie_obs_size,sites,nsites,/*charpos*/0);
#endif
	if (nsites > MAX_SITES_ALLOCATED) {
	  FREE(sites);
	}
      }

      FREE(coords);
      break;

    case ACCEPTOR:
      nsites = 0;
      coords = IIT_get_lows_for_high(&ncoords,splicing_iit,divno,splicesites[j] + INTRON_HIGH_TO_LOW - chroffset);

      j1 = j - 1;
      k = 0;
      while (j1 >= 0 && k < ncoords) {
	if (splicetypes[j1] == DONOR && splicesites[j1] == coords[k] + chroffset) {
	  nsites += List_length(splicestrings[j1]);
	  k++;
	}
	j1--;
      }
      /* assertion may not hold if last few introns are invalidated */
      /* assert(k == ncoords); */

      if (nsites == 0) {
#if 0
	(*trieoffsets_obs)[j] = NULL_POINTER;
#endif

      } else {
	if (nsites > MAX_SITES_ALLOCATED) {
	  sites = (Splicestring_T *) CALLOC(nsites,sizeof(Splicestring_T));
	} else {
	  sites = &(sites_allocated[0]);
	}
	nsites = 0;
	j1 = j - 1;
	k = 0;
	while (j1 >= 0 && k < ncoords) {
	  if (splicetypes[j1] == DONOR && splicesites[j1] == coords[k] + chroffset) {
	    for (p = splicestrings[j1]; p != NULL; p = List_next(p)) {
	      sites[nsites++] = (Splicestring_T) List_head(p);
	    }
	    k++;
	  }
	  j1--;
	}

	qsort(sites,nsites,sizeof(Splicestring_T),Splicestring_cmp);
#if 0
	triecontents_obs_list = Trie_output_list(&((*trieoffsets_obs)[j]),&nprinted_obs,triecontents_obs_list,
						sites,nsites,/*charpos*/0);
#else
	*trie_obs_size = Trie_output_size(*trie_obs_size,sites,nsites,/*charpos*/0);
#endif
	if (nsites > MAX_SITES_ALLOCATED) {
	  FREE(sites);
	}
      }

      FREE(coords);
      break;

    case ANTIDONOR:
      nsites = 0;
      coords = IIT_get_lows_for_high(&ncoords,splicing_iit,divno,splicesites[j] + INTRON_HIGH_TO_LOW - chroffset);

      j1 = j - 1;
      k = 0;
      while (j1 >= 0 && k < ncoords) {
	if (splicetypes[j1] == ANTIACCEPTOR && splicesites[j1] == coords[k] + chroffset) {
	  nsites += List_length(splicestrings[j1]);
	  k++;
	}
	j1--;
      }
      /* assertion may not hold if last few introns are invalidated */
      /* assert(k == ncoords); */

      if (nsites == 0) {
#if 0
	(*trieoffsets_obs)[j] = NULL_POINTER;
#endif

      } else {
	if (nsites > MAX_SITES_ALLOCATED) {
	  sites = (Splicestring_T *) CALLOC(nsites,sizeof(Splicestring_T));
	} else {
	  sites = &(sites_allocated[0]);
	}
	nsites = 0;
	j1 = j - 1;
	k = 0;
	while (j1 >= 0 && k < ncoords) {
	  if (splicetypes[j1] == ANTIACCEPTOR && splicesites[j1] == coords[k] + chroffset) {
	    for (p = splicestrings[j1]; p != NULL; p = List_next(p)) {
	      sites[nsites++] = (Splicestring_T) List_head(p);
	    }
	    k++;
	  }
	  j1--;
	}

	qsort(sites,nsites,sizeof(Splicestring_T),Splicestring_cmp);
#if 0
	triecontents_obs_list = Trie_output_list(&((*trieoffsets_obs)[j]),&nprinted_obs,triecontents_obs_list,
						sites,nsites,/*charpos*/0);
#else
	*trie_obs_size = Trie_output_size(*trie_obs_size,sites,nsites,/*charpos*/0);
#endif
	if (nsites > MAX_SITES_ALLOCATED) {
	  FREE(sites);
	}
      }

      FREE(coords);
      break;

    case ANTIACCEPTOR:
      nsites = 0;
      coords = IIT_get_highs_for_low(&ncoords,splicing_iit,divno,splicesites[j]-chroffset);

      j1 = j + 1;
      k = 0;
      while (j1 < nsplicesites && k < ncoords) {
	if (splicetypes[j1] == ANTIDONOR && splicesites[j1] == coords[k] - INTRON_HIGH_TO_LOW + chroffset) {
	  nsites += List_length(splicestrings[j1]);
	  k++;
	}
	j1++;
      }
      /* assertion may not hold if last few introns are invalidated */
      /* assert(k == ncoords); */

      if (nsites == 0) {
#if 0
	(*trieoffsets_obs)[j] = NULL_POINTER;
#endif

      } else {
	if (nsites > MAX_SITES_ALLOCATED) {
	  sites = (Splicestring_T *) CALLOC(nsites,sizeof(Splicestring_T));
	} else {
	  sites = &(sites_allocated[0]);
	}
	nsites = 0;
	j1 = j + 1;
	k = 0;
	while (j1 < nsplicesites && k < ncoords) {
	  if (splicetypes[j1] == ANTIDONOR && splicesites[j1] == coords[k] - INTRON_HIGH_TO_LOW + chroffset) {
	    for (p = splicestrings[j1]; p != NULL; p = List_next(p)) {
	      sites[nsites++] = (Splicestring_T) List_head(p);
	    }
	    k++;
	  }
	  j1++;
	}

	qsort(sites,nsites,sizeof(Splicestring_T),Splicestring_cmp);
#if 0
	triecontents_obs_list = Trie_output_list(&((*trieoffsets_obs)[j]),&nprinted_obs,triecontents_obs_list,
						sites,nsites,/*charpos*/0);
#else
	*trie_obs_size = Trie_output_size(*trie_obs_size,sites,nsites,/*charpos*/0);
#endif
	if (nsites > MAX_SITES_ALLOCATED) {
	  FREE(sites);
	}
      }

      FREE(coords);
      break;

    default: 
      fprintf(stderr,"Unexpected splicetype %d\n",splicetypes[j]);
      abort();
    }
  }

  return;
}


void
Splicetrie_build_via_introns (Triecontent_T **triecontents_obs, Trieoffset_T **trieoffsets_obs,
			      Univcoord_T *splicesites, Splicetype_T *splicetypes,
			      List_T *splicestrings, int nsplicesites,
			      Univ_IIT_T chromosome_iit, IIT_T splicing_iit, int *splicing_divint_crosstable) {
#ifdef USE_LIST
  Uintlist_T triecontents_obs_list = NULL;
#else
  unsigned int *triecontents_obs_ptr;
#endif
  List_T p;
  int nsites, j, j1;
  Splicestring_T *sites, sites_allocated[MAX_SITES_ALLOCATED];
  int nprinted_obs = 0;
  int trie_obs_size = 0;

  Univcoord_T chroffset, chrhigh;
  Chrpos_T chrlength;
  Chrpos_T *coords;
  int ncoords, k;
  Chrnum_T chrnum;
  int divno;


  *trieoffsets_obs = (Trieoffset_T *) CALLOC(nsplicesites,sizeof(Trieoffset_T));

  build_via_introns_size(&trie_obs_size,splicesites,splicetypes,splicestrings,nsplicesites,
			 chromosome_iit,splicing_iit,splicing_divint_crosstable);
#ifndef USE_LIST
  triecontents_obs_ptr = *triecontents_obs = (Triecontent_T *) MALLOC(trie_obs_size * sizeof(Triecontent_T));
#endif


  chrhigh = 0U;
  for (j = 0; j < nsplicesites; j++) {
    if (splicesites[j] > chrhigh) {
      chrnum = Univ_IIT_get_one(chromosome_iit,splicesites[j],splicesites[j]);
      Univ_IIT_interval_bounds(&chroffset,&chrhigh,&chrlength,chromosome_iit,chrnum,/*circular_typeint*/-1);
      /* chrhigh += 1U; */

      divno = splicing_divint_crosstable[chrnum];
      assert(divno > 0);
    }

    switch (splicetypes[j]) {
    case DONOR:
      nsites = 0;
      coords = IIT_get_highs_for_low(&ncoords,splicing_iit,divno,splicesites[j]-chroffset);

      j1 = j + 1;
      k = 0;
      while (j1 < nsplicesites && k < ncoords) {
	if (splicetypes[j1] == ACCEPTOR && splicesites[j1] == coords[k] - INTRON_HIGH_TO_LOW + chroffset) {
	  nsites += List_length(splicestrings[j1]);
	  k++;
	}
	j1++;
      }
      /* assertion may not hold if last few introns are invalidated */
      /* assert(k == ncoords); */

      if (nsites == 0) {
	(*trieoffsets_obs)[j] = NULL_POINTER;

      } else {
	if (nsites > MAX_SITES_ALLOCATED) {
	  sites = (Splicestring_T *) CALLOC(nsites,sizeof(Splicestring_T));
	} else {
	  sites = &(sites_allocated[0]);
	}
	nsites = 0;
	j1 = j + 1;
	k = 0;
	while (j1 < nsplicesites && k < ncoords) {
	  if (splicetypes[j1] == ACCEPTOR && splicesites[j1] == coords[k] - INTRON_HIGH_TO_LOW + chroffset) {
	    for (p = splicestrings[j1]; p != NULL; p = List_next(p)) {
	      sites[nsites++] = (Splicestring_T) List_head(p);
	    }
	    k++;
	  }
	  j1++;
	}

	qsort(sites,nsites,sizeof(Splicestring_T),Splicestring_cmp);
#ifdef USE_LIST
	triecontents_obs_list = Trie_output_list(&((*trieoffsets_obs)[j]),&nprinted_obs,triecontents_obs_list,
					      sites,nsites,/*charpos*/0);
#else
	triecontents_obs_ptr = Trie_output_array(&((*trieoffsets_obs)[j]),&nprinted_obs,triecontents_obs_ptr,
					      sites,nsites,/*charpos*/0);
#endif
	if (nsites > MAX_SITES_ALLOCATED) {
	  FREE(sites);
	}
      }

      FREE(coords);
      break;

    case ACCEPTOR:
      nsites = 0;
      coords = IIT_get_lows_for_high(&ncoords,splicing_iit,divno,splicesites[j] + INTRON_HIGH_TO_LOW - chroffset);

      j1 = j - 1;
      k = 0;
      while (j1 >= 0 && k < ncoords) {
	if (splicetypes[j1] == DONOR && splicesites[j1] == coords[k] + chroffset) {
	  nsites += List_length(splicestrings[j1]);
	  k++;
	}
	j1--;
      }
      /* assertion may not hold if last few introns are invalidated */
      /* assert(k == ncoords); */

      if (nsites == 0) {
	(*trieoffsets_obs)[j] = NULL_POINTER;

      } else {
	if (nsites > MAX_SITES_ALLOCATED) {
	  sites = (Splicestring_T *) CALLOC(nsites,sizeof(Splicestring_T));
	} else {
	  sites = &(sites_allocated[0]);
	}
	nsites = 0;
	j1 = j - 1;
	k = 0;
	while (j1 >= 0 && k < ncoords) {
	  if (splicetypes[j1] == DONOR && splicesites[j1] == coords[k] + chroffset) {
	    for (p = splicestrings[j1]; p != NULL; p = List_next(p)) {
	      sites[nsites++] = (Splicestring_T) List_head(p);
	    }
	    k++;
	  }
	  j1--;
	}

	qsort(sites,nsites,sizeof(Splicestring_T),Splicestring_cmp);
#ifdef USE_LIST
	triecontents_obs_list = Trie_output_list(&((*trieoffsets_obs)[j]),&nprinted_obs,triecontents_obs_list,
						 sites,nsites,/*charpos*/0);
#else
	triecontents_obs_ptr = Trie_output_array(&((*trieoffsets_obs)[j]),&nprinted_obs,triecontents_obs_ptr,
						 sites,nsites,/*charpos*/0);
#endif
	if (nsites > MAX_SITES_ALLOCATED) {
	  FREE(sites);
	}
      }

      FREE(coords);
      break;

    case ANTIDONOR:
      nsites = 0;
      coords = IIT_get_lows_for_high(&ncoords,splicing_iit,divno,splicesites[j] + INTRON_HIGH_TO_LOW - chroffset);

      j1 = j - 1;
      k = 0;
      while (j1 >= 0 && k < ncoords) {
	if (splicetypes[j1] == ANTIACCEPTOR && splicesites[j1] == coords[k] + chroffset) {
	  nsites += List_length(splicestrings[j1]);
	  k++;
	}
	j1--;
      }
      /* assertion may not hold if last few introns are invalidated */
      /* assert(k == ncoords); */

      if (nsites == 0) {
	(*trieoffsets_obs)[j] = NULL_POINTER;

      } else {
	if (nsites > MAX_SITES_ALLOCATED) {
	  sites = (Splicestring_T *) CALLOC(nsites,sizeof(Splicestring_T));
	} else {
	  sites = &(sites_allocated[0]);
	}
	nsites = 0;
	j1 = j - 1;
	k = 0;
	while (j1 >= 0 && k < ncoords) {
	  if (splicetypes[j1] == ANTIACCEPTOR && splicesites[j1] == coords[k] + chroffset) {
	    for (p = splicestrings[j1]; p != NULL; p = List_next(p)) {
	      sites[nsites++] = (Splicestring_T) List_head(p);
	    }
	    k++;
	  }
	  j1--;
	}

	qsort(sites,nsites,sizeof(Splicestring_T),Splicestring_cmp);
#ifdef USE_LIST
	triecontents_obs_list = Trie_output_list(&((*trieoffsets_obs)[j]),&nprinted_obs,triecontents_obs_list,
						sites,nsites,/*charpos*/0);
#else
	triecontents_obs_ptr = Trie_output_array(&((*trieoffsets_obs)[j]),&nprinted_obs,triecontents_obs_ptr,
                                        	 sites,nsites,/*charpos*/0);
#endif
	if (nsites > MAX_SITES_ALLOCATED) {
	  FREE(sites);
	}
      }

      FREE(coords);
      break;

    case ANTIACCEPTOR:
      nsites = 0;
      coords = IIT_get_highs_for_low(&ncoords,splicing_iit,divno,splicesites[j]-chroffset);

      j1 = j + 1;
      k = 0;
      while (j1 < nsplicesites && k < ncoords) {
	if (splicetypes[j1] == ANTIDONOR && splicesites[j1] == coords[k] - INTRON_HIGH_TO_LOW + chroffset) {
	  nsites += List_length(splicestrings[j1]);
	  k++;
	}
	j1++;
      }
      /* assertion may not hold if last few introns are invalidated */
      /* assert(k == ncoords); */

      if (nsites == 0) {
	(*trieoffsets_obs)[j] = NULL_POINTER;

      } else {
	if (nsites > MAX_SITES_ALLOCATED) {
	  sites = (Splicestring_T *) CALLOC(nsites,sizeof(Splicestring_T));
	} else {
	  sites = &(sites_allocated[0]);
	}
	nsites = 0;
	j1 = j + 1;
	k = 0;
	while (j1 < nsplicesites && k < ncoords) {
	  if (splicetypes[j1] == ANTIDONOR && splicesites[j1] == coords[k] - INTRON_HIGH_TO_LOW + chroffset) {
	    for (p = splicestrings[j1]; p != NULL; p = List_next(p)) {
	      sites[nsites++] = (Splicestring_T) List_head(p);
	    }
	    k++;
	  }
	  j1++;
	}

	qsort(sites,nsites,sizeof(Splicestring_T),Splicestring_cmp);
#ifdef USE_LIST
	triecontents_obs_list = Trie_output_list(&((*trieoffsets_obs)[j]),&nprinted_obs,triecontents_obs_list,
						 sites,nsites,/*charpos*/0);
#else
	triecontents_obs_ptr = Trie_output_array(&((*trieoffsets_obs)[j]),&nprinted_obs,triecontents_obs_ptr,
						 sites,nsites,/*charpos*/0);
#endif
	if (nsites > MAX_SITES_ALLOCATED) {
	  FREE(sites);
	}
      }

      FREE(coords);
      break;

    default: 
      fprintf(stderr,"Unexpected splicetype %d\n",splicetypes[j]);
      abort();
    }
  }

  assert(nprinted_obs == trie_obs_size);

  fprintf(stderr,"splicetrie_obs has %d entries...",nprinted_obs);
#ifdef USE_LIST
  triecontents_obs_list = Uintlist_reverse(triecontents_obs_list);
  *triecontents_obs = Uintlist_to_array(&nprinted_obs,triecontents_obs_list);
  Uintlist_free(&triecontents_obs_list);
#endif

  return;
}

