static char rcsid[] = "$Id: snpindex.c 168395 2015-06-26 17:13:13Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#ifndef HAVE_MEMCPY
# define memcpy(d,s,n) bcopy((s),(d),(n))
#endif
#ifndef HAVE_MEMMOVE
# define memmove(d,s,n) bcopy((s),(d),(n))
#endif

#ifdef WORDS_BIGENDIAN
#include "bigendian.h"
#else
#include "littleendian.h"
#endif

#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>		/* For memset */
#include <ctype.h>		/* For toupper */
#include <sys/mman.h>		/* For munmap */
#include <math.h>		/* For qsort */
#ifdef HAVE_UNISTD_H
#include <unistd.h>		/* For lseek and close */
#endif
#ifdef HAVE_SYS_TYPES_H
#include <sys/types.h>		/* For off_t */
#endif
#if HAVE_DIRENT_H
# include <dirent.h>
# define NAMLEN(dirent) strlen((dirent)->d_name)
#else
# define dirent direct
# define NAMLEN(dirent) (dirent)->d_namlen
# if HAVE_SYS_NDIR_H
#  include <sys/ndir.h>
# endif
# if HAVE_SYS_DIR_H
#  include <sys/dir.h>
# endif
# if HAVE_NDIR_H
#  include <ndir.h>
# endif
#endif

#include "mem.h"
#include "fopen.h"
#include "access.h"
#include "types.h"

#include "compress-write.h"
#include "compress.h"
#include "interval.h"
#include "complement.h"

#include "bool.h"
#include "genomicpos.h"
#include "iitdef.h"
#include "uintlist.h"
#include "chrnum.h"
#include "genome.h"
#include "datadir.h"
#include "iit-read-univ.h"
#include "iit-read.h"
#include "indexdb.h"
#include "indexdb-write.h"
#include "bitpack64-write.h"
#include "getopt.h"

#define POSITIONS8_HIGH_SHIFT 32
#define POSITIONS8_LOW_MASK 0xFFFFFFFF

#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif

/* SNP blocks and writing of positions */
#ifdef DEBUG1
#define debug1(x) x
#else
#define debug1(x)
#endif

#ifdef CHECK
#define check(x) x
#else
#define check(x)
#endif



static char *user_sourcedir = NULL;
static char *user_destdir = NULL;
static char *dbroot = NULL;
static char *fileroot = NULL;
static int index1part = 15;
static int required_index1part = 0;
static int index1interval = 3;
static int required_interval = 0;

static char *dbversion = NULL;
static int circular_typeint = -1;

static char *snps_root = NULL;
static bool show_warnings_p = true;
static int max_warnings = -1;


static struct option long_options[] = {
  /* Input options */
  {"sourcedir", required_argument, 0, 'D'},	/* user_sourcedir */
  {"db", required_argument, 0, 'd'}, /* dbroot */
  {"kmer", required_argument, 0, 'k'},	   /* required_index1part */
  {"sampling", required_argument, 0, 'q'},  /* required_interval */
  {"destdir", required_argument, 0, 'V'},	/* user_destdir */
  {"snpsdb", required_argument, 0, 'v'}, /* snps_root */

  /* Output options */
  {"max-warnings", required_argument, 0, 'w'}, /* max_warnings */

  /* Help options */
  {"version", no_argument, 0, 0}, /* print_program_version */
  {"help", no_argument, 0, 0}, /* print_program_usage */
  {0, 0, 0, 0}
};


static void
print_program_version () {
  fprintf(stdout,"\n");
  fprintf(stdout,"SNPINDEX: Builds GMAP index files for known SNPs\n");
  fprintf(stdout,"Part of GMAP package, version %s\n",PACKAGE_VERSION);
  fprintf(stdout,"Default gmap directory: %s\n",GMAPDB);
  fprintf(stdout,"Thomas D. Wu, Genentech, Inc.\n");
  fprintf(stdout,"Contact: twu@gene.com\n");
  fprintf(stdout,"\n");
  return;
}

static void
print_program_usage ();



static Oligospace_T
power (int base, int exponent) {
  Oligospace_T result = 1UL;
  int i;

  for (i = 0; i < exponent; i++) {
    result *= base;
  }
  return result;
}


/*                87654321 */
#define RIGHT_A 0x00000000
#define RIGHT_C 0x00000001
#define RIGHT_G 0x00000002
#define RIGHT_T 0x00000003

/*                      87654321 */
#define LOW_TWO_BITS  0x00000003

static char *
shortoligo_nt (Storedoligomer_T oligo, int oligosize) {
  char *nt;
  int i, j;
  Storedoligomer_T lowbits;

  nt = (char *) CALLOC(oligosize+1,sizeof(char));
  j = oligosize-1;
  for (i = 0; i < oligosize; i++) {
    lowbits = oligo & LOW_TWO_BITS;
    switch (lowbits) {
    case RIGHT_A: nt[j] = 'A'; break;
    case RIGHT_C: nt[j] = 'C'; break;
    case RIGHT_G: nt[j] = 'G'; break;
    case RIGHT_T: nt[j] = 'T'; break;
    }
    oligo >>= 2;
    j--;
  }

  return nt;
}

static char
check_acgt (char nt, IIT_T snps_iit, int divno, char *divstring,
	    int origindex, Interval_T interval) {
  char *label;
  bool allocp;

  if (nt == 'N') {
    return nt;
  } else if (nt != 'A' && nt != 'C' && nt != 'G' && nt != 'T') {
    label = IIT_label(snps_iit,origindex,&allocp);
    fprintf(stderr,"\nFor %s at %s:%u, alternate allele %c is not ACGT, so using 'N' as alternate allele",
	    label,divstring,Interval_low(interval),nt);
    if (allocp) {
      FREE(label);
    }
    return 'N';
  } else {
    return nt;
  }
}


static Storedoligomer_T
revise_oligomer (Storedoligomer_T oligo, char nt) {
  switch (nt) {
  case 'A': return (oligo << 2);
  case 'C': return (oligo << 2) | 1U;
  case 'G': return (oligo << 2) | 2U;
  case 'T': return (oligo << 2) | 3U;
  }
  fprintf(stderr,"altstring char is %c\n",nt);
  exit(9);
  return 0;
}


#if 0
static bool
samep (char *string1, char *string2, int length) {
  int i;

  for (i = 0; i < length; i++) {
    if (string1[i] != string2[i]) {
      return false;
    }
  }
  return true;
}
#endif


typedef struct Labeled_interval_T *Labeled_interval_T;
struct Labeled_interval_T {
  int origindex;
  Interval_T interval;
};

static Labeled_interval_T
Labeled_interval_new (int origindex, Interval_T interval) {
  Labeled_interval_T new = (Labeled_interval_T) MALLOC(sizeof(*new));

  new->origindex = origindex;
  new->interval = interval;
  return new;
}

static void
Labeled_interval_free (Labeled_interval_T *old, bool free_interval_p) {
  if (free_interval_p == true) {
    Interval_free(&((*old)->interval));
  }
  FREE(*old);
  return;
}

static int
Labeled_interval_cmp (const void *a, const void *b) {
  Labeled_interval_T x = * (Labeled_interval_T *) a;
  Labeled_interval_T y = * (Labeled_interval_T *) b;

  return Interval_cmp((const void *) &(x->interval),(const void *) &(y->interval));
}



static int
process_snp_block (int *nwarnings, Positionsptr_T *offsets, UINT4 *positions4, UINT8 *positions8,
		   Labeled_interval_T *intervals, int nintervals,
		   Univcoord_T chroffset, Genome_T genome,
		   Genomecomp_T *snp_blocks, int divno, char *divstring, int intervali,
		   IIT_T snps_iit, Univ_IIT_T chromosome_iit, int index1part,
		   bool coord_values_8p) {
  int nerrors = 0;
  bool *snpp;
  char *refstring;
  char *altstring;
#ifdef DEBUG1
  char *nt;
#endif

  int index, length;
  int nsnps, stringi, starti, shift, i, k;
  char *snptype, *label, refnt, altnt;
  Univcoord_T ptr;
  Univcoord_T snpposition, startposition, endposition, first_snppos, last_snppos, position;
  Chrpos_T chrpos;
  Chrnum_T chrnum;
  Interval_T interval;
  int nunknowns;
  bool badcharp, allocp;
  
  Uintlist_T oligomers, newoligomers, p;
  Storedoligomer_T oligo;
#ifdef WORDS_BIGENDIAN
  Genomecomp_T high, low, flags;
#endif


  /* Subtract 1 because snps_iit is 1-based */
  first_snppos = Interval_low(intervals[0]->interval) - 1U;
  last_snppos = Interval_low(intervals[nintervals-1]->interval) - 1U;
  debug1(
	 if (nintervals == 1) {
	   printf("Processing snp at chrpos %s:%u\n",divstring,first_snppos+1U);
	 } else {
	   printf("Processing block of %d snps from chrpos %s:%u to %u\n",
		  nintervals,divstring,first_snppos+1U,last_snppos+1U);
	 }
	 );

  if (first_snppos < (index1part - 1)) {
    startposition = chroffset;
  } else {
    startposition = chroffset + first_snppos - (index1part - 1);
  }
  endposition = chroffset + last_snppos + (index1part - 1);
  length = endposition - startposition + 1;

  snpp = (bool *) CALLOC(length,sizeof(bool));
  refstring = (char *) CALLOC(length + 1,sizeof(char));
  altstring = (char *) CALLOC(length + 1,sizeof(char));

#if 0
  /* Doesn't work with circular chromosomes */
  Genome_fill_buffer(&chrnum,&nunknowns,genome,/*left*/startposition,length,
		     refstring,chromosome_iit);
#else
  Genome_fill_buffer_simple(genome,/*left*/startposition,length,refstring);
#endif

  for (i = 0; i < nintervals; i++) {
    interval = intervals[i]->interval;
    snpposition = chroffset + Interval_low(interval) - 1U; /* Subtract 1 because snps_iit is 1-based */

    debug1(label = IIT_label(snps_iit,intervals[i]->origindex,&allocp));
    debug1(printf("  Neighbor %s at %s:%u\n",label,divstring,Interval_low(interval)));

    stringi = snpposition - startposition;

    snptype = IIT_typestring(snps_iit,Interval_type(interval));
    if (strlen(snptype) != 2) {
      fprintf(stderr,"Unrecognized snptype %s\n",snptype);
      abort();
    } else {
      refnt = refstring[stringi];
      if (refnt == snptype[0]) {
	if (altstring[stringi] != '\0' && altstring[stringi] != snptype[1]) {
	  nerrors++;
	  if (show_warnings_p == true) {
	    label = IIT_label(snps_iit,intervals[i]->origindex,&allocp);
	    fprintf(stderr,"\nFor %s at %s:%u, saw two different alternate alleles %c and %c, so using N as alternate allele.",
		    label,divstring,Interval_low(interval),altstring[stringi],snptype[1]);
	    if (allocp == true) {
	      FREE(label);
	    }
	    if (++(*nwarnings) == max_warnings) {
	      fprintf(stderr,"\nMaximum of %d warnings reached.  No more warnings will be shown\n",max_warnings);
	      show_warnings_p = false;
	    }
	  }
	  altnt = 'N';
	} else {
	  altnt = check_acgt(snptype[1],snps_iit,divno,divstring,intervals[i]->origindex,intervals[i]->interval);
	}
      } else if (refnt == snptype[1]) {
	if (altstring[stringi] != '\0' && altstring[stringi] != snptype[0]) {
	  nerrors++;
	  if (show_warnings_p == true) {
	    label = IIT_label(snps_iit,intervals[i]->origindex,&allocp);
	    fprintf(stderr,"\nFor %s at %s:%u, saw two different alternate alleles %c and %c, so using N as alternate allele.",
		    label,divstring,Interval_low(interval),altstring[stringi],snptype[0]);
	    if (allocp == true) {
	      FREE(label);
	    }
	    if (++(*nwarnings) == max_warnings) {
	      fprintf(stderr,"\nMaximum of %d warnings reached.  No more warnings will be shown\n",max_warnings);
	      show_warnings_p = false;
	    }
	  }
	  altnt = 'N';
	} else {
	  altnt = check_acgt(snptype[0],snps_iit,divno,divstring,intervals[i]->origindex,intervals[i]->interval);
	}
      } else {
	nerrors++;
	if (show_warnings_p == true) {
	  label = IIT_label(snps_iit,intervals[i]->origindex,&allocp);
	  fprintf(stderr,"\nFor %s at %s:%u, snptype %s not consistent with reference allele %c, so ignoring.",
		  label,divstring,Interval_low(interval),snptype,refstring[stringi]);
	  if (allocp == true) {
	    FREE(label);
	  }
	  if (++(*nwarnings) == max_warnings) {
	    fprintf(stderr,"\nMaximum of %d warnings reached.  No more warnings will be shown\n",max_warnings);
	    show_warnings_p = false;
	  }
	}
	altnt = '\0';		/* Ignoring */
      }

      if (altnt == '\0') {
	/* Skip */
      } else if (altnt == refnt) {
	fprintf(stderr,"\nAt %s:%u, alternate allele %c is same as reference allele\n",
		divstring,Interval_low(interval),altnt);
      } else {
	altstring[stringi] = altnt;
	snpp[stringi] = true;
	if (snp_blocks != NULL) {
	  /* revising genome */
	  ptr = snpposition/32U*3;
	  shift = snpposition % 32U;

#ifdef WORDS_BIGENDIAN
	  flags = Bigendian_convert_uint(snp_blocks[ptr+2]);
	  flags |= (1 << shift);
	  snp_blocks[ptr+2] = Bigendian_convert_uint(flags);
#else
	  snp_blocks[ptr+2] |= (1 << shift);	/* Flags.  Change even for 'N'. */
#endif
	  if (altnt == 'N') {
	    /* refnt + flag indicates 'N' */
	    altnt = refnt;	/* Change back to refnt, if necessary. */
	  }

	  if (shift >= 16) {
	    /* high */
	    shift -= 16;
#ifdef WORDS_BIGENDIAN
	    high = Bigendian_convert_uint(snp_blocks[ptr]);
	    high &= ~(0x3 << (2*shift)); /* clear bits */
	    switch (altnt) {	/* set bits */
	    case 'A': break;
	    case 'C': high |= (0x1 << (2*shift)); break;
	    case 'G': high |= (0x2 << (2*shift)); break;
	    case 'T': high |= (0x3 << (2*shift)); break;
	    default: abort();
	    }
	    snp_blocks[ptr] = Bigendian_convert_uint(high);
#else
	    snp_blocks[ptr] &= ~(0x3 << (2*shift)); /* clear bits */
	    switch (altnt) {	/* set bits */
	    case 'A': break;
	    case 'C': snp_blocks[ptr] |= (0x1 << (2*shift)); break;
	    case 'G': snp_blocks[ptr] |= (0x2 << (2*shift)); break;
	    case 'T': snp_blocks[ptr] |= (0x3 << (2*shift)); break;
	    default: abort();
	    }
#endif
	  } else {
	    /* low */
#ifdef WORDS_BIGENDIAN
	    low = Bigendian_convert_uint(snp_blocks[ptr+1]);
	    low &= ~(0x3 << (2*shift)); /* clear bits */
	    switch (altnt) {	/* set bits */
	    case 'A': break;
	    case 'C': low |= (0x1 << (2*shift)); break;
	    case 'G': low |= (0x2 << (2*shift)); break;
	    case 'T': low |= (0x3 << (2*shift)); break;
	    default: abort();
	    }
	    snp_blocks[ptr+1] = Bigendian_convert_uint(low);
#else
	    snp_blocks[ptr+1] &= ~(0x3 << (2*shift)); /* clear bits */
	    switch (altnt) {	/* set bits */
	    case 'A': break;
	    case 'C': snp_blocks[ptr+1] |= (0x1 << (2*shift)); break;
	    case 'G': snp_blocks[ptr+1] |= (0x2 << (2*shift)); break;
	    case 'T': snp_blocks[ptr+1] |= (0x3 << (2*shift)); break;
	    default: abort();
	    }
#endif
	  }
	}
      }
    }
  }
  
  for (starti = 0, position = startposition, chrpos = startposition - chroffset; 
       position <= endposition-index1part+1U; starti++, position++, chrpos++) {
    if (chrpos % index1interval == 0) {
      /* chrpos % 3 == 0 is same as the condition in indexdb.c for storing a position */
      nsnps = 0;
      badcharp = false;
      for (k = starti; k < starti + index1part; k++) {
	if (snpp[k] == true) {
	  nsnps++;
	}
	if (refstring[k] != 'A' && refstring[k] != 'C' && refstring[k] != 'G' && refstring[k] != 'T') {
	  badcharp = true;
	}
      }
      if (nsnps == 0) {
	/* no snps */
	/* fprintf(stderr,"\nNo snps at position %llu, %s:%u",(unsigned long long) position,divstring,Interval_low(interval)); */
      } else if (nsnps > 4) {
	/* too many snps */
      } else if (badcharp == true) {
	/* bad reference char */
      } else {
	oligomers = Uintlist_push(NULL,0U);
	for (k = starti, i = 0; k < starti + index1part; k++, i++) {
	  newoligomers = NULL;
	  if (snpp[k] == false) {
	    for (p = oligomers; p != NULL; p = Uintlist_next(p)) {
	      newoligomers = Uintlist_push(newoligomers,revise_oligomer(Uintlist_head(p),refstring[k]));
	    }
	    
	  } else if (altstring[k] != 'N') {
	    for (p = oligomers; p != NULL; p = Uintlist_next(p)) {
	      newoligomers = Uintlist_push(newoligomers,revise_oligomer(Uintlist_head(p),refstring[k]));
	      newoligomers = Uintlist_push(newoligomers,revise_oligomer(Uintlist_head(p),altstring[k]));
	    }

	  } else {
	    for (p = oligomers; p != NULL; p = Uintlist_next(p)) {
	      newoligomers = Uintlist_push(newoligomers,revise_oligomer(Uintlist_head(p),refstring[k]));
	      if (refstring[k] != 'A') {
		newoligomers = Uintlist_push(newoligomers,revise_oligomer(Uintlist_head(p),'A'));
	      }
	      if (refstring[k] != 'C') {
		newoligomers = Uintlist_push(newoligomers,revise_oligomer(Uintlist_head(p),'C'));
	      }
	      if (refstring[k] != 'G') {
		newoligomers = Uintlist_push(newoligomers,revise_oligomer(Uintlist_head(p),'G'));
	      }
	      if (refstring[k] != 'T') {
		newoligomers = Uintlist_push(newoligomers,revise_oligomer(Uintlist_head(p),'T'));
	      }
	    }
	  }

	  Uintlist_free(&oligomers);
	  oligomers = Uintlist_reverse(newoligomers);
	}

#ifdef CHECK
	for (p = Uintlist_next(oligomers); p != NULL; p = Uintlist_next(p)) {
	  oligo = Uintlist_head(p);
	  nt = shortoligo_nt(oligo,index1part);
	  if (samep(nt,&(refstring[starti]),index1part) == true) {
	    fprintf(stderr,"Storing oligomer %s that is the same as the reference at %llu (%s:%u)\n",
		    nt,(unsigned long long) position,divstring,chrpos+1U);
	    abort();
	  }
	  FREE(nt);
	}
#endif	

	/* Ignore the first element in oligomers, which is all reference */
	if (positions4 == NULL && positions8 == NULL) {
	  /* writing offsets */
	  for (p = Uintlist_next(oligomers); p != NULL; p = Uintlist_next(p)) {
	    oligo = Uintlist_head(p);
	    offsets[oligo + 1U] += 1;
	    debug1(nt = shortoligo_nt(oligo,index1part);
		   printf("Storing %s at %llu (%s:%u)\n",nt,(unsigned long long) position,divstring,chrpos+1U);
		   FREE(nt));
	  }

	} else {
	  /* writing positions */
	  if (coord_values_8p == true) {
	    for (p = Uintlist_next(oligomers); p != NULL; p = Uintlist_next(p)) {
	      oligo = Uintlist_head(p);
	      positions8[offsets[oligo]++] = position;
	      debug1(nt = shortoligo_nt(oligo,index1part);
		     printf("Storing %s at %llu (%s:%u)\n",nt,(unsigned long long) position,divstring,chrpos+1U);
		     FREE(nt));
	    }
	    
	  } else {
	    for (p = Uintlist_next(oligomers); p != NULL; p = Uintlist_next(p)) {
	      oligo = Uintlist_head(p);
	      positions4[offsets[oligo]++] = (UINT4) position;
	      debug1(nt = shortoligo_nt(oligo,index1part);
		     printf("Storing %s at %llu (%s:%u)\n",nt,(unsigned long long) position,divstring,chrpos+1U);
		     FREE(nt));
	    }
	  }
	}
	Uintlist_free(&oligomers);
      }
    }
  }

  FREE(altstring);
  FREE(refstring);
  FREE(snpp);

  return nerrors;
}


static Positionsptr_T *
compute_offsets (IIT_T snps_iit, Univ_IIT_T chromosome_iit, Genome_T genome, Genomecomp_T *snp_blocks,
		 Oligospace_T oligospace, int index1part) {
  Positionsptr_T *offsets;

  Labeled_interval_T *intervals;
  int origindex;
  Interval_T interval, copy;
  int nintervals, nintervals_alias, nerrors, divno, i, j;
  Oligospace_T oligoi;
  char *divstring;
  Chrnum_T chrnum;
  Univcoord_T chroffset;
  Chrpos_T chrlength;
  int nwarnings = 0;

  offsets = (Positionsptr_T *) CALLOC(oligospace+1,sizeof(Positionsptr_T));

  for (divno = 1; divno < snps_iit->ndivs; divno++) {
    divstring = IIT_divstring(snps_iit,divno);
    fprintf(stderr,"Processing offsets for chromosome %s...",divstring);
    if ((chrnum = Univ_IIT_find_one(chromosome_iit,divstring)) <= 0) {
      fprintf(stderr,"not found in chromosome iit\n");
    } else {
      nerrors = 0;
      fprintf(stderr,"has %d snps...",IIT_nintervals(snps_iit,divno));
      chroffset = Univ_IIT_interval_low(chromosome_iit,chrnum);
      chrlength = Univ_IIT_interval_length(chromosome_iit,chrnum);
      nintervals = IIT_nintervals(snps_iit,divno);

      if (Univ_IIT_interval_type(chromosome_iit,chrnum) == circular_typeint) {
	fprintf(stderr,"and is circular...");
	nintervals_alias = 2*nintervals;
	intervals = (Labeled_interval_T *) CALLOC(nintervals_alias,sizeof(Labeled_interval_T));
	for (i = 0; i < nintervals; i++) {
	  origindex = IIT_index(snps_iit,divno,i);
	  intervals[i] = Labeled_interval_new(origindex,/*interval*/&(snps_iit->intervals[divno][i]));
	}
	for (j = 0; i < nintervals_alias; i++, j++) {
	  origindex = IIT_index(snps_iit,divno,j);
	  interval = intervals[j]->interval;
	  copy = Interval_new(/*low*/Interval_low(interval) + chrlength,
			      /*high*/Interval_high(interval) + chrlength,
			      Interval_type(interval));
	  intervals[i] = Labeled_interval_new(origindex,/*interval*/copy);
	}

      } else {
	nintervals_alias = nintervals;
	intervals = (Labeled_interval_T *) CALLOC(nintervals,sizeof(Labeled_interval_T));
	for (i = 0; i < nintervals; i++) {
	  origindex = IIT_index(snps_iit,divno,i);
	  intervals[i] = Labeled_interval_new(origindex,/*interval*/&(snps_iit->intervals[divno][i]));
	}
      }

      qsort(intervals,nintervals,sizeof(Labeled_interval_T),Labeled_interval_cmp);

      i = 0;
      while (i < nintervals_alias) {
	j = i + 1;
	while (j < nintervals_alias && Interval_low(intervals[j]->interval) < Interval_low(intervals[j-1]->interval) + index1part) {
	  j++;
	}
	nerrors += process_snp_block(&nwarnings,offsets,/*positions4*/NULL,/*positions8*/NULL,
				     &(intervals[i]),/*nintervals*/j-i,
				     chroffset,genome,snp_blocks,
				     divno,divstring,/*intervali*/i,snps_iit,chromosome_iit,index1part,
				     /*coord_values_8p (irrelevant)*/false);
	i = j;
      }

      for (i = nintervals; i < nintervals_alias; i++) {
	Labeled_interval_free(&(intervals[i]),/*free_interval_p*/true);
      }
      for (i = 0; i < nintervals; i++) {
	Labeled_interval_free(&(intervals[i]),/*free_interval_p*/false);
      }

      FREE(intervals);
      fprintf(stderr,"done (%d snps inconsistent with reference genome)\n",nerrors);
    }
  }

  /* Do not add extra for sentinel */

  for (oligoi = 1; oligoi <= oligospace; oligoi++) {
    offsets[oligoi] = offsets[oligoi] + offsets[oligoi-1];
    debug(if (offsets[oligoi] != offsets[oligoi-1]) {
	    printf("Offset for %lX: %u\n",oligoi,offsets[oligoi]);
	  });
  }

  return offsets;
}


static void *
compute_positions (Positionsptr_T *offsets, IIT_T snps_iit, Univ_IIT_T chromosome_iit,
		   Genome_T genome, Oligospace_T oligospace, bool coord_values_8p) {
  UINT4 *positions4 = NULL;
  UINT8 *positions8 = NULL;

  Labeled_interval_T *intervals;
  int origindex;
  Interval_T interval, copy;
  int nintervals, nintervals_alias, divno, i, j;
  Oligospace_T oligoi;
  char *divstring;
  Chrnum_T chrnum;
  Univcoord_T chroffset;
  Chrpos_T chrlength;
  Positionsptr_T *pointers, totalcounts, block_start, block_end, npositions;
  int nwarnings = 0;

  totalcounts = offsets[oligospace];
  if (totalcounts == 0) {
    fprintf(stderr,"Something is wrong with the offsets.  Total counts is zero.\n");
    fprintf(stderr,"Do the chromosomes in the IIT file match those in the genome?\n");
    fprintf(stderr,"Here are known chromosomes in the genome: ");
    Univ_IIT_dump_labels(stderr,chromosome_iit);
    fprintf(stderr,"Here are chromosomes in the SNPs IIT file: ");
    IIT_dump_divstrings(stderr,snps_iit);
    exit(9);
  } else if (coord_values_8p == true) {
    fprintf(stderr,"Trying to allocate %u*%d bytes of memory...",totalcounts,(int) sizeof(UINT8));
    positions4 = (UINT4 *) NULL;
    positions8 = (UINT8 *) CALLOC_NO_EXCEPTION(totalcounts,sizeof(UINT8));
    if (positions8 == NULL) {
      fprintf(stderr,"failed.  Need a computer with sufficient memory.\n");
      exit(9);
    } else {
      fprintf(stderr,"done\n");
    }
  } else {
    fprintf(stderr,"Trying to allocate %u*%d bytes of memory...",totalcounts,(int) sizeof(UINT4));
    positions8 = (UINT8 *) NULL;
    positions4 = (UINT4 *) CALLOC_NO_EXCEPTION(totalcounts,sizeof(UINT4));
    if (positions4 == NULL) {
      fprintf(stderr,"failed.  Need a computer with sufficient memory.\n");
      exit(9);
    } else {
      fprintf(stderr,"done\n");
    }
  }

  /* Copy offsets */
  pointers = (Positionsptr_T *) CALLOC(oligospace+1,sizeof(Positionsptr_T));
  for (i = 0; i <= oligospace; i++) {
    pointers[i] = offsets[i];
  }

  for (divno = 1; divno < snps_iit->ndivs; divno++) {
    divstring = IIT_divstring(snps_iit,divno);
    fprintf(stderr,"Processing positions for chromosome %s...",divstring);
    if ((chrnum = Univ_IIT_find_one(chromosome_iit,divstring)) <= 0) {
      fprintf(stderr,"not found in chromosome iit\n");
    } else {
      fprintf(stderr,"has %d snps...",IIT_nintervals(snps_iit,divno));
      chroffset = Univ_IIT_interval_low(chromosome_iit,chrnum);
      chrlength = Univ_IIT_interval_length(chromosome_iit,chrnum);
      nintervals = IIT_nintervals(snps_iit,divno);

      if (Univ_IIT_interval_type(chromosome_iit,chrnum) == circular_typeint) {
	fprintf(stderr,"and is circular...");
	nintervals_alias = 2*nintervals;
	intervals = (Labeled_interval_T *) CALLOC(nintervals_alias,sizeof(Labeled_interval_T));
	for (i = 0; i < nintervals; i++) {
	  origindex = IIT_index(snps_iit,divno,i);
	  intervals[i] = Labeled_interval_new(origindex,/*interval*/&(snps_iit->intervals[divno][i]));
	}
	for (j = 0; i < nintervals_alias; i++, j++) {
	  origindex = IIT_index(snps_iit,divno,j);
	  interval = intervals[j]->interval;
	  copy = Interval_new(/*low*/Interval_low(interval) + chrlength,
			      /*high*/Interval_high(interval) + chrlength,
			      Interval_type(interval));
	  intervals[i] = Labeled_interval_new(origindex,/*interval*/copy);
	}

      } else {
	nintervals_alias = nintervals;
	intervals = (Labeled_interval_T *) CALLOC(nintervals,sizeof(Labeled_interval_T));
	for (i = 0; i < nintervals; i++) {
	  origindex = IIT_index(snps_iit,divno,i);
	  intervals[i] = Labeled_interval_new(origindex,/*interval*/&(snps_iit->intervals[divno][i]));
	}
      }

      qsort(intervals,nintervals,sizeof(Labeled_interval_T),Labeled_interval_cmp);

      i = 0;
      while (i < nintervals_alias) {
	j = i + 1;
	while (j < nintervals_alias && Interval_low(intervals[j]->interval) < Interval_low(intervals[j-1]->interval) + index1part) {
	  j++;
	}
	process_snp_block(&nwarnings,/*offsets*/pointers,positions4,positions8,
			  &(intervals[i]),/*nintervals*/j-i,
			  chroffset,genome,/*snp_blocks*/NULL,
			  divno,divstring,/*intervali*/i,snps_iit,chromosome_iit,index1part,
			  coord_values_8p);
	i = j;
      }

      for (i = nintervals; i < nintervals_alias; i++) {
	Labeled_interval_free(&(intervals[i]),/*free_interval_p*/true);
      }
      for (i = 0; i < nintervals; i++) {
	Labeled_interval_free(&(intervals[i]),/*free_interval_p*/false);
      }

      FREE(intervals);
      fprintf(stderr,"done\n");
    }
  }

  FREE(pointers);

  /* Sort positions in each block */
  if (coord_values_8p == true) {
    for (oligoi = 0; oligoi < oligospace; oligoi++) {
      block_start = offsets[oligoi];
      block_end = offsets[oligoi+1];
      if ((npositions = block_end - block_start) > 1) {
	qsort(&(positions8[block_start]),npositions,sizeof(UINT8),UINT8_compare);
      }
    }
    return positions8;
  } else {
    for (oligoi = 0; oligoi < oligospace; oligoi++) {
      block_start = offsets[oligoi];
      block_end = offsets[oligoi+1];
      if ((npositions = block_end - block_start) > 1) {
	qsort(&(positions4[block_start]),npositions,sizeof(UINT4),UINT4_compare);
      }
    }
    return positions4;
  }
}


static void
merge_positions8 (FILE *positions_high_fp, FILE *positions_low_fp,
		  UINT8 *start1, UINT8 *end1, UINT8 *start2, UINT8 *end2, 
		  Storedoligomer_T oligo, int index1part) {
  unsigned char position8_high;
  UINT4 position8_low;
  UINT8 *ptr1 = start1, *ptr2 = start2;
  char *nt;
#ifdef WORDS_BIGENDIAN
  UINT8 position1, position2;
#endif

  while (ptr1 < end1 && ptr2 < end2) {
#ifdef WORDS_BIGENDIAN
    position1 = Bigendian_convert_uint8(*ptr1);
    position2 = Bigendian_convert_uint8(*ptr2);
    if (position1 < position2) {
      position8_high = position1 >> POSITIONS8_HIGH_SHIFT;
      position8_low = position1 & POSITIONS8_LOW_MASK;
      FWRITE_CHAR(position8_high,positions_high_fp);
      FWRITE_UINT(position8_low,positions_low_fp);
      ptr1++;
    } else if (position2 < position1) {
      position8_high = position2 >> POSITIONS8_HIGH_SHIFT;
      position8_low = position2 & POSITIONS8_LOW_MASK;
      FWRITE_CHAR(position8_high,positions_high_fp);
      FWRITE_UINT(position8_low,positions_low_fp);
      ptr2++;
    } else {
      nt = shortoligo_nt(oligo,index1part);
      fprintf(stderr,"Problem: saw duplicate positions %u in oligo %s\n",position1,nt);
      FREE(nt);
      abort();
      /*
      FWRITE_UINT8(*ptr1,positions_fp);
      ptr1++;
      ptr2++;
      */
    }

#else
    if (*ptr1 < *ptr2) {
      position8_high = *ptr1 >> POSITIONS8_HIGH_SHIFT;
      position8_low = *ptr1 & POSITIONS8_LOW_MASK;
      FWRITE_CHAR(position8_high,positions_high_fp);
      FWRITE_UINT(position8_low,positions_low_fp);
      ptr1++;
    } else if (*ptr2 < *ptr1) {
      position8_high = *ptr2 >> POSITIONS8_HIGH_SHIFT;
      position8_low = *ptr2 & POSITIONS8_LOW_MASK;
      FWRITE_CHAR(position8_high,positions_high_fp);
      FWRITE_UINT(position8_low,positions_low_fp);
      ptr2++;
    } else {
      nt = shortoligo_nt(oligo,index1part);
      fprintf(stderr,"Problem: saw duplicate positions %llu in oligo %s\n",(unsigned long long) *ptr1,nt);
      FREE(nt);
      abort();
      /*
      FWRITE_UINT8(*ptr1,positions_fp);
      ptr1++;
      ptr2++;
      */
    }
#endif
  }

  while (ptr1 < end1) {
#ifdef WORDS_BIGENDIAN
    position1 = Bigendian_convert_uint8(*ptr1);
    position8_high = position1 >> POSITIONS8_HIGH_SHIFT;
    position8_low = position1 & POSITIONS8_LOW_MASK;
#else
    position8_high = *ptr1 >> POSITIONS8_HIGH_SHIFT;
    position8_low = *ptr1 & POSITIONS8_LOW_MASK;
#endif
    FWRITE_CHAR(position8_high,positions_high_fp);
    FWRITE_UINT(position8_low,positions_low_fp);
    ptr1++;
  }

  while (ptr2 < end2) {
#ifdef WORDS_BIGENDIAN
    position2 = Bigendian_convert_uint8(*ptr2);
    position8_high = position2 >> POSITIONS8_HIGH_SHIFT;
    position8_low = position2 & POSITIONS8_LOW_MASK;
#else
    position8_high = *ptr2 >> POSITIONS8_HIGH_SHIFT;
    position8_low = *ptr2 & POSITIONS8_LOW_MASK;
#endif
    FWRITE_CHAR(position8_high,positions_high_fp);
    FWRITE_UINT(position8_low,positions_low_fp);
    ptr2++;
  }

  return;
}


static void
merge_positions4 (FILE *positions_fp, UINT4 *start1, UINT4 *end1,
		  UINT4 *start2, UINT4 *end2, Storedoligomer_T oligo, int index1part) {
  UINT4 *ptr1 = start1, *ptr2 = start2;
  char *nt;
#ifdef WORDS_BIGENDIAN
  UINT4 position1, position2;
#endif

  while (ptr1 < end1 && ptr2 < end2) {
#ifdef WORDS_BIGENDIAN
    position1 = Bigendian_convert_uint(*ptr1);
    position2 = Bigendian_convert_uint(*ptr2);
    if (position1 < position2) {
      FWRITE_UINT(position1,positions_fp);
      ptr1++;
    } else if (position2 < position1) {
      FWRITE_UINT(position2,positions_fp);
      ptr2++;
    } else {
      nt = shortoligo_nt(oligo,index1part);
      fprintf(stderr,"Problem: saw duplicate positions %u in oligo %s\n",*ptr1,nt);
      FREE(nt);
      abort();
      /*
      FWRITE_UINT(*ptr1,positions_fp);
      ptr1++;
      ptr2++;
      */
    }

#else

    if (*ptr1 < *ptr2) {
      FWRITE_UINT(*ptr1,positions_fp);
      ptr1++;
    } else if (*ptr2 < *ptr1) {
      FWRITE_UINT(*ptr2,positions_fp);
      ptr2++;
    } else {
      nt = shortoligo_nt(oligo,index1part);
      fprintf(stderr,"Problem: saw duplicate positions %u in oligo %s\n",*ptr1,nt);
      FREE(nt);
      abort();
      /*
      FWRITE_UINT(*ptr1,positions_fp);
      ptr1++;
      ptr2++;
      */
    }
#endif
  }

  while (ptr1 < end1) {
#ifdef WORDS_BIGENDIAN
    position1 = Bigendian_convert_uint(*ptr1);
    FWRITE_UINT(position1,positions_fp);
#else    
    FWRITE_UINT(*ptr1,positions_fp);
#endif
    ptr1++;
  }

  while (ptr2 < end2) {
#ifdef WORDS_BIGENDIAN
    position2 = Bigendian_convert_uint(*ptr2);
    FWRITE_UINT(position2,positions_fp);
#else
    FWRITE_UINT(*ptr2,positions_fp);
#endif
    ptr2++;
  }

  return;
}


/* Usage: snpindex -d <genome> -V <destdir> -v <snps_root> */


/* Note: Depends on having gmapindex sampled on mod 3 bounds */
int
main (int argc, char *argv[]) {
  char *sourcedir = NULL, *destdir = NULL, *mapdir = NULL;
  int compression_type;
  Univ_IIT_T chromosome_iit;
  IIT_T snps_iit;
  Genome_T genome;
  Positionsptr_T *offsets, *snp_offsets, *ref_offsets;
  size_t totalcounts, i;
#ifdef EXTRA_ALLOCATION
  Positionsptr_T npositions;
#endif

  int shmid;
  unsigned char *ref_positions8_high;
  UINT4 *ref_positions8_low;
  UINT8 *snp_positions8, *ref_positions8;
  UINT4 *snp_positions4, *ref_positions4;
  Univcoord_T nblocks;
  Genomecomp_T *snp_blocks;
  Oligospace_T oligospace, oligoi;
#ifndef HAVE_MMAP
  double seconds;
#endif

  bool coord_values_8p;
  Filenames_T filenames;
  char *filename, *filename1, *filename2;
  FILE *genome_fp, *genomebits_fp, *positions_high_fp, *positions_low_fp;
  int ref_positions_high_fd, ref_positions_low_fd;
  size_t ref_positions_high_len, ref_positions_low_len;
#ifdef WORDS_BIGENDIAN
  Positionsptr_T offset1, offset2;
#endif

  int opt;
  extern int optind;
  extern char *optarg;
  int long_option_index = 0;
  const char *long_name;

  while ((opt = getopt_long(argc,argv,"D:d:k:q:V:v:w:",
			    long_options,&long_option_index)) != -1) {
    switch (opt) {
    case 0: 
      long_name = long_options[long_option_index].name;
      if (!strcmp(long_name,"version")) {
	print_program_version();
	exit(0);
      } else if (!strcmp(long_name,"help")) {
	print_program_usage();
	exit(0);
      } else {
	/* Shouldn't reach here */
	fprintf(stderr,"Don't recognize option %s.  For usage, run 'snpindex --help'",long_name);
	exit(9);
      }
      break;

    case 'D': user_sourcedir = optarg; break;
    case 'd': dbroot = optarg; break;
    case 'k': required_index1part = atoi(optarg); break;
    case 'q': required_interval = atoi(optarg); break;
    case 'V': user_destdir = optarg; break;
    case 'v': snps_root = optarg; break;
    case 'w': max_warnings = atoi(optarg); break;
    }
  }
  argc -= (optind - 1);
  argv += (optind - 1);

  if (dbroot == NULL) {
    fprintf(stderr,"Missing name of genome database.  Must specify with -d flag.\n");
    print_program_usage();
    exit(9);
  } else if (snps_root == NULL) {
    fprintf(stderr,"Missing name of SNP database.  Must specify with -v flag.\n");
    print_program_usage();
    exit(9);
  } else {
    sourcedir = Datadir_find_genomesubdir(&fileroot,&dbversion,user_sourcedir,dbroot);
    fprintf(stderr,"Reading source files from %s\n",sourcedir);
  }

  if (argc > 1) {
    /* IIT file provided as an argument */
    if (Access_file_exists_p(argv[1]) == false) {
      fprintf(stderr,"SNP IIT file %s not found\n",argv[1]);
      exit(9);
    } else {
      fprintf(stderr,"Reading SNPs IIT file %s...",argv[1]);
      if ((snps_iit = IIT_read(argv[1],/*name*/NULL,/*readonlyp*/true,
			       /*divread*/READ_ALL,/*divstring*/NULL,/*add_iit_p*/false,
			       /*labels_read_p*/true)) == NULL) {
	fprintf(stderr,"SNP IIT file %s is not valid\n",argv[1]);
	exit(9);
      }
      fprintf(stderr,"done\n");
    } 

  } else {
    mapdir = Datadir_find_mapdir(user_sourcedir,sourcedir,fileroot);
    filename = (char *) CALLOC(strlen(mapdir)+strlen("/")+
			      strlen(snps_root)+strlen(".iit")+1,sizeof(char));
    sprintf(filename,"%s/%s.iit",mapdir,snps_root);
    if (Access_file_exists_p(filename) == false) {
      fprintf(stderr,"Map file %s.iit not found in %s.  Available files:\n",snps_root,mapdir);
      Datadir_list_directory(stderr,mapdir);
      fprintf(stderr,"Either install file %s.iit or specify IIT file as an argument\n",snps_root);
      exit(9);
    } else {
      if ((snps_iit = IIT_read(filename,/*name*/NULL,/*readonlyp*/true,
			       /*divread*/READ_ALL,/*divstring*/NULL,/*add_iit_p*/false,
			       /*labels_read_p*/true)) == NULL) {
	fprintf(stderr,"SNP IIT file %s is not valid\n",filename);
	exit(9);
      }
      fprintf(stderr,"done\n");
    }
    FREE(filename);
    FREE(mapdir);
  }


  /* Chromosome IIT file */
  filename = (char *) CALLOC(strlen(sourcedir)+strlen("/")+
			    strlen(fileroot)+strlen(".chromosome.iit")+1,sizeof(char));
  sprintf(filename,"%s/%s.chromosome.iit",sourcedir,fileroot);
  if ((chromosome_iit = Univ_IIT_read(filename,/*readonlyp*/true,/*add_iit_p*/false)) == NULL) {
    fprintf(stderr,"IIT file %s is not valid\n",filename);
    exit(9);
  } else {
    circular_typeint = Univ_IIT_typeint(chromosome_iit,"circular");
  }
  FREE(filename);
  coord_values_8p = Univ_IIT_coord_values_8p(chromosome_iit);

  fprintf(stderr,"Chromosomes in the genome: ");
  Univ_IIT_dump_labels(stderr,chromosome_iit);
  fprintf(stderr,"Chromosomes in the SNPs IIT file: ");
  IIT_dump_divstrings(stderr,snps_iit);

  genome = Genome_new(sourcedir,fileroot,/*snps_root*/NULL,/*genometype*/GENOME_OLIGOS,
		      /*uncompressedp*/false,/*access*/USE_MMAP_ONLY,/*sharedp*/false);

  /* Copy genome */
  nblocks = Genome_totallength(genome)/32U;
  snp_blocks = (Genomecomp_T *) CALLOC(nblocks*3,sizeof(Genomecomp_T));
  fprintf(stderr,"Allocating %llu*3*%d bytes for compressed genome\n",
	  (unsigned long long) nblocks,(int) sizeof(Genomecomp_T));
  memcpy(snp_blocks,Genome_blocks(genome),nblocks*3*sizeof(Genomecomp_T));

  /* Prepare for write */
  if (user_destdir == NULL) {
    destdir = sourcedir;
  } else {
    destdir = user_destdir;
  }
  fprintf(stderr,"Writing snpindex files to %s\n",destdir);

  if (max_warnings == 0) {
    show_warnings_p = false;
  }


  filenames = Indexdb_get_filenames(&compression_type,&index1part,&index1interval,
				    sourcedir,fileroot,IDX_FILESUFFIX,/*snps_root*/NULL,
				    required_index1part,required_interval,
				    /*offsets_only_p*/false);

  /* Compute offsets and write genomecomp */
  oligospace = power(4,index1part);
  snp_offsets = compute_offsets(snps_iit,chromosome_iit,genome,snp_blocks,oligospace,index1part);
  fprintf(stderr,"last offset = %u\n",snp_offsets[oligospace]);
  filename = (char *) CALLOC(strlen(destdir)+strlen("/")+strlen(fileroot)+
			     strlen(".genomecomp.")+strlen(snps_root)+1,sizeof(char));
  sprintf(filename,"%s/%s.genomecomp.%s",destdir,fileroot,snps_root);
  if ((genome_fp = FOPEN_WRITE_BINARY(filename)) == NULL) {
    fprintf(stderr,"Can't open file %s for writing genome\n",filename);
    exit(9);
  }
  fprintf(stderr,"Writing filename %s...",filename);
  FWRITE_UINTS(snp_blocks,nblocks*3,genome_fp);
  fclose(genome_fp);
  FREE(snp_blocks);
  fprintf(stderr,"done\n");


  /* Write genomebits */
  if ((genome_fp = FOPEN_READ_BINARY(filename)) == NULL) {
    fprintf(stderr,"Can't open file %s for genome\n",filename);
    exit(9);
  }
  FREE(filename);

  filename = (char *) CALLOC(strlen(destdir)+strlen("/")+strlen(fileroot)+
			     strlen(".genomebits128.")+strlen(snps_root)+1,sizeof(char));
  sprintf(filename,"%s/%s.genomebits128.%s",destdir,fileroot,snps_root);
  if ((genomebits_fp = FOPEN_WRITE_BINARY(filename)) == NULL) {
    fprintf(stderr,"Can't open file %s for writing genome\n",filename);
    exit(9);
  }
  fprintf(stderr,"Writing filename %s...",filename);
  Compress_unshuffle_bits128(/*out*/genomebits_fp,/*in*/genome_fp);
  fclose(genomebits_fp);
  fclose(genome_fp);
  FREE(filename);


  /* Compute positions */
  show_warnings_p = false;	/* Already shown in compute_offsets */
  if (coord_values_8p == true) {
    snp_positions8 = compute_positions(snp_offsets,snps_iit,chromosome_iit,genome,oligospace,/*coord_values_8p*/true);
  } else {
    snp_positions4 = compute_positions(snp_offsets,snps_iit,chromosome_iit,genome,oligospace,/*coord_values_8p*/false);
  }


  /* Read reference offsets and update */
#ifdef EXTRA_ALLOCATION
  if (compression_type == BITPACK64_COMPRESSION) {
    ref_offsets = Indexdb_offsets_from_bitpack(filenames->pointers_filename,filenames->offsets_filename,index1part);
  } else {
    abort();
  }
  offsets = (Positionsptr_T *) CALLOC(oligospace+1,sizeof(Positionsptr_T));
  offsets[0] = 0U;
  for (oligoi = 1; oligoi <= oligospace; oligoi++) {
#ifdef WORDS_BIGENDIAN
    npositions = (Bigendian_convert_uint(ref_offsets[oligoi]) - Bigendian_convert_uint(ref_offsets[oligoi-1])) + 
      (snp_offsets[oligoi] - snp_offsets[oligoi-1]);
#else
    npositions = (ref_offsets[oligoi] - ref_offsets[oligoi-1]) + (snp_offsets[oligoi] - snp_offsets[oligoi-1]);
#endif
    offsets[oligoi] = offsets[oligoi-1] + npositions;
  }

#else
  /* This version saves on one allocation of oligospace * sizeof(Positionsptr_T) */
  if (compression_type == BITPACK64_COMPRESSION) {
    offsets = Indexdb_offsets_from_bitpack(filenames->pointers_filename,filenames->offsets_filename,index1part);
  } else {
    abort();
  }
  for (oligoi = oligospace; oligoi >= 1; oligoi--) {
#ifdef WORDS_BIGENDIAN
    offsets[oligoi] = Bigendian_convert_uint(offsets[oligoi]) - Bigendian_convert_uint(offsets[oligoi-1]);
#else
    offsets[oligoi] = offsets[oligoi] - offsets[oligoi-1];
#endif
  }

  for (oligoi = 1; oligoi <= oligospace; oligoi++) {
    offsets[oligoi] = offsets[oligoi-1] + offsets[oligoi] + (snp_offsets[oligoi] - snp_offsets[oligoi-1]);
  }
#endif


  /* Write offsets */
  if (filenames->pointers_basename_ptr == NULL) {
    filename1 = (char *) NULL;
  } else {
    filename1 = (char *) CALLOC(strlen(destdir)+strlen("/")+strlen(filenames->pointers_basename_ptr)+
				strlen(".")+strlen(snps_root)+1,sizeof(char));
    sprintf(filename1,"%s/%s.%s",destdir,filenames->pointers_basename_ptr,snps_root);
  }

  filename2 = (char *) CALLOC(strlen(destdir)+strlen("/")+strlen(filenames->offsets_basename_ptr)+
			      strlen(".")+strlen(snps_root)+1,sizeof(char));
  sprintf(filename2,"%s/%s.%s",destdir,filenames->offsets_basename_ptr,snps_root);


  fprintf(stderr,"Writing %llu offsets with %llu total positions\n",
	  (unsigned long long) oligospace+1,(unsigned long long) offsets[oligospace]);
  if (compression_type == BITPACK64_COMPRESSION) {
    Bitpack64_write_differential(/*ptrsfile*/filename1,/*compfile*/filename2,offsets,oligospace);
  } else {
    abort();
  }
  FREE(filename2);
  if (filename1 != NULL) {
    FREE(filename1);
  }
  FREE(offsets);


  /* Read reference positions and merge */
#if 0
  if ((positions_low_fp = FOPEN_READ_BINARY(filenames->positions_low_filename)) == NULL) {
    fprintf(stderr,"Can't open file %s\n",filenames->positions_low_filename);
    exit(9);
  } else {
    fclose(positions_low_fp);
  }
#endif

  if (coord_values_8p == true) {
#ifdef HAVE_MMAP
    ref_positions8_high = (unsigned char *) Access_mmap(&ref_positions_high_fd,&ref_positions_high_len,
							filenames->positions_high_filename,sizeof(unsigned char),/*randomp*/false);
    ref_positions8_low = (UINT4 *) Access_mmap(&ref_positions_low_fd,&ref_positions_low_len,
					       filenames->positions_low_filename,sizeof(UINT4),/*randomp*/false);
#else
    ref_positions8_high = (unsigned char *) Access_allocate(&shmid,&ref_positions_high_len,&seconds,
							    filenames->positions_high_filename,sizeof(unsigned char),
							    /*sharedp*/false);
    ref_positions8_low = (UINT4 *) Access_allocate(&shmid,&ref_positions_low_len,&seconds,
						   filenames->positions_low_filename,sizeof(UINT4),/*sharedp*/false);
#endif
    /* Unpack */
    totalcounts = ref_positions_high_len/sizeof(unsigned char);
    if (totalcounts > 4294967295) {
      fprintf(stderr,"Program not yet designed to handle huge genomes\n");
      abort();
    }
    ref_positions8 = (UINT8 *) CALLOC_NO_EXCEPTION(totalcounts,sizeof(UINT8));
    for (i = 0; i < totalcounts; i++) {
      ref_positions8[i] = ((UINT8) ref_positions8_high[i] << 32) + ref_positions8_low[i];
    }
#ifdef HAVE_MMAP
    munmap((void *) ref_positions8_high,ref_positions_high_len);
    munmap((void *) ref_positions8_low,ref_positions_low_len);
    close(ref_positions_high_fd);
    close(ref_positions_low_fd);
#else
    FREE(ref_positions8_high);
    FREE(ref_positions8_low);
#endif

  } else {
#ifdef HAVE_MMAP
    ref_positions4 = (UINT4 *) Access_mmap(&ref_positions_low_fd,&ref_positions_low_len,
					   filenames->positions_low_filename,sizeof(UINT4),/*randomp*/false);
#else
    ref_positions4 = (UINT4 *) Access_allocate(&shmid,&ref_positions_low_len,&seconds,
					       filenames->positions_low_filename,sizeof(UINT4),/*sharedp*/false);
#endif
  }


  if (coord_values_8p == true) {
    filename = (char *) CALLOC(strlen(destdir)+strlen("/")+strlen(filenames->positions_high_basename_ptr)+
			       strlen(".")+strlen(snps_root)+1,sizeof(char));
    sprintf(filename,"%s/%s.%s",destdir,filenames->positions_high_basename_ptr,snps_root);
    if ((positions_high_fp = FOPEN_WRITE_BINARY(filename)) == NULL) {
      fprintf(stderr,"Can't open file %s for writing\n",filename);
      exit(9);
    }
    FREE(filename);
  }

  filename = (char *) CALLOC(strlen(destdir)+strlen("/")+strlen(filenames->positions_low_basename_ptr)+
			     strlen(".")+strlen(snps_root)+1,sizeof(char));
  sprintf(filename,"%s/%s.%s",destdir,filenames->positions_low_basename_ptr,snps_root);
  if ((positions_low_fp = FOPEN_WRITE_BINARY(filename)) == NULL) {
    fprintf(stderr,"Can't open file %s for writing\n",filename);
    exit(9);
  }
  FREE(filename);


  fprintf(stderr,"Merging snp positions with reference positions\n");
#ifndef EXTRA_ALLOCATION
  if (compression_type == BITPACK64_COMPRESSION) {
    ref_offsets = Indexdb_offsets_from_bitpack(filenames->pointers_filename,filenames->offsets_filename,index1part);
  } else {
    abort();
  }
#endif

  if (coord_values_8p == true) {
    for (oligoi = 0; oligoi < oligospace; oligoi++) {
#ifdef WORDS_BIGENDIAN
      offset1 = Bigendian_convert_uint(ref_offsets[oligoi]);
      offset2 = Bigendian_convert_uint(ref_offsets[oligoi+1]);
      merge_positions8(positions_high_fp,positions_low_fp,
		       &(snp_positions8[snp_offsets[oligoi]]),&(snp_positions8[snp_offsets[oligoi+1]]),
		       &(ref_positions8[offset1]),&(ref_positions8[offset2]),oligoi,index1part);
#else
      merge_positions8(positions_high_fp,positions_low_fp,
		       &(snp_positions8[snp_offsets[oligoi]]),&(snp_positions8[snp_offsets[oligoi+1]]),
		       &(ref_positions8[ref_offsets[oligoi]]),&(ref_positions8[ref_offsets[oligoi+1]]),oligoi,index1part);
#endif
    }
  } else {
    for (oligoi = 0; oligoi < oligospace; oligoi++) {
#ifdef WORDS_BIGENDIAN
      offset1 = Bigendian_convert_uint(ref_offsets[oligoi]);
      offset2 = Bigendian_convert_uint(ref_offsets[oligoi+1]);
      merge_positions4(positions_low_fp,&(snp_positions4[snp_offsets[oligoi]]),&(snp_positions4[snp_offsets[oligoi+1]]),
		       &(ref_positions4[offset1]),&(ref_positions4[offset2]),oligoi,index1part);
#else
      merge_positions4(positions_low_fp,&(snp_positions4[snp_offsets[oligoi]]),&(snp_positions4[snp_offsets[oligoi+1]]),
		       &(ref_positions4[ref_offsets[oligoi]]),&(ref_positions4[ref_offsets[oligoi+1]]),oligoi,index1part);
#endif
    }
  }


  FREE(ref_offsets);
  fclose(positions_low_fp);
  if (coord_values_8p == true) {
    fclose(positions_high_fp);
  }


  /* Clean up */
  if (coord_values_8p == true) {
    /* For both mmap and allocated, since we have already combined positions_high and positions_low */
    FREE(ref_positions8);
  } else {
#ifdef HAVE_MMAP
    munmap((void *) ref_positions4,ref_positions_low_len);
    close(ref_positions_low_fd);
#else
    FREE(ref_positions4);
#endif
  }


  if (coord_values_8p == true) {
    FREE(snp_positions8);
  } else {
    FREE(snp_positions4);
  }
  FREE(snp_offsets);

  Genome_free(&genome);
  Univ_IIT_free(&chromosome_iit);
  IIT_free(&snps_iit);

  Filenames_free(&filenames);


  if (argc <= 1) {
    /* Program called using -v flag only.  No need to install. */
    /* fprintf(stderr,"IIT file already present in .maps directory\n"); */
#if 0
    /* Old code used for copying IIT file to .maps directory */
    /* To use this code, cannot free mapdir earlier */
    fprintf(stderr,"Now copying IIT file from %s/%s.iit to %s...",mapdir,snps_root,filename);
    filename1 = (char *) CALLOC(strlen(mapdir)+strlen("/")+
				strlen(snps_root)+strlen(".iit")+1,sizeof(char));
    sprintf(filename1,"%s/%s.iit",mapdir,snps_root);
    Access_file_copy(/*dest*/filename,/*source*/filename1);
    fprintf(stderr,"done\n");
    FREE(filename1);
    FREE(mapdir);
#endif

  } else {
    /* Install IIT file */
    if (!strcmp(destdir,sourcedir)) {
      filename = (char *) CALLOC(strlen(destdir)+strlen("/")+ strlen(fileroot) + strlen(".maps/") +
				 strlen(snps_root)+strlen(".iit")+1,sizeof(char));
      sprintf(filename,"%s/%s.maps/%s.iit",destdir,fileroot,snps_root);
    } else {
      filename = (char *) CALLOC(strlen(destdir)+strlen("/")+
				 strlen(snps_root)+strlen(".iit")+1,sizeof(char));
      sprintf(filename,"%s/%s.iit",destdir,snps_root);
    }

    fprintf(stderr,"SNP genome indices created.\n");
    if (Access_file_exists_p(filename) == true) {
      if (Access_file_equal(filename,argv[1]) == false) {
	fprintf(stderr,"IIT file already present as %s, but it is different from the given file %s\n",
		filename,argv[1]);
	fprintf(stderr,"Please copy file %s as %s\n",argv[1],filename);
      } else {
	fprintf(stderr,"IIT file already present as %s, and it is the same as given file %s\n",
		filename,argv[1]);
      }

    } else {
      fprintf(stderr,"Now installing IIT file %s as %s...",
	      argv[1],filename);
      Access_file_copy(/*dest*/filename,/*source*/argv[1]);
      fprintf(stderr,"done\n");
    }

    FREE(filename);
  }

  FREE(dbversion);
  FREE(fileroot);
  FREE(sourcedir);

  return 0;
}


static void
print_program_usage () {
  fprintf(stdout,"\
Usage: snpindex [OPTIONS...] -d <genome> -v <snpsdb> [<iitfile>]\n\
\n\
If iitfile is provided as a non-flag argument, then use that iitfile and create SNP database\n\
as named by -v flag.  Otherwise, try to find iit file named <snpsdb>.iit in GMAP index files\n\
for <genome>.\n\
\n\
");

  /* Input options */
  fprintf(stdout,"Options (must include -d)\n");
  fprintf(stdout,"\
  -D, --sourcedir=directory      Directory where to read genome index files (default is\n\
                                   GMAP genome directory specified at compile time)\n\
  -d, --db=STRING                Genome database\n\
  -k, --kmer=INT                 kmer size to use in genome database (allowed values: 16 or less).\n\
                                   If not specified, the program will find the highest available\n\
                                   kmer size in the genome database\n\
  -q, --sampling=INT             Sampling to use in genome database.  If not specified, the program\n\
                                   will find the smallest available sampling value in the genome database\n\
                                   within selected k-mer size\n\
  -V, --destdir=directory        Directory where to write SNP index files (default is\n\
                                   GMAP genome directory specified at compile time)\n\
  -v, --snpsdb=STRING            Name of SNP database\n\
  -w, --max-warnings=INT         Maximum number of warnings to print to stderr about\n\
                                   inconsistencies relative to the reference genome.\n\
                                   A value of 0 turns off all warnings.  A negative value\n\
                                   prints all warnings.  (default -1, meaning no limit)\n\
\n\
  --version                      Show version\n\
  --help                         Show this help message\n\
");
  return;
}

