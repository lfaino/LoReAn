static char rcsid[] = "$Id: sarray-write.c 170326 2015-07-22 17:49:55Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "sarray-write.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/mman.h>		/* For munmap */
#include <math.h>		/* For rint */

#include "bool.h"
#include "access.h"
#include "mem.h"
#include "genomicpos.h"
#include "assert.h"
#include "compress.h"
#include "bitpack64-write.h"
#include "bitpack64-read.h"
#include "bitpack64-access.h"	/* For Sarray_plcp_compare */
#include "bytecoding.h"
#include "fopen.h"
#include "saca-k.h"
#include "genome128_hr.h"
#include "uintlist.h"
#include "intlist.h"


#ifdef WORDS_BIGENDIAN
#include "bigendian.h"
#else
#include "littleendian.h"
#endif


/* make_index */
#ifdef DEBUG1
#define debug1(x) x
#else
#define debug1(x)
#endif

/* Sarray_compute_child */
#ifdef DEBUG2
#define debug2(x) x
#else
#define debug2(x)
#endif

/* Sarray_discriminating_chars */
#ifdef DEBUG4
#define debug4(x) x
#else
#define debug4(x)
#endif

/* correctness of using Genome_consecutive_matches_pair */
#ifdef DEBUG14
#define debug14(x) x
#else
#define debug14(x)
#endif

/* correctness of make_index_incremental */
#ifdef DEBUG15
#define debug15(x) x
#else
#define debug15(x)
#endif



/* For computing LCP.  Comment out because mmap is faster than fread */
/* #define READ_SA_FROM_FILE 1 */

#define MONITOR_INTERVAL 100000000 /* 100 million nt */
#define RW_BATCH  10000000	/* 10 million elements */

/* For standard genome */
void
Sarray_write_array (char *sarrayfile, Genome_T genomecomp, UINT4 genomelength) {
  UINT4 *SA;
  UINT4 n = genomelength, ii;
  unsigned char *gbuffer;
  FILE *fp;
  void *p;
  
  SA = (UINT4 *) MALLOC((n+1)*sizeof(UINT4));
  gbuffer = (unsigned char *) CALLOC(n+1,sizeof(unsigned char));
  Genome_fill_buffer_int_string(genomecomp,/*left*/0,/*length*/n,gbuffer,/*conversion*/NULL);
  gbuffer[n] = 0;		       /* Tried N/X, but SACA_K fails */
  SACA_K(gbuffer,SA,n+/*virtual sentinel*/1,/*K, alphabet_size*/5,/*m*/n+1,/*level*/0);

  if ((fp = FOPEN_WRITE_BINARY(sarrayfile)) == NULL) {
    fprintf(stderr,"Can't write to file %s\n",sarrayfile);
    exit(9);
  } else {
#if 0
    FWRITE_UINTS(SA,n+1,fp);
#else
    for (ii = 0; ii + RW_BATCH <= n; ii += RW_BATCH) {
      p = (void *) &(SA[ii]);
      FWRITE_UINTS(p,RW_BATCH,fp);
    }
    if (ii <= n) {
      p = (void *) &(SA[ii]);
      FWRITE_UINTS(p,n - ii + 1,fp);
    }
#endif
    fclose(fp);
  }

  FREE(gbuffer);
  FREE(SA);

  return;
}


void
Sarray_write_array_from_genome (char *sarrayfile, unsigned char *gbuffer, UINT4 genomelength) {
  UINT4 *SA;
  UINT4 n = genomelength, ii;
  FILE *fp;
  void *p;


  SA = (UINT4 *) MALLOC((n+1)*sizeof(UINT4));
  SACA_K(gbuffer,SA,n+/*virtual sentinel*/1,/*K, alphabet_size*/5,/*m*/n+1,/*level*/0);

  if ((fp = FOPEN_WRITE_BINARY(sarrayfile)) == NULL) {
    fprintf(stderr,"Can't write to file %s\n",sarrayfile);
    exit(9);
  } else {
#if 0
    FWRITE_UINTS(SA,n+1,fp);
#else
    for (ii = 0; ii + RW_BATCH <= n; ii += RW_BATCH) {
      p = (void *) &(SA[ii]);
      FWRITE_UINTS(p,RW_BATCH,fp);
    }
    if (ii <= n) {
      p = (void *) &(SA[ii]);
      FWRITE_UINTS(p,n - ii + 1,fp);
    }
#endif
    fclose(fp);
  }

  FREE(SA);

  return;
}


#define MIN_INDEXSIZE 12
#define MAX_INDEXSIZE 12
#define INDEX_MONITOR_INTERVAL 100000

static UINT4
power (int base, int exponent) {
  UINT4 result = 1U;
  int i;

  for (i = 0; i < exponent; i++) {
    result *= base;
  }
  return result;
}


static void
sarray_search_char (Sarrayptr_T *initptr, Sarrayptr_T *finalptr, char desired_char,
		    Genome_T genomecomp, UINT4 *SA, int n, char *chartable) {
  Sarrayptr_T low, high, mid;
  Univcoord_T pos;
  char c;


  low = 1;
  high = n + 1;

  while (low < high) {
    /* Compute mid for unsigned ints.  Want floor((low+high)/2). */
    mid = low/2 + high/2;
    if (low % 2 == 1 && high % 2 == 1) {
      mid += 1;
    }
#ifdef WORDS_BIGENDIAN
    pos = Bigendian_convert_uint(SA[mid]);
#else
    pos = SA[mid];
#endif
    c = Genome_get_char_lex(genomecomp,pos,n,chartable);
    if (desired_char > c) {
      low = mid + 1;
    } else {
      high = mid;
    }
  }

  *initptr = low;

  low--;
  high = n;
  while (low < high) {
    /* Compute mid for unsigned ints.  Want ceil((low+high)/2). */
    mid = low/2 + high/2;
    if (low % 2 == 1 || high % 2 == 1) {
      mid += 1;
    }
#ifdef WORDS_BIGENDIAN
    pos = Bigendian_convert_uint(SA[mid]);
#else
    pos = SA[mid];
#endif
    c = Genome_get_char_lex(genomecomp,pos,n,chartable);
    if (desired_char >= c) {
      low = mid;
    } else {
      high = mid - 1;
    }
  }

  *finalptr = high;
  return;
}


#define LOW_TWO_BITS 0x3
#define RIGHT_A 0
#define RIGHT_C 1
#define RIGHT_G 2
#define RIGHT_T 3


static void
oligo_nt (char *nt, UINT4 oligo, int oligosize) {
  int i, j;
  UINT4 lowbits;

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

  return;
}


/* Taken from Johannes Fischer, Advanced Text Indexing Techniques, Algorithm 1 */
/* Does not use LCP, so time is O(m * log(n)) */
/* Result should be within [i..j] */
static void
sarray_search_simple (Sarrayptr_T *initptr, Sarrayptr_T *finalptr, char *query,
		      int querylength, Genome_T genomecomp, UINT4 *SA,
		      UINT4 i, UINT4 j, UINT4 n, char *chartable) {
  Sarrayptr_T low, high, mid;
  Univcoord_T pos;
  int nmatches;
  char c;


  low = i;
  high = j+1;

  while (low < high) {
    /* Compute mid for unsigned ints.  Want floor((low+high)/2). */
    mid = low/2 + high/2;
    if (low % 2 == 1 && high % 2 == 1) {
      mid += 1;
    }

    nmatches = 0;
#ifdef WORDS_BIGENDIAN
    pos = Bigendian_convert_uint(SA[mid]);
#else
    pos = SA[mid];
#endif

    while (nmatches < querylength && (c = Genome_get_char_lex(genomecomp,pos,n,chartable)) == query[nmatches]) {
      nmatches++;
      pos++;
    }
    if (nmatches == querylength || c > query[nmatches]) {
      high = mid;
    } else {
      low = mid + 1;
    }
  }

  *initptr = low;

  low--;
  high = j;

  while (low < high) {
    /* Compute mid for unsigned ints.  Want ceil((low+high)/2). */
    mid = low/2 + high/2;
    if (low % 2 == 1 || high % 2 == 1) {
      mid += 1;
    }

    nmatches = 0;
#ifdef WORDS_BIGENDIAN
    pos = Bigendian_convert_uint(SA[mid]);
#else
    pos = SA[mid];
#endif

    while (nmatches < querylength && (c = Genome_get_char_lex(genomecomp,pos,n,chartable)) == query[nmatches]) {
      nmatches++;
      pos++;
    }
    if (nmatches == querylength || c < query[nmatches]) {
      low = mid;
    } else {
      high = mid - 1;
    }
  }

  *finalptr = high;

  return;
}


/* Separate representation: Two arrays, one for start (i) and one for
   end (j).  The end index is inclusive, so suffix array entries are in
   [indexi[oligo],indexj[oligo]]. */

/* Interleaved representation.  Alternate start (i) and end (j).  The
   end index is exclusive, so suffix array entries are in
   [indexi[oligo],indexj[oligo]-1].  Need to store j+1 to maintain
   monotonicity, since an empty interval has j < i (or j == i - 1). */


static UINT4
make_index_separate (Sarrayptr_T *saindexi, Sarrayptr_T *saindexj,
		     UINT4 oligospace, int querylength, Genome_T genomecomp, UINT4 *SA, UINT4 n,
		     char chartable[]) {
  UINT4 noccupied = 0;
  char *queryuc_ptr;
  UINT4 oligo;
  char *comma;


  queryuc_ptr = (char *) CALLOC(querylength+1,sizeof(char));

  for (oligo = 0; oligo < oligospace; oligo++) {
    oligo_nt(queryuc_ptr,oligo,querylength);
    sarray_search_simple(&(saindexi[oligo]),&(saindexj[oligo]),queryuc_ptr,querylength,
			 genomecomp,SA,/*i*/1,/*j*/n,n,chartable);
    if (saindexi[oligo] <= saindexj[oligo]) {
      debug1(printf("%u\t%s\t%u\t%u\n",oligo,queryuc_ptr,saindexi[oligo],saindexj[oligo]));
      noccupied++;
    }

#if 0
    if (oligo % INDEX_MONITOR_INTERVAL == 0) {
      comma = Genomicpos_commafmt(oligo);
      fprintf(stderr,"Computing index %s\n",comma);
      FREE(comma);
    }
#endif
  }
  FREE(queryuc_ptr);

  return noccupied;
}


static UINT4
make_index_interleaved (Sarrayptr_T *saindex,
			UINT4 oligospace, int querylength, Genome_T genomecomp, UINT4 *SA, UINT4 n,
			char chartable[]) {
  UINT4 noccupied = 0;
  char *queryuc_ptr;
  UINT4 oligo;
  char *comma;
  Sarrayptr_T indexi, indexj;

  queryuc_ptr = (char *) CALLOC(querylength+1,sizeof(char));

  for (oligo = 0; oligo < oligospace; oligo++) {
    oligo_nt(queryuc_ptr,oligo,querylength);
    sarray_search_simple(&indexi,&indexj,queryuc_ptr,querylength,
			 genomecomp,SA,/*i*/1,/*j*/n,n,chartable);
#if 0
    if (indexj < indexi) {
      fprintf(stderr,"Warning: j %u < i %u\n",indexj,indexi);
    }
#endif

    /* Need to add 1 to indexj, because an empty lcp-interval has j == i - 1 */
    saindex[2*oligo] = indexi;
    saindex[2*oligo+1] = indexj+1;
    assert(indexi <= indexj+1);

    if (indexi <= indexj) {
      debug1(printf("%u\t%s\t%u\t%u\n",oligo,queryuc_ptr,indexi,indexj));
      noccupied++;
    }

#if 0
    if (oligo % INDEX_MONITOR_INTERVAL == 0) {
      comma = Genomicpos_commafmt(oligo);
      fprintf(stderr,"Computing index %s\n",comma);
      FREE(comma);
    }
#endif
  }
  FREE(queryuc_ptr);

  return noccupied;
}



#if 0
/* oligo is based on old indexsize. indexsize is new indexsize. */
static int
sarray_search_incremental (Sarrayptr_T *initptrs, Sarrayptr_T *finalptrs,
			   Oligospace_T oligo, char *query, Genome_T genomecomp, UINT4 *SA,
			   Sarrayptr_T *saindexi, Sarrayptr_T *saindexj,
			   int indexsize, Oligospace_T prev_oligospace, UINT4 n) {
  int noccupied = 0, k;
  Sarrayptr_T i, j;
  Oligospace_T prev_oligo;

  prev_oligo = oligo/4;
  if (saindexj[prev_oligo] < saindexi[prev_oligo]) {
    /* Oligo from 0..(indexsize)-1 does not match, so lengthening also will not match */
    initptrs[0] = initptrs[1] = initptrs[2] = initptrs[3] = saindexi[prev_oligo];
    finalptrs[0] = finalptrs[1] = finalptrs[2] = finalptrs[3] = saindexj[prev_oligo];

  } else {
    if (oligo == 0) {
      i = 1;
      j = saindexj[prev_oligo + 1];
    } else if (prev_oligo + 1 == prev_oligospace) {
      i = saindexi[prev_oligo - 1];
      j = n;
    } else {
      i = saindexi[prev_oligo - 1];
      j = saindexj[prev_oligo + 1];
    }

    for (k = 0; k < 4; k++) {
      oligo_nt(query,oligo+k,indexsize);
      sarray_search_simple(&(initptrs[k]),&(finalptrs[k]),query,/*querylength*/indexsize,
			   genomecomp,SA,i,j,n,chartable);
      if (initptrs[k] <= finalptrs[k]) {
	noccupied++;
      }
    }
  }

  return noccupied;
}
#endif


#if 0
static UINT4
make_index_incremental (Sarrayptr_T *initptrs, Sarrayptr_T *finalptrs,
			Genome_T genomecomp, UINT4 *SA, Sarrayptr_T *saindexi, Sarrayptr_T *saindexj, 
			int indexsize, Oligospace_T oligospace, Oligospace_T prev_oligospace, UINT4 n) {
  int noccupied = 0;
  char *query;
  Oligospace_T oligo;
  char *comma;

  query = (char *) CALLOC(indexsize+1,sizeof(char));

  for (oligo = 0; oligo < oligospace; oligo += 4) {
    noccupied += sarray_search_incremental(&(initptrs[oligo]),&(finalptrs[oligo]),oligo,query,
					   genomecomp,SA,saindexi,saindexj,indexsize,prev_oligospace,n,
					   chartable);
#if 0
    if (oligo % INDEX_MONITOR_INTERVAL == 0) {
      comma = Genomicpos_commafmt(oligo);
      fprintf(stderr,"Computing index %s\n",comma);
      FREE(comma);
    }
#endif
  }

  FREE(query);

  return noccupied;
}
#endif


void
Sarray_write_index_separate (char *indexiptrsfile, char *indexicompfile, char *indexjptrsfile,char *indexjcompfile,
			     char *sarrayfile, Genome_T genomecomp, UINT4 genomelength, bool compressp,
			     char chartable[]) {
  UINT4 n = genomelength;
  Oligospace_T oligospace, prev_oligospace, noccupied;
  /* Oligospace_T prev_noccupied; */
  Sarrayptr_T *saindexi_new, *saindexj_new, *saindexi_old, *saindexj_old;
  UINT4 *SA;
  int sa_fd;
  size_t sa_len;
  int indexsize;
  FILE *fp;


  SA = (UINT4 *) Access_mmap(&sa_fd,&sa_len,sarrayfile,sizeof(UINT4),/*randomp*/true);

  indexsize = MIN_INDEXSIZE;
  oligospace = power(4,/*querylength*/indexsize);
  saindexi_old = (Sarrayptr_T *) CALLOC(oligospace,sizeof(Sarrayptr_T));
  saindexj_old = (Sarrayptr_T *) CALLOC(oligospace,sizeof(Sarrayptr_T));
  /* prev_noccupied = 0; */
  noccupied = make_index_separate(saindexi_old,saindexj_old,
				  oligospace,/*querylength*/indexsize,genomecomp,SA,n,chartable);
  fprintf(stderr,"For indexsize %d, occupied %u/%u\n",indexsize,noccupied,oligospace);

#if 0
  indexsize++;
  while (indexsize <= MAX_INDEXSIZE && noccupied > prev_noccupied) {
    prev_noccupied = noccupied;
    prev_oligospace = oligospace;
    oligospace = power(4,/*querylength*/indexsize);
    saindexi_new = (Sarrayptr_T *) CALLOC(oligospace,sizeof(Sarrayptr_T));
    saindexj_new = (Sarrayptr_T *) CALLOC(oligospace,sizeof(Sarrayptr_T));
    noccupied = make_index_incremental(saindexi_new,saindexj_new,genomecomp,SA,
				       saindexi_old,saindexj_old,indexsize,
				       oligospace,prev_oligospace,n);
    fprintf(stderr,"For indexsize %d, occupied %u/%u\n",indexsize,noccupied,oligospace);
    if (noccupied > prev_noccupied) {
      FREE(saindexj_old);
      FREE(saindexi_old);
      saindexi_old = saindexi_new;
      saindexj_old = saindexj_new;
      indexsize++;
    } else {
      FREE(saindexj_new);
      FREE(saindexi_new);
    }
  }
  indexsize--;
#endif

  oligospace = power(4,/*querylength*/indexsize);
  fprintf(stderr,"Optimal indexsize = %d\n",indexsize);

#ifdef DEBUG15
  /* For comparison */
  fprintf(stderr,"Checking...");
  query = (char *) CALLOC(indexsize+1,sizeof(char));
  saindexi_new = (Sarrayptr_T *) CALLOC(oligospace,sizeof(Sarrayptr_T));
  saindexj_new = (Sarrayptr_T *) CALLOC(oligospace,sizeof(Sarrayptr_T));
  noccupied = make_index(saindexi_new,saindexj_new,oligospace,
			 /*querylength*/indexsize,genomecomp,SA,n,chartable);
  
  for (oligo = 0; oligo < oligospace; oligo++) {
    if (saindexi_old[oligo] != saindexi_new[oligo]) {
      oligo_nt(query,oligo,indexsize);
      printf("%u\t%s\t%u\t%u\t%u\t%u\n",oligo,query,saindexi_old[oligo],saindexj_old[oligo],saindexi_new[oligo],saindexj_new[oligo]);
      abort();
    } else if (saindexj_old[oligo] != saindexj_new[oligo]) {
      oligo_nt(query,oligo,indexsize);
      printf("%u\t%s\t%u\t%u\t%u\t%u\n",oligo,query,saindexi_old[oligo],saindexj_old[oligo],saindexi_new[oligo],saindexj_new[oligo]);
      abort();
    }
  }
  FREE(query);
  FREE(saindexj_new);
  FREE(saindexi_new);
  fprintf(stderr,"done\n");
#endif

  if (compressp == false) {
    fp = fopen(indexicompfile,"w");
    FWRITE_UINTS(saindexi_old,oligospace,fp);
    fclose(fp);

    fp = fopen(indexjcompfile,"w");
    FWRITE_UINTS(saindexj_old,oligospace,fp);
    fclose(fp);
    
  } else {
    Bitpack64_write_differential(indexiptrsfile,indexicompfile,saindexi_old,oligospace-1);
    Bitpack64_write_differential(indexjptrsfile,indexjcompfile,saindexj_old,oligospace-1);
  }

  FREE(saindexj_old);
  FREE(saindexi_old);

  munmap((void *) SA,sa_len);
  close(sa_fd);

  return;
}


void
Sarray_write_index_interleaved (char *indexptrsfile, char *indexcompfile,
				char *sarrayfile, Genome_T genomecomp, UINT4 genomelength, bool compressp,
				char chartable[]) {
  UINT4 n = genomelength;
  Oligospace_T oligospace, prev_oligospace, noccupied;
  /* Oligospace_T prev_noccupied; */
  Sarrayptr_T *saindex_new, *saindex_old;
  UINT4 *SA;
  int sa_fd;
  size_t sa_len;
  int indexsize;
  FILE *fp;


  SA = (UINT4 *) Access_mmap(&sa_fd,&sa_len,sarrayfile,sizeof(UINT4),/*randomp*/true);

  indexsize = MIN_INDEXSIZE;
  oligospace = power(4,/*querylength*/indexsize);
  saindex_old = (Sarrayptr_T *) CALLOC(2*oligospace,sizeof(Sarrayptr_T));
  /* prev_noccupied = 0; */
  noccupied = make_index_interleaved(saindex_old,
				     oligospace,/*querylength*/indexsize,genomecomp,SA,n,chartable);
  fprintf(stderr,"For indexsize %d, occupied %u/%u\n",indexsize,noccupied,oligospace);

#if 0
  indexsize++;
  while (indexsize <= MAX_INDEXSIZE && noccupied > prev_noccupied) {
    prev_noccupied = noccupied;
    prev_oligospace = oligospace;
    oligospace = power(4,/*querylength*/indexsize);
    saindexi_new = (Sarrayptr_T *) CALLOC(oligospace,sizeof(Sarrayptr_T));
    saindexj_new = (Sarrayptr_T *) CALLOC(oligospace,sizeof(Sarrayptr_T));
    noccupied = make_index_incremental(saindexi_new,saindexj_new,genomecomp,SA,
				       saindexi_old,saindexj_old,indexsize,
				       oligospace,prev_oligospace,n);
    fprintf(stderr,"For indexsize %d, occupied %u/%u\n",indexsize,noccupied,oligospace);
    if (noccupied > prev_noccupied) {
      FREE(saindexj_old);
      FREE(saindexi_old);
      saindexi_old = saindexi_new;
      saindexj_old = saindexj_new;
      indexsize++;
    } else {
      FREE(saindexj_new);
      FREE(saindexi_new);
    }
  }
  indexsize--;
#endif

  oligospace = power(4,/*querylength*/indexsize);
  fprintf(stderr,"Optimal indexsize = %d\n",indexsize);

#ifdef DEBUG15
  /* For comparison */
  fprintf(stderr,"Checking...");
  query = (char *) CALLOC(indexsize+1,sizeof(char));
  saindexi_new = (Sarrayptr_T *) CALLOC(oligospace,sizeof(Sarrayptr_T));
  saindexj_new = (Sarrayptr_T *) CALLOC(oligospace,sizeof(Sarrayptr_T));
  noccupied = make_index(saindexi_new,saindexj_new,oligospace,
			 /*querylength*/indexsize,genomecomp,SA,n);
  
  for (oligo = 0; oligo < oligospace; oligo++) {
    if (saindexi_old[oligo] != saindexi_new[oligo]) {
      oligo_nt(query,oligo,indexsize);
      printf("%u\t%s\t%u\t%u\t%u\t%u\n",oligo,query,saindexi_old[oligo],saindexj_old[oligo],saindexi_new[oligo],saindexj_new[oligo]);
      abort();
    } else if (saindexj_old[oligo] != saindexj_new[oligo]) {
      oligo_nt(query,oligo,indexsize);
      printf("%u\t%s\t%u\t%u\t%u\t%u\n",oligo,query,saindexi_old[oligo],saindexj_old[oligo],saindexi_new[oligo],saindexj_new[oligo]);
      abort();
    }
  }
  FREE(query);
  FREE(saindexj_new);
  FREE(saindexi_new);
  fprintf(stderr,"done\n");
#endif

  if (compressp == false) {
    fp = fopen(indexcompfile,"w");
    FWRITE_UINTS(saindex_old,2*oligospace,fp);
    fclose(fp);
  } else {
    Bitpack64_write_differential(indexptrsfile,indexcompfile,saindex_old,2*oligospace-1);
  }

  FREE(saindex_old);

  munmap((void *) SA,sa_len);
  close(sa_fd);

  return;
}


/* phi is the successor array: [0..n] */
void
Sarray_write_csa (char **csaptrfiles, char **csacompfiles, char *sasampleqfile, char *sasamplesfile, char *saindex0file,
		  char *sarrayfile, char *rankfile, Genome_T genomecomp, UINT4 genomelength, char chartable[]) {
  UINT4 *CSA, *SA, *SA_inv, sa_i;
  FILE *fp, *sa_fp, *samples_fp;
  /* FILE *csa_fp; */
  UINT4 n = genomelength, n_plus_one, ii, i, b;
  int chari, k;
  Sarrayptr_T saindexi[5], saindexj[5], saindexn, indexX;
  int sa_fd, rank_fd;
  size_t sa_len, rank_len;
  int indexsize;
  UINT4 *read_buffer, *write_buffer, ignore;
  char *queryuc_ptr;
  int csa_sampling;

  /* Write SA sampling interval */
  fp = fopen(sasampleqfile,"wb");
  csa_sampling = rint(log((double) genomelength)/log(2.0));
  fprintf(stderr,"CSA sampling: %d\n",csa_sampling);
  FWRITE_INT(csa_sampling,fp);
  fclose(fp);


  /* Determine sizes of each csa */
  SA = (UINT4 *) Access_mmap(&sa_fd,&sa_len,sarrayfile,sizeof(UINT4),/*randomp*/true);
  queryuc_ptr = (char *) CALLOC(/*querylength*/1+1,sizeof(char));

  /* A */
  oligo_nt(queryuc_ptr,/*oligo*/0,/*querylength*/1);
  sarray_search_simple(&(saindexi[0]),&(saindexj[0]),queryuc_ptr,/*querylength*/1,
		       genomecomp,SA,/*i*/1,/*j*/n,n,chartable);
  printf("A: %u..%u\n",saindexi[0],saindexj[0]);

  /* C */
  oligo_nt(queryuc_ptr,/*oligo*/1,/*querylength*/1);
  sarray_search_simple(&(saindexi[1]),&(saindexj[1]),queryuc_ptr,/*querylength*/1,
		       genomecomp,SA,/*i*/1,/*j*/n,n,chartable);
  printf("C: %u..%u\n",saindexi[1],saindexj[1]);

  /* G */
  oligo_nt(queryuc_ptr,/*oligo*/2,/*querylength*/1);
  sarray_search_simple(&(saindexi)[2],&(saindexj[2]),queryuc_ptr,/*querylength*/1,
		       genomecomp,SA,/*i*/1,/*j*/n,n,chartable);
  printf("G: %u..%u\n",saindexi[2],saindexj[2]);

  /* T */
  oligo_nt(queryuc_ptr,/*oligo*/3,/*querylength*/1);
  sarray_search_simple(&(saindexi[3]),&(saindexj[3]),queryuc_ptr,/*querylength*/1,
		       genomecomp,SA,/*i*/1,/*j*/n,n,chartable);
  printf("T: %u..%u\n",saindexi[3],saindexj[3]);

  /* X */
  saindexi[4] = saindexj[3] + 1;
  saindexj[4] = genomelength;
  printf("X: %u..%u\n",saindexi[4],saindexj[4]);

  munmap((void *) SA,sa_len);
  close(sa_fd);

  fp = fopen(saindex0file,"wb");
  FWRITE_UINT(saindexi[0],fp);
  FWRITE_UINT(saindexi[1],fp);
  FWRITE_UINT(saindexi[2],fp);
  FWRITE_UINT(saindexi[3],fp);
  FWRITE_UINT(saindexi[4],fp);

  n_plus_one = genomelength + 1;
  FWRITE_UINT(n_plus_one,fp);	/* Needed by sarray-read to find genomiclength */

  fclose(fp);


  /* Process suffix array */
  read_buffer = (UINT4 *) MALLOC(RW_BATCH * sizeof(UINT4));

  SA_inv = (UINT4 *) Access_mmap(&rank_fd,&rank_len,rankfile,sizeof(UINT4),/*randomp*/true);
  /* csa_fp = fopen(csafile,"wb");*/
  sa_fp = fopen(sarrayfile,"rb");
  samples_fp = fopen(sasamplesfile,"wb");

  CSA = (UINT4 *) MALLOC((n+1)*sizeof(UINT4));

  /* Ignore csa[0] which corresponds to end-of-string terminator */
  FREAD_UINT(&sa_i,sa_fp);
  FWRITE_UINT(sa_i,samples_fp);
  CSA[0] = genomelength;
  /* FWRITE_UINT(CSA[0],csa_fp); */

  ii = 1;
  while (ii + RW_BATCH <= n) {
    FREAD_UINTS(read_buffer,RW_BATCH,sa_fp);
    for (b = 0, i = ii; b < RW_BATCH; b++, i++) {
      if ((i % csa_sampling) == 0) {
	FWRITE_UINT(read_buffer[b],samples_fp);
      }
      CSA[i] = SA_inv[read_buffer[b] + 1];
    }
    /* FWRITE_UINTS(&(CSA[ii]),RW_BATCH,csa_fp); */
    ii += RW_BATCH;
  }
  
  /* Final partial batch */
  for (i = ii; i <= n; i++) {	/* final partial batch */
    FREAD_UINT(&sa_i,sa_fp);
    if ((i % csa_sampling) == 0) {
      FWRITE_UINT(sa_i,samples_fp);
    }
    CSA[i] = SA_inv[sa_i + 1];
    /* FWRITE_UINT(CSA[i],csa_fp);*/
  }
  /* fclose(csa_fp); */

  fclose(samples_fp);
  fclose(sa_fp);
  munmap((void *) SA_inv,rank_len);
  close(rank_fd);
  FREE(read_buffer);

  for (chari = 0; chari < 5; chari++) {
    if (saindexj[chari] < saindexi[chari]) {
      fp = fopen(csaptrfiles[chari],"wb");
      fclose(fp);
      fp = fopen(csacompfiles[chari],"wb");
      fclose(fp);
    } else {
      saindexn = saindexj[chari] - saindexi[chari] + 1;
      /* Provide (n-1) to write values [0..n] */
      Bitpack64_write_differential(csaptrfiles[chari],csacompfiles[chari],
				   &(CSA[saindexi[chari]]),saindexn-1);
    }
  }

  fprintf(stderr,"done\n");

  FREE(CSA);

  return;
}



#if 0
UINT4 *
Sarray_compute_lcp_kasai (UINT4 *SA, UINT4 n) {
  UINT4 *lcp;
  UINT4 *rank, h;
  UINT4 i, j;
  char *comma;
#ifdef DEBUG14
  UINT4 horig;
#endif

  lcp = (UINT4 *) MALLOC((n+1)*sizeof(UINT4));

  rank = (UINT4 *) CALLOC(n+1,sizeof(UINT4));
  for (i = 0; i <= n; i++) {
    rank[SA[i]] = i;
  }

#if 0
  /* Used for comparison with Manzini */
  for (i = 0; i <= n; i++) {
    printf("%u %u\n",i,rank[i]);
  }
  printf("End of Kasai\n\n");
#endif

  lcp[0] = 0;			/* -1 ? */
  h = 0;
  for (i = 0; i <= n; i++) {
    if (rank[i] > 0) {
      j = SA[rank[i] - 1];
#ifdef DEBUG14
      horig = h;
      while (i + h < n && j + h < n && s[i+h] == s[j+h]) {
	h++;
      }
      if ((h - horig) != Genome_consecutive_matches_pair(i+horig,j+horig,/*genomelength*/n)) {
	abort();
      }
#else
      h += Genome_consecutive_matches_pair(i+h,j+h,/*genomelength*/n);
#endif
      lcp[rank[i]] = h;
      if (h > 0) {
	h--;
      }
    }
    if (i % MONITOR_INTERVAL == 0) {
      comma = Genomicpos_commafmt(i);
      fprintf(stderr,"Computing lcp index %s\n",comma);
      FREE(comma);
    }
  }

  FREE(rank);

  return lcp;
}
#endif


#if 0
/* Puts rank in file, to save on memory */
/* Rank file contains the inverse suffix array, needed to compute the compressed suffix array */
UINT4 *
Sarray_compute_lcp (char *rankfile, UINT4 *SA, UINT4 n) {
  UINT4 *lcp;
  UINT4 *rank, rank_i, h;
  UINT4 i, j;
  char *comma;

  FILE *fp;
#ifdef DEBUG14
  UINT4 horig;
#endif

  /* Compute rank and store in temporary file */
  rank = (UINT4 *) MALLOC((n+1)*sizeof(UINT4));
  for (i = 0; i <= n; i++) {
    rank[SA[i]] = i;
  }

  fprintf(stderr,"Writing temporary file %s...",rankfile);
  fp = fopen(rankfile,"w");
  for (i = 0; i + FWRITE_BATCH <= n; i += FWRITE_BATCH) {
    fwrite((void *) &(rank[i]),sizeof(UINT4),FWRITE_BATCH,fp);
  }

  if (i <= n) {
    fwrite((void *) &(rank[i]),sizeof(UINT4),n - i + 1,fp);
  }
  fclose(fp);
  FREE(rank);
  fprintf(stderr,"done\n");

  
  /* Now allocate memory for lcp */
  fp = fopen(rankfile,"r");

  lcp = (UINT4 *) MALLOC((n+1)*sizeof(UINT4));

  lcp[0] = 0;			/* -1 ? */
  h = 0;
  for (i = 0; i <= n; i++) {
    FREAD_UINT(&rank_i,fp);
    if (rank_i > 0) {
      j = SA[rank_i - 1];

      h += Genome_consecutive_matches_pair(i+h,j+h,/*genomelength*/n);

      lcp[rank_i] = h;
      if (h > 0) {
	h--;
      }
    }
    if (i % MONITOR_INTERVAL == 0) {
      comma = Genomicpos_commafmt(i);
      fprintf(stderr,"Computing lcp index %s\n",comma);
      FREE(comma);
    }
  }
  fclose(fp);

  remove(rankfile);

  return lcp;
}
#endif


/* Puts rank and permuted suffix array in file, to save on memory even further */
/* Rank file is the same as the inverted suffix array, needed to compute the compressed suffix array */
UINT4 *
Sarray_compute_lcp (char *rankfile, char *permuted_sarray_file, char *sarrayfile, UINT4 n) {
  UINT4 *lcp;
  UINT4 *SA, SA_i, zero = 0;
  UINT4 *rank, rank_i, h;
  UINT4 i, ii, b, j;
  char *comma;
  UINT4 *read_buffer_1, *read_buffer_2, *write_buffer;
  void *p;

  int sa_fd;
  size_t sa_len;
  FILE *fp, *permsa_fp;


  read_buffer_1 = (UINT4 *) MALLOC(RW_BATCH * sizeof(UINT4));

  /* Compute rank */
  fp = fopen(sarrayfile,"rb");
  rank = (UINT4 *) MALLOC((n+1)*sizeof(UINT4));

  for (ii = 0; ii + RW_BATCH <= n; ii += RW_BATCH) {
    FREAD_UINTS(read_buffer_1,RW_BATCH,fp);
    for (b = 0, i = ii; b < RW_BATCH; b++, i++) {
      rank[read_buffer_1[b]] = i;      /* rank[SA_i] = i; */
    }
    if (ii % MONITOR_INTERVAL == 0) {
      comma = Genomicpos_commafmt(ii);
      fprintf(stderr,"Computing rank %s\n",comma);
      FREE(comma);
    }
  }
  for (i = ii; i <= n; i++) {	/* final partial batch */
    FREAD_UINT(&SA_i,fp);
    rank[SA_i] = i;
  }

  fclose(fp);			/* sarrayfile */


  /* Store rank in temporary file */
  fprintf(stderr,"Writing temporary file for rank...");
  fp = fopen(rankfile,"wb");
  for (ii = 0; ii + RW_BATCH <= n; ii += RW_BATCH) {
    p = (void *) &(rank[ii]);
    FWRITE_UINTS(p,RW_BATCH,fp);
  }
  if (ii <= n) {
    p = (void *) &(rank[ii]);
    FWRITE_UINTS(p,n - ii + 1,fp);
  }
  fclose(fp);			/* rankfile */
  FREE(rank);
  fprintf(stderr,"done\n");

  
  /* Write permuted sarray */
  fprintf(stderr,"Writing temporary file for permuted sarray...");
  write_buffer = (UINT4 *) MALLOC(RW_BATCH * sizeof(UINT4));
  fp = fopen(rankfile,"rb");
  permsa_fp = fopen(permuted_sarray_file,"wb");
  SA = (UINT4 *) Access_mmap(&sa_fd,&sa_len,sarrayfile,sizeof(UINT4),/*randomp*/false);

  for (ii = 0; ii + RW_BATCH <= n; ii += RW_BATCH) {
    FREAD_UINTS(read_buffer_1,RW_BATCH,fp);
    for (b = 0, i = ii; b < RW_BATCH; b++, i++) {
      rank_i = read_buffer_1[b];
      if (rank_i > 0) {
#ifdef WORDS_BIGENDIAN
	write_buffer[b] = Bigendian_convert_uint(SA[rank_i - 1]);
#else
	write_buffer[b] = SA[rank_i - 1];
#endif
      } else {
	write_buffer[b] = 0;	/* Will be ignored */
      }
    }
    FWRITE_UINTS(write_buffer,RW_BATCH,permsa_fp);
  }
  for (i = ii; i <= n; i++) {	/* final partial batch */
    FREAD_UINT(&rank_i,fp);
    if (rank_i > 0) {
#ifdef WORDS_BIGENDIAN
      FWRITE_UINT(Bigendian_convert_uint(SA[rank_i - 1]),permsa_fp);
#else
      FWRITE_UINT(SA[rank_i - 1],permsa_fp);
#endif
    } else {
      FWRITE_UINT(zero,permsa_fp); /* Will be ignored */
    }
  }

  munmap((void *) SA,sa_len);
  close(sa_fd);
  fclose(permsa_fp);		/* permuted_sarray_file */
  fclose(fp);			/* rankfile */
  FREE(write_buffer);
  fprintf(stderr,"done\n");


  /* Now allocate memory for lcp and compute */
  read_buffer_2 = (UINT4 *) MALLOC(RW_BATCH * sizeof(UINT4));
  fp = fopen(rankfile,"rb");
  permsa_fp = fopen(permuted_sarray_file,"rb");

  lcp = (UINT4 *) MALLOC((n+1)*sizeof(UINT4));

  lcp[0] = 0;			/* -1 ? */
  h = 0;

  for (ii = 0; ii + RW_BATCH <= n; ii += RW_BATCH) {
    FREAD_UINTS(read_buffer_1,RW_BATCH,fp);
    FREAD_UINTS(read_buffer_2,RW_BATCH,permsa_fp);
    for (b = 0, i = ii; b < RW_BATCH; b++, i++) {
      rank_i = read_buffer_1[b];
      j = read_buffer_2[b];	/* j = SA[rank_i - 1] */
      if (rank_i > 0) {
	h += Genome_consecutive_matches_pair(i+h,j+h,/*genomelength*/n);
	lcp[rank_i] = h;
	if (h > 0) {
	  h--;
	}
      }
    }

    if (ii % MONITOR_INTERVAL == 0) {
      comma = Genomicpos_commafmt(ii);
      fprintf(stderr,"Computing lcp index %s\n",comma);
      FREE(comma);
    }
  }

  for (i = ii; i <= n; i++) {	/* final partial batch */
    FREAD_UINT(&rank_i,fp);
    FREAD_UINT(&j,permsa_fp);  /* j = SA[rank_i - 1] */
    if (rank_i > 0) {
      h += Genome_consecutive_matches_pair(i+h,j+h,/*genomelength*/n);
      lcp[rank_i] = h;
      if (h > 0) {
	h--;
      }
    }
  }

  fclose(permsa_fp);		/* permuted_sarray_file */
  fclose(fp);			/* rankfile */
  FREE(read_buffer_2);
  FREE(read_buffer_1);

  remove(permuted_sarray_file);
#ifndef USE_CSA
  remove(rankfile);
#endif

  return lcp;
}


#if 0
/* Based on Manzini, 2004 */
/* eos_pos: end-of-string pos? */
static UINT4
compute_next_rank (UINT4 *next_rank, UINT4 *SA, Genome_T genomecomp, UINT4 n, char *chartable) {
  UINT4 eos_pos;
  Univcoord_T nACGT, na, nc, ng, nt;
  UINT4 i, j;
  UINT4 count[5];
  int numeric[128];
  unsigned char c;
  char *comma;

  nACGT = Genome_ntcounts(&na,&nc,&ng,&nt,genomecomp,/*left*/0,/*length*/n);
  fprintf(stderr,"Genome content: A %u, C %u, G %u, T %u, other %u\n",na,nc,ng,nt,n - nACGT);
  count[0] = 0;
  count[1] = na;
  count[2] = na + nc;
  count[3] = na + nc + ng;
  count[4] = na + nc + ng + nt;

  for (c = 0; c < 128; c++) {
    numeric[c] = 4;
  }
  numeric['A'] = 0;
  numeric['C'] = 1;
  numeric['G'] = 2;
  numeric['T'] = 3;

  c = Genome_get_char_lex(genomecomp,/*pos*/n,n,chartable);
  j = ++count[numeric[c]];
  next_rank[j] = 0;

  for (i = 1; i <= n; i++) {
    if (SA[i] == 1) {
      eos_pos = i;
    } else {
      c = Genome_get_char_lex(genomecomp,/*pos*/SA[i] - 1,n,chartable);
      j = ++count[numeric[c]];
      next_rank[j] = i;
    }

    if (i % MONITOR_INTERVAL == 0) {
      comma = Genomicpos_commafmt(i);
      fprintf(stderr,"Computing rank %s\n",comma);
      FREE(comma);
    }
  }

  for (i = 0; i <= n; i++) {
    printf("%u %u\n",i,next_rank[i]);
  }
  printf("Returning %u\n",eos_pos);
  exit(0);


  return eos_pos;
}


/* Use rank_i instead of k, and next_rank_i instead of nextk */
/* Buggy: Result differs from that obtained by Kasai procedure */
UINT4 *
Sarray_compute_lcp_manzini (UINT4 *SA, Genome_T genomecomp, UINT4 n, char *chartable) {
  UINT4 *lcp;
  UINT4 h, i, j, rank_i, next_rank_i;
  char *comma;

  lcp = (UINT4 *) MALLOC((n+1)*sizeof(UINT4));
  rank_i = compute_next_rank(lcp,SA,genomecomp,n,chartable); /* Re-use lcp for next_rank */
  h = 0;

  for (i = 1; i <= n; i++) {
    next_rank_i = lcp[rank_i];
    if (rank_i > 0) {
      j = SA[rank_i - 1];
      h += Genome_consecutive_matches_pair(i+h,j+h,/*genomelength*/n);
      lcp[rank_i] = h;
      if (h > 0) {
	h--;
      }
    }
    rank_i = next_rank_i;

    if (i % MONITOR_INTERVAL == 0) {
      comma = Genomicpos_commafmt(i);
      fprintf(stderr,"Computing lcp index %s\n",comma);
      FREE(comma);
    }
  }

  return lcp;
}
#endif


/* Used by cmetindex and atoiindex */
UINT4 *
Sarray_compute_lcp_from_genome (UINT4 *SA, unsigned char *gbuffer, UINT4 n) {
  UINT4 *lcp;
  UINT4 *rank, h;
  UINT4 i, j;
  char *comma;
  /* UINT4 horig; */

  lcp = (UINT4 *) MALLOC((n+1)*sizeof(UINT4));

  rank = (UINT4 *) CALLOC(n+1,sizeof(UINT4));
  for (i = 0; i <= n; i++) {
    rank[SA[i]] = i;
  }

  lcp[0] = 0;			/* -1 ? */
  h = 0;
  for (i = 0; i <= n; i++) {
    if (rank[i] > 0) {
      j = SA[rank[i] - 1];
      /* horig = h; */
      while (i + h < n && j + h < n && gbuffer[i+h] == gbuffer[j+h]) {
	h++;
      }
      lcp[rank[i]] = h;
      if (h > 0) {
	h--;
      }
    }
    if (i % MONITOR_INTERVAL == 0) {
      comma = Genomicpos_commafmt(i);
      fprintf(stderr,"Computing lcp index %s\n",comma);
      FREE(comma);
    }
  }

  FREE(rank);

  return lcp;
}


/* Computes permuted lcp (Karkkainen, CPM 2009) */
/* Two methods for storing plcp as a cumulative sum: (1) Cum value at k is
   sum_{i=1}^k (plcp[i] + 1 - plcp[i-1]), or more simply (2) plcp[k] + k */
static void
compute_plcp (UINT4 *plcp, UINT4 *SA, UINT4 n) {
  UINT4 *phi, h;
  UINT4 i, j;
  char *comma;

  phi = plcp;			/* Use space allocated for plcp */

  fprintf(stderr,"Inverting suffix array (via mmap)...");
  for (i = 1; i <= n; i++) {
    phi[SA[i]] = SA[i-1];
  }
  fprintf(stderr,"done\n");
  /* Note that phi[n] is not assigned, because SA[i] == n for i == 0, and we don't look up i == 0 */


  h = 0;
  for (i = 0; i < n; i++) {
    j = phi[i];			/* To be overwritten by plcp[i] */
    h += Genome_consecutive_matches_pair(i+h,j+h,/*genomelength*/n);
    plcp[i] = h;				 /* overwrites phi[i] */
    if (h > 0) {
      h--;
    }

    if (i % MONITOR_INTERVAL == 0) {
      comma = Genomicpos_commafmt(i);
      fprintf(stderr,"Computing permuted lcp index %s\n",comma);
      FREE(comma);
    }
  }

#if 0
  /* This makes lcp[0] = -1, because lcp[0] = plcp[SA[0]] = plcp[n] = -1 */
  plcp[n] = -1;
#else
  /* This makes lcp[0] = 0, because lcp[0] = plcp[SA[0]] = plcp[n] = 0 */
  plcp[n] = 0;
#endif

  return;
}



#if 0
static int
get_all_children (bool *filledp, Sarrayptr_T *l, Sarrayptr_T *r, Sarrayptr_T i, Sarrayptr_T j,
		  UINT4 *child, UINT8 *nextp, Genome_T genomecomp, UINT4 *SA, UINT4 *plcpptrs, UINT4 *plcpcomp,
		  int indexsize, UINT4 n) {
  int noccupied = 0;
  UINT4 up, nextl;
  Sarrayptr_T sa_nextl;
  UINT4 lcp_whole;
  UINT4 pos;
  char c;

  /* Test for child[j] being up: lcp[j] > lcp[j+1] */
  up = child[j];		/* childtab[j+1].up */
  if (i < up && up <= j) {
    nextl = up;
  } else {
    nextl = child[i];	/* down */
  }
  sa_nextl = SA[nextl];
  lcp_whole = Bitpack64_read_one(sa_nextl,plcpptrs,plcpcomp) - sa_nextl;

  if (lcp_whole != (UINT4) (indexsize - 1)) {
    /* Not at desired level, so exit this procedure */
    return 0;
  } else {
    debug1(printf("Filling children for lcp-interval %u..%u with lcp_whole %d\n",i,j,lcp_whole));
  }

  pos = SA[i] + lcp_whole;
  c = Genome_get_char_lex(genomecomp,pos,n);

  debug1(printf("For char %c, creating interval %u..%u\n",c,i,nextl-1));
  switch (c) {
  case 'A': l[0] = i; r[0] = nextl - 1; filledp[0] = true; noccupied++; break;
  case 'C': l[1] = i; r[1] = nextl - 1; filledp[1] = true; noccupied++; break;
  case 'G': l[2] = i; r[2] = nextl - 1; filledp[2] = true; noccupied++; break;
  case 'T': l[3] = i; r[3] = nextl - 1; filledp[3] = true; noccupied++; break;
  }

  
  /* Test for child[i] being down: lcp[child[i]] > lcp[i] */
  /* Test for child[i] being next_lindex: lcp[child[i]] == lcp[i] */
  while (get_bit(nextp,nextl) != 0) {
    pos = SA[nextl] + lcp_whole;
    c = Genome_get_char_lex(genomecomp,pos,n,chartable);

    debug1(printf("For char %c, creating interval %u..%u\n",c,nextl,child[nextl]-1));
    switch (c) {
    case 'A': l[0] = nextl; r[0] = child[nextl] - 1; filledp[0] = true; noccupied++; break;
    case 'C': l[1] = nextl; r[1] = child[nextl] - 1; filledp[1] = true; noccupied++; break;
    case 'G': l[2] = nextl; r[2] = child[nextl] - 1; filledp[2] = true; noccupied++; break;
    case 'T': l[3] = nextl; r[3] = child[nextl] - 1; filledp[3] = true; noccupied++; break;
    }
    
    nextl = child[nextl];
  }

  pos = SA[nextl] + lcp_whole;
  c = Genome_get_char_lex(genomecomp,pos,n,chartable);

  debug1(printf("For char %c, creating interval %u..%u\n",c,nextl,j));
  switch (c) {
  case 'A': l[0] = nextl; r[0] = j; filledp[0] = true; noccupied++; break;
  case 'C': l[1] = nextl; r[1] = j; filledp[1] = true; noccupied++; break;
  case 'G': l[2] = nextl; r[2] = j; filledp[2] = true; noccupied++; break;
  case 'T': l[3] = nextl; r[3] = j; filledp[3] = true; noccupied++; break;
  }

  return noccupied;
}
#endif



/* Uses permuted lcp for speed and reduced memory usage */
void
Sarray_write_plcp (char *plcpptrsfile, char *plcpcompfile, UINT4 *SA, UINT4 genomelength) {
  UINT4 *plcp;
  UINT4 *ramp, *p;

  UINT4 n = genomelength, i;
  UINT4 ii;
  FILE *fp;

  plcp = (UINT4 *) MALLOC((n+1)*sizeof(UINT4));
  ramp = plcp;

  compute_plcp(plcp,SA,n);

#ifdef USE_CUMDEV
  /* Compute deviations from downward ramp */
  for (i = n; i >= 1; i--) {
    dev[i] = (plcp[i] + 1) - plcp[i-1];
  }
  /* dev[0] = plcp[0]; */
#else
  for (i = 0; i <= n; i++) {
    ramp[i] = plcp[i] + i;
  }
#endif

  fprintf(stderr,"Writing permuted lcp file...");
  /* Provide n to write values [0..n] */

#if 0
  /* Print plcp as an array */
  fp = fopen("plcp","wb");
  for (ii = 0; ii + RW_BATCH <= n; ii += RW_BATCH) {
    p = (void *) &(ramp[ii]);
    FWRITE_UINTS(p,RW_BATCH,fp);
  }
  if (ii <= n) {
    p = (void *) &(ramp[ii]);
    FWRITE_UINTS(p,n - ii + 1,fp);
  }
  fclose(fp);
#else
  Bitpack64_write_differential(plcpptrsfile,plcpcompfile,ramp,n);
#endif

  fprintf(stderr,"done\n");

  FREE(plcp);

  return;
}



/* Without encoding, would be child[index] = value */
#define encode_up(child,index,value) child[index-1] = (index - 1) - (value)
#define encode_down(child,index,value) child[index] = (value) - 1 - (index)
#define encode_next(child,index,value) child[index] = (value) - 1 - (index)

/* For child[index+1].up, just calling child[index] */
#define decode_up(child_i,index) index - child_i
#define decode_down(child_i,index) child_i + index + 1
#define decode_next(child_i,index) child_i + index + 1


#if 0
static UINT4 *
make_child_twopass (UINT8 **nextp, UINT4 *nbytes, UINT4 *SA, UINT4 *plcpptrs, UINT4 *plcpcomp, UINT4 n) {
  UINT4 *child;
  UINT4 lastindex, i;
  Uintlist_T lcpstack, indexstack;
  UINT4 sa_i, lcp_i, lcp_lastindex;
  char *comma;

  *nbytes = ((n+1) + WORDSIZE-1)/WORDSIZE;
  *nextp = (UINT8 *) CALLOC(*nbytes,sizeof(UINT8));

#if 0
  child = (UINT4 *) MALLOC((n+1) * sizeof(UINT4));
  for (i = 0; i <= n; i++) {
    child[i] = -1U;
  }
#else
  child = (UINT4 *) CALLOC(n+1,sizeof(UINT4));
#endif

  
  /* Because we sort suffixes with $ < rest of alphabet, we never use
     the entry at 0, where SA[0] = n and lcp[0] = 0 */

   /* Algorithm 6.2: Compute up and down values */
   lastindex = 0;

   fprintf(stderr,"Computing child up/down index 0\n");
   i = 1;
   indexstack = Uintlist_push(NULL,i);

   sa_i = SA[i];
   lcp_i = Bitpack64_read_one(sa_i,plcpptrs,plcpcomp) - sa_i;
   lcpstack = Uintlist_push(NULL,lcp_i);

   for (i = 2; i <= n; i++) {
     sa_i = SA[i];
     lcp_i = Bitpack64_read_one(sa_i,plcpptrs,plcpcomp) - sa_i;
     while (lcp_i < Uintlist_head(lcpstack)) {
       lcpstack = Uintlist_pop(lcpstack,&lcp_lastindex);
       indexstack = Uintlist_pop(indexstack,&lastindex);
       if (lcp_i <= Uintlist_head(lcpstack) && Uintlist_head(lcpstack) != lcp_lastindex) {
 #ifdef NO_ENCODING
	 child[Uintlist_head(indexstack)] = lastindex; /* down */
 #else
	 encode_down(child,Uintlist_head(indexstack),lastindex);
 #endif
       }
     }
     /* Now lcp[i] >= lcp[stack->first] holds */
     if (lastindex != 0) {
 #ifdef NO_ENCODING
       child[i-1] = lastindex;	/* up */
 #else
       encode_up(child,i,lastindex);
 #endif
       lastindex = 0;
     }
     indexstack = Uintlist_push(indexstack,i);
     lcpstack = Uintlist_push(lcpstack,lcp_i);

     if (i % MONITOR_INTERVAL == 0) {
       comma = Genomicpos_commafmt(i);
       fprintf(stderr,"Computing child up/down index %s\n",comma);
       FREE(comma);
     }
   }

   /* Handle end of suffix array */
   lcp_i = 0;
   while (lcp_i < Uintlist_head(lcpstack)) {
     lcpstack = Uintlist_pop(lcpstack,&lcp_lastindex);
     indexstack = Uintlist_pop(indexstack,&lastindex);
     if (lcp_i <= Uintlist_head(lcpstack) && Uintlist_head(lcpstack) != lcp_lastindex) {
 #ifdef NO_ENCODING
       child[Uintlist_head(indexstack)] = lastindex; /* down */
 #else
       encode_down(child,Uintlist_head(indexstack),lastindex);
 #endif
     }
   }
   if (lastindex != 0) {
 #ifdef NO_ENCODING
     child[i-1] = lastindex;	/* up */
 #else
     encode_up(child,i,lastindex);
 #endif
     /* lastindex = 0; */
   }

   Uintlist_free(&lcpstack);
   Uintlist_free(&indexstack);


   /* Algorithm 6.5: Compute next l-index values */
   fprintf(stderr,"Computing child next index 0\n");
   i = 1;
   indexstack = Uintlist_push(NULL,i);

   sa_i = SA[i];
   lcp_i = Bitpack64_read_one(sa_i,plcpptrs,plcpcomp) - sa_i;
   lcpstack = Uintlist_push(NULL,lcp_i);

   for (i = 2; i <= n; i++) {
     sa_i = SA[i];
     lcp_i = Bitpack64_read_one(sa_i,plcpptrs,plcpcomp) - sa_i;
     while (lcp_i < Uintlist_head(lcpstack)) {
       lcpstack = Uintlist_pop(lcpstack,&lcp_lastindex);
       indexstack = Uintlist_pop(indexstack,&lastindex);
     }
     if (lcp_i == Uintlist_head(lcpstack)) {
       lcpstack = Uintlist_pop(lcpstack,&lcp_lastindex);
       indexstack = Uintlist_pop(indexstack,&lastindex);
 #ifdef NO_ENCODING
       child[lastindex] = i;
 #else
       encode_next(child,lastindex,i);
 #endif
       set_bit(*nextp,lastindex);
     }
     indexstack = Uintlist_push(indexstack,i);
     lcpstack = Uintlist_push(lcpstack,lcp_i);

     if (i % MONITOR_INTERVAL == 0) {
       comma = Genomicpos_commafmt(i);
       fprintf(stderr,"Computing child next index %s\n",comma);
       FREE(comma);
     }
   }

   Uintlist_free(&lcpstack);
   Uintlist_free(&indexstack);

   return child;
}
#endif


#if 0
/* Reads SA one element at a time */
/* For adjoining chars, need to store 15 possibilities, or one nibble */
/* Possibilities: $a, $c, $g, $t, $x, ac, ag, at, ax, cg, ct, cx, gt, gx, tx */
unsigned char *
Sarray_discriminating_chars (UINT4 *nbytes, char *sarrayfile, Genome_T genome,
			     unsigned char *lcp_bytes, UINT4 *lcp_guide, UINT4 *lcp_exceptions, int guide_interval,
			     UINT4 n, char *chartable) {
  unsigned char *discrim_chars;
  char char_before, char_at;
  UINT4 i;
  UINT4 lcp_i;

  FILE *fp;
  UINT4 SA_i_minus_1, SA_i;


  *nbytes = ((n+1) + 1)/2;
  discrim_chars = (unsigned char *) CALLOC(*nbytes,sizeof(unsigned char));

  fp = fopen(sarrayfile,"r");
  FREAD_UINT(&SA_i_minus_1,fp);

  for (i = 1; i <= n; i++) {
    FREAD_UINT(&SA_i,fp);

    lcp_i = Bytecoding_read_wguide(i,lcp_bytes,lcp_guide,lcp_exceptions,/*lcp_guide_interval*/1024);
    char_before = Genome_get_char_lex(genome,/*left: SA[i-1]*/SA_i_minus_1 + lcp_i,/*genomelength*/n,chartable);
    char_at = Genome_get_char_lex(genome,/*left: SA[i]*/SA_i + lcp_i,/*genomelength*/n,chartable);
    debug4(printf("i = %u, SA = %u and %u, and lcp_i = %u => %c %c\n",
		  i,SA_i_minus_1,SA_i,lcp_i,char_before == 0 ? '$' : char_before,char_at));

    if (i % 2 == 0) {
      /* Even, put into low nibble of byte */
      switch (char_before) {
      case 0:
	switch (char_at) {
	case 'A': discrim_chars[i/2] |= 0x01; break;
	case 'C': discrim_chars[i/2] |= 0x02; break;
	case 'G': discrim_chars[i/2] |= 0x03; break;
	case 'T': discrim_chars[i/2] |= 0x04; break;
	case 'X': discrim_chars[i/2] |= 0x05; break;
	default: abort();
	}
	break;

      case 'A':
	switch (char_at) {
	case 'C': discrim_chars[i/2] |= 0x06; break;
	case 'G': discrim_chars[i/2] |= 0x07; break;
	case 'T': discrim_chars[i/2] |= 0x08; break;
	case 'X': discrim_chars[i/2] |= 0x09; break;
	default: abort();
	}
	break;

      case 'C':
	switch (char_at) {
	case 'G': discrim_chars[i/2] |= 0x0A; break;
	case 'T': discrim_chars[i/2] |= 0x0B; break;
	case 'X': discrim_chars[i/2] |= 0x0C; break;
	default: abort();
	}
	break;

      case 'G':
	switch (char_at) {
	case 'T': discrim_chars[i/2] |= 0x0D; break;
	case 'X': discrim_chars[i/2] |= 0x0E; break;
	default: abort();
	}
	break;

      case 'T':
	switch (char_at) {
	case 'X': discrim_chars[i/2] |= 0x0F; break;
	default: abort();
	}
	break;
      }

    } else {
      /* Odd, put into high nibble of byte */
      switch (char_before) {
      case 0:
	switch (char_at) {
	case 'A': discrim_chars[i/2] |= 0x10; break;
	case 'C': discrim_chars[i/2] |= 0x20; break;
	case 'G': discrim_chars[i/2] |= 0x30; break;
	case 'T': discrim_chars[i/2] |= 0x40; break;
	case 'X': discrim_chars[i/2] |= 0x50; break;
	default: abort();
	}
	break;

      case 'A':
	switch (char_at) {
	case 'C': discrim_chars[i/2] |= 0x60; break;
	case 'G': discrim_chars[i/2] |= 0x70; break;
	case 'T': discrim_chars[i/2] |= 0x80; break;
	case 'X': discrim_chars[i/2] |= 0x90; break;
	default: abort();
	}
	break;

      case 'C':
	switch (char_at) {
	case 'G': discrim_chars[i/2] |= 0xA0; break;
	case 'T': discrim_chars[i/2] |= 0xB0; break;
	case 'X': discrim_chars[i/2] |= 0xC0; break;
	default: abort();
	}
	break;

      case 'G':
	switch (char_at) {
	case 'T': discrim_chars[i/2] |= 0xD0; break;
	case 'X': discrim_chars[i/2] |= 0xE0; break;
	default: abort();
	}
	break;

      case 'T':
	switch (char_at) {
	case 'X': discrim_chars[i/2] |= 0xF0; break;
	default: abort();
	}
	break;
      }
    }

    SA_i_minus_1 = SA_i;
  }

  fclose(fp);

  return discrim_chars;
}
#endif


/* Reads SA in batches for faster I/O */
/* For adjoining chars, need to store 15 possibilities, or one nibble */
/* Possibilities: $a, $c, $g, $t, $x, ac, ag, at, ax, cg, ct, cx, gt, gx, tx */
unsigned char *
Sarray_discriminating_chars (UINT4 *nbytes, char *sarrayfile, Genome_T genome,
			     unsigned char *lcp_bytes, UINT4 *lcp_guide, UINT4 *lcp_exceptions, int guide_interval,
			     UINT4 n, char *chartable) {
  unsigned char *discrim_chars;
  char char_before, char_at;
  UINT4 i, ii, b;
  UINT4 lcp_i;
  char *comma;

  FILE *fp;
  UINT4 *read_buffer;
  UINT4 SA_i_minus_1, SA_i;


  *nbytes = ((n+1) + 1)/2;
  discrim_chars = (unsigned char *) CALLOC(*nbytes,sizeof(unsigned char));


  read_buffer = (UINT4 *) MALLOC(RW_BATCH * sizeof(UINT4));

  fp = fopen(sarrayfile,"rb");
  FREAD_UINT(&SA_i_minus_1,fp);	/* Initializes SA[0] */

  for (ii = 1; ii + RW_BATCH <= n; ii += RW_BATCH) {
    FREAD_UINTS(read_buffer,RW_BATCH,fp);

    for (b = 0, i = ii; b < RW_BATCH; b++, i++) {
      SA_i = read_buffer[b];

      lcp_i = Bytecoding_read_wguide(i,lcp_bytes,lcp_guide,lcp_exceptions,/*lcp_guide_interval*/1024);
      char_before = Genome_get_char_lex(genome,/*left: SA[i-1]*/SA_i_minus_1 + lcp_i,/*genomelength*/n,chartable);
      char_at = Genome_get_char_lex(genome,/*left: SA[i]*/SA_i + lcp_i,/*genomelength*/n,chartable);
      debug4(printf("i = %u, SA = %u and %u, and lcp_i = %u => %c %c\n",
		    i,SA_i_minus_1,SA_i,lcp_i,char_before == 0 ? '$' : char_before,char_at));

      if (i % 2 == 0) {
	/* Even, put into low nibble of byte */
	switch (char_before) {
	case 0:
	  switch (char_at) {
	  case 'A': discrim_chars[i/2] |= 0x01; break;
	  case 'C': discrim_chars[i/2] |= 0x02; break;
	  case 'G': discrim_chars[i/2] |= 0x03; break;
	  case 'T': discrim_chars[i/2] |= 0x04; break;
	  case 'X': discrim_chars[i/2] |= 0x05; break;
	  default: abort();
	  }
	  break;

	case 'A':
	  switch (char_at) {
	  case 'C': discrim_chars[i/2] |= 0x06; break;
	  case 'G': discrim_chars[i/2] |= 0x07; break;
	  case 'T': discrim_chars[i/2] |= 0x08; break;
	  case 'X': discrim_chars[i/2] |= 0x09; break;
	  default: abort();
	  }
	  break;

	case 'C':
	  switch (char_at) {
	  case 'G': discrim_chars[i/2] |= 0x0A; break;
	  case 'T': discrim_chars[i/2] |= 0x0B; break;
	  case 'X': discrim_chars[i/2] |= 0x0C; break;
	  default: abort();
	  }
	  break;

	case 'G':
	  switch (char_at) {
	  case 'T': discrim_chars[i/2] |= 0x0D; break;
	  case 'X': discrim_chars[i/2] |= 0x0E; break;
	  default: abort();
	  }
	  break;

	case 'T':
	  switch (char_at) {
	  case 'X': discrim_chars[i/2] |= 0x0F; break;
	  default: abort();
	  }
	  break;
	}

      } else {
	/* Odd, put into high nibble of byte */
	switch (char_before) {
	case 0:
	  switch (char_at) {
	  case 'A': discrim_chars[i/2] |= 0x10; break;
	  case 'C': discrim_chars[i/2] |= 0x20; break;
	  case 'G': discrim_chars[i/2] |= 0x30; break;
	  case 'T': discrim_chars[i/2] |= 0x40; break;
	  case 'X': discrim_chars[i/2] |= 0x50; break;
	  default: abort();
	  }
	  break;

	case 'A':
	  switch (char_at) {
	  case 'C': discrim_chars[i/2] |= 0x60; break;
	  case 'G': discrim_chars[i/2] |= 0x70; break;
	  case 'T': discrim_chars[i/2] |= 0x80; break;
	  case 'X': discrim_chars[i/2] |= 0x90; break;
	  default: abort();
	  }
	  break;

	case 'C':
	  switch (char_at) {
	  case 'G': discrim_chars[i/2] |= 0xA0; break;
	  case 'T': discrim_chars[i/2] |= 0xB0; break;
	  case 'X': discrim_chars[i/2] |= 0xC0; break;
	  default: abort();
	  }
	  break;

	case 'G':
	  switch (char_at) {
	  case 'T': discrim_chars[i/2] |= 0xD0; break;
	  case 'X': discrim_chars[i/2] |= 0xE0; break;
	  default: abort();
	  }
	  break;

	case 'T':
	  switch (char_at) {
	  case 'X': discrim_chars[i/2] |= 0xF0; break;
	  default: abort();
	  }
	  break;
	}
      }

      SA_i_minus_1 = SA_i;
    }

    /* Need (ii - 1) because we start with ii = 1 */
    if ((ii - 1) % MONITOR_INTERVAL == 0) {
      comma = Genomicpos_commafmt(ii-1);
      fprintf(stderr,"Computing DC array %s\n",comma);
      FREE(comma);
    }
  }

  for (i = ii; i <= n; i++) {
    FREAD_UINT(&SA_i,fp);

    lcp_i = Bytecoding_read_wguide(i,lcp_bytes,lcp_guide,lcp_exceptions,/*lcp_guide_interval*/1024);
    char_before = Genome_get_char_lex(genome,/*left: SA[i-1]*/SA_i_minus_1 + lcp_i,/*genomelength*/n,chartable);
    char_at = Genome_get_char_lex(genome,/*left: SA[i]*/SA_i + lcp_i,/*genomelength*/n,chartable);
    debug4(printf("i = %u, SA = %u and %u, and lcp_i = %u => %c %c\n",
		  i,SA_i_minus_1,SA_i,lcp_i,char_before == 0 ? '$' : char_before,char_at));

    if (i % 2 == 0) {
      /* Even, put into low nibble of byte */
      switch (char_before) {
      case 0:
	switch (char_at) {
	case 'A': discrim_chars[i/2] |= 0x01; break;
	case 'C': discrim_chars[i/2] |= 0x02; break;
	case 'G': discrim_chars[i/2] |= 0x03; break;
	case 'T': discrim_chars[i/2] |= 0x04; break;
	case 'X': discrim_chars[i/2] |= 0x05; break;
	default: abort();
	}
	break;

      case 'A':
	switch (char_at) {
	case 'C': discrim_chars[i/2] |= 0x06; break;
	case 'G': discrim_chars[i/2] |= 0x07; break;
	case 'T': discrim_chars[i/2] |= 0x08; break;
	case 'X': discrim_chars[i/2] |= 0x09; break;
	default: abort();
	}
	break;

      case 'C':
	switch (char_at) {
	case 'G': discrim_chars[i/2] |= 0x0A; break;
	case 'T': discrim_chars[i/2] |= 0x0B; break;
	case 'X': discrim_chars[i/2] |= 0x0C; break;
	default: abort();
	}
	break;

      case 'G':
	switch (char_at) {
	case 'T': discrim_chars[i/2] |= 0x0D; break;
	case 'X': discrim_chars[i/2] |= 0x0E; break;
	default: abort();
	}
	break;

      case 'T':
	switch (char_at) {
	case 'X': discrim_chars[i/2] |= 0x0F; break;
	default: abort();
	}
	break;
      }

    } else {
      /* Odd, put into high nibble of byte */
      switch (char_before) {
      case 0:
	switch (char_at) {
	case 'A': discrim_chars[i/2] |= 0x10; break;
	case 'C': discrim_chars[i/2] |= 0x20; break;
	case 'G': discrim_chars[i/2] |= 0x30; break;
	case 'T': discrim_chars[i/2] |= 0x40; break;
	case 'X': discrim_chars[i/2] |= 0x50; break;
	default: abort();
	}
	break;

      case 'A':
	switch (char_at) {
	case 'C': discrim_chars[i/2] |= 0x60; break;
	case 'G': discrim_chars[i/2] |= 0x70; break;
	case 'T': discrim_chars[i/2] |= 0x80; break;
	case 'X': discrim_chars[i/2] |= 0x90; break;
	default: abort();
	}
	break;

      case 'C':
	switch (char_at) {
	case 'G': discrim_chars[i/2] |= 0xA0; break;
	case 'T': discrim_chars[i/2] |= 0xB0; break;
	case 'X': discrim_chars[i/2] |= 0xC0; break;
	default: abort();
	}
	break;

      case 'G':
	switch (char_at) {
	case 'T': discrim_chars[i/2] |= 0xD0; break;
	case 'X': discrim_chars[i/2] |= 0xE0; break;
	default: abort();
	}
	break;

      case 'T':
	switch (char_at) {
	case 'X': discrim_chars[i/2] |= 0xF0; break;
	default: abort();
	}
	break;
      }
    }

    SA_i_minus_1 = SA_i;
  }

  fclose(fp);

  FREE(read_buffer);

  return discrim_chars;
}


/* Onepass method */
UINT4 *
Sarray_compute_child (unsigned char *lcp_bytes, UINT4 *lcp_guide, UINT4 *lcp_exceptions, UINT4 n) {
  UINT4 *child;
  UINT4 lastindex, i;
  Uintlist_T indexstack;
  Uintlist_T lcpstack;
  UINT4 lcp_i, lcp_lastindex;
  char *comma;

  child = (UINT4 *) CALLOC(n+1,sizeof(UINT4));

  /* Allocate based on two nibbles per byte */
  /* We consider 0 and (n+1) to be at lcp -1, and 1 at lcp 0 */

  i = 1;
  lastindex = 1;
  lcp_i = 0;
  indexstack = Uintlist_push(NULL,i);
  lcpstack = Uintlist_push(NULL,lcp_i);

  for (i = 2; i <= n; i++) {
    lcp_i = Bytecoding_read_wguide(i,lcp_bytes,lcp_guide,lcp_exceptions,/*lcp_guide_interval*/1024);
    debug2(printf("index %d, lcp %u\n",i,lcp_i));

    while (lcp_i < Uintlist_head(lcpstack)) {
      indexstack = Uintlist_pop(indexstack,&lastindex);
      lcpstack = Uintlist_pop(lcpstack,&lcp_lastindex);
      /* Mark as a right child.  Either down or next, which are
	 treated the same.  Therefore, the next conditional is not
	 really necessary. */
      if (lcp_lastindex == Uintlist_head(lcpstack)) {
#if 0
	debug2(printf("Encoding next %u at %u\n",lastindex,Uintlist_head(indexstack)));
#ifdef NO_ENCODING
	child[Uintlist_head(indexstack)] = lastindex;
#else
	encode_next(child,Uintlist_head(indexstack),lastindex);
#endif
#endif
      } else {
	debug2(printf("Encoding down %u at %u\n",lastindex,Uintlist_head(indexstack)));
#ifdef NO_ENCODING
	child[Uintlist_head(indexstack)] = lastindex;
#else
	encode_down(child,Uintlist_head(indexstack),lastindex);
#endif
      }
    }

    if (lastindex != 0) {
      debug2(printf("Encoding up %u at %u\n",lastindex,i));
#ifdef NO_ENCODING
      child[i-1] = lastindex;	/* up */
#else
      encode_up(child,i,lastindex);
#endif
      lastindex = 0;
    }

    /* This check is still necessary, even if we handle "next" within the previous loop */
    if (lcp_i == Uintlist_head(lcpstack)) {
      /* This is a right sibling, so mark previous index as having a right sibling */
      debug2(printf("Encoding next %u at %u\n",i,Uintlist_head(indexstack)));
#ifdef NO_ENCODING
      child[Uintlist_head(indexstack)] = i;
#else
      encode_next(child,Uintlist_head(indexstack),i);
#endif
    }

    indexstack = Uintlist_push(indexstack,i);
    lcpstack = Uintlist_push(lcpstack,lcp_i);

    if (i % MONITOR_INTERVAL == 0) {
      comma = Genomicpos_commafmt(i);
      fprintf(stderr,"Computing child index %s\n",comma);
      FREE(comma);
    }
    debug2(printf("\n"));
  }

  /* Previously, thought there was no need to clean out stack, because
     all of the next links have been written.  However, skipping the
     section below gave rise to an incorrect child array in the T section */
  debug2(printf("stack still has %d entries\n",Uintlist_length(lcpstack)));
#if 1
  lcp_i = 0;
  while (lcp_i < Uintlist_head(lcpstack)) {
    indexstack = Uintlist_pop(indexstack,&lastindex);
    lcpstack = Uintlist_pop(lcpstack,&lcp_lastindex);
    /* Mark as a right child.  Either down or next, which are
       treated the same.  Therefore, the next conditional is not
       really necessary. */
    if (lcp_lastindex == Uintlist_head(lcpstack)) {
#if 0
      debug2(printf("Final: Encoding next %u at %u\n",lastindex,Uintlist_head(indexstack)));
#ifdef NO_ENCODING
      child[Uintlist_head(indexstack)] = lastindex; /* next */
#else
      encode_next(child,Uintlist_head(indexstack),lastindex);
#endif
#endif
    } else {
      debug2(printf("Final: Encoding down %u at %u\n",lastindex,Uintlist_head(indexstack)));
#ifdef NO_ENCODING
      child[Uintlist_head(indexstack)] = lastindex; /* down */
#else
      encode_down(child,Uintlist_head(indexstack),lastindex);
#endif
    }
  }
#endif


  /* Two choices for final value for up at n+1 (stored at n).  If we
     want it to fail the test (1 < up && up <= n), we could use 0, 1,
     or n+1.  The value n+1 causes -1 to be encoded.  If we want it to
     succeed, then the correct value should be already stored as
     child[1] as a next value. */

#ifdef NO_ENCODING
  debug2(printf("Final: Encoding up %u at %u\n",child[1],n+1));
  child[n] = child[1];
#else
  debug2(printf("Final: Encoding up %u at %u\n",child[1]+1+1,n+1));
  encode_up(child,n+1,child[1] + /*index*/1 + 1);
#endif


#if 0
  if (lastindex != 0) {
    debug2(printf("Final: Encoding up %u at %u\n",lastindex,i));
#ifdef NO_ENCODING
    child[/*i-1*/n] = 1;	/* up: start with node 1 */
#else
    encode_up(child,/*i*/n+1,1);
#endif
    lastindex = 0;
  }

  /* Top-level siblings, at lcp_i == 0.  This conditional should always be true. */
  /* This check is still necessary, even if we handle "next" within the previous loop */
  i = n+1;
  if (lcp_i == Uintlist_head(lcpstack)) {
    /* This is a right sibling, so mark previous index as having a right sibling */
    debug2(printf("Final: Encoding next %u at %u\n",i,Uintlist_head(indexstack)));
#ifdef NO_ENCODING
    child[Uintlist_head(indexstack)] = i;
#else
    encode_next(child,Uintlist_head(indexstack),i);
#endif
  }
#endif

      
  Uintlist_free(&lcpstack);
  Uintlist_free(&indexstack);

#ifdef DEBUG2
  for (i = 0; i <= n; i++) {
    printf("%d %u\n",i,child[i]);
  }
#endif

  return child;
}


#if 0
/* Modeled after Kasai */
UINT4 *
Sarray_traverse_lcp_intervals (unsigned char *lcp_bytes, UINT4 *lcp_exceptions, int n_lcp_exceptions, UINT4 n) {
  UINT4 *child;
  UINT4 lastindex, i;
  Uintlist_T lcpstack, indexstack, lbstack; /* leftbound stack */
  UINT4 lcp_i, lcp_lastindex, leftbound;
  char *comma;

  child = (UINT4 *) CALLOC(n+1,sizeof(UINT4));

  i = 1;
  lastindex = 1;
  lcp_i = 0;			/* Just as good as -1 and avoids comparison with -1U */
  indexstack = Uintlist_push(NULL,i);
  lcpstack = Uintlist_push(NULL,lcp_i);
  lbstack = Uintlist_push(NULL,i-1);
  printf("Pushing %u\n",i-1);

  for (i = 2; i <= n; i++) {
    lcp_i = Bytecoding_read(i,lcp_bytes,lcp_exceptions,n_lcp_exceptions);
    printf("index %d, lcp %u, stack %s\n",i-1,lcp_i,Uintlist_to_string(indexstack));

    leftbound = i-1;
    while (lcp_i < Uintlist_head(lcpstack)) {
      indexstack = Uintlist_pop(indexstack,&lastindex);
      lcpstack = Uintlist_pop(lcpstack,&lcp_lastindex);
      lbstack = Uintlist_pop(lbstack,&leftbound);
      /* rightbound = i-1; */
      printf("lcp interval %u..%u\n",leftbound-1,(i-1)-1);

      if (lcp_i <= Uintlist_head(lcpstack)) {
	printf("  child of %u\n",Uintlist_head(indexstack));
      }

    }

    if (lcp_i > Uintlist_head(lcpstack)) {
      indexstack = Uintlist_push(indexstack,i);
      lcpstack = Uintlist_push(lcpstack,lcp_i);
      lbstack = Uintlist_push(lbstack,leftbound);
    }

    if (i % MONITOR_INTERVAL == 0) {
      comma = Genomicpos_commafmt(i);
      fprintf(stderr,"Computing child table %s\n",comma);
      FREE(comma);
    }
    printf("\n");
  }

  /* i == n + 1.  lcp[n+1] = -1 */
  lcp_i = 0;
  while (lcp_i < Uintlist_head(lcpstack)) {
    indexstack = Uintlist_pop(indexstack,&lastindex);
    lcpstack = Uintlist_pop(lcpstack,&lcp_lastindex);
    lbstack = Uintlist_pop(lbstack,&leftbound);
    printf("lcp interval %u..%u\n",leftbound-1,n-1);
  }
      
  Uintlist_free(&lbstack);
  Uintlist_free(&lcpstack);
  Uintlist_free(&indexstack);

  return child;
}
#endif


void
Sarray_array_uncompress (Genome_T genomecomp, char *sarrayfile, char *plcpptrsfile, char *plcpcompfile,
			 UINT4 genomelength, UINT4 start, UINT4 end) {
  UINT4 n = genomelength, pos, match, h;
  unsigned char *gbuffer;

  int shmid;
  UINT4 *SA, *plcpptrs, *plcpcomp;

  int sa_fd, plcpcomp_fd;
  size_t sa_len, plcpptrs_len, plcpcomp_len;

  double seconds;
  UINT4 sa_i, lcp_i, sa_nexti, lcp_nexti;


  gbuffer = (unsigned char *) CALLOC(n+1,sizeof(unsigned char));
  Genome_fill_buffer_simple(genomecomp,/*left*/0,/*length*/n,gbuffer);
  gbuffer[n] = 0;		       /* '\0', terminator */

  if (end == 0) {
    start = 0;
    end = n;
  }

  SA = (UINT4 *) Access_mmap(&sa_fd,&sa_len,sarrayfile,sizeof(UINT4),/*randomp*/false);
  plcpptrs = (UINT4 *) Access_allocate(&shmid,&plcpptrs_len,&seconds,plcpptrsfile,sizeof(UINT4),/*sharedp*/false);
  plcpcomp = (UINT4 *) Access_mmap(&plcpcomp_fd,&plcpcomp_len,plcpcompfile,sizeof(UINT4),
				  /*randomp*/true);
  plcpcomp = (UINT4 *) Access_mmap(&plcpcomp_fd,&plcpcomp_len,plcpcompfile,sizeof(UINT4),
				  /*randomp*/true);

  /* Bitpack64_read_setup(); */


  printf("i\tSA\tLCP\n");

  pos = start;
#ifdef WORDS_BIGENDIAN
  sa_i = Bigendian_convert_uint(SA[pos]);
#else
  sa_i = SA[pos];
#endif
  lcp_i = Bitpack64_read_one(sa_i,plcpptrs,plcpcomp) - sa_i;

#ifdef WORDS_BIGENDIAN
  sa_nexti = Bigendian_convert_uint(SA[pos+1]);
#else
  sa_nexti = SA[pos+1];
#endif
  lcp_nexti = Bitpack64_read_one(sa_nexti,plcpptrs,plcpcomp) - sa_nexti;

  if (pos == 0) {
    /* lcp_i is -1 */
    printf("%u\t%u\t-1\t",pos,sa_i);
    printf("R");
    printf("\n");
  }

  for (pos = start + 1; pos < end; pos++) {
    sa_i = sa_nexti;
    lcp_i = lcp_nexti;

#ifdef WORDS_BIGENDIAN
    sa_nexti = Bigendian_convert_uint(SA[pos+1]);
#else
    sa_nexti = SA[pos+1];
#endif
    lcp_nexti = Bitpack64_read_one(sa_nexti,plcpptrs,plcpcomp) - sa_nexti;

    printf("%u\t%u\t%u\t",pos,sa_i,lcp_i);
    if (lcp_i > lcp_nexti) {
      printf("L");
    } else {
      printf("R");
    }

    printf("\t%c",gbuffer[sa_i]);
    for (h = 1; (h <= lcp_i || h <= lcp_nexti) && h < 25; h++) {
      if (gbuffer[sa_i+h] == '\0') {
	printf("$");
      } else {
	printf("%c",gbuffer[sa_i+h]);
      }
    }
    if (h <= lcp_i || h <= lcp_nexti) {
      printf("...");
    }
    
    printf("\n");
  }

  sa_i = sa_nexti;
  lcp_i = lcp_nexti;
  if (pos == n) {
    printf("%u\t%u\t%u\t",pos,sa_i,lcp_i);
    printf("L");

    printf("\t%c",gbuffer[sa_i]);
    for (h = 1; h <= lcp_i && h < 25; h++) {
      if (gbuffer[sa_i+h] == '\0') {
	printf("$");
      } else {
	printf("%c",gbuffer[sa_i+h]);
      }
    }
    if (h <= lcp_i) {
      printf("...");
    }

    printf("\n");
  }

  FREE(plcpptrs);
  munmap((void *) plcpcomp,plcpcomp_len);
  close(plcpcomp_fd);
  munmap((void *) SA,sa_len);
  close(sa_fd);

  return;
}


/* May not be possible to uncompress in linear form.  May need to uncompress by lcp-intervals */
void
Sarray_child_uncompress (Genome_T genomecomp, unsigned char *lcpchilddc, UINT4 *lcp_guide, UINT4 *lcp_exceptions,
			 int n_lcp_exceptions, UINT4 *child_guide, UINT4 *child_exceptions, int n_child_exceptions,
			 UINT4 *SA, UINT4 genomelength, UINT4 start, UINT4 end) {
  UINT4 n = genomelength, pos, h;
  unsigned char *gbuffer;
  UINT4 sa_i, lcp_i, child_i, sa_nexti, lcp_nexti, child_nexti;
  char c1, c2;


  gbuffer = (unsigned char *) CALLOC(n+1,sizeof(unsigned char));
  Genome_fill_buffer_simple(genomecomp,/*left*/0,/*length*/n,gbuffer);
  gbuffer[n] = 0;		       /* '\0', terminator */

  if (end == 0) {
    start = 0;
    end = n;
  }

  printf("i\tSA\tLCP\tDC\n");

  pos = start;

  for (pos = start; pos <= end; pos++) {
#ifdef WORDS_BIGENDIAN
    sa_i = Bigendian_convert_uint(SA[pos]);
#else
    sa_i = SA[pos];
#endif
    lcp_i = Bytecoding_lcpchilddc_lcp(pos,lcpchilddc,lcp_exceptions,n_lcp_exceptions); /* lcp(i,j) */
    c2 = Bytecoding_lcpchilddc_dc(&c1,pos,lcpchilddc);

    printf("%u\t%u\t%u\t%c%c",pos,sa_i,lcp_i,c1,c2);

    printf("\t%c",gbuffer[sa_i]);
    for (h = 1; (h <= lcp_i || h <= lcp_nexti) && h < 25; h++) {
      if (gbuffer[sa_i+h] == '\0') {
	printf("$");
      } else {
	printf("%c",gbuffer[sa_i+h]);
      }
    }
    if (h <= lcp_i || h <= lcp_nexti) {
      printf("...");
    }
    
    printf("\n");
  }

  return;
}




