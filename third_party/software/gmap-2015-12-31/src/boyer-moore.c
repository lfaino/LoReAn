static char rcsid[] = "$Id: boyer-moore.c 145990 2014-08-25 21:47:32Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifndef STANDALONE
#include "boyer-moore.h"
#endif
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#ifdef STANDALONE
#define CALLOC calloc
#define FREE free
#else
#include "mem.h"
#endif
#include "bool.h"
#include "complement.h"
#include "genome.h"


#define ASIZE 5			/* In genomic sequence: A, C, G, T, other */
#define MAX(a,b) ((a) > (b) ? (a) : (b))

#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif

/* Successes */
#ifdef DEBUG1
#define debug1(x) x
#else
#define debug1(x)
#endif

/* PMAP good suffix calculation */
#ifdef DEBUG2
#define debug2(x) x
#else
#define debug2(x)
#endif

#ifdef PMAP
/* Same as in dynprog.c */
/* Handle only cases in iupac table in dynprog.c */
static bool matchtable[26][26] = 
/*  A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q,R,S,T,U,V,W,X,Y,Z */
  {{1,0,0,0,0,0,0,1,0,0,0,0,1,1,0,0,0,1,0,0,0,0,1,0,0,0}, /* A */
   {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}, /* B */
   {0,0,1,0,0,0,0,1,0,0,0,0,1,1,0,0,0,0,1,0,0,0,0,0,1,0}, /* C */
   {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}, /* D */
   {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}, /* E */
   {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}, /* F */
   {0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,1,1,0,0,0,0,0,0,0}, /* G */
   {1,0,1,0,0,0,0,1,0,0,0,0,1,1,0,0,0,1,1,1,0,0,1,0,1,0}, /* H = [ACT] */
   {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}, /* I */
   {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}, /* J */
   {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}, /* K */
   {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}, /* L */
   {1,0,1,0,0,0,0,1,0,0,0,0,1,1,0,0,0,1,1,0,0,0,1,0,1,0}, /* M = [AC] */
   {1,0,1,0,0,0,1,1,0,0,0,0,1,1,0,0,0,1,1,1,0,0,1,0,1,0}, /* N = [ACGT] */
   {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}, /* O */
   {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}, /* P */
   {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}, /* Q */
   {1,0,0,0,0,0,1,1,0,0,0,0,1,1,0,0,0,1,1,0,0,0,1,0,0,0}, /* R = [AG] */
   {0,0,1,0,0,0,1,1,0,0,0,0,1,1,0,0,0,1,1,0,0,0,0,0,1,0}, /* S = [CG] */
   {0,0,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,1,0,0,1,0,1,0}, /* T */
   {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}, /* U */
   {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}, /* V */
   {1,0,0,0,0,0,0,1,0,0,0,0,1,1,0,0,0,1,0,1,0,0,1,0,1,0}, /* W = [AT] */
   {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}, /* X */
   {0,0,1,0,0,0,0,1,0,0,0,0,1,1,0,0,0,0,1,1,0,0,1,0,1,0}, /* Y = [CT] */
   {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}}; /* Z */

static bool inverse_matchtable[26][26] = 
/*  A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q,R,S,T,U,V,W,X,Y,Z */
  {{0,0,1,0,0,0,1,1,0,0,0,0,1,1,0,0,0,1,1,1,0,0,1,0,1,0}, /* not A = [CGT] */
   {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}, /* B */
   {1,0,0,0,0,0,1,1,0,0,0,0,1,1,0,0,0,1,1,1,0,0,1,0,1,0}, /* not C = [AGT] */
   {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}, /* D */
   {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}, /* E */
   {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}, /* F */
   {1,0,1,0,0,0,0,1,0,0,0,0,1,1,0,0,0,1,1,1,0,0,1,0,1,0}, /* not G = [ACT] */
   {0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,1,1,0,0,0,0,0,0,0}, /* not H = [G] */
   {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}, /* I */
   {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}, /* J */
   {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}, /* K */
   {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}, /* L */
   {0,0,0,0,0,0,1,1,0,0,0,0,0,1,0,0,0,1,1,1,0,0,1,0,1,0}, /* not M = [GT] */
   {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}, /* not N = [] */
   {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}, /* O */
   {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}, /* P */
   {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}, /* Q */
   {0,0,1,0,0,0,0,1,0,0,0,0,1,1,0,0,0,0,1,1,0,0,1,0,1,0}, /* not R = [CT] */
   {1,0,0,0,0,0,0,1,0,0,0,0,1,1,0,0,0,1,0,1,0,0,1,0,1,0}, /* not S = [AT] */
   {1,0,1,0,0,0,1,1,0,0,0,0,1,1,0,0,0,1,1,0,0,0,1,0,1,0}, /* not T = [ACG] */
   {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}, /* U */
   {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}, /* V */
   {0,0,1,0,0,0,1,1,0,0,0,0,1,1,0,0,0,1,1,0,0,0,0,0,1,0}, /* not W = [CG] */
   {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}, /* X */
   {1,0,0,0,0,0,1,1,0,0,0,0,1,1,0,0,0,1,1,0,0,0,1,0,0,0}, /* not Y = [AG] */
   {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}}; /* Z */
#endif


static int
na_index (char c) {
  switch (c) {
  case 'A': return 0;
  case 'C': return 1;
  case 'G': return 2;
  case 'T': return 3;
  default: return 4;
  }
}

static void
precompute_bad_char_shift (int *bad_char_shift, char *query, int querylen) {
  int i;
  char *p;

  for (i = 0; i < ASIZE; i++) {
    bad_char_shift[i] = querylen;
  }
  for (i = 0, p = query; i < querylen - 1; i++, p++) {
#ifdef PMAP
    if (matchtable[*p-'A']['A'-'A'] == true) {
      bad_char_shift[na_index('A')] = querylen - i - 1;
    }
    if (matchtable[*p-'A']['C'-'A'] == true) {
      bad_char_shift[na_index('C')] = querylen - i - 1;
    }
    if (matchtable[*p-'A']['G'-'A'] == true) {
      bad_char_shift[na_index('G')] = querylen - i - 1;
    }
    if (matchtable[*p-'A']['T'-'A'] == true) {
      bad_char_shift[na_index('T')] = querylen - i - 1;
    }
#else
    bad_char_shift[na_index(*p)] = querylen - i - 1;
#endif
  }

  return;
}


#ifdef PMAP
static bool
suffixp (char *query, int querylen, int subseqlen, int j) {
  int i;

  debug2(printf("Entering suffixp at subseqlen = %d, j = %d\n",subseqlen,j));
  for (i = querylen - 1; i >= querylen - 1 - subseqlen + 1 && j >= 0; i--, j--) {
    if (matchtable[query[i]-'A'][query[j]-'A'] == false) {
      debug2(printf("Returning false at i = %d\n",i));
      return false;
    }
  }
  if (j >= 0) {
    if (inverse_matchtable[query[i]-'A'][query[j]-'A'] == false) {
      debug2(printf("Returning false at i = %d\n",i));
      return false;
    }
  }
  debug2(printf("Returning true\n"));
  return true;
}

static int
suffix_jump (char *query, int querylen, int subseqlen) {
  int j;

  for (j = querylen - 1; j >= 0; j--) {
    if (suffixp(query,querylen,subseqlen,j) == true) {
      return querylen - 1 - j;
    }
  }
  return querylen;
}

static void
precompute_good_suffix_shift (int *good_suffix_shift, char *query, int querylen) {
  int subseqlen;

  for (subseqlen = 0; subseqlen < querylen; subseqlen++) {
    good_suffix_shift[querylen - 1 - subseqlen] = suffix_jump(query,querylen,subseqlen);
    debug2(printf("Assigning %d to position %d\n",good_suffix_shift[querylen - 1 - subseqlen],querylen-1-subseqlen));
  }

  return;
}

#else

static void
precompute_suffix (int *suffix, char *query, int querylen) {
  /* Note: initialization of f not necessary, since i < g in first iteration */
  int f, g, i;

  suffix[querylen - 1] = querylen;
  g = querylen - 1;
  for (i = querylen - 2; i >= 0; i--) {
    if (i > g && suffix[i + querylen - 1 - f] < i - g) {
      suffix[i] = suffix[i + querylen - 1 - f];
    } else {
      if (i < g) {
	g = i;
      }
      f = i;
      while (g >= 0 && query[g] == query[g + querylen - 1 - f]) {
	g--;
      }
      suffix[i] = f - g;
    }
  }
  debug(
	printf("suffix:\n");
	for (i = 0; i < querylen; i++) {
	  printf("%d %d\n",i,suffix[i]);
	}
	printf("\n");
	);

  return;
}
  
static void
precompute_good_suffix_shift (int *good_suffix_shift, char *query, int querylen) {
  int i, j, *suffix;

  suffix = (int *) MALLOCA(querylen * sizeof(int));
  precompute_suffix(suffix,query,querylen);
  
  for (i = 0; i < querylen; i++) {
    good_suffix_shift[i] = querylen;
  }
  j = 0;
  for (i = querylen - 1; i >= -1; i--) {
    if (i == -1 || suffix[i] == i + 1) {
      for ( ; j < querylen - 1 - i; j++) {
	if (good_suffix_shift[j] == querylen) {
	  good_suffix_shift[j] = querylen - 1 - i;
	}
      }
    }
  }
  for (i = 0; i <= querylen - 2; i++) {
    good_suffix_shift[querylen - 1 - suffix[i]] = querylen - 1 - i;
  }

  FREEA(suffix);

  return;
}
#endif

static bool
query_okay (char *query, int querylen) {
  int i;
  char *p, c;

#ifdef PMAP
  for (i = 0, p = query; i < querylen; i++, p++) {
    c = *p;
    if (c < 'A' || c > 'Y') {
      return false;
    }
  }
#else
  for (i = 0, p = query; i < querylen; i++, p++) {
    c = *p;
    if (c != 'A' && c != 'C' && c != 'G' && c != 'T') {
      return false;
    }
  }
#endif
  return true;
}

#ifdef STANDALONE
static void
#else
Intlist_T
#endif
BoyerMoore (char *query, int querylen, char *text, int textlen) {
#ifndef STANDALONE
  Intlist_T hits = NULL;
#endif
  int i, j, *good_suffix_shift;
  int bad_char_shift[ASIZE];

  if (query_okay(query,querylen)) {
    good_suffix_shift = (int *) MALLOCA(querylen * sizeof(int));
    precompute_good_suffix_shift(good_suffix_shift,query,querylen);
    precompute_bad_char_shift(bad_char_shift,query,querylen);

    debug(
	  printf("bad_char_shift:\n");
	  for (i = 0; i < ASIZE; i++) {
	    printf("%d %d\n",i,bad_char_shift[i]);
	  }
	  printf("\n");
	  printf("good_suffix_shift:\n");
	  for (i = 0; i < querylen; i++) {
	    printf("%d %d\n",i,good_suffix_shift[i]);
	  }
	  );

    j = 0;
    while (j <= textlen - querylen) {
#ifdef PMAP
      for (i = querylen - 1; i >= 0 && matchtable[query[i]-'A'][text[i+j]-'A'] == true; i--) ;
#else
      for (i = querylen - 1; i >= 0 && query[i] == text[i+j]; i--) ;
#endif
      if (i < 0) {
#ifndef STANDALONE
	hits = Intlist_push(hits,j);
#endif
	
	debug1(printf("Success at %d\n",j));
	debug(printf("Shift by %d (Gs[0])\n",good_suffix_shift[0]));
	j += good_suffix_shift[0];
      } else {
	debug(
	      if (good_suffix_shift[i] > 
		  bad_char_shift[na_index(text[i+j])] - querylen + 1 + i) {
		printf("Shift by %d (Gs[%d])\n",
		       good_suffix_shift[i],i);
	      } else {
		printf("Shift by %d (Gs[%d] == Bc[%c] - %d + %d)\n",
		       bad_char_shift[na_index(text[i+j])] - querylen + 1 + i,
		       i,text[i+j],querylen,i+1);
	      }
	      );
	j += MAX(good_suffix_shift[i],
		 bad_char_shift[na_index(text[i+j])] - querylen + 1 + i);
      }
    }

    FREEA(good_suffix_shift);
  }

#ifndef STANDALONE
  return hits;
#endif
}


#if 0
static char complCode[128] = COMPLEMENT_LC;

static char
get_genomic_nt (char *g_alt, int genomicpos,
		Univcoord_T chroffset, Univcoord_T chrhigh, bool watsonp) {
  char c2, c2_alt;

  if (watsonp) {
#if 0
    if (genome) {
      return Genome_get_char(genome,chroffset + genomicpos);
    } else {
#endif
      return Genome_get_char_blocks(&(*g_alt),chroffset + genomicpos);
#if 0
    }
#endif

  } else {
#if 0
    if (genome) {
      c2 = Genome_get_char(genome,chrhigh - genomicpos);
    } else {
#endif
      c2 = Genome_get_char_blocks(&c2_alt,chrhigh - genomicpos);
#if 0
    }
#endif
    *g_alt = complCode[(int) c2_alt];
    return complCode[(int) c2];
  }
}
#endif


Intlist_T
BoyerMoore_nt (char *query, int querylen, int textoffset, int textlen,
	       Univcoord_T chroffset, Univcoord_T chrhigh, bool watsonp) {
  Intlist_T hits = NULL;
  int i, j, *good_suffix_shift;
  int bad_char_shift[ASIZE];
  char *text, *text_alt;

  if (query_okay(query,querylen)) {
    good_suffix_shift = (int *) MALLOCA(querylen * sizeof(int));
    text = (char *) MALLOC((textlen+querylen+1) * sizeof(char)); /* alloca could cause stack overflow */
    text_alt = (char *) MALLOC((textlen+querylen+1) * sizeof(char)); /* alloca could cause stack overflow */

    precompute_good_suffix_shift(good_suffix_shift,query,querylen);
    precompute_bad_char_shift(bad_char_shift,query,querylen);
    if (watsonp) {
      Genome_get_segment_blocks_right(text,text_alt,/*left*/chroffset+textoffset,/*length*/textlen+querylen,chrhigh,/*revcomp*/false);
    } else {
      Genome_get_segment_blocks_left(text,text_alt,/*right*/chrhigh-textoffset+1,/*length*/textlen+querylen,chroffset,/*revcomp*/true);
    }
    if (text[0] == '\0') {
      FREE(text_alt);
      FREE(text);
      FREEA(good_suffix_shift);
      return hits;
    }

    /* This makes text[i+j] == get_genomic_nt(&g_alt,textoffset+i+j,chroffset,chrhigh,watsonp) */

    debug(
	  printf("bad_char_shift:\n");
	  for (i = 0; i < ASIZE; i++) {
	    printf("%d %d\n",i,bad_char_shift[i]);
	  }
	  printf("\n");
	  printf("good_suffix_shift:\n");
	  for (i = 0; i < querylen; i++) {
	    printf("%d %d\n",i,good_suffix_shift[i]);
	  }
	  );

    j = 0;
    while (j <= textlen - querylen) {
#ifdef PMAP
      for (i = querylen - 1; i >= 0 && matchtable[query[i]-'A'][text[i+j]-'A'] == true; i--) ;
#else
      for (i = querylen - 1; i >= 0 && query[i] == text[i+j]; i--) ;
#endif
      if (i < 0) {
	hits = Intlist_push(hits,j);
	
	debug1(printf("Success at %d\n",j));
	debug(printf("Shift by %d (Gs[0])\n",good_suffix_shift[0]));
	j += good_suffix_shift[0];
      } else {
	debug(
	      if (good_suffix_shift[i] > 
		  bad_char_shift[na_index(text[i+j])] - querylen + 1 + i) {
		printf("Shift by %d (Gs[%d])\n",
		       good_suffix_shift[i],i);
	      } else {
		printf("Shift by %d (Gs[%d] == Bc[%c] - %d + %d)\n",
		       bad_char_shift[na_index(text[i+j])] - querylen + 1 + i,
		       i,text[i+j],querylen,i+1);
	      }
	      );
	j += MAX(good_suffix_shift[i],
		 bad_char_shift[na_index(text[i+j])] - querylen + 1 + i);
      }
    }

    FREE(text_alt);
    FREE(text);
    FREEA(good_suffix_shift);
  }

  return hits;
}


void
BoyerMoore_bad_char_shift (int *bad_char_shift, char *query, int querylen) {
#ifdef DEBUG
  int i;
#endif

  if (query_okay(query,querylen)) {
    precompute_bad_char_shift(bad_char_shift,query,querylen);
  } else {
    fprintf(stderr,"Query cannot have bad characters\n");
    abort();
  }

  debug(
	printf("bad_char_shift:\n");
	for (i = 0; i < ASIZE; i++) {
	  printf("%d %d\n",i,bad_char_shift[i]);
	}
	printf("\n");
	);

  return;
}


int
BoyerMoore_maxprefix (char *query, int querylen, char *text, int textlen,
		      int *bad_char_shift) {
  int maxprefix = 0;
  int i, j;

  debug(printf("Query: %s\n",query));
  debug(printf("Text:  %s\n",text));

  j = -querylen + 1;
  while (j <= textlen - querylen) {
#ifdef PMAP
    for (i = querylen - 1; i >= 0 && matchtable[query[i]-'A'][text[i+j]-'A'] == true; i--) ;
#else
    for (i = querylen - 1; i+j >= 0 && query[i] == text[i+j]; i--) ;
#endif
    if (i < 0 || i+j < 0) {
      maxprefix = querylen+j;
      debug1(printf("Success at %d (length %d)\n",j,querylen+j));
      debug(printf("Shift by 1\n"));
      j += 1;
    } else {
      debug1(printf("Failure at j=%d, i=%d\n",j,i));
      debug(
	    printf("Shift by %d (Gs[%d] == Bc[%c] - %d + %d)\n",
		   bad_char_shift[na_index(text[i+j])],
		   i,text[i+j],querylen,i+1);
	    );
      j += bad_char_shift[na_index(text[i+j])];
    }
  }

  return maxprefix;
}



#if 0
static char *
string_reverse (char *original, int length) {
  char *reverse;
  int i, j;

  reverse = (char *) CALLOC(length,sizeof(char));

  for (i = 0, j = length-1; i < length; i++, j--) {
    reverse[i] = original[j];
  }
  return reverse;
}
#endif


#ifdef STANDALONE
int
main (int argc, char *argv[]) {
  char *query, *text;
  char *query_rev, *text_rev;
  int querylen, textlen;
  int bad_char_shift[ASIZE];

#ifdef PMAP
  int i, j;
  for (i = 0; i < 26; i++) {
    for (j = 0; j < 26; j++) {
      if (matchtable[i][j] != matchtable[j][i]) {
	fprintf(stderr,"Assymetry problem at %d,%d\n",i,j);
      }
    }
  }
#endif

  query = argv[1];
  text = argv[2];
  querylen = strlen(query);
  textlen = strlen(text);

#if 0
  BoyerMoore(query,querylen,text,textlen);
#else
  /* Test of maxprefix */
  query_rev = string_reverse(query,querylen);
  text_rev = string_reverse(text,textlen);

  printf("Query rev: %s\n",query_rev);
  printf("Text rev: %s\n",text_rev);
  BoyerMoore_bad_char_shift(bad_char_shift,query_rev,querylen);
  BoyerMoore_maxprefix(query_rev,querylen,text_rev,textlen,bad_char_shift);
#endif

  return 0;
}
#endif
