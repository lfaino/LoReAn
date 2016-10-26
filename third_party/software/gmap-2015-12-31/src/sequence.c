static char rcsid[] = "$Id: sequence.c 170023 2015-07-17 16:47:21Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#ifndef HAVE_MEMCPY
# define memcpy(d,s,n) bcopy((s),(d),(n))
#endif
#ifndef HAVE_MEMMOVE
# define memmove(d,s,n) bcopy((s),(d),(n))
#endif

#include "sequence.h"

#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>		/* For rindex */
#include <ctype.h>		/* For iscntrl, isspace, and toupper */

#ifdef HAVE_ZLIB
#include <zlib.h>
#define GZBUFFER_SIZE 131072
#endif


#include "assert.h"
#include "mem.h"
#include "complement.h"
#include "intlist.h"
#include "fopen.h"
#include "md5.h"
#ifdef PMAP
#include "dynprog.h"
#endif

/* Before setting DEBUG, may want to reduce MAXSEQLEN in sequence.h */
#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif

/* Pointers for first half and second half */
#ifdef DEBUG1
#define debug1(x) x
#else
#define debug1(x)
#endif



/***********************************************************************
 *    Definitions:
 *
 *   TTTTTT ACGT ...... ACGT AAAAAA
 *          <- trimlength ->
 *   <-------- fulllength -------->
 *          ^trimstart
 *   ^contents
 *
 ************************************************************************/


#define T Sequence_T
struct T {
  char *acc;			/* Accession */
  char *restofheader;		/* Rest of header */
  char *contents;		/* Original sequence, ends with '\0' */
  char *contents_alloc;		/* Allocation */
  int fulllength;		/* Full length (not including chopped sequence) */

  int trimstart;		/* Start of trim */
  int trimend;			/* End of trim */
#ifdef PMAP
  int fulllength_given;		/* Full length minus implicit stop codon at end */
#endif
  int subseq_offset;		/* Used only for subsequences */
  int skiplength;		/* Used only for sequences longer than MAXSEQLEN */
  /* bool free_contents_p; */

#ifndef PMAP
  char *quality;		/* For Illumina short reads read via extended FASTA */
  char *quality_alloc;		/* Allocation */
#endif

  bool firstp;			/* First end indicated by '>', second end by '<' in extended FASTA file */
};

bool
Sequence_firstp (T this) {
  return this->firstp;
}

char *
Sequence_accession (T this) {
  return this->acc;
}

char *
Sequence_fullpointer (T this) {
  return this->contents;
}

char *
Sequence_trimpointer (T this) {
  return &(this->contents[this->trimstart]);
}

#ifndef PMAP
char *
Sequence_quality_string (T this) {
  return this->quality;
}
#endif


int
Sequence_ntlength (T this) {
#ifdef PMAP
  return 3*this->fulllength;
#else
  return this->fulllength;
#endif
}

int
Sequence_fulllength (T this) {
  return this->fulllength;
}

int
Sequence_fulllength_given (T this) {
#ifdef PMAP
  return this->fulllength_given;
#else
  return this->fulllength;
#endif
}

char *
Sequence_subseq_pointer (T this, int querystart) {
  return &(this->contents[querystart]);
}

int
Sequence_subseq_length (T this, int querystart) {
  return this->fulllength - querystart;
}



int
Sequence_trimlength (T this) {
  return this->trimend - this->trimstart;
}

void
Sequence_trim (T this, int trim_start, int trim_end) {
  this->trimstart = trim_start;
  this->trimend = trim_end;
  return;
}

int
Sequence_trim_start (T this) {
  return this->trimstart;
}

int
Sequence_trim_end (T this) {
  return this->trimend;
}

int
Sequence_subseq_offset (T this) {
  return this->subseq_offset;
}

int
Sequence_skiplength (T this) {
  return this->skiplength;
}

void
Sequence_free (T *old) {
  if (*old) {
    if ((*old)->restofheader != NULL) {
      FREE_IN((*old)->restofheader);
    }
    if ((*old)->acc != NULL) {
      FREE_IN((*old)->acc);
    }

#ifndef PMAP
    if ((*old)->quality_alloc != NULL) {
      FREE_IN((*old)->quality_alloc);
    }
#endif

    FREE_IN((*old)->contents_alloc);

    FREE_IN(*old);
  }
  return;
}

#ifdef PMAP
static int aa_index_table[128] =
  { -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,

  /*     *                                  */
    -1, 20, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1,

  /* A,  B,  C,  D,  E,  F,  G,  H,  I,  J, */
     0, -1,  1,  2,  3,  4,  5,  6,  7, -1,

  /* K,  L,  M,  N,  O,  P,  Q,  R,  S,  T  */
     8,  9, 10, 11, -1, 12, 13, 14, 15, 16,

  /* U,  V,  W,  X,  Y,  Z  */
    21, 17, 18, -1, 19, -1,

    -1, -1, -1, -1, -1, -1,

  /* a,  b,  c,  d,  e,  f,  g,  h,  i,  j, */
     0, -1,  1,  2,  3,  4,  5,  6,  7, -1,

  /* k,  l,  m,  n,  o,  p,  q,  r,  s,  t  */
     8,  9, 10, 11, -1, 12, 13, 14, 15, 16,

  /* u,  v,  w,  x,  y,  z  */
    21, 17, 18, -1, 19, -1,

    -1, -1, -1, -1, -1};


static char *iupac_table[21+1] =
  {"GCN",			/* A */
   "TGY",			/* C */
   "GAY",			/* D */
   "GAR",			/* E */
   "TTY",			/* F */
   "GGN",			/* G */
   "CAY",			/* H */
   "ATH",			/* I */
   "AAR",			/* K */
   "YTN",			/* L */
   "ATG",			/* M */
   "AAY",			/* N */
   "CCN",			/* P */
   "CAR",			/* Q */
   "MGN",			/* R */
   "WSN",			/* S */
   "ACN",			/* T */
   "GTN",			/* V */
   "TGG",			/* W */
   "TAY",			/* Y */
   "TRR",			/* STOP */
   "TGA"};			/* U */


static char aa_table[21+1] = "ACDEFGHIKLMNPQRSTVWY*U";

char
Sequence_codon_char (char aa, int codonpos) {
  int index;
  char *codon;

  if ((index = aa_index_table[(int) aa]) < 0) {
    return 'N';
  } else {
    codon = iupac_table[index];
    return codon[codonpos];
  }
}



static char *
instantiate_codons (char *aasequence, int ntlength) {
  char *newntsequence;
  int aapos, ntpos;
  int index, frame;
  char *codon, aa;

  newntsequence = (char *) CALLOC(ntlength+1,sizeof(char));

  for (ntpos = 0; ntpos < ntlength; ntpos++) {
    aapos = (/*offset+*/ntpos)/3;
	
    aa = aasequence[aapos];
    index = aa_index_table[(int) aa];
    if (index < 0) {
      newntsequence[ntpos] = 'N';
    } else {
      codon = iupac_table[index];

      frame = (/*offset+*/ntpos) % 3;
      newntsequence[ntpos] = codon[frame];
    }
  }

#if 0
  printf("New sequence: %s\n",newntsequence);
#endif
  return newntsequence;
}


T
Sequence_convert_to_nucleotides (T this) {
  T new = (T) MALLOC_IN(sizeof(*new));
  int i;

  new->acc = (char *) NULL;
  new->restofheader = (char *) NULL;
  new->fulllength = this->fulllength*3;
  new->fulllength_given = this->fulllength_given*3;
  new->contents = new->contents_alloc =
    instantiate_codons(/*aasequence*/this->contents,/*ntlength*/this->fulllength*3);
#if 0
  for (i = 0; i < new->fulllength; i++) {
    new->contents[i] = '?';
  }
#endif

  new->trimstart = 0;
  new->trimend = new->fulllength_given;
  /* new->free_contents_p = true; */
  new->subseq_offset = 0;
  new->skiplength = this->skiplength;
  new->firstp = this->firstp;

  return new;
}
#endif


#if 0
int
Sequence_count_bad (T this, int pos, int max, int direction) {
  int nbad = 0;

  if (direction > 0) {
    while (--max >= 0 && pos < this->fulllength) {
      if (this->contents[pos] == 'X') {
	nbad++;
      }
      pos++;
    }
  } else {
    while (--max >= 0 && pos >= 0) {
      if (this->contents[pos] == 'X') {
	nbad++;
      }
      pos--;
    }
  }

  return nbad;
}
#endif


#define HEADERLEN 512
#define DISCARDLEN 8192

static char Header[HEADERLEN];
static char Discard[DISCARDLEN];

static char Sequence[1+MAXSEQLEN+1]; /* Used by Sequence_read_unlimited */
static char Sequence1[HALFLEN+1]; /* 1 at end for '\0' */
static char Sequence2[HALFLEN+3]; /* 1 at end for '\0' and 2 extra in cyclic part for '\n' and '\0' */

static char *Firsthalf;
static char *Secondhalf;
/* static int Initc = '\0'; */


/* The first element of Sequence is always the null character, to mark
   the end of the string */

/* Skipping of dashes might still be buggy */
/*
#define DASH '-'
*/


/* Returns '@', '>', or '<' if FASTA file, first sequence char if not */
int
Sequence_input_init (FILE *fp) {
  int c;
  bool okayp = false;

  Header[0] = '\0';
  Sequence[0] = '\0';
  Firsthalf = &(Sequence1[0]);
  Secondhalf = &(Sequence2[0]);

  while (okayp == false && (c = fgetc(fp)) != EOF) {
    debug(printf("Read character %c\n",c));
    if (iscntrl(c)) {
#ifdef DASH
    } else if (c == DASH) {
#endif
    } else if (isspace(c)) {
    } else {
      okayp = true;
    }
  }

  debug(printf("Returning initial character %c\n",c));
  return c;
}

#ifdef HAVE_ZLIB
/* Returns '@', '>', or '<' if FASTA file, first sequence char if not */
int
Sequence_input_init_gzip (gzFile fp) {
  int c;
  bool okayp = false;

  Header[0] = '\0';
  Sequence[0] = '\0';
  Firsthalf = &(Sequence1[0]);
  Secondhalf = &(Sequence2[0]);

  while (okayp == false && (c = gzgetc(fp)) != EOF) {
    debug(printf("Read character %c\n",c));
    if (iscntrl(c)) {
#ifdef DASH
    } else if (c == DASH) {
#endif
    } else if (isspace(c)) {
    } else {
      okayp = true;
    }
  }

  debug(printf("Returning initial character %c\n",c));
  return c;
}
#endif

#ifdef HAVE_BZLIB
/* Returns '@', '>', or '<' if FASTA file, first sequence char if not */
int
Sequence_input_init_bzip2 (Bzip2_T fp) {
  int c;
  bool okayp = false;

  Header[0] = '\0';
  Sequence[0] = '\0';
  Firsthalf = &(Sequence1[0]);
  Secondhalf = &(Sequence2[0]);

  while (okayp == false && (c = bzgetc(fp)) != EOF) {
    debug(printf("Read character %c\n",c));
    if (iscntrl(c)) {
#ifdef DASH
    } else if (c == DASH) {
#endif
    } else if (isspace(c)) {
    } else {
      okayp = true;
    }
  }

  debug(printf("Returning initial character %c\n",c));
  return c;
}
#endif


static void
blank_header (T this) {
  this->acc = (char *) CALLOC_IN(strlen("NO_HEADER")+1,sizeof(char));
  strcpy(this->acc,"NO_HEADER");
  this->restofheader = (char *) CALLOC_IN(1,sizeof(char));
  this->restofheader[0] = '\0';
  return;
}

static char *
input_header (FILE *fp, T this) {
  char *p;
  size_t length;

  if (feof(fp)) {
    return NULL;
  } else if (fgets(&(Header[0]),HEADERLEN,fp) == NULL) {
    /* File must terminate after > */
    return NULL;
  }

  if (Header[0] == '\n') {
    Header[0] = '\0';
  } else if ((p = rindex(&(Header[0]),'\n')) != NULL) {
    if (p[-1] == '\r') {
      p--;
    }
    *p = '\0';
  } else {
    /* Eliminate rest of header from input */
    while (fgets(&(Discard[0]),DISCARDLEN,fp) != NULL &&
	   rindex(&(Discard[0]),'\n') == NULL) ;
  }

  p = &(Header[0]);
  while (*p != '\0' && !isspace((int) *p)) {
    p++;
  }
  if (*p == '\0') {
    /* Accession only */
    length = (p - &(Header[0]))/sizeof(char);
    this->acc = (char *) CALLOC_IN(length+1,sizeof(char));
    strcpy(this->acc,Header);
    this->restofheader = (char *) CALLOC_IN(1,sizeof(char));
    this->restofheader[0] = '\0';
  } else {
    *p = '\0';
    length = (p - &(Header[0]))/sizeof(char);
    this->acc = (char *) CALLOC_IN(length+1,sizeof(char));
    strcpy(this->acc,Header);
    p++;
    this->restofheader = (char *) CALLOC_IN(strlen(p)+1,sizeof(char));
    strcpy(this->restofheader,p);
  }

  return this->acc;
} 

static bool
skip_header (FILE *fp) {

  if (feof(fp)) {
    return false;
  } else if (fgets(&(Header[0]),HEADERLEN,fp) == NULL) {
    /* File must terminate after > */
    return false;
  }

  if (rindex(&(Header[0]),'\n') == NULL) {
    /* Eliminate rest of header from input */
    while (fgets(&(Discard[0]),DISCARDLEN,fp) != NULL &&
	   rindex(&(Discard[0]),'\n') == NULL) ;
  }

  return true;
} 

#ifdef DEBUG
static void
print_contents (char *p, int length) {
  int i;
  FILE *fp = stdout;
	
  fprintf(fp,"\"");
  for (i = 0; i < length; i++) {
    if (*p == '\0') {
      fprintf(fp,"_");
    } else {
      fprintf(fp,"%c",*p);
    }
    p++;
  }
  fprintf(fp,"\"\n");
  return;
}
#endif


#define CONTROLM 13		/* From PC */
#define SPACE 32

#if 0
static char *
find_bad_char (char *line) {
  char *first, *p1, *p2;
#ifdef DASH
  char *p3;
#endif

  /* p1 = index(line,CONTROLM); */
  p2 = index(line,SPACE);
  /* p3 = index(line,DASH); */

  if (/* p1 == NULL && */ p2 == NULL /* && p3 == NULL*/) {
    return NULL;
  } else {
#if 0
    if (p1) {
      first = p1;
    }
    if (p2) {
      first = p2;
    }
    /*
    if (p3) {
      first = p3;
    }
    */
    if (p1 && p1 < first) {
      first = p1;
    }
    if (p2 && p2 < first) {
      first = p2;
    }
    /*
    if (p3 && p3 < first) {
      first = p3;
    }
    */
    return first;
#else
    return p2;
#endif
  }
}
#endif


static int
read_first_half (int *nextchar, bool *eolnp, FILE *fp, bool possible_fasta_header_p) {
  int remainder, strlenp;
  char *ptr, *p = NULL;
  int c;
  bool init_char_p = true;

  ptr = &(Firsthalf[0]);
  if (possible_fasta_header_p == false || (*nextchar != '@' && *nextchar != '>' && *nextchar != '<' && *nextchar != '+')) {
    *ptr++ = (char) *nextchar;
  }
  remainder = (&(Firsthalf[HALFLEN]) - ptr)/sizeof(char);

  while (1) {
    if (remainder <= 0) {
      debug(printf("remainder <= 0.  Returning false\n"));
      *nextchar = EOF;
      debug1(printf("read_first_half returning length1 of %d\n",(ptr - &(Firsthalf[0]))/sizeof(char)));
      return (ptr - &(Firsthalf[0]))/sizeof(char);

    } else if (feof(fp)) {
      /* EOF in middle of line */
      debug(printf("EOF.  Returning true\n"));
      *nextchar = EOF;
      debug1(printf("read_first_half returning length1 of %d\n",(ptr - &(Firsthalf[0]))/sizeof(char)));
      return (ptr - &(Firsthalf[0]))/sizeof(char);

    } else if (*eolnp == true) {
      /* Peek at character after eoln */
      if ((c = fgetc(fp)) == EOF || (init_char_p == false && (c == '@' || c == '>' || c == '<' || c == '+'))) {
	debug(printf("c == EOF or > or < or +.   Returning true\n"));
	*nextchar = c;
	return (ptr - &(Firsthalf[0]))/sizeof(char);
      } else if (iscntrl(c)) {
	debug(printf("c == control char.  Continuing\n"));
#ifdef DASH
      } else if (c == DASH) {
	debug(printf("c == dash.  Continuing\n"));
#endif
      } else if (isspace(c)) {
	*eolnp = true;
	debug(printf("c == NULL.  Continuing\n"));
      } else {
	*ptr++ = (char) c;
	remainder--;
	*eolnp = false;
	p = NULL;
	debug(printf("c == sth.  Continuing\n"));
      }

    } else {
      debug(printf("Trying to read remainder of %d\n",remainder));
      if (p != NULL) {
	strlenp = strlen(p);
	memmove(ptr,p,strlenp);
        ptr[strlenp] = '\0';
	debug(printf("Did memmove of %d chars at %p to %p\n",strlenp,p,ptr));
      } else {
	p = fgets(ptr,remainder+1,fp);
      }
      if (p == NULL) {
	debug(printf("line == NULL.  Returning true\n"));
	*nextchar = EOF;
	debug1(printf("read_first_half returning length1 of %d\n",(ptr - &(Firsthalf[0]))/sizeof(char)));
	return (ptr - &(Firsthalf[0]))/sizeof(char);
      } else {
	debug(printf("Read %s.\n",ptr));
	/* Was a call to find_bad_char */
	while ((p = index(ptr,SPACE)) != NULL) {
	  ptr = p++;
	  strlenp = strlen(p);
	  memmove(ptr,p,strlenp);
	  ptr[strlenp] = '\0';
	  debug(printf("Found space.  Did memmove of %d chars at %p to %p\n",strlenp,p,ptr));
	}
	if (*ptr == '\n') {
	  *eolnp = true;
	  debug(printf("line == EOLN.  Continuing\n"));
	} else if ((p = index(ptr,'\n')) != NULL) {
	  if (p[-1] == '\r') {
	    p--;
	  }
	  ptr = p;
	  *eolnp = true;
	  debug(printf("line == EOLN.  Continuing\n"));
	} else {
	  ptr += strlen(ptr);
	  *eolnp = false;
	  p = NULL;
	  debug(printf("line != EOLN.  Continuing\n"));
	}
	remainder = (&(Firsthalf[HALFLEN]) - ptr)/sizeof(char);
      }
    }

    debug(print_contents(&(Firsthalf[0]),HALFLEN+1));
    init_char_p = false;
  }
}

/* returns skip length */
static int
read_second_half (int *nextchar, char **pointer2a, int *length2a, char **pointer2b, int *length2b,
		  bool eolnp, FILE *fp) {
  int skiplength, ncycles = 0, remainder, terminator, strlenp;
  char *ptr;
  char *p = NULL;
  int c;
  bool init_char_p = true;
  
  ptr = &(Secondhalf[0]);
  remainder = (&(Secondhalf[HALFLEN+2]) - ptr)/sizeof(char);

  while (1) {
    debug(printf("\nEnd: %d\n",remainder));

    if (feof(fp)) {
      debug(printf("EOF.  Returning\n"));
      *nextchar = EOF;
      break;

    } else if (remainder <= 0) {
      ptr = &(Secondhalf[0]);
      remainder = (&(Secondhalf[HALFLEN+2]) - ptr)/sizeof(char);
      ncycles++;
      debug(printf("remainder <= 0.  Cycling\n"));

    } else if (eolnp == true) {
      /* Peek at character after eoln */
      if ((c = fgetc(fp)) == EOF || (init_char_p == false && (c == '@' || c == '>' || c == '<' || c == '+'))) {
	debug(printf("c == EOF or > or < or +.  Returning\n"));
	*nextchar = c;
	break;
      } else if (iscntrl(c)) {
	debug(printf("c == control char.  Continuing\n"));
#ifdef DASH
      } else if (c == DASH) {
	debug(printf("c == dash.  Continuing\n"));
#endif
      } else if (isspace(c)) {
	debug(printf("c == NULL.  Continuing\n"));
      } else {
	*ptr++ = (char) c;
	remainder--;
	eolnp = false;
	p = NULL;
	debug(printf("c == sth.  Continuing\n"));
      }
      
    } else {
      if (p != NULL) {
        strlenp = strlen(p);
	memmove(ptr,p,strlenp);
        ptr[strlenp] = '\0';
	debug(printf("Did memmove of %d chars at %p to %p\n",strlenp,p,ptr));
      } else {
	p = fgets(ptr,remainder+1,fp);
      }
      if (p == NULL) {
	debug(printf("line == NULL.  Returning\n"));
	*nextchar = EOF;
	break;
      } else {
	debug(printf("Read %s.\n",ptr));
	/* Was a call to find_bad_char */
	while ((p = index(ptr,SPACE)) != NULL) {
	  ptr = p++;
	  strlenp = strlen(p);
	  memmove(ptr,p,strlenp);
	  ptr[strlenp] = '\0';
	  debug(printf("Found space.  Did memmove of %d chars at %p to %p\n",strlenp,p,ptr));
	} 
	if (*ptr == '\n') {
	  eolnp = true;
	  debug(printf("line == EOLN.  Continuing\n"));
	} else if ((p = index(ptr,'\n')) != NULL) {
	  if (p[-1] == '\r') {
	    p--;
	  }
	  ptr = p;
	  eolnp = true;
	  debug(printf("line == EOLN.  Continuing\n"));
	} else {
	  ptr += strlen(ptr);
	  eolnp = false;
	  p = NULL;
	  debug(printf("line != EOLN.  Continuing\n"));
	}
	remainder = (&(Secondhalf[HALFLEN+2]) - ptr)/sizeof(char);
      }
    }

    debug(print_contents(&(Secondhalf[0]),HALFLEN+3));
    init_char_p = false;
  }

  terminator = (ptr - &(Secondhalf[0]))/sizeof(char);
  debug(printf("ncycles = %d, terminator is %d\n",ncycles,terminator));
  if (ncycles == 0) {
    *length2a = 0;
    if (terminator < HALFLEN) {
      skiplength = 0;
    } else {
      skiplength = terminator-HALFLEN;
    }
  } else {
    *length2a = HALFLEN-terminator;
    skiplength = ncycles*(HALFLEN+2) + terminator-HALFLEN;
  }
  if (*length2a <= 0) {
    *length2a = 0;
    *pointer2a = (char *) NULL;
  } else {
    *pointer2a = &(Secondhalf[HALFLEN+2-(*length2a)]);
  }
  if (terminator == 0) {
    *length2b = 0;
    *pointer2b = (char *) NULL;
  } else if (terminator > HALFLEN) {
    *length2b = HALFLEN;
    *pointer2b = &(ptr[-(*length2b)]);
  } else {
    *length2b = terminator;
    *pointer2b = &(Secondhalf[0]);
  }

  return skiplength;
}


/* Returns sequence length */
static int
input_sequence (int *nextchar, char **pointer1, int *length1, char **pointer2a, int *length2a,
		char **pointer2b, int *length2b, int *skiplength, FILE *fp, bool possible_fasta_header_p) {
  bool eolnp = true;

  *pointer1 = &(Firsthalf[0]);
  *pointer2a = (char *) NULL;
  *length2a = 0;
  *pointer2b = (char *) NULL;
  *length2b = 0;

  /* printf("Beginning input_sequence with nextchar = %c\n",*nextchar); */
  if ((*length1 = read_first_half(&(*nextchar),&eolnp,fp,possible_fasta_header_p)) == 0) {
    *pointer1 = (char *) NULL;
    *skiplength = 0;
  } else if (*length1 < HALFLEN) {
    *skiplength = 0;
  } else {
    *skiplength = read_second_half(&(*nextchar),&(*pointer2a),&(*length2a),
				   &(*pointer2b),&(*length2b),eolnp,fp);
    debug1(printf("read_second_half returns skiplength of %d, length2a=%d, length2b=%d\n",
		  *skiplength,*length2a,*length2b));
  }

  debug1(printf("length1 = %d, length2a = %d, length2b = %d\n",
		*length1,*length2a,*length2b));

  return (*length1) + (*length2a) + (*length2b);
}


/* Used only by extern procedures (outside of this file).  Internal
   procedures have their own specialized creators. */
T
Sequence_genomic_new (char *contents, int length, bool copyp) {
  T new = (T) MALLOC(sizeof(*new));
  char *copy;

  new->acc = (char *) NULL;
  new->restofheader = (char *) NULL;

  new->trimstart = 0;
  new->trimend = new->fulllength = length;
#ifdef PMAP
  new->fulllength_given = length;
#endif

  if (copyp == true) {
    copy = (char *) CALLOC(length+1,sizeof(char));
    strncpy(copy,contents,length);
    new->contents = copy;
    new->contents_alloc = copy;
  } else {
    /* new->free_contents_p = false; */
    new->contents = contents;
    new->contents_alloc = contents;
    /* new->contents_uc_alloc = (char *) NULL; -- only for GSNAP */
  }

#ifndef PMAP
  new->quality = new->quality_alloc = (char *) NULL;
#endif

  new->subseq_offset = 0;
  new->skiplength = 0;
  new->firstp = true;

  return new;
}


static char complCode[128] = COMPLEMENT_LC;

static char *
make_complement (char *sequence, unsigned int length) {
  char *complement;
  int i, j;

  complement = (char *) CALLOC_IN(length+1,sizeof(char));
  for (i = length-1, j = 0; i >= 0; i--, j++) {
    complement[j] = complCode[(int) sequence[i]];
  }
  complement[length] = '\0';
  return complement;
}

static char *
make_reverse (char *sequence, unsigned int length) {
  char *reverse;
  int i, j;

  if (sequence == NULL) {
    return (char *) NULL;
  } else {
    reverse = (char *) CALLOC_IN(length+1,sizeof(char));
    for (i = length-1, j = 0; i >= 0; i--, j++) {
      reverse[j] = sequence[i];
    }
    reverse[length] = '\0';
    return reverse;
  }
}


#if 0
static char *
make_complement_uppercase (char *sequence, unsigned int length) {
  char *complement;
  char uppercaseCode[128] = UPPERCASE_U2T;
  int i, j;

  complement = (char *) CALLOC_IN(length+1,sizeof(char));
  for (i = length-1, j = 0; i >= 0; i--, j++) {
    complement[j] = uppercaseCode[(int) complCode[(int) sequence[i]]];
  }
  complement[length] = '\0';
  return complement;
}
#endif


static void
make_complement_buffered (char *complement, char *sequence, unsigned int length) {
  int i, j;

  /* complement = (char *) CALLOC_IN(length+1,sizeof(char)); */
  for (i = length-1, j = 0; i >= 0; i--, j++) {
    complement[j] = complCode[(int) sequence[i]];
  }
  complement[length] = '\0';
  return;
}


static void
make_complement_inplace (char *sequence, unsigned int length) {
  char temp;
  unsigned int i, j;

  for (i = 0, j = length-1; i < length/2; i++, j--) {
    temp = complCode[(int) sequence[i]];
    sequence[i] = complCode[(int) sequence[j]];
    sequence[j] = temp;
  }
  if (i == j) {
    sequence[i] = complCode[(int) sequence[i]];
  }

  return;
}


/************************************************************************
 *  Original:
 *   TTTTTT ACGT ...... ACGT AAAAAA
 *          ^trimstart     ^trimend
 *   ^contents
 ************************************************************************
 *  Subsequence:
 *       ^start                ^end
 *          ^trimstart     ^trimend
 *       ^contents
 ************************************************************************/

T
Sequence_subsequence (T this, int start, int end) {
  T new;

#ifdef PMAP
  start /= 3;
  end /= 3;
#endif

  if (start < 0) {
    start = 0;
  }
  if (end > this->fulllength) {
    end = this->fulllength;
  }

  if (end <= start) {
    return NULL;
  } else {
    new = (T) MALLOC_IN(sizeof(*new));

    new->acc = (char *) NULL;
    new->restofheader = (char *) NULL;
    new->contents = &(this->contents[start]); 

    new->fulllength = end - start;
#ifdef PMAP
    new->fulllength_given = new->fulllength;
#endif
    if ((new->trimstart = this->trimstart - start) < 0) {
      new->trimstart = 0;
    }
    if ((new->trimend = this->trimend - start) > new->fulllength) {
      new->trimend = new->fulllength;
    }

    /* new->free_contents_p = false; */
    new->contents_alloc = (char *) NULL;
    /* new->contents_uc_alloc = (char *) NULL; -- only for GSNAP */

#ifndef PMAP
    if (this->quality == NULL) {
      new->quality = (char *) NULL;
    } else {
      new->quality = &(this->quality[start]);
    }
    new->quality_alloc = (char *) NULL;
#endif

#ifdef PMAP
    new->subseq_offset = 3*start;
#else
    new->subseq_offset = start;
#endif
    new->skiplength = this->skiplength;

    new->firstp = this->firstp;

    return new;
  }
}


T
Sequence_revcomp (T this) {
  T new = (T) MALLOC_IN(sizeof(*new));

  new->acc = (char *) NULL;
  new->restofheader = (char *) NULL;
  new->contents = new->contents_alloc = make_complement(this->contents,this->fulllength);

#ifndef PMAP
  new->quality = new->quality_alloc = make_reverse(this->quality,this->fulllength);
#endif

  new->fulllength = this->fulllength;
#ifdef PMAP
  new->fulllength_given = this->fulllength_given;
#endif
  new->trimstart = this->trimstart;
  new->trimend = this->trimend;
  /* new->free_contents_p = true; */
  new->subseq_offset = 0;	/* Not sure if this is right */
  new->skiplength = this->skiplength;
  new->firstp = this->firstp;
  return new;
}


static char *
make_uppercase (char *sequence, unsigned int length) {
  char *uppercase;
#ifdef PMAP
  char uppercaseCode[128] = UPPERCASE_STD;
#else
  char uppercaseCode[128] = UPPERCASE_U2T;
#endif
  unsigned int i;

  uppercase = (char *) CALLOC_IN(length+1,sizeof(char));
  for (i = 0; i < length; i++) {
    uppercase[i] = uppercaseCode[(int) sequence[i]];
  }
  uppercase[length] = '\0';
  return uppercase;
}


T
Sequence_uppercase (T this) {
  T new = (T) MALLOC_IN(sizeof(*new));

  new->acc = (char *) NULL;
  new->restofheader = (char *) NULL;
  new->contents = new->contents_alloc = make_uppercase(this->contents,this->fulllength);

#ifndef PMAP
  if (this->quality_alloc == NULL) {
    new->quality = new->quality_alloc = (char *) NULL;
  } else {
    new->quality = new->quality_alloc =(char *) CALLOC_IN(this->fulllength+1,sizeof(char));
    strcpy(new->quality,this->quality);
  }
#endif

  new->fulllength = this->fulllength;
#ifdef PMAP
  new->fulllength_given = this->fulllength_given;
#endif
  new->trimstart = this->trimstart;
  new->trimend = this->trimend;
  /* new->free_contents_p = true; */
  new->subseq_offset = this->subseq_offset;
  new->skiplength = this->skiplength;
  new->firstp = this->firstp;
  return new;
}


T
Sequence_alias (T this) {
  T new = (T) MALLOC_IN(sizeof(*new));

  new->acc = (char *) NULL;
  new->restofheader = (char *) NULL;
  new->contents = this->contents;

#ifndef PMAP
  new->quality = this->quality;
  new->quality_alloc = (char *) NULL;
#endif

  new->fulllength = this->fulllength;
#ifdef PMAP
  new->fulllength_given = this->fulllength_given;
#endif
  new->trimstart = this->trimstart;
  new->trimend = this->trimend;

  /* new->free_contents_p = false; */
  new->contents_alloc = (char *) NULL;
  /* new->contents_uc_alloc = (char *) NULL; -- only for GSNAP */


  new->subseq_offset = this->subseq_offset;
  new->firstp = this->firstp;
  return new;
}


/*
void
Sequence_endstream () {
  Initc = '\0';
  return;
}
*/


T
Sequence_read (int *nextchar, FILE *input) {
  T new;
  int fulllength, skiplength;
  char *pointer1, *pointer2a, *pointer2b;
  int length1, length2a, length2b;
#ifdef PMAP
  char lastchar = '*';
#else
  int quality_length;
#endif

  if (feof(input)) {
    *nextchar = EOF;
    return NULL;
  }

  if (*nextchar == '\0') {
    if ((*nextchar = Sequence_input_init(input)) == EOF) {
      *nextchar = EOF;
      return NULL;
    }
  }

  new = (T) MALLOC_IN(sizeof(*new));

  if (*nextchar != '@' && *nextchar != '>' && *nextchar != '<') {
    new->firstp = true;		/* by default */
    blank_header(new);
  } else if (input_header(input,new) == NULL) {
    /* File ends after >.  Don't process. */
    *nextchar = EOF;
    FREE_IN(new);
    return NULL;
  } else if (*nextchar == '@') {
    new->firstp = true;		/* by default */
  } else if (*nextchar == '>') {
    new->firstp = true;
  } else if (*nextchar == '<') {
    new->firstp = false;
  } else {
    abort();
  }

  if ((fulllength = input_sequence(&(*nextchar),&pointer1,&length1,&pointer2a,&length2a,
				   &pointer2b,&length2b,&skiplength,input,/*possible_fasta_header_p*/true)) == 0) {
    /* File ends during header.  Continue with a sequence of length 0. */
    /* fprintf(stderr,"File ends after header\n"); */
  }

  if (skiplength > 0) {
    fprintf(stderr,"Warning: cDNA sequence length of %d exceeds maximum length of %d.  Truncating %d chars in middle.\n",
	    fulllength+skiplength,MAXSEQLEN,skiplength);
    fprintf(stderr,"  (For long sequences, perhaps you want maponly mode, by providing the '-1' flag.)\n");
  }

#ifdef PMAP
  if (length1 > 0) {
    lastchar = pointer1[length1-1];
    if (length2a > 0) {
      lastchar = pointer2a[length2a-1];
    }
    if (length2b > 0) {
      lastchar = pointer2b[length2b-1];
    }
  }

  new->fulllength_given = fulllength;
  if (lastchar != '*') {
    debug(printf("Sequence does not end with *, so adding it\n"));
    fulllength++;
  }
#endif

  debug(printf("fulllength = %d\n",fulllength));
  new->fulllength = fulllength;
  new->skiplength = skiplength;

  new->trimstart = 0;
#ifdef PMAP
  new->trimend = new->fulllength_given;
#else
  new->trimend = fulllength;
#endif

  new->contents = new->contents_alloc = (char *) CALLOC_IN(fulllength+1,sizeof(char));
  if (length1 > 0) {
    strncpy(new->contents,pointer1,length1);
    if (length2a > 0) {
      strncpy(&(new->contents[length1]),pointer2a,length2a);
    }
    if (length2b > 0) {
      strncpy(&(new->contents[length1+length2a]),pointer2b,length2b);
    }
  }

#ifdef PMAP
  if (lastchar != '*') {
    new->contents[fulllength-1] = '*';
  }
#endif
  /* new->free_contents_p = true; */
  new->subseq_offset = 0;

#ifndef PMAP
  /* Quality string */
  new->quality = new->quality_alloc = (char *) NULL;
  if (*nextchar == '+') {
    skip_header(input);
    *nextchar = fgetc(input);
    quality_length = input_sequence(&(*nextchar),&pointer1,&length1,&pointer2a,&length2a,
				    &pointer2b,&length2b,&skiplength,input,/*possible_fasta_header_p*/false);
    if (quality_length != fulllength) {
      fprintf(stderr,"Length %d of quality score differs from length %d of nucleotides in sequence %s\n",
	      quality_length,fulllength,new->acc);
      exit(9);
    } else {
      new->quality = new->quality_alloc = (char *) CALLOC_IN(fulllength+1,sizeof(char));
      if (length1 > 0) {
	strncpy(new->quality,pointer1,length1);
	if (length2a > 0) {
	  strncpy(&(new->quality[length1]),pointer2a,length2a);
	}
	if (length2b > 0) {
	  strncpy(&(new->quality[length1+length2a]),pointer2b,length2b);
	}
      }
    }
  }
#endif

  debug(printf("Final query sequence is:\n"));
  debug(Sequence_print(stdout,new,/*uppercasep*/false,/*wraplength*/60,/*trimmedp*/false));
  return new;
}


T
Sequence_read_multifile (int *nextchar, FILE **input, char ***files, int *nfiles) {
  T queryseq;

  while (1) {
    if (*input == NULL || feof(*input)) {
      if (*input != NULL) {
	fclose(*input);
	*input = NULL;
      }

      if (*nfiles == 0) {
	*nextchar = EOF;
	return NULL;
      } else {
	while (*nfiles > 0 && (*input = FOPEN_READ_TEXT((*files)[0])) == NULL) {
	  fprintf(stderr,"Can't open file %s => skipping it.\n",(*files)[0]);
	  (*files)++;
	  (*nfiles)--;
	}
	if (*input == NULL) {
	  *nextchar = EOF;
	  return NULL;
	} else {
	  (*files)++;
	  (*nfiles)--;
	  *nextchar = '\0';
	}
      }
    }
    if ((queryseq = Sequence_read(&(*nextchar),*input)) != NULL) {
      return queryseq;
    }
  }
}


/* This is intended for user genomicseg, which will start with
   standard FASTA ">" */
T
Sequence_read_unlimited (int *nextchar, FILE *input) {
  T new;
  Intlist_T intlist = NULL;
  char *p;
  int length, startpos = 1, maxseqlen = MAXSEQLEN;
  bool eolnp;

  if (feof(input)) {
    *nextchar = EOF;
    return NULL;
  }

  if (*nextchar == '\0') {
    if ((*nextchar = Sequence_input_init(input)) == EOF) {
      *nextchar = EOF;
      return NULL;
    }
  }

  new = (T) MALLOC_IN(sizeof(*new));

  if (*nextchar != '>' /* && *nextchar != '<' */) {
    new->firstp = true;		/* by default */
    blank_header(new);
    Sequence[startpos++] = (char) *nextchar;
    maxseqlen--;
  } else if (input_header(input,new) == NULL) {
    /* File ends after >.  Don't process. */
    FREE_IN(new);
    return NULL;
  } else if (*nextchar == '>') {
    new->firstp = true;
#if 0
  } else if (*nextchar == '<') {
    new->firstp = false;
#endif
  } else {
    abort();
  }

  /* Don't touch Sequence[0], because subsequent calls to
     Sequence_read depend on it being '\0'. */
  eolnp = true;
  while (fgets(&(Sequence[startpos]),maxseqlen,input) != NULL &&
	 (eolnp == false || (Sequence[1] != '>' /* && Sequence[1] != '<' */))) {
    for (p = &(Sequence[1]); *p != '\r' && *p != '\n' && *p != '\0'; p++) {
      if (!iscntrl((int) *p)
#ifdef DASH
	  && *p != DASH
#endif
	  ) {
	intlist = Intlist_push(intlist,(int) *p);
      }
    }
    if (*p == '\r' || *p == '\n') {
      eolnp = true;
    } else {
      eolnp = false;
    }
    startpos = 1;
    maxseqlen = MAXSEQLEN;
  }
  intlist = Intlist_reverse(intlist);
  new->contents = new->contents_alloc = Intlist_to_char_array(&length,intlist);

  Intlist_free(&intlist);

  if (length == 0) {
    return NULL;
  } else {
    new->fulllength = new->trimend = length;
#ifdef PMAP
    new->fulllength_given = length;
#endif
    new->trimstart = 0;

#ifndef PMAP
    /* user genomic segment should not have quality */
    new->quality = new->quality_alloc = (char *) NULL;
#endif

    /* new->free_contents_p = true; */
    new->subseq_offset = 0;
    new->skiplength = 0;

    /* Important to initialize for subsequent cDNA reads */
    *nextchar = '\0';

    return new;
  }
}


void
Sequence_print_digest (Filestring_T fp, T this) {
  unsigned char *digest;

  digest = MD5_compute((unsigned char *) this->contents,this->fulllength);
  MD5_print(fp,digest);
  FREE(digest);
  return;
}

/* Calling procedure needs to print the initial ">", if desired */
void
Sequence_print_header (Filestring_T fp, T this, bool checksump) {

  if (this->acc == NULL) {
    FPRINTF(fp,"NO_HEADER");
  } else {
    if (this->restofheader == NULL || this->restofheader[0] == '\0') {
      FPRINTF(fp,"%s",this->acc);
    } else {
      FPRINTF(fp,"%s %s",this->acc,this->restofheader);
    }

    if (checksump == true) {
      FPRINTF(fp," md5:");
      Sequence_print_digest(fp,this);
    }
  }

  FPRINTF(fp,"\n");

  return;
}

#if 0
/* Used by revcomp.c */
void
Sequence_stdout_header_revcomp (T this) {
  if (this->restofheader == NULL || this->restofheader[0] == '\0') {
    printf(">%s",this->acc);
  } else {
    printf(">%s %s",this->acc,this->restofheader);
  }
  printf(" REVCOMP");
  printf("\n");
  return;
}
#endif


void
Sequence_print (Filestring_T fp, T this, bool uppercasep, int wraplength, bool trimmedp) {
  int i = 0, pos, start, end;
  char uppercaseCode[128] = UPPERCASE_STD;

  if (trimmedp == true) {
    start = this->trimstart;
    end = this->trimend;
  } else {
    start = 0;
    end = this->fulllength;
  }

  if (uppercasep == true) {
    for (pos = start; pos < end; pos++, i++) {
      PUTC(uppercaseCode[(int) this->contents[i]],fp);
      if ((i+1) % wraplength == 0) {
	PUTC('\n',fp);
      }
    }
  } else {
    for (pos = start; pos < end; pos++, i++) {
      PUTC(this->contents[i],fp);
      if ((i+1) % wraplength == 0) {
	PUTC('\n',fp);
      }
    }
  }
  if (i % wraplength != 0) {
    PUTC('\n',fp);
  }

  return;
}


void
Sequence_stdout (T this, bool uppercasep, int wraplength, bool trimmedp) {
  int i = 0, pos, start, end;
  char uppercaseCode[128] = UPPERCASE_STD;

  if (trimmedp == true) {
    start = this->trimstart;
    end = this->trimend;
  } else {
    start = 0;
    end = this->fulllength;
  }

  if (uppercasep == true) {
    for (pos = start; pos < end; pos++, i++) {
      putchar(uppercaseCode[(int) this->contents[i]]);
      if ((i+1) % wraplength == 0) {
	putchar('\n');
      }
    }
  } else {
    for (pos = start; pos < end; pos++, i++) {
      putchar(this->contents[i]);
      if ((i+1) % wraplength == 0) {
	putchar('\n');
      }
    }
  }
  if (i % wraplength != 0) {
    putchar('\n');
  }

  return;
}


void
Sequence_stdout_alt (T ref, T alt, T snp, bool uppercasep, int wraplength) {
  int i = 0, pos, start, end;
  char uppercaseCode[128] = UPPERCASE_STD;

  start = 0;
  end = alt->fulllength;

  if (uppercasep == true) {
    for (pos = start; pos < end; pos++, i++) {
      if (snp->contents[i] == ' ') {
	/* Not a SNP, so print reference */
	printf("%c",uppercaseCode[(int) ref->contents[i]]);
      } else if (uppercaseCode[(int) alt->contents[i]] == uppercaseCode[(int) ref->contents[i]]) {
	/* Wildcard SNP */
	printf("N");
      } else {
	printf("%c",uppercaseCode[(int) alt->contents[i]]);
      }
      if ((i+1) % wraplength == 0) {
	printf("\n");
      }
    }
  } else {
    for (pos = start; pos < end; pos++, i++) {
      if (snp->contents[i] == ' ') {
	/* Not a SNP, so print reference */
	printf("%c",ref->contents[i]);
      } else if (uppercaseCode[(int) alt->contents[i]] == uppercaseCode[(int) ref->contents[i]]) {
	/* Wildcard SNP */
	printf("N");
      } else {
	printf("%c",alt->contents[i]);
      }
      if ((i+1) % wraplength == 0) {
	printf("\n");
      }
    }
  }
  if (i % wraplength != 0) {
    printf("\n");
  }
  return;
}


void
Sequence_stdout_two (T ref, T alt, bool uppercasep, int wraplength) {
  int i = 0, pos, pos2, startpos, end;
  char uppercaseCode[128] = UPPERCASE_STD;

  end = ref->fulllength;

  pos = 0;
  i = 0;
  if (uppercasep == true) {
    printf("ref\t");
    startpos = pos;
    while (pos < end) {
      printf("%c",uppercaseCode[(int) ref->contents[pos]]);
      if (++i % wraplength == 0) {
	printf("\n");
	printf("alt\t");
	for (pos2 = startpos, i = 0; i < wraplength; pos2++, i++) {
	  if (uppercaseCode[(int) alt->contents[pos2]] == uppercaseCode[(int) ref->contents[pos2]]) {
	    /* Wildcard SNP */
	    printf("N");
	  } else {
	    printf("%c",uppercaseCode[(int) alt->contents[pos2]]);
	  }
	}
	printf("\n\n");
	if (pos+1 < end) {
	  printf("ref\t");
	}
	startpos = pos+1;
      }
      pos++;
    }
    if (i % wraplength != 0) {
      printf("\n");
      printf("alt\t");
      for (pos2 = startpos; pos2 < end; pos2++) {
	if (uppercaseCode[(int) alt->contents[pos2]] == uppercaseCode[(int) ref->contents[pos2]]) {
	  /* Wildcard SNP */
	  printf("N");
	} else {
	  printf("%c",uppercaseCode[(int) alt->contents[pos2]]);
	}
      }
      printf("\n\n");
    }

  } else {
    printf("ref\t");
    startpos = pos;
    while (pos < end) {
      printf("%c",ref->contents[pos]);
      if (++i % wraplength == 0) {
	printf("\n");
	printf("alt\t");
	for (pos2 = startpos, i = 0; i < wraplength; pos2++, i++) {
	  if (uppercaseCode[(int) alt->contents[pos2]] == uppercaseCode[(int) ref->contents[pos2]]) {
	    /* Wildcard SNP */
	    printf("N");
	  } else {
	    printf("%c",alt->contents[pos2]);
	  }
	}
	printf("\n\n");
	if (pos+1 < end) {
	  printf("ref\t");
	}
	startpos = pos+1;
      }
      pos++;
    }
    if (i % wraplength != 0) {
      printf("\n");
      printf("alt\t");
      for (pos2 = startpos; pos2 < end; pos2++) {
	if (uppercaseCode[(int) alt->contents[pos2]] == uppercaseCode[(int) ref->contents[pos2]]) {
	  /* Wildcard SNP */
	  printf("N");
	} else {
	  printf("%c",alt->contents[pos2]);
	}
      }
      printf("\n\n");
    }
  }

  return;
}


void
Sequence_stdout_raw (T this) {
  int i = 0, pos, start, end;

  start = 0;
  end = this->fulllength;

  for (pos = start; pos < end; pos++, i++) {
    printf("%d\n",(int) this->contents[i]);
  }
  return;
}

void
Sequence_stdout_stream_chars (T this) {
  int i = 0, pos, start, end;

  start = 0;
  end = this->fulllength;

  for (pos = start; pos < end; pos++, i++) {
    switch (toupper(this->contents[i])) {
    case 'A': putchar('A'); break;
    case 'C': putchar('C'); break;
    case 'G': putchar('G'); break;
    case 'T': putchar('T'); break;
    default: putchar('X');
    }
  }
  return;
}

void
Sequence_stdout_stream_ints (T this) {
  int i = 0, pos, start, end;

  start = 0;
  end = this->fulllength;

  for (pos = start; pos < end; pos++, i++) {
    switch (toupper(this->contents[i])) {
    case 'A': putchar(0); break;
    case 'C': putchar(1); break;
    case 'G': putchar(2); break;
    case 'T': putchar(3); break;
    default: putchar(4);
    }
  }
  return;
}


T
Sequence_substring (T usersegment, unsigned int left, unsigned int length, 
		    bool revcomp) {
  char *gbuffer;

  gbuffer = (char *) CALLOC(length+1,sizeof(char));

  memcpy(gbuffer,&(usersegment->contents[left]),length*sizeof(char));
  gbuffer[length] = '\0';

  if (revcomp == true) {
    /* make_complement_buffered(gbuffer2,gbuffer1,length); */
    make_complement_inplace(gbuffer,length);
    debug(fprintf(stderr,"Got sequence at %u with length %u, revcomp\n",left,length));
    return Sequence_genomic_new(gbuffer,length,/*copyp*/false);
  } else {
    debug(fprintf(stderr,"Got sequence at %u with length %u, forward\n",left,length));
    return Sequence_genomic_new(gbuffer,length,/*copyp*/false);
  }
}


