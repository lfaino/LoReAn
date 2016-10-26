static char rcsid[] = "$Id: genomicpos.c 155282 2014-12-12 19:42:54Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#ifndef HAVE_MEMCPY
# define memcpy(d,s,n) bcopy((s),(d),(n))
#endif
#ifndef HAVE_MEMMOVE
# define memmove(d,s,n) bcopy((s),(d),(n))
#endif

#include "genomicpos.h"
#include <string.h>
#include "mem.h"

#define T Genomicpos_T

/* Adapted from public domain code by Bob Stout and Mark Kamradt */

#define BUFSIZE 100

/* Needs to handle large files */
char *
Genomicpos_commafmt (
#ifdef HAVE_64_BIT
		     UINT8 N
#else
		     UINT4 N
#endif
		     ) {
  char *string, *buffer;
  int len, posn = 1;
  char *ptr, *start;

  buffer = (char *) CALLOC(BUFSIZE+1,sizeof(char));
  start = ptr = &(buffer[BUFSIZE]);
  buffer[BUFSIZE] = '\0';

  if (N == 0UL) {
    *--ptr = '0';
  } else {
    while (N > 0UL) {
      *--ptr = (char)((N % 10UL) + '0');
      N /= 10UL;
      if (N > 0UL) {
	if ((posn % 3) == 0) {
	  *--ptr = ',';
	}
      }
      posn++;
    }
  }

  len = start - ptr;		/* Not including terminal '\0'. */
  string = (char *) CALLOC(len+1,sizeof(char));
  memcpy(string,ptr,len+1);
  FREE(buffer);
  return string;
}


#ifdef MEMUSAGE
/* Does not allocate memory.  Used for reporting MEMUSAGE results. */
void
Genomicpos_commafmt_fill (char *string,
#ifdef HAVE_64_BIT
		     UINT8 N
#else
		     UINT4 N
#endif
		     ) {
  char *buffer;
  int len, posn = 1;
  char *ptr, *start;

  buffer = (char *) CALLOC(BUFSIZE+1,sizeof(char));
  start = ptr = &(buffer[BUFSIZE]);
  buffer[BUFSIZE] = '\0';

  if (N == 0UL) {
    *--ptr = '0';
  } else {
    while (N > 0UL) {
      *--ptr = (char)((N % 10UL) + '0');
      N /= 10UL;
      if (N > 0UL) {
	if ((posn % 3) == 0) {
	  *--ptr = ',';
	}
      }
      posn++;
    }
  }

  len = start - ptr;		/* Not including terminal '\0'. */
  memcpy(string,ptr,len+1);
  FREE(buffer);
  return;
}
#endif


int
UINT8_compare (const void *a, const void *b) {
  UINT8 x = * (UINT8 *) a;
  UINT8 y = * (UINT8 *) b;

  if (x < y) {
    return -1;
  } else if (y < x) {
    return 1;
  } else {
    return 0;
  }
}

int
UINT4_compare (const void *a, const void *b) {
  UINT4 x = * (UINT4 *) a;
  UINT4 y = * (UINT4 *) b;

  if (x < y) {
    return -1;
  } else if (y < x) {
    return 1;
  } else {
    return 0;
  }
}


int
Univcoord_compare (const void *a, const void *b) {
  Univcoord_T x = * (Univcoord_T *) a;
  Univcoord_T y = * (Univcoord_T *) b;

  if (x < y) {
    return -1;
  } else if (y < x) {
    return 1;
  } else {
    return 0;
  }
}

int
Chrpos_compare (const void *a, const void *b) {
  Chrpos_T x = * (Chrpos_T *) a;
  Chrpos_T y = * (Chrpos_T *) b;

  if (x < y) {
    return -1;
  } else if (y < x) {
    return 1;
  } else {
    return 0;
  }
}

