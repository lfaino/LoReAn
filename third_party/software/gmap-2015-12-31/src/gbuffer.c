static char rcsid[] = "$Id: gbuffer.c 40271 2011-05-28 02:29:18Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "gbuffer.h"

#include <stddef.h>
#include <stdlib.h>
#include "mem.h"
#include "pair.h"
#include "sequence.h"


#define T Gbuffer_T
struct T {
  int gbufferlen;		/* Actual gbufferlen */

  char *chars1;			/* For retrieval of genomic sequence after stage 1 */
  char *chars2;
  char *chars3;
};


int
Gbuffer_gbufferlen (T this) {
  return this->gbufferlen;
}

char *
Gbuffer_chars1 (T this) {
  return this->chars1;
}

char *
Gbuffer_chars2 (T this) {
  return this->chars2;
}

char *
Gbuffer_chars3 (T this) {
  return this->chars3;
}

void
Gbuffer_free_contents (T this) {
  FREE(this->chars3);
  FREE(this->chars2);
  FREE(this->chars1);
  return;
}

void
Gbuffer_alloc_contents (T this, int gbufferlen) {

  if (this->chars1 != NULL) {
    FREE(this->chars3);
    FREE(this->chars2);
    FREE(this->chars1);
  }

  this->chars1 = (char *) CALLOC(gbufferlen+1,sizeof(char));
  this->chars2 = (char *) CALLOC(gbufferlen+1,sizeof(char));
  this->chars3 = (char *) CALLOC(gbufferlen+1,sizeof(char));
  
  this->gbufferlen = gbufferlen;

  return;
}


void
Gbuffer_free (T *old) {
  if (*old) {
    Gbuffer_free_contents(*old);
    FREE(*old);
  }
  return;
}

T
Gbuffer_new () {
  T new = (T) MALLOC(sizeof(*new));

  new->chars1 = (char *) NULL;
  new->chars2 = (char *) NULL;
  new->chars3 = (char *) NULL;

  return new;
}

