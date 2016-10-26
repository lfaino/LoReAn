static char rcsid[] = "$I$";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "univdiag.h"
#include "univdiagdef.h"
#include "mem.h"
#include <stdio.h>
#include <stdlib.h>


#define T Univdiag_T


T
Univdiag_new (int querystart, int queryend, Univcoord_T univdiagonal) {
  T new = (T) MALLOC(sizeof(*new));

  new->univdiagonal = univdiagonal;
  new->querystart = querystart;
  new->queryend = queryend;
  new->nconsecutive = queryend - querystart + 1;
  new->nmismatches_known_p = true;

  new->intscore = 0;
  new->nlinked = 0;
  new->prev = (Univdiag_T) NULL;

  return new;
}


T
Univdiag_new_fillin (int querystart, int queryend, int indexsize, Univcoord_T univdiagonal) {
  T new = (T) MALLOC(sizeof(*new));

  new->univdiagonal = univdiagonal;
  new->querystart = querystart;
  new->queryend = queryend + indexsize - 1;
  new->nconsecutive = new->queryend - querystart + 1;
  new->nmismatches_known_p = false;

  new->intscore = 0;
  new->nlinked = 0;
  new->prev = (Univdiag_T) NULL;

  return new;
}


void
Univdiag_free (T *old) {
  FREE(*old);
  return;
}

void
Univdiag_gc (List_T *list) {
  T univdiagonal;
  List_T p;

  for (p = *list; p != NULL; p = List_next(p)) {
    univdiagonal = (T) List_head(p);
    FREE(univdiagonal);
  }
  List_free(&(*list));
  return;
}


int
Univdiag_ascending_cmp (const void *a, const void *b) {
  T x = * (T *) a;
  T y = * (T *) b;

  if (x->querystart < y->querystart) {
    return -1;
  } else if (y->querystart < x->querystart) {
    return +1;
  } else if (x->queryend < y->queryend) {
    return -1;
  } else if (y->queryend < x->queryend) {
    return +1;
  } else if (x->univdiagonal < y->univdiagonal) {
    return -1;
  } else if (y->univdiagonal < x->univdiagonal) {
    return +1;
  } else {
    return 0;
  }
}


int
Univdiag_descending_cmp (const void *a, const void *b) {
  T x = * (T *) a;
  T y = * (T *) b;

  if (x->querystart > y->querystart) {
    return -1;
  } else if (y->querystart > x->querystart) {
    return +1;
  } else if (x->queryend > y->queryend) {
    return -1;
  } else if (y->queryend > x->queryend) {
    return +1;
  } else if (x->univdiagonal > y->univdiagonal) {
    return -1;
  } else if (y->univdiagonal > x->univdiagonal) {
    return +1;
  } else {
    return 0;
  }
}


