static char rcsid[] = "$Id: diagpool.c 166641 2015-05-29 21:13:04Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "diagpool.h"
#include <stdio.h>
#include <stdlib.h>
#include "mem.h"
#include "comp.h"
#include "diagdef.h"
#include "listdef.h"

#define CHUNKSIZE 20000

#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif

/* For mechanics of memory allocation and deallocation */
#ifdef DEBUG1
#define debug1(x) x
#else
#define debug1(x)
#endif

/* For popping */
#ifdef DEBUG2
#define debug2(x) x
#else
#define debug2(x)
#endif


#define T Diagpool_T
struct T {
  int ndiags;
  int diagctr;
  struct Diag_T *diagptr;
  List_T diagchunks;

  int nlistcells;
  int listcellctr;
  struct List_T *listcellptr;
  List_T listcellchunks;
};

void
Diagpool_free (T *old) {
  List_T p;
  struct Diag_T *diagptr;
  struct List_T *listcellptr;

  if (*old) {
    for (p = (*old)->diagchunks; p != NULL; p = List_next(p)) {
      diagptr = (struct Diag_T *) List_head(p);
      FREE(diagptr);
    }
    List_free(&(*old)->diagchunks);
    for (p = (*old)->listcellchunks; p != NULL; p = List_next(p)) {
      listcellptr = (struct List_T *) List_head(p);
      FREE(listcellptr);
    }
    List_free(&(*old)->listcellchunks);
    FREE(*old);
  }
  return;
}

void
Diagpool_free_memory (T this) {
  List_T p;
  struct Diag_T *diagptr;
  struct List_T *listcellptr;

  for (p = this->diagchunks; p != NULL; p = List_next(p)) {
    diagptr = (struct Diag_T *) List_head(p);
    FREE_KEEP(diagptr);
  }
  List_free_keep(&this->diagchunks);
  for (p = this->listcellchunks; p != NULL; p = List_next(p)) {
    listcellptr = (struct List_T *) List_head(p);
    FREE_KEEP(listcellptr);
  }
  List_free_keep(&this->listcellchunks);

  this->ndiags = 0;
  this->diagctr = 0;
  this->diagchunks = NULL;
  /* this->diagptr = add_new_diagchunk(this); */

  this->nlistcells = 0;
  this->listcellctr = 0;
  this->listcellchunks = NULL;
  /* this->listcellptr = add_new_listcellchunk(this); */

  return;
}

void
Diagpool_report_memory (T this) {
  printf("Diagpool has %d diagchunks and %d listcellchunks\n",
	 List_length(this->diagchunks),List_length(this->listcellchunks));
  return;
}

static struct Diag_T *
add_new_diagchunk (T this) {
  struct Diag_T *chunk;

  chunk = (struct Diag_T *) MALLOC_KEEP(CHUNKSIZE*sizeof(struct Diag_T));
  this->diagchunks = List_push_keep(this->diagchunks,(void *) chunk);
  debug1(printf("Adding a new chunk of diags.  Ptr for diag %d is %p\n",
		this->ndiags,chunk));

  this->ndiags += CHUNKSIZE;

  return chunk;
}

static struct List_T *
add_new_listcellchunk (T this) {
  struct List_T *chunk;

  chunk = (struct List_T *) MALLOC_KEEP(CHUNKSIZE*sizeof(struct List_T));
  this->listcellchunks = List_push_keep(this->listcellchunks,(void *) chunk);
  debug1(printf("Adding a new chunk of listcells.  Ptr for listcell %d is %p\n",
	       this->nlistcells,chunk));

  this->nlistcells += CHUNKSIZE;

  return chunk;
}

T
Diagpool_new (void) {
  T new = (T) MALLOC(sizeof(*new));

  new->ndiags = 0;
  new->diagctr = 0;
  new->diagchunks = NULL;
  /* new->diagptr = add_new_diagchunk(new); */

  new->nlistcells = 0;
  new->listcellctr = 0;
  new->listcellchunks = NULL;
  /* new->listcellptr = add_new_listcellchunk(new); */

  return new;
}

void
Diagpool_reset (T this) {
  this->diagctr = 0;
  this->listcellctr = 0;
  return;
}

List_T
Diagpool_push (List_T list, T this, int diagonal, int querystart, int queryend, int nconsecutive) {
  List_T listcell;
  Diag_T diag;
  List_T p;
  int n;

  if (this->diagctr >= this->ndiags) {
    this->diagptr = add_new_diagchunk(this);
  } else if ((this->diagctr % CHUNKSIZE) == 0) {
    for (n = this->ndiags - CHUNKSIZE, p = this->diagchunks;
	 n > this->diagctr; p = p->rest, n -= CHUNKSIZE) ;
    this->diagptr = (struct Diag_T *) p->first;
    debug1(printf("Located diag %d at %p\n",this->diagctr,this->diagptr));
  }    
  diag = this->diagptr++;
  this->diagctr++;

  diag->diagonal = diagonal;
  diag->querystart = querystart;
  diag->queryend = queryend;
  diag->nconsecutive = nconsecutive;
  diag->dominatedp = false;
  diag->score = 0.0;

  debug(printf("Creating %p: %d %d..%d\n",diag,diag->diagonal,diag->querystart,diag->queryend));

  if (this->listcellctr >= this->nlistcells) {
    this->listcellptr = add_new_listcellchunk(this);
  } else if ((this->listcellctr % CHUNKSIZE) == 0) {
    for (n = this->nlistcells - CHUNKSIZE, p = this->listcellchunks;
	 n > this->listcellctr; p = p->rest, n -= CHUNKSIZE) ;
    this->listcellptr = (struct List_T *) p->first;
    debug1(printf("Located listcell %d at %p\n",this->listcellctr,this->listcellptr));
  }
  listcell = this->listcellptr++;
  this->listcellctr++;

  listcell->first = (void *) diag;
  listcell->rest = list;

  return listcell;
}


List_T
Diagpool_pop (List_T list, Diag_T *x) {
  List_T head;

  if (list != NULL) {
    head = list->rest;
    *x = (Diag_T) list->first;
    return head;
  } else {
    return list;
  }
}


List_T
Diagpool_push_existing (List_T list, T this, Diag_T diag) {
  List_T listcell;
  List_T p;
  int n;

  if (this->listcellctr >= this->nlistcells) {
    this->listcellptr = add_new_listcellchunk(this);
  } else if ((this->listcellctr % CHUNKSIZE) == 0) {
    for (n = this->nlistcells - CHUNKSIZE, p = this->listcellchunks;
	 n > this->listcellctr; p = p->rest, n -= CHUNKSIZE) ;
    this->listcellptr = (struct List_T *) p->first;
    debug1(printf("Located listcell %d at %p\n",this->listcellctr,this->listcellptr));
  }
  listcell = this->listcellptr++;
  this->listcellctr++;

  listcell->first = (void *) diag;
  listcell->rest = list;

  return listcell;
}
