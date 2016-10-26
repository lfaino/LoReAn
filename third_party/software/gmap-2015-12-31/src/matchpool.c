static char rcsid[] = "$Id: matchpool.c 99737 2013-06-27 19:33:03Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "matchpool.h"
#include <stdio.h>
#include <stdlib.h>
#include "mem.h"
#include "comp.h"
#include "matchdef.h"
#include "listdef.h"
#include "univinterval.h"


#define CHUNKSIZE 16384

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


#define T Matchpool_T
struct T {
  int nmatches;
  int matchctr;
  struct Match_T *matchptr;
  List_T matchchunks;

  int matchctr_save;
  struct Match_T *matchptr_save;

  int nlistcells;
  int listcellctr;
  struct List_T *listcellptr;
  List_T listcellchunks;

  int listcellctr_save;
  struct List_T *listcellptr_save;
};

void
Matchpool_free (T *old) {
  List_T p;
  struct Match_T *matchptr;
  struct List_T *listcellptr;

  if (*old) {
    for (p = (*old)->matchchunks; p != NULL; p = List_next(p)) {
      matchptr = (struct Match_T *) List_head(p);
      FREE(matchptr);
    }
    List_free(&(*old)->matchchunks);
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
Matchpool_free_memory (T this) {
  List_T p;
  struct Match_T *matchptr;
  struct List_T *listcellptr;

  for (p = this->matchchunks; p != NULL; p = List_next(p)) {
    matchptr = (struct Match_T *) List_head(p);
    FREE(matchptr);
  }
  List_free(&this->matchchunks);
  for (p = this->listcellchunks; p != NULL; p = List_next(p)) {
    listcellptr = (struct List_T *) List_head(p);
    FREE(listcellptr);
  }
  List_free(&this->listcellchunks);

  this->nmatches = 0;
  this->matchctr = 0;
  this->matchchunks = NULL;
  /* this->matchptr = add_new_matchchunk(this); */

  this->nlistcells = 0;
  this->listcellctr = 0;
  this->listcellchunks = NULL;
  /* this->listcellptr = add_new_listcellchunk(this); */

  return;
}


static struct Match_T *
add_new_matchchunk (T this) {
  struct Match_T *chunk;

  chunk = (struct Match_T *) MALLOC(CHUNKSIZE*sizeof(struct Match_T));
  this->matchchunks = List_push(this->matchchunks,(void *) chunk);
  debug1(printf("Adding a new chunk of matches.  Ptr for match %d is %p\n",
		this->nmatches,chunk));

  this->nmatches += CHUNKSIZE;
  return chunk;
}

static struct List_T *
add_new_listcellchunk (T this) {
  struct List_T *chunk;

  chunk = (struct List_T *) MALLOC(CHUNKSIZE*sizeof(struct List_T));
  this->listcellchunks = List_push(this->listcellchunks,(void *) chunk);
  debug1(printf("Adding a new chunk of listcells.  Ptr for listcell %d is %p\n",
	       this->nlistcells,chunk));

  this->nlistcells += CHUNKSIZE;
  return chunk;
}

T
Matchpool_new (void) {
  T new = (T) MALLOC(sizeof(*new));

  new->nmatches = 0;
  new->matchctr = 0;
  new->matchchunks = NULL;
  /* new->matchptr = add_new_matchchunk(new); */

  new->nlistcells = 0;
  new->listcellctr = 0;
  new->listcellchunks = NULL;
  /* new->listcellptr = add_new_listcellchunk(new); */

  return new;
}

void
Matchpool_reset (T this) {
  this->matchctr = 0;
  this->listcellctr = 0;
  return;
}

void
Matchpool_save (T this) {
  this->matchctr_save = this->matchctr;
  this->matchptr_save = this->matchptr;
  this->listcellctr_save = this->listcellctr;
  this->listcellptr_save = this->listcellptr;
  return;
}

void
Matchpool_restore (T this) {
  this->matchctr = this->matchctr_save;
  this->matchptr = this->matchptr_save;
  this->listcellctr = this->listcellctr_save;
  this->listcellptr = this->listcellptr_save;
  return;
}

List_T
Matchpool_push (List_T list, T this, int querypos, int querylength, bool forwardp, bool fivep,
		Univcoord_T diagonal, Univ_IIT_T chromosome_iit) {
  List_T listcell;
  Match_T match;
  List_T p;
  int n;
  int index;

  if (this->matchctr >= this->nmatches) {
    this->matchptr = add_new_matchchunk(this);
  } else if ((this->matchctr % CHUNKSIZE) == 0) {
    for (n = this->nmatches - CHUNKSIZE, p = this->matchchunks;
	 n > this->matchctr; p = p->rest, n -= CHUNKSIZE) ;
    this->matchptr = (struct Match_T *) p->first;
    debug1(printf("Located match %d at %p\n",this->matchctr,this->matchptr));
  }    
  match = this->matchptr++;
  this->matchctr++;

  match->querypos = querypos;
  match->weight = 0.0;		/* Will be entered later */
  if ((match->forwardp = forwardp) == true) {
    match->position = diagonal + querypos - querylength;
  } else {
    match->position = diagonal - querypos;
  }
  match->fivep = fivep;
  match->npairings = 0;

  if (chromosome_iit == NULL) {
    match->chrnum = 0;
    match->chrpos = match->position;
  } else {
    index = Univ_IIT_get_one(chromosome_iit,match->position,match->position);
    match->chrpos = match->position - Univinterval_low(Univ_IIT_interval(chromosome_iit,index));
    match->chrnum = index;
  }

  debug(
	printf("Creating %d: %d %d %d %u\n",
	       this->matchctr-1,match->querypos,match->forwardp,match->fivep,match->position);
	);
  
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

  listcell->first = (void *) match;
  listcell->rest = list;

  return listcell;
}

List_T
Matchpool_push_existing (List_T list, T this, Match_T match) {
  List_T listcell;
  List_T p;
  int n;

  debug(
	printf("Pushing: %d %d %d %u\n",
	       match->querypos,match->forwardp,match->fivep,match->position);
	);
  
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

  listcell->first = (void *) match;
  listcell->rest = list;

  return listcell;
}


List_T
Matchpool_pop (List_T list, Match_T *x) {
  List_T head;

  if (list != NULL) {
    head = list->rest;
    *x = (Match_T) list->first;
    return head;
  } else {
    return list;
  }
}


List_T
Matchpool_transfer (List_T dest, List_T source) {
  List_T p, next;
#ifdef DEBUG
  Match_T match;
#endif

  for (p = source; p != NULL; p = next) {
    debug(
	  match = List_head(p);
	  printf("Transferring: %d %d %d %u\n",
		 match->querypos,match->forwardp,match->fivep,match->position);
	  );
    next = p->rest;
    p->rest = dest;
    dest = p;
  }
  debug(printf("Done with transfer\n"));
  
  return dest;
}

