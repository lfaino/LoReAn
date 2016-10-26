static char rcsid[] = "$Id: pairpool.c 147823 2014-09-15 23:13:11Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "pairpool.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>		/* For memcpy */
#include "assert.h"
#include "mem.h"
#include "comp.h"
#include "pairdef.h"
#include "listdef.h"
#include "intron.h"
#include "complement.h"


#define CHUNKSIZE 20000
#define MICROINTRON_LENGTH 9	/* For Pairpool_add_genomeskip */

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

/* Getting genomic nt */
#ifdef DEBUG8
#define debug8(x) x
#else
#define debug8(x)
#endif

/* joining ends */
#ifdef DEBUG15
#define debug15(x) x
#else 
#define debug15(x)
#endif

/* clean_join */
#ifdef DEBUG16
#define debug16(x) x
#else 
#define debug16(x)
#endif


#define T Pairpool_T
struct T {
  int npairs;
  int pairctr;
  struct Pair_T *pairptr;
  List_T pairchunks;

  int nlistcells;
  int listcellctr;
  struct List_T *listcellptr;
  List_T listcellchunks;
};

void
Pairpool_free (T *old) {
  List_T p;
  struct Pair_T *pairptr;
  struct List_T *listcellptr;

  if (*old) {
    for (p = (*old)->pairchunks; p != NULL; p = List_next(p)) {
      pairptr = (struct Pair_T *) List_head(p);
      FREE(pairptr);
    }
    List_free(&(*old)->pairchunks);
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
Pairpool_free_memory (T this) {
  List_T p;
  struct Pair_T *pairptr;
  struct List_T *listcellptr;

  for (p = this->pairchunks; p != NULL; p = List_next(p)) {
    pairptr = (struct Pair_T *) List_head(p);
    FREE_KEEP(pairptr);
  }
  List_free_keep(&this->pairchunks);
  for (p = this->listcellchunks; p != NULL; p = List_next(p)) {
    listcellptr = (struct List_T *) List_head(p);
    FREE_KEEP(listcellptr);
  }
  List_free_keep(&this->listcellchunks);

  this->npairs = 0;
  this->pairctr = 0;
  this->pairchunks = NULL;
  /* this->pairptr = add_new_pairchunk(this); */

  this->nlistcells = 0;
  this->listcellctr = 0;
  this->listcellchunks = NULL;
  /* this->listcellptr = add_new_listcellchunk(this); */

  return;
}


void
Pairpool_report_memory (T this) {
  printf("Pairpool has %d pairchunks and %d listcellchunks\n",
	 List_length(this->pairchunks),List_length(this->listcellchunks));
  return;
}


static struct Pair_T *
add_new_pairchunk (T this) {
  struct Pair_T *chunk;

  chunk = (struct Pair_T *) MALLOC_KEEP(CHUNKSIZE*sizeof(struct Pair_T));
  this->pairchunks = List_push_keep(this->pairchunks,(void *) chunk);
  debug1(printf("Adding a new chunk of pairs.  Ptr for pair %d is %p\n",
		this->npairs,chunk));

  this->npairs += CHUNKSIZE;

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
Pairpool_new (void) {
  T new = (T) MALLOC(sizeof(*new));

  new->npairs = 0;
  new->pairctr = 0;
  new->pairchunks = NULL;
  /* new->pairptr = add_new_pairchunk(new); */

  new->nlistcells = 0;
  new->listcellctr = 0;
  new->listcellchunks = NULL;
  /* new->listcellptr = add_new_listcellchunk(new); */

  return new;
}

void
Pairpool_reset (T this) {
  this->pairctr = 0;
  this->listcellctr = 0;
  return;
}

/* gapp should be false for the following comps: MATCH_COMP,
   DYNPROG_MATCH_COMP, AMBIGUOUS_COMP, MISMATCH_COMP, INDEL_COMP,
   SHORTGAP_COMP */

List_T
Pairpool_push (List_T list, T this, int querypos, int genomepos, char cdna, char comp,
	       char genome, char genomealt, int dynprogindex) {
  List_T listcell;
  Pair_T pair;
  List_T p;
  int n;

  assert(querypos >= 0);

  if (this->pairctr >= this->npairs) {
    this->pairptr = add_new_pairchunk(this);
  } else if ((this->pairctr % CHUNKSIZE) == 0) {
    for (n = this->npairs - CHUNKSIZE, p = this->pairchunks;
	 n > this->pairctr; p = p->rest, n -= CHUNKSIZE) ;
    this->pairptr = (struct Pair_T *) p->first;
    debug1(printf("Located pair %d at %p\n",this->pairctr,this->pairptr));
  }    
  pair = this->pairptr++;
  this->pairctr++;

  pair->querypos = querypos;
  pair->genomepos = genomepos;
  pair->aapos = 0;
  pair->aaphase_g = -1;
  pair->aaphase_e = -1;
  pair->cdna = cdna;
  pair->comp = comp;
  pair->genome = genome;
  pair->genomealt = genomealt;
  pair->dynprogindex = dynprogindex;

  pair->aa_g = ' ';
  pair->aa_e = ' ';
  pair->shortexonp = false;
  pair->gapp = false;
  pair->knowngapp = false;
  pair->introntype = NONINTRON;
  if (comp == EXTRAEXON_COMP) {
    pair->extraexonp = true;
  } else {
    pair->extraexonp = false;
  }
  
  pair->queryjump = 0;
  pair->genomejump = 0;

  pair->state = GOOD;
  pair->protectedp = false;
  pair->disallowedp = false;
  pair->donor_prob = 0.0;
  pair->acceptor_prob = 0.0;
  pair->end_intron_p = false;

  debug(
	printf("Creating %p: %d %d %c %c %c\n",
	       pair,pair->querypos,pair->genomepos,pair->cdna,pair->comp,pair->genome);
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

  listcell->first = (void *) pair;
  listcell->rest = list;

  return listcell;
}


List_T
Pairpool_push_copy (List_T list, T this, Pair_T orig) {
  List_T listcell;
  Pair_T pair;
  List_T p;
  int n;

  if (this->pairctr >= this->npairs) {
    this->pairptr = add_new_pairchunk(this);
  } else if ((this->pairctr % CHUNKSIZE) == 0) {
    for (n = this->npairs - CHUNKSIZE, p = this->pairchunks;
	 n > this->pairctr; p = p->rest, n -= CHUNKSIZE) ;
    this->pairptr = (struct Pair_T *) p->first;
    debug1(printf("Located pair %d at %p\n",this->pairctr,this->pairptr));
  }    
  pair = this->pairptr++;
  this->pairctr++;

  memcpy(pair,orig,sizeof(struct Pair_T));

  debug(
	printf("Copying %p <= %p: %d %d %c %c %c\n",
	       pair,orig,pair->querypos,pair->genomepos,pair->cdna,pair->comp,pair->genome);
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

  listcell->first = (void *) pair;
  listcell->rest = list;

  return listcell;
}


List_T
Pairpool_push_gapalign (List_T list, T this, int querypos, int genomepos, char cdna, char comp,
			int introntype, char genome, char genomealt, bool extraexonp) {
  List_T listcell;
  Pair_T pair;
  List_T p;
  int n;

  if (this->pairctr >= this->npairs) {
    this->pairptr = add_new_pairchunk(this);
  } else if ((this->pairctr % CHUNKSIZE) == 0) {
    for (n = this->npairs - CHUNKSIZE, p = this->pairchunks;
	 n > this->pairctr; p = p->rest, n -= CHUNKSIZE) ;
    this->pairptr = (struct Pair_T *) p->first;
    debug1(printf("Located pair %d at %p\n",this->pairctr,this->pairptr));
  }    
  pair = this->pairptr++;
  this->pairctr++;

  pair->querypos = querypos;
  pair->genomepos = genomepos;
  pair->aapos = 0;
  pair->aaphase_g = -1;
  pair->aaphase_e = -1;
  pair->cdna = cdna;
  pair->comp = comp;
  pair->genome = genome;
  pair->genomealt = genomealt;
  pair->dynprogindex = 0;

  pair->aa_g = ' ';
  pair->aa_e = ' ';
  pair->shortexonp = false;
  pair->gapp = true;
  pair->knowngapp = false;
  pair->introntype = introntype;
  pair->extraexonp = extraexonp;
  
  pair->queryjump = 0;
  pair->genomejump = 0;

  pair->state = GOOD;
  pair->protectedp = false;
  pair->disallowedp = false;
  pair->donor_prob = 0.0;
  pair->acceptor_prob = 0.0;
  pair->end_intron_p = false;

  debug(
	printf("Creating %p: %d %d %c %c %c introntype %d\n",
	       pair,pair->querypos,pair->genomepos,pair->cdna,pair->comp,pair->genome,pair->introntype);
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

  listcell->first = (void *) pair;
  listcell->rest = list;

  return listcell;
}

List_T
Pairpool_push_gapholder (List_T list, T this, int queryjump, int genomejump,
			 Pair_T leftpair, Pair_T rightpair, bool knownp) {
  List_T listcell;
  Pair_T pair;
  List_T p;
  int n;

  if (this->pairctr >= this->npairs) {
    this->pairptr = add_new_pairchunk(this);
  } else if ((this->pairctr % CHUNKSIZE) == 0) {
    for (n = this->npairs - CHUNKSIZE, p = this->pairchunks;
	 n > this->pairctr; p = p->rest, n -= CHUNKSIZE) ;
    this->pairptr = (struct Pair_T *) p->first;
    debug1(printf("Located pair %d at %p\n",this->pairctr,this->pairptr));
  }    
  pair = this->pairptr++;
  this->pairctr++;

  pair->querypos = -1;
  pair->genomepos = -1;

  pair->aapos = 0;
  pair->aaphase_g = -1;
  pair->aaphase_e = -1;
  pair->cdna = ' ';
  pair->comp = ' ';
  pair->genome = ' ';
  pair->genomealt = ' ';
  pair->dynprogindex = 0;

  pair->aa_g = ' ';
  pair->aa_e = ' ';
  pair->shortexonp = false;
  pair->gapp = true;
  if (knownp == true) {
    pair->knowngapp = true;
    pair->donor_prob = 2.0;
    pair->acceptor_prob = 2.0;
  } else {
    pair->knowngapp = false;
    pair->donor_prob = 0.0;
    pair->acceptor_prob = 0.0;
  }
  pair->introntype = NONINTRON;
  pair->extraexonp = false;

  if (leftpair && rightpair) {
    queryjump = rightpair->querypos - leftpair->querypos - 1;
    genomejump = rightpair->genomepos - leftpair->genomepos - 1;
    if (leftpair->cdna == ' ') queryjump++;
    if (leftpair->genome == ' ') genomejump++;
  }

  pair->queryjump = queryjump;
  pair->genomejump = genomejump;

  pair->state = GOOD;
  pair->protectedp = false;
  pair->disallowedp = false;
  pair->end_intron_p = false;

  debug(printf("Creating gap %p, queryjump=%d, genomejump=%d\n",pair,queryjump,genomejump));

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

  listcell->first = (void *) pair;
  listcell->rest = list;

  return listcell;
}

List_T
Pairpool_push_existing (List_T list, T this, Pair_T pair) {
  List_T listcell;
  List_T p;
  int n;

  debug(
	Pair_T head;
	if (pair->gapp == true) {
	  printf("Pushing gap %p: queryjump=%d, genomejump=%d onto ",
		 pair,pair->queryjump,pair->genomejump);
	} else {
	  printf("Pushing %p: %d %d %c %c %c onto ",
		 pair,pair->querypos,pair->genomepos,pair->cdna,pair->comp,pair->genome);
	}
	if (list == NULL) {
	  printf("NULL\n");
	} else {
	  head = list->first;
	  if (head->gapp == true) {
	    printf("gap %p: queryjump=%d, genomejump=%d\n",
		   head,head->queryjump,head->genomejump);
	  } else {
	    printf("%p: %d %d %c %c %c\n",
		   head,head->querypos,head->genomepos,head->cdna,head->comp,head->genome);
	  }
	}
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

  listcell->first = (void *) pair;
  listcell->rest = list;

  return listcell;
}


/* Note: this does not free the list cell */
List_T
Pairpool_pop (List_T list, Pair_T *x) {
  List_T head;

  if (list != NULL) {
    head = list->rest;
    *x = (Pair_T) list->first;
    debug2(
	   if ((*x)->gapp == true) {
	     printf("Popping gap: queryjump=%d, genomejump=%d\n",
		    (*x)->queryjump,(*x)->genomejump);
	   } else {
	     printf("Popping: %d %d %c %c %c\n",
		    (*x)->querypos,(*x)->genomepos,(*x)->cdna,(*x)->comp,(*x)->genome);
	   }
	   );
    return head;
  } else {
    return list;
  }
}


List_T
Pairpool_transfer (List_T dest, List_T source) {
  List_T p, next;
#ifdef DEBUG
  Pair_T pair;
#endif

  for (p = source; p != NULL; p = next) {
    debug(
	  pair = List_head(p);
	  if (pair->cdna == '\0' || pair->genome == '\0') {
	    abort();
	  }
	  if (pair->gapp) {
	    printf("Transferring gap %p: queryjump=%d, genomejump=%d\n",
		   pair,pair->queryjump,pair->genomejump);
	  } else {
	    printf("Transferring %p: %d %d %c %c %c\n",
		   pair,pair->querypos,pair->genomepos,pair->cdna,pair->comp,pair->genome);
	  }
	  );
    next = p->rest;
    p->rest = dest;
    dest = p;
  }
  return dest;
}

List_T
Pairpool_transfer_n (List_T dest, List_T source, int n) {
  List_T p, next;
#ifdef DEBUG
  Pair_T pair;
#endif

  for (p = source; p != NULL && --n >= 0; p = next) {
    debug(
	  pair = List_head(p);
	  if (pair->cdna == '\0' || pair->genome == '\0') {
	    abort();
	  }
	  if (pair->gapp) {
	    printf("Transferring gap %p: queryjump=%d, genomejump=%d\n",
		   pair,pair->queryjump,pair->genomejump);
	  } else {
	    printf("Transferring %p: %d %d %c %c %c\n",
		   pair,pair->querypos,pair->genomepos,pair->cdna,pair->comp,pair->genome);
	  }
	  );
    next = p->rest;
    p->rest = dest;
    dest = p;
  }
  return dest;
}

int
Pairpool_count_bounded (int *nstart, List_T source, int minpos, int maxpos) {
  int npairs = 0;
  List_T p, next;
  Pair_T pair;

  *nstart = 0;
  for (p = source; p != NULL; p = next) {
    pair = List_head(p);
    next = p->rest;
    if (pair->querypos < minpos) {
      *nstart += 1;
    } else if (pair->querypos < maxpos) {
      npairs++;
    } else {
      p = NULL;			/* Terminate transfer */
    }
  }
  return npairs;
}


#if 0
List_T
Pairpool_transfer_bounded (List_T dest, List_T source, int minpos, int maxpos) {
  List_T p, next;
  Pair_T pair;

  for (p = source; p != NULL; p = next) {
    debug(
	  pair = (Pair_T) List_head(p);
	  if (pair->cdna == '\0' || pair->genome == '\0') {
	    abort();
	  }
	  printf("Transferring %p: %d %d %c %c %c\n",
		 pair,pair->querypos,pair->genomepos,pair->cdna,pair->comp,pair->genome);
	  );
    pair = (Pair_T) List_head(p);
    next = p->rest;
    if (pair->querypos == minpos) {
      if (dest != NULL) {
	/* Pop last querypos off the stack, because we want only one of them */
	dest = dest->rest;
      }
      p->rest = dest;
      dest = p;
    } else if (pair->querypos == maxpos) {
      p->rest = dest;
      dest = p;
      p = NULL;			/* Terminate transfer */
    } else if (pair->querypos > minpos && pair->querypos < maxpos) {
      p->rest = dest;
      dest = p;
    }
  }

  return dest;
}
#endif


/* Originally prohibited copying of gaps */
List_T
Pairpool_copy (List_T source, T this) {
  List_T dest = NULL;

  while (source != NULL) {
    dest = Pairpool_push_copy(dest,this,/*orig*/source->first);
    source = source->rest;
  }
  return List_reverse(dest);
}


struct Pair_T *
Pairpool_copy_array (struct Pair_T *source, int npairs) {
  struct Pair_T *dest;

  dest = (struct Pair_T *) MALLOC_OUT(npairs * sizeof(struct Pair_T));
  memcpy(dest,source,npairs*sizeof(struct Pair_T));
  return dest;
}


void
Pairpool_clean_join (List_T *left_path, List_T *right_pairs) {
  Pair_T leftpair, rightpair;
  int queryjump, genomejump;


  debug16(printf("Entered clean_join\n"));
  debug16(printf("left path:\n"));
  debug16(Pair_dump_list(*left_path,true));
  debug16(printf("right pairs:\n"));
  debug16(Pair_dump_list(*right_pairs,true));


  while (*left_path != NULL && ((Pair_T) (*left_path)->first)->gapp == true) {
    debug16(printf("Clearing gap on left\n"));
    *left_path = Pairpool_pop(*left_path,&leftpair);
  }
  while (*right_pairs != NULL && ((Pair_T) (*right_pairs)->first)->gapp == true) {
    debug16(printf("Clearing gap on right\n"));
    *right_pairs = Pairpool_pop(*right_pairs,&rightpair);
  }
  
  if (*left_path != NULL && *right_pairs != NULL) {
    leftpair = (Pair_T) (*left_path)->first;
    rightpair = (Pair_T) (*right_pairs)->first;
    queryjump = rightpair->querypos - leftpair->querypos - 1;
    genomejump = rightpair->genomepos - leftpair->genomepos - 1;
    debug16(printf("queryjump %d, genomejump %d\n",queryjump,genomejump));

    /* Fix overlap */
    while (*left_path != NULL && *right_pairs != NULL && (queryjump < 0 || genomejump < 0)) {
      while (*left_path != NULL && ((Pair_T) (*left_path)->first)->gapp == true) {
	debug16(printf("Clearing gap on left\n"));
	*left_path = Pairpool_pop(*left_path,&leftpair);
      }
      while (*right_pairs != NULL && ((Pair_T) (*right_pairs)->first)->gapp == true) {
	debug16(printf("Clearing gap on right\n"));
	*right_pairs = Pairpool_pop(*right_pairs,&rightpair);
      }
      *left_path = Pairpool_pop(*left_path,&leftpair);
      *right_pairs = Pairpool_pop(*right_pairs,&rightpair);
      queryjump = rightpair->querypos - leftpair->querypos - 1;
      genomejump = rightpair->genomepos - leftpair->genomepos - 1;
      debug16(printf("Revising queryjump to be %d = %d - %d - 1\n",queryjump,rightpair->querypos,leftpair->querypos));
    }
  
    while (*left_path != NULL && ((Pair_T) (*left_path)->first)->gapp == true) {
      debug16(printf("Clearing gap on left\n"));
      *left_path = Pairpool_pop(*left_path,&leftpair);
    }
    while (*right_pairs != NULL && ((Pair_T) (*right_pairs)->first)->gapp == true) {
      debug16(printf("Clearing gap on right\n"));
      *right_pairs = Pairpool_pop(*right_pairs,&rightpair);
    }
  }

  return;
}


List_T
Pairpool_remove_gapholders (List_T pairs) {
  List_T path = NULL;
  Pair_T pair;

  while (pairs != NULL) {
    /* pairptr = pairs; */
    /* pairs = Pairpool_pop(pairs,&pair); */
    pair = (Pair_T) pairs->first;
    if (pair->gapp == true) {
      debug(printf("Removing a gap with queryjump = %d, genomejump = %d\n",
		   pair->queryjump,pair->genomejump));
      pairs = Pairpool_pop(pairs,&pair);
    } else {
#ifdef WASTE
      path = Pairpool_push_existing(path,pairpool,pair);
#else
      path = List_transfer_one(path,&pairs);
#endif
    }
  }

  return List_reverse(path);
}


#ifdef CHECK_ASSERTIONS
static List_T
check_for_gapholders (List_T pairs) {
  List_T path = NULL;
  Pair_T pair;

  while (pairs != NULL) {
    /* pairptr = pairs; */
    /* pairs = Pairpool_pop(pairs,&pair); */
    pair = (Pair_T) pairs->first;
    if (pair->gapp == true) {
      abort();
    } else {
#ifdef WASTE
      path = Pairpool_push_existing(path,pairpool,pair);
#else
      path = List_transfer_one(path,&pairs);
#endif
    }
  }

  return List_reverse(path);
}
#endif


List_T
Pairpool_join_end3 (List_T path_orig, List_T end3_pairs_orig, Pairpool_T pairpool,
		    bool copy_end_p) {
  List_T path, end3_pairs;
  Pair_T pair, leftpair;
  int queryjump = -1, genomejump = -1;
  
#ifdef CHECK_ASSERTIONS
  path_orig = check_for_gapholders(path_orig);
  end3_pairs_orig = check_for_gapholders(end3_pairs_orig);
#endif

  path = Pairpool_copy(path_orig,pairpool);
  if (copy_end_p == true) {
    end3_pairs = Pairpool_copy(end3_pairs_orig,pairpool);
  } else {
    end3_pairs = end3_pairs_orig;
  }

  if (path == NULL && end3_pairs == NULL) {
    return (List_T) NULL;
  } else if (path == NULL) {
    return List_reverse(end3_pairs);
  } else if (end3_pairs == NULL) {
    return path;
  }

  debug15(printf("Entered join_end3\n"));
  debug15(printf("path:\n"));
  debug15(Pair_dump_list(path,true));
  debug15(printf("end3_pairs:\n"));
  debug15(Pair_dump_list(end3_pairs,true));


  leftpair = (Pair_T) path->first;
  pair = (Pair_T) end3_pairs->first;
  queryjump = pair->querypos - leftpair->querypos - 1;
  genomejump = pair->genomepos - leftpair->genomepos - 1;
  debug15(printf("queryjump %d, genomejump %d\n",queryjump,genomejump));

  if (queryjump == 0 && genomejump == 0) {
    /* Do nothing, although this is unexpected */
  } else if (queryjump >= 0 && genomejump >= 0) {
    /* Insert a gapholder */
    path = Pairpool_push_gapholder(path,pairpool,queryjump,genomejump,
				   /*leftpair*/NULL,/*rightpair*/NULL,/*knownp*/false);
  } else {
    /* Fix overlap */
    while (path != NULL && end3_pairs != NULL && (queryjump < 0 || genomejump < 0)) {
      pair = (Pair_T) end3_pairs->first;

      if (path != NULL) {
	path = Pairpool_pop(path,&leftpair);
      }
      if (end3_pairs != NULL) {
	end3_pairs = Pairpool_pop(end3_pairs,&pair);
      }
      queryjump = pair->querypos - leftpair->querypos - 1;
      genomejump = pair->genomepos - leftpair->genomepos - 1;
      debug15(printf("Revising queryjump to be %d = %d - %d - 1\n",queryjump,pair->querypos,leftpair->querypos));
    }

    path = Pairpool_push_existing(path,pairpool,leftpair);
    if (queryjump == 0 && genomejump == 0) {
      /* No gapholder needed */
    } else {
      path = Pairpool_push_gapholder(path,pairpool,queryjump,genomejump,
				     /*leftpair*/NULL,/*rightpair*/NULL,/*knownp*/false);
    }
    path = Pairpool_push_existing(path,pairpool,pair);
  }
  
  while (end3_pairs != NULL) {
    path = List_transfer_one(path,&end3_pairs);
  }
    
  debug15(printf("joined path:\n"));
  debug15(Pair_dump_list(path,true));
  debug15(printf("\n"));

  return path;
}


List_T
Pairpool_join_end5 (List_T pairs_orig, List_T end5_path_orig, Pairpool_T pairpool,
		    bool copy_end_p) {
  List_T pairs, end5_path;
  Pair_T pair, rightpair;
  int queryjump = -1, genomejump = -1;
  
#ifdef CHECK_ASSERTIONS
  pairs_orig = check_for_gapholders(pairs_orig);
  end5_path_orig = check_for_gapholders(end5_path_orig);
#endif

  pairs = Pairpool_copy(pairs_orig,pairpool);
  if (copy_end_p == true) {
    end5_path = Pairpool_copy(end5_path_orig,pairpool);
  } else {
    end5_path = end5_path_orig;
  }

  if (pairs == NULL && end5_path == NULL) {
    return (List_T) NULL;
  } else if (pairs == NULL) {
    return List_reverse(end5_path);
  } else if (end5_path == NULL) {
    return pairs;
  }

  debug15(printf("Entered join_end5\n"));
  debug15(printf("pairs:\n"));
  debug15(Pair_dump_list(pairs,true));
  debug15(printf("end5_path:\n"));
  debug15(Pair_dump_list(end5_path,true));


  rightpair = (Pair_T) pairs->first;
  pair = (Pair_T) end5_path->first;
  queryjump = rightpair->querypos - pair->querypos - 1;
  genomejump = rightpair->genomepos - pair->genomepos - 1;
  debug15(printf("queryjump %d, genomejump %d\n",queryjump,genomejump));

  if (queryjump == 0 && genomejump == 0) {
    /* Do nothing, although this is unexpected */
  } else if (queryjump >= 0 && genomejump >= 0) {
    /* Insert a gapholder */
    pairs = Pairpool_push_gapholder(pairs,pairpool,queryjump,genomejump,
				    /*leftpair*/NULL,/*rightpair*/NULL,/*knownp*/false);
  } else {
    /* Fix overlap */
    while (pairs != NULL && end5_path != NULL && (queryjump < 0 || genomejump < 0)) {
      pair = (Pair_T) end5_path->first;

      if (pairs != NULL) {
	pairs = Pairpool_pop(pairs,&rightpair);
      }
      if (end5_path != NULL) {
	end5_path = Pairpool_pop(end5_path,&pair);
      }
      queryjump = rightpair->querypos - pair->querypos - 1;
      genomejump = rightpair->genomepos - pair->genomepos - 1;
      debug15(printf("Revising queryjump to be %d = %d - %d - 1\n",queryjump,pair->querypos,rightpair->querypos));
    }

    pairs = Pairpool_push_existing(pairs,pairpool,rightpair);
    if (queryjump == 0 && genomejump == 0) {
      /* No gapholder needed */
    } else {
      pairs = Pairpool_push_gapholder(pairs,pairpool,queryjump,genomejump,
				      /*leftpair*/NULL,/*rightpair*/NULL,/*knownp*/false);
    }
    pairs = Pairpool_push_existing(pairs,pairpool,pair);
  }
  
  while (end5_path != NULL) {
    pairs = List_transfer_one(pairs,&end5_path);
  }
    
  debug15(printf("joined pairs:\n"));
  debug15(Pair_dump_list(pairs,true));
  debug15(printf("\n"));

  return pairs;
}


List_T
Pairpool_add_queryskip (List_T pairs, int r, int c, int dist, char *querysequence,
			int queryoffset, int genomeoffset, Pairpool_T pairpool, bool revp, int dynprogindex) {
  int j;
  char c1;
  int querycoord, genomecoord, step;

  querycoord = r-1;
  genomecoord = c-1;

  if (revp == true) {
    querycoord = -querycoord;
    genomecoord = -genomecoord;
    step = +1;
  } else {
    /* Advance to next genomepos */
    /* genomecoord++;  -- No, this leads to bugs */
    step = -1;
  }
  debug(printf("Entered Pairpool_add_genomeskip with r %d, c %d, revp %d => querycoord %d\n",r,c,revp,querycoord));

  for (j = 0; j < dist; j++) {
    c1 = querysequence[querycoord];
    debug(printf("Pushing %d,%d [%d,%d] (%c,-), ",r,c,queryoffset+querycoord,genomeoffset+genomecoord,c1));
    pairs = Pairpool_push(pairs,pairpool,queryoffset+querycoord,genomeoffset+genomecoord,
			  c1,INDEL_COMP,/*genome*/' ',/*genomealt*/' ',dynprogindex);
    debug(r--);
    querycoord += step;
  }

  return pairs;
}


static char complCode[128] = COMPLEMENT_LC;

static char
get_genomic_nt (char *g_alt, int genomicpos, Univcoord_T chroffset, Univcoord_T chrhigh,
		bool watsonp) {
  char c2, c2_alt;
  Univcoord_T pos;

#if 0
  /* If the read has a deletion, then we will extend beyond 0 or genomiclength, so do not restrict. */
  if (genomicpos < 0) {
    return '*';

  } else if (genomicpos >= genomiclength) {
    return '*';

  }
#endif

  if (watsonp) {
    if ((pos = chroffset + genomicpos) < chroffset) { /* Must be <, and not <=, or dynamic programming will fail */
      *g_alt = '*';
      return '*';

    } else if (pos >= chrhigh) {
      *g_alt = '*';
      return '*';

#if 0
    } else if (genome) {
      /* Not necessary, because Genome_get_char_blocks should work */
      debug8(printf("At %u, genomicnt is %c\n",
		    genomicpos,Genome_get_char(genome,pos)));
      return Genome_get_char(genome,pos);
#endif

    } else {
      /* GMAP with user-supplied genomic segment */
      debug8(printf("At %u, genomicnt is %c\n",
		    genomicpos,Genome_get_char_blocks(&(*g_alt),pos)));
      return Genome_get_char_blocks(&(*g_alt),pos);
    }

  } else {
    if ((pos = chrhigh - genomicpos) < chroffset) { /* Must be <, and not <=, or dynamic programming will fail */
      *g_alt = '*';
      return '*';

    } else if (pos >= chrhigh) {
      *g_alt = '*';
      return '*';

#if 0
    } else if (genome) {
      /* Not necessary, because Genome_get_char_blocks should work */
      c2 = Genome_get_char(genome,pos);
#endif

    } else {
      /* GMAP with user-supplied genomic segment */
      c2 = Genome_get_char_blocks(&c2_alt,pos);
    }
    debug8(printf("At %u, genomicnt is %c\n",genomicpos,complCode[(int) c2]));
    *g_alt = complCode[(int) c2_alt];
    return complCode[(int) c2];
  }
}


List_T
Pairpool_add_genomeskip (bool *add_dashes_p, List_T pairs, int r, int c, int dist,
			 char *genomesequence, char *genomesequenceuc,
			 int queryoffset, int genomeoffset, Pairpool_T pairpool, bool revp,
			 Univcoord_T chroffset, Univcoord_T chrhigh,
			 int cdna_direction, bool watsonp, int dynprogindex, bool use_genomicseg_p) {
  int j;
  char left1, left2, right2, right1, left1_alt, left2_alt, right2_alt, right1_alt, c2, c2_alt;
  int introntype;
  int querycoord, leftgenomecoord, rightgenomecoord, genomecoord, temp, step;

  querycoord = r-1;
  leftgenomecoord = c-dist;
  rightgenomecoord = c-1;

  if (revp == true) {
    querycoord = -querycoord;
    temp = leftgenomecoord;
    leftgenomecoord = -rightgenomecoord;
    rightgenomecoord = -temp;
    step = +1;
  } else {
    /* Advance to next querypos */
    /* querycoord++; -- No, this leads to bugs */
    step = -1;
  }
  debug(printf("Entered Pairpool_add_genomeskip with r %d, c %d, revp %d => querycoord %d\n",r,c,revp,querycoord));

  if (dist < MICROINTRON_LENGTH) {
    *add_dashes_p = true;
  } else {
    /* Check for intron */
    if (use_genomicseg_p) {
      left1 = left1_alt = genomesequenceuc[leftgenomecoord];
      left2 = left2_alt = genomesequenceuc[leftgenomecoord+1];
      right2 = right2_alt = genomesequenceuc[rightgenomecoord-1];
      right1 = right1_alt = genomesequenceuc[rightgenomecoord];
    } else {
      left1 = get_genomic_nt(&left1_alt,genomeoffset+leftgenomecoord,chroffset,chrhigh,watsonp);
      left2 = get_genomic_nt(&left2_alt,genomeoffset+leftgenomecoord+1,chroffset,chrhigh,watsonp);
      right2 = get_genomic_nt(&right2_alt,genomeoffset+rightgenomecoord-1,chroffset,chrhigh,watsonp);
      right1 = get_genomic_nt(&right1_alt,genomeoffset+rightgenomecoord,chroffset,chrhigh,watsonp);
    }
#ifdef EXTRACT_GENOMICSEG
    assert(left1 == genomesequenceuc[leftgenomecoord]);
    assert(left2 == genomesequenceuc[leftgenomecoord+1]);
    assert(right2 == genomesequenceuc[rightgenomecoord-1]);
    assert(right1 == genomesequenceuc[rightgenomecoord]);
#endif
	
#ifdef PMAP
    introntype = Intron_type(left1,left2,right2,right1,
			     left1_alt,left2_alt,right2_alt,right1_alt,
			     /*cdna_direction*/+1);
#else
    introntype = Intron_type(left1,left2,right2,right1,
			     left1_alt,left2_alt,right2_alt,right1_alt,
			     cdna_direction);
#endif

#if 0
    if (introntype == NONINTRON) {
      *add_dashes_p = true;
    } else {
      *add_dashes_p = false;
    }
#endif
    *add_dashes_p = false;
  }

  if (*add_dashes_p == true) {
    if (revp == true) {
      genomecoord = leftgenomecoord;
    } else {
      genomecoord = rightgenomecoord;
    }
    for (j = 0; j < dist; j++) {
      if (use_genomicseg_p) {
	c2 = c2_alt = genomesequence[genomecoord];
      } else {
	c2 = get_genomic_nt(&c2_alt,genomeoffset+genomecoord,chroffset,chrhigh,watsonp);
      }
#ifdef EXTRACT_GENOMICSEG
      assert(c2 == genomesequence[genomecoord]);
#endif

      debug(printf("Pushing %d,%d [%d,%d] (-,%c),",r,c,queryoffset+querycoord,genomeoffset+genomecoord,c2));
      pairs = Pairpool_push(pairs,pairpool,queryoffset+querycoord,genomeoffset+genomecoord,
			    ' ',INDEL_COMP,c2,c2_alt,dynprogindex);
      debug(c--);
      genomecoord += step;
    }
  } else {
    debug(printf("Large gap %c%c..%c%c.  Adding gap of type %d.\n",
		 left1,left2,right2,right1,introntype));
#ifndef NOGAPHOLDER
    pairs = Pairpool_push_gapholder(pairs,pairpool,/*queryjump*/0,/*genomejump*/dist,
				    /*leftpair*/NULL,/*rightpair*/NULL,/*knownp*/false);
#endif
  }
  
  return pairs;
}


/* Used for compacting pairs from working pairpool to final pairpool
   (dest) at end of stage2 and stage3, so the working pairpool can be
   re-used for the next gregion */
List_T
Pairpool_compact_copy (List_T list, T dest) {
  List_T newlist = NULL, listcell;
  Pair_T pair, orig;
  List_T p, q;
  int n;

  for (q = list; q != NULL; q = q->rest) {
    orig = (Pair_T) q->first;
    
    /* Following code taken from Pairpool_push_copy */
    if (dest->pairctr >= dest->npairs) {
      dest->pairptr = add_new_pairchunk(dest);
    } else if ((dest->pairctr % CHUNKSIZE) == 0) {
      for (n = dest->npairs - CHUNKSIZE, p = dest->pairchunks;
	   n > dest->pairctr; p = p->rest, n -= CHUNKSIZE) ;
      dest->pairptr = (struct Pair_T *) p->first;
      debug1(printf("Located pair %d at %p\n",dest->pairctr,dest->pairptr));
    }    
    pair = dest->pairptr++;
    dest->pairctr++;

    memcpy(pair,orig,sizeof(struct Pair_T));

    debug(
	  printf("Copying %p: %d %d %c %c %c\n",
		 pair,pair->querypos,pair->genomepos,pair->cdna,pair->comp,pair->genome);
	  );

    if (dest->listcellctr >= dest->nlistcells) {
      dest->listcellptr = add_new_listcellchunk(dest);
    } else if ((dest->listcellctr % CHUNKSIZE) == 0) {
      for (n = dest->nlistcells - CHUNKSIZE, p = dest->listcellchunks;
	   n > dest->listcellctr; p = p->rest, n -= CHUNKSIZE) ;
      dest->listcellptr = (struct List_T *) p->first;
      debug1(printf("Located listcell %d at %p\n",dest->listcellctr,dest->listcellptr));
    }
    listcell = dest->listcellptr++;
    dest->listcellctr++;

    listcell->first = (void *) pair;
    listcell->rest = newlist;
    newlist = listcell;
  }

  return List_reverse(newlist);
}

