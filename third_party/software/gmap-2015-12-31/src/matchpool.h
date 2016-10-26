/* $Id: matchpool.h 157221 2015-01-22 18:38:57Z twu $ */
#ifndef MATCHPOOL_INCLUDED
#define MATCHPOOL_INCLUDED

#include "bool.h"
#include "iit-read-univ.h"
#include "genomicpos.h"
#include "types.h"
#include "match.h"
#include "list.h"

#define T Matchpool_T
typedef struct T *T;

extern void
Matchpool_free (T *old);
extern void
Matchpool_free_memory (T this);
extern T
Matchpool_new (void);
extern void
Matchpool_reset (T this);
extern void
Matchpool_save (T this);
extern void
Matchpool_restore (T this);
extern List_T
Matchpool_push (List_T list, T this, int querypos, int querylength, bool forwardp, bool fivep,
		Univcoord_T diagonal, Univ_IIT_T chromosome_iit);
extern List_T
Matchpool_push_existing (List_T list, T this, Match_T match);
extern List_T
Matchpool_pop (List_T list, Match_T *x);
extern List_T
Matchpool_transfer (List_T dest, List_T source);

#undef T
#endif
