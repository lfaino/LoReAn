/* $Id: cellpool.h 135447 2014-05-07 22:25:45Z twu $ */
#ifndef CELLPOOL_INCLUDED
#define CELLPOOL_INCLUDED

#include "bool.h"
#include "list.h"


typedef struct Cell_T *Cell_T;
struct Cell_T {
  int rootposition;
  int querypos;
  int hit;
  bool fwdp;
  int score;
};


#define T Cellpool_T
typedef struct T *T;

extern void
Cellpool_free (T *old);
extern void
Cellpool_free_memory (T this);
extern void
Cellpool_report_memory (T this);
extern T
Cellpool_new (void);
extern void
Cellpool_reset (T this);
extern List_T
Cellpool_push (List_T list, T this, int rootposition, int querypos, int hit, bool fwdp, int score);
extern List_T
Cellpool_pop (List_T list, Cell_T *x);

#undef T
#endif


