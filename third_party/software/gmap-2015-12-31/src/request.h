/* $Id: request.h 155282 2014-12-12 19:42:54Z twu $ */
#ifndef REQUEST_INCLUDED
#define REQUEST_INCLUDED

#ifdef GSNAP
#include "shortread.h"
#else
#include "sequence.h"
#endif

#define T Request_T
typedef struct T *T;

extern int
Request_id (T this);

#ifdef GSNAP

extern Shortread_T
Request_queryseq1 (T this);
extern Shortread_T
Request_queryseq2 (T this);
extern T
Request_new (int id, Shortread_T queryseq1, Shortread_T queryseq2);

#else

extern Sequence_T
Request_queryseq (T this);
extern T
Request_new (int id, Sequence_T queryseq);

#endif


extern void
Request_free (T *old);

#undef T
#endif
