/* $Id: junction.h 166641 2015-05-29 21:13:04Z twu $ */
#ifndef JUNCTION_INCLUDED
#define JUNCTION_INCLUDED

typedef enum {NO_JUNCTION, INS_JUNCTION, DEL_JUNCTION, SPLICE_JUNCTION,
	      CHIMERA_JUNCTION, AMB_JUNCTION, END_JUNCTION} Junctiontype_T;

#include "types.h"
#include "genomicpos.h"
#include "bool.h"
#include "genome.h"
#include "list.h"


#define T Junction_T
typedef struct T *T;

extern void
Junction_print (T this);
extern void
Junction_free (T *old);
extern void
Junction_gc (List_T *list);

extern T
Junction_new_insertion (int nindels);
extern T
Junction_new_deletion (int nindels, Univcoord_T deletionpos);
extern T
Junction_new_splice (Chrpos_T splice_distance, int sensedir, double donor_prob, double acceptor_prob);

extern T
Junction_new_chimera (int sensedir, double donor_prob, double acceptor_prob);

extern T
Junction_copy (T old);


extern Junctiontype_T
Junction_type (T this);
extern int
Junction_sensedir (T this);
extern double
Junction_prob (T this);
extern double
Junction_donor_prob (T this);
extern double
Junction_acceptor_prob (T this);

extern int
Junction_nindels (T this);
extern int
Junction_adj (T this);
extern char *
Junction_deletion_string (T this, Genome_T genome, bool plusp);
extern Chrpos_T
Junction_splice_distance (T this);
extern void
Junction_set_unambiguous (T this, Chrpos_T distance, double donor_prob, double acceptor_prob);

#undef T
#endif

