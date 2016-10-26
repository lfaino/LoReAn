/* $Id: oligo.h 157221 2015-01-22 18:38:57Z twu $ */
#ifndef OLIGO_INCLUDED
#define OLIGO_INCLUDED

#include "bool.h"
#include "genomicpos.h"
#include "indexdb.h"
#include "reader.h"

typedef enum {INIT, DONE, INVALID, VALID} Oligostate_T;

extern void
Oligo_setup (int index1part);

extern char *
Oligo_one_nt (Storedoligomer_T oligo, int oligosize);

extern int
Oligo_lookup (Univcoord_T **positions, Indexdb_T indexdb, Storedoligomer_T storedoligo);

extern Oligostate_T
Oligo_next (Oligostate_T last_state, int *querypos, Storedoligomer_T *forward, 
	    Storedoligomer_T *revcomp, Reader_T reader, cDNAEnd_T cdnaend);
extern Oligostate_T
Oligo_skip (Oligostate_T last_state, int *querypos, Storedoligomer_T *forward,
	    Storedoligomer_T *revcomp, Reader_T reader, cDNAEnd_T cdnaend, int nskip);

extern bool
Oligo_repetitive_p (Storedoligomer_T oligo);

#endif
