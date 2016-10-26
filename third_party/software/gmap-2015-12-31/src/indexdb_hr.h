/* $Id: indexdb_hr.h 157221 2015-01-22 18:38:57Z twu $ */
#ifndef INDEXDB_HR_INCLUDED
#define INDEXDB_HR_INCLUDED

#include "bool.h"
#include "genomicpos.h"
#include "reader.h"
#include "indexdb.h"


typedef struct Compoundpos_T *Compoundpos_T;


#define T Indexdb_T

extern void
Indexdb_hr_setup (int index1part_in);
extern void
Compoundpos_init_positions_free (bool positions_fileio_p);
extern void
Compoundpos_set (Compoundpos_T compoundpos);
extern void
Compoundpos_reset (Compoundpos_T compoundpos);
extern void
Compoundpos_print_sizes (Compoundpos_T compoundpos);
extern void
Compoundpos_dump (Compoundpos_T compoundpos, int diagterm);
extern void
Compoundpos_free (Compoundpos_T *old);
extern void
Compoundpos_heap_init (Compoundpos_T compoundpos, int querylength, int diagterm);
extern bool
Compoundpos_find (bool *emptyp, Compoundpos_T compoundpos, Univcoord_T local_goal);
extern int
Compoundpos_search (Univcoord_T *value, Compoundpos_T compoundpos, Univcoord_T local_goal);
extern Compoundpos_T
Indexdb_compoundpos_left_subst_2 (T this, Storedoligomer_T oligo);
extern Compoundpos_T
Indexdb_compoundpos_left_subst_1 (T this, Storedoligomer_T oligo);
extern Compoundpos_T
Indexdb_compoundpos_right_subst_2 (T this, Storedoligomer_T oligo);
extern Compoundpos_T
Indexdb_compoundpos_right_subst_1 (T this, Storedoligomer_T oligo);

extern Univcoord_T *
Indexdb_merge_compoundpos (int *nmerged, Compoundpos_T compoundpos, int diagterm);

extern Univcoord_T *
Indexdb_read_left_subst_2 (int *npositions, T this, Storedoligomer_T oligo);
extern Univcoord_T *
Indexdb_read_left_subst_1 (int *npositions, T this, Storedoligomer_T oligo);
extern Univcoord_T *
Indexdb_read_right_subst_2 (int *npositions, T this, Storedoligomer_T oligo);
extern Univcoord_T *
Indexdb_read_right_subst_1 (int *npositions, T this, Storedoligomer_T oligo);

extern int
Indexdb_count_no_subst (T this, Storedoligomer_T oligo);
extern int
Indexdb_count_left_subst_2 (T this, Storedoligomer_T oligo);
extern int
Indexdb_count_left_subst_1 (T this, Storedoligomer_T oligo);
extern int
Indexdb_count_right_subst_2 (T this, Storedoligomer_T oligo);
extern int
Indexdb_count_right_subst_1 (T this, Storedoligomer_T oligo);

#undef T
#endif

