/* $Id: spanningelt.h 146624 2014-09-02 21:32:50Z twu $ */
#ifndef SPANNINGELT_INCLUDED
#define SPANNINGELT_INCLUDED

typedef struct Spanningelt_T *Spanningelt_T;

#include "bool.h"
#include "genomicpos.h"
#include "types.h"
#include "indexdb_hr.h"
#include "indexdb.h"
#include "list.h"


#define T Spanningelt_T
struct T {
  bool partnerp;
  int querylength;
  int partner_querypos;		/* for debugging */
  int querypos;			/* for debugging */

  /* Intersectionr results are in native format, not littleendian */
  Univcoord_T *intersection_diagonals;
  int intersection_ndiagonals;

#ifdef LARGE_GENOMES
  unsigned char *partner_positions_high;
  UINT4 *partner_positions_low;
#else
  Univcoord_T *partner_positions;
#endif
  int partner_diagterm;
  int partner_npositions;

  Compoundpos_T compoundpos;
  int compoundpos_diagterm;

#ifdef LARGE_GENOMES
  unsigned char *positions_high_allocated;
  unsigned char *positions_high;
  UINT4 *positions_low_allocated;
  UINT4 *positions_low;
#else
  Univcoord_T *positions_allocated;
  Univcoord_T *positions;
#endif
  int diagterm;
  int npositions;
  
  int candidates_score;		/* score used for generating candidates */
  int pruning_score;		 /* score used for pruning */
  int miss_querypos5; /* If partnerp is true, this is the overlap of the two partners */
  int miss_querypos3;

  /* Reset values */
  Univcoord_T *intersection_diagonals_reset;
  int intersection_ndiagonals_reset;
#ifdef LARGE_GENOMES
  unsigned char *partner_positions_high_reset;
  UINT4 *partner_positions_low_reset;
  unsigned char *positions_high_reset;
  UINT4 *positions_low_reset;
#else
  Univcoord_T *partner_positions_reset;
  Univcoord_T *positions_reset;
#endif
  int partner_npositions_reset;
  int npositions_reset;
};

extern Width_T
Spanningelt_setup (Width_T index1part_in, Width_T index1interval_in);
extern void
Spanningelt_init_positions_free (bool positions_fileio_p);
extern void
Spanningelt_gc (T old);
extern T
Spanningelt_reset (T this);
extern void
Spanningelt_print (T this);
extern void
Spanningelt_print_array (Spanningelt_T *array, int nelts);
extern void
Spanningelt_print_set (List_T spanningset);

extern void
Spanningelt_set (int *minscore, int *nelts, Spanningelt_T *array,
		 Storedoligomer_T *stage1_oligos, bool **stage1_retrievedp,
#ifdef LARGE_GENOMES
		 unsigned char ***stage1_positions_high, UINT4 ***stage1_positions_low,
#else
		 Univcoord_T ***stage1_positions,
#endif
		 int **stage1_npositions, Indexdb_T indexdb, int query_lastpos, int querylength,
		 int mod, bool plusp);
extern int
Spanningelt_candidates_cmp (const void *a, const void *b);
extern int
Spanningelt_pruning_cmp (const void *a, const void *b);

extern Univcoord_T *
Spanningelt_diagonals (int *ndiagonals, T this, int *miss_querypos5, int *miss_querypos3);

#undef T
#endif

