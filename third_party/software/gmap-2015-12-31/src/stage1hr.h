/* $Id: stage1hr.h 173896 2015-09-12 00:11:40Z twu $ */
#ifndef STAGE1HR_INCLUDED
#define STAGE1HR_INCLUDED


#include "bool.h"
#include "types.h"
#include "mode.h"
#include "genomicpos.h"
#include "indexdb.h"
#include "shortread.h"
#include "iit-read-univ.h"
#include "genome.h"
#include "splicetrie.h"
#include "resulthr.h"		/* For Pairtype_T */
#include "stage2.h"
#include "stage3hr.h"

#ifdef PMAP
#include "oligoindex_pmap.h"
#else
#include "oligoindex_hr.h"
#endif
#include "pairpool.h"
#include "diagpool.h"
#include "cellpool.h"
#include "dynprog.h"


#if 0
/* Obsolete */
typedef enum {MASK_NONE, MASK_FREQUENT, MASK_REPETITIVE, MASK_GREEDY_FREQUENT, MASK_GREEDY_REPETITIVE} Masktype_T;
#endif

#define GMAP_IMPROVEMENT 1
#define GMAP_TERMINAL 2
#define GMAP_INDEL_KNOWNSPLICE 4
#define GMAP_PAIRSEARCH 8


typedef struct Floors_T *Floors_T;

extern void
Floors_free (Floors_T *old);
extern void
Floors_free_keep (Floors_T *old);


#define T Stage1_T
typedef struct T *T;

extern void
Stage1_init_positions_free (bool positions_fileio_p);
extern void
Stage1_free (T *old, int querylength);


extern Stage3end_T *
Stage1_single_read (int *npaths, int *first_absmq, int *second_absmq,
		    Shortread_T queryseq, Indexdb_T indexdb_fwd, Indexdb_T indexdb_rev,
		    int indexdb_size_threshold, Floors_T *floors_array,
		    double user_maxlevel_float, double user_mincoverage_float,
		    int indel_penalty_middle, int indel_penalty_end,
		    bool allow_end_indels_p, int max_end_insertions, int max_end_deletions, int min_indel_end_matches,
		    int localsplicing_penalty, int distantsplicing_penalty, int min_shortend,
		    Oligoindex_array_T oligoindices_major, Oligoindex_array_T oligoindices_minor,
		    Pairpool_T pairpool, Diagpool_T diagpool, Cellpool_T cellpool,
		    Dynprog_T dynprogL, Dynprog_T dynprogM, Dynprog_T dynprogR,
		    bool keep_floors_p);

extern Stage3pair_T *
Stage1_paired_read (int *npaths, int *first_absmq, int *second_absmq, Pairtype_T *final_pairtype,
		    Stage3end_T **stage3array5, int *nhits5, int *first_absmq5, int *second_absmq5,
		    Stage3end_T **stage3array3, int *nhits3, int *first_absmq3, int *second_absmq3,
		    Shortread_T queryseq5, Shortread_T queryseq3,
		    Indexdb_T indexdb_fwd, Indexdb_T indexdb_rev, int indexdb_size_threshold,
		    Floors_T *floors_array,
		    double usermax_level_float, double user_mincoverage_float,
		    int indel_penalty_middle, int indel_penalty_end,
		    bool allow_end_indels_p, int max_end_insertions, int max_end_deletions, int min_indel_end_matches,
		    int localsplicing_penalty, int distantsplicing_penalty, int min_shortend,
		    Oligoindex_array_T oligoindices_major, Oligoindex_array_T oligoindices_minor,
		    Pairpool_T pairpool, Diagpool_T diagpool, Cellpool_T cellpool,
		    Dynprog_T dynprogL, Dynprog_T dynprogM, Dynprog_T dynprogR,
		    Chrpos_T pairmax, bool keep_floors_p);

extern void
Stage1hr_cleanup ();

extern void
Stage1hr_setup (bool use_sarray_p_in, bool use_only_sarray_p_in, int index1part_in, int index1interval_in,
		int spansize_in, Univ_IIT_T chromosome_iit_in, int nchromosomes_in,
		Genome_T genome_in, Genome_T genomealt, Mode_T mode_in, int maxpaths_search_in,

		Univcoord_T *splicesites_in, Splicetype_T *splicetypes_in,
		Chrpos_T *splicedists_in, int nsplicesites_in,

		bool novelsplicingp_in, bool knownsplicingp_in, bool find_dna_chimeras_p_in,
		bool distances_observed_p_in, int subopt_levels_in,
		Chrpos_T max_middle_insertions_in, Chrpos_T max_middle_deletions_in,
		Chrpos_T shortsplicedist_in, Chrpos_T shortsplicedist_known_in, Chrpos_T shortsplicedist_novelend_in,
		Chrpos_T min_intronlength_in,

		int min_distantsplicing_end_matches_in, int min_distantsplicing_identity_in,

		int nullgap_in, int maxpeelback_in, int maxpeelback_distalmedial_in,
		int extramaterial_end_in, int extramaterial_paired_in,
		int gmap_mode, int trigger_score_for_gmap_in, int gmap_allowance_in,
		int max_gmap_pairsearch_in, int max_gmap_terminal_in,
		int max_gmap_improvement_in, int antistranded_penalty_in);


#undef T
#endif

