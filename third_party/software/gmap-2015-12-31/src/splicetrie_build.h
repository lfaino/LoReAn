#ifndef SPLICETRIE_BUILD_INCLUDED
#define SPLICETRIE_BUILD_INCLUDED

typedef enum {DONOR, ANTIDONOR, ACCEPTOR, ANTIACCEPTOR} Splicetype_T;

#include "bool.h"
#include "genomicpos.h"
#include "types.h"
#include "iit-read-univ.h"
#include "iit-read.h"
#include "genome.h"
#include "list.h"
#include "splicestringpool.h"


/* For offsets */
/* #define USE_2BYTE_RELOFFSETS 1 */
#ifdef USE_2BYTE_RELOFFSETS
#define NULL_POINTER 65535
#else
#define NULL_POINTER -1U    /* Note: 0 does not work */
#endif

#define MAX_DUPLICATES 1000
#define DUPLICATE_NODE -1000U	/* Needs to be -MAX_DUPLICATES */

#define INTERNAL_NODE -1U

#define single_leaf_p(x) (x) < DUPLICATE_NODE
#define multiple_leaf_p(x) (x) < INTERNAL_NODE

extern char *
Splicetype_string (Splicetype_T splicetype);

extern Univcoord_T *
Splicetrie_retrieve_via_splicesites (bool *distances_observed_p,
#ifdef GSNAP
				     Genomecomp_T **splicecomp,
#endif
				     Splicetype_T **splicetypes, Chrpos_T **splicedists,
				     List_T **splicestrings, Genomecomp_T **splicefrags_ref, Genomecomp_T **splicefrags_alt,
				     int *nsplicesites, IIT_T splicing_iit, int *splicing_divint_crosstable,
				     int donor_typeint, int acceptor_typeint, Univ_IIT_T chromosome_iit,
				     Genome_T genome, Genome_T genomealt, Chrpos_T shortsplicedist,
				     Splicestringpool_T splicestringpool);
extern Univcoord_T *
Splicetrie_retrieve_via_introns (
#ifdef GSNAP
				 Genomecomp_T **splicecomp,
#endif
				 Splicetype_T **splicetypes, Chrpos_T **splicedists,
				 List_T **splicestrings, Genomecomp_T **splicefrags_ref, Genomecomp_T **splicefrags_alt,
				 int *nsplicesites, IIT_T splicing_iit, int *splicing_divint_crosstable,
				 Univ_IIT_T chromosome_iit, Genome_T genome, Genome_T genomealt,
				 Splicestringpool_T splicestringpool);
extern void
Splicetrie_npartners (int **nsplicepartners_skip, int **nsplicepartners_obs, int **nsplicepartners_max,
		      Univcoord_T *splicesites, Splicetype_T *splicetypes,
		      Chrpos_T *splicedists, List_T *splicestrings, int nsplicesites,
		      Univ_IIT_T chromosome_iit, Chrpos_T max_distance,
		      bool distances_observed_p);

extern void
Splicetrie_build_via_splicesites (Triecontent_T **triecontents_obs, Trieoffset_T **trieoffsets_obs,
				  Triecontent_T **triecontents_max, Trieoffset_T **trieoffsets_max,
				  int *nsplicepartners_skip, int *nsplicepartners_obs, int *nsplicepartners_max,
				  Splicetype_T *splicetypes, List_T *splicestrings, int nsplicesites);

extern void
Splicetrie_build_via_introns (Triecontent_T **triecontents_obs, Trieoffset_T **trieoffsets_obs,
			      Univcoord_T *splicesites, Splicetype_T *splicetypes,
			      List_T *splicestrings, int nsplicesites,
			      Univ_IIT_T chromosome_iit, IIT_T splicing_iit, int *splicing_divint_crosstable);

#endif

