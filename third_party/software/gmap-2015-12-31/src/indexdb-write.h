/* $Id: indexdb-write.h 165969 2015-05-20 00:18:07Z twu $ */
#ifndef INDEXDB_WRITE_INCLUDED
#define INDEXDB_WRITE_INCLUDED
#ifdef HAVE_CONFIG_H
#include <config.h>		/* For HAVE_64_BIT */
#endif

#include "bool.h"
#include "types.h"		/* For Oligospace_T */
#include "iit-read-univ.h"

#ifdef PMAP
#include "alphabet.h"
#define FWD_FILESUFFIX "pf"
#define REV_FILESUFFIX "pr"

#else
#define IDX_FILESUFFIX "ref"
#endif

#define OFFSETS_FILESUFFIX "offsets"
#define POSITIONS_HIGH_FILESUFFIX "positionsh"
#define POSITIONS_LOW_FILESUFFIX "positions"

#ifdef HAVE_64_BIT
extern UINT8
Indexdb_count_offsets (FILE *sequence_fp, Univ_IIT_T chromosome_iit,
#ifdef PMAP
		       Width_T index1part_aa, bool watsonp,
#else
		       Width_T index1part,
#endif
		       Width_T index1interval, bool genome_lc_p, char *fileroot, bool mask_lowercase_p);
#endif

extern void
Indexdb_write_offsets (char *destdir, char interval_char, FILE *sequence_fp, Univ_IIT_T chromosome_iit,
#ifdef PMAP
		       Alphabet_T alphabet, Width_T index1part_aa, bool watsonp,
#else
		       Width_T index1part,
#endif
		       Width_T index1interval, bool genome_lc_p, char *fileroot, bool mask_lowercase_p,
		       int compression_types);
#ifdef HAVE_64_BIT
extern void
Indexdb_write_offsets_huge (char *destdir, char interval_char, FILE *sequence_fp, Univ_IIT_T chromosome_iit,
#ifdef PMAP
			    Alphabet_T alphabet, Width_T index1part_aa, bool watsonp,
#else
			    Width_T index1part,
#endif
			    Width_T index1interval, bool genome_lc_p, char *fileroot, bool mask_lowercase_p,
			    int compression_types);
#endif


extern void
Indexdb_write_positions (char *positionsfile_high, char *positionsfile_low, char *pointersfile, char *offsetsfile,
			 FILE *sequence_fp, Univ_IIT_T chromosome_iit,
#ifdef PMAP
			 Alphabet_T alphabet, Width_T index1part_aa, bool watsonp,
#else
			 Width_T index1part,
#endif
			 Width_T index1interval, Univcoord_T genomelength, bool genome_lc_p, bool writefilep,
			 char *fileroot, bool mask_lowercase_p, int compression_type,
			 bool coord_values_8p);


#ifdef HAVE_64_BIT
extern void
Indexdb_write_positions_huge (char *positionsfile_high, char *positionsfile_low, char *pagesfile, char *pointersfile, char *offsetsfile,
			      FILE *sequence_fp, Univ_IIT_T chromosome_iit,
#ifdef PMAP
			      Alphabet_T alphabet, int index1part_aa, bool watsonp,
#else
			      int index1part,
#endif
			      int index1interval, Univcoord_T genomelength, bool genome_lc_p, bool writefilep,
			      char *fileroot, bool mask_lowercase_p, int compression_type,
			      bool coord_values_8p);
#endif

#endif

