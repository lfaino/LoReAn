/* $Id: genome128_hr.h 166739 2015-06-02 01:23:18Z twu $ */
#ifndef GENOME128_HR_INCLUDED
#define GENOME128_HR_INCLUDED
#include "types.h"
#include "mode.h"
#include "genomicpos.h"
#include "compress.h"

extern void
Genome_hr_setup (Genomecomp_T *ref_blocks_in, Genomecomp_T *snp_blocks_in,
		 bool query_unk_mismatch_p_in, bool genome_unk_mismatch_p_in,
		 Mode_T mode_in);

extern void
Genome_hr_user_setup (UINT4 *ref_blocks_in,
		      bool query_unk_mismatch_p_in, bool genome_unk_mismatch_p_in,
		      Mode_T mode_in);

extern int
Genome_consecutive_matches_rightward (Compress_T query_compress, Univcoord_T left, int pos5, int pos3,
				      bool plusp, int genestrand, bool first_read_p);
extern int
Genome_consecutive_matches_leftward (Compress_T query_compress, Univcoord_T left, int pos5, int pos3,
				     bool plusp, int genestrand, bool first_read_p);
extern int
Genome_consecutive_matches_pair (UINT4 lefta, UINT4 leftb, UINT4 genomelength);

extern int
Genome_count_mismatches_limit (Compress_T query_compress, Univcoord_T left, int pos5, int pos3,
			       int max_mismatches, bool plusp, int genestrand, bool first_read_p);
extern int
Genome_count_mismatches_substring_ref (Compress_T query_compress, Univcoord_T left, int pos5, int pos3,
				       bool plusp, int genestrand, bool first_read_p);
extern int
Genome_count_mismatches_substring (Compress_T query_compress, Univcoord_T left, int pos5, int pos3,
				   bool plusp, int genestrand, bool first_read_p);

extern int
Genome_count_mismatches_fragment_left (Compress_T query_compress, int pos5, int pos3,
				       Genomecomp_T ref_fragment, Genomecomp_T alt_fragment);
extern int
Genome_count_mismatches_fragment_right (Compress_T query_compress, int pos5, int pos3,
					Genomecomp_T ref_fragment, Genomecomp_T alt_fragment);

extern int
Genome_mismatches_left (int *mismatch_positions, int max_mismatches, Compress_T query_compress,
			Univcoord_T left, int pos5, int pos3, bool plusp, int genestrand, bool first_read_p);
extern int
Genome_mismatches_left_trim (int *mismatch_positions, int max_mismatches, Compress_T query_compress,
			     Univcoord_T left, int pos5, int pos3, bool plusp, int genestrand, bool first_read_p);
extern int
Genome_mismatches_right (int *mismatch_positions, int max_mismatches, Compress_T query_compress,
			 Univcoord_T left, int pos5, int pos3, bool plusp, int genestrand, bool first_read_p);
extern int
Genome_mismatches_right_trim (int *mismatch_positions, int max_mismatches, Compress_T query_compress,
			      Univcoord_T left, int pos5, int pos3, bool plusp, int genestrand, bool first_read_p);

extern int
Genome_mark_mismatches_ref (char *genomic, int querylength, Compress_T query_compress,
			    Univcoord_T left, int pos5, int pos3,
			    bool plusp, int genestrand, bool first_read_p);
extern int
Genome_mark_mismatches (char *genomic, int querylength, Compress_T query_compress,
			Univcoord_T left, int pos5, int pos3,
			bool plusp, int genestrand, bool first_read_p);

extern int
Genome_trim_left (Compress_T query_compress, Univcoord_T left, int pos5, int pos3,
		  bool plusp, int genestrand, bool first_read_p);

extern int
Genome_trim_right (Compress_T query_compress, Univcoord_T left, int pos5, int pos3,
		   bool plusp, int genestrand, bool first_read_p);

#endif

