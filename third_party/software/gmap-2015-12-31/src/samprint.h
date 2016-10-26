/* $Id: samprint.h 183725 2016-02-04 00:40:15Z twu $ */
#ifndef SAMPRINT_INCLUDED
#define SAMPRINT_INCLUDED

#include <stdio.h>
#include "iit-read-univ.h"
#include "iit-read.h"
#include "genomicpos.h"
#include "types.h"
#include "substring.h"
#include "bool.h"
#include "intlist.h"
#include "filestring.h"

#ifdef GSNAP
#include "shortread.h"
#include "stage3hr.h"
#include "resulthr.h"
#include "genome.h"
#endif


#ifdef GSNAP
extern void
SAM_setup (bool add_paired_nomappers_p_in, bool paired_flag_means_concordant_p_in,
	   bool quiet_if_excessive_p_in, int maxpaths_report_in,
	   char *failedinput_root_in, bool fastq_format_p_in, bool hide_soft_clips_p_in,
	   bool clip_overlap_p_in, bool merge_overlap_p_in, bool sam_multiple_primaries_p_in,
	   bool force_xs_direction_p_in, bool md_lowercase_variant_p_in, IIT_T snps_iit_in,
	   Univ_IIT_T chromosome_iit_in, Genome_T genome_in);

extern Chrpos_T
SAM_compute_chrpos (int hardclip_low, int hardclip_high, Stage3end_T this, int querylength,
		    bool first_read_p);

extern unsigned int
SAM_compute_flag (bool plusp, Stage3end_T mate, Resulttype_T resulttype,
		  bool first_read_p, int pathnum, int npaths, bool artificial_mate_p, int npaths_mate,
		  int absmq_score, int first_absmq, bool invertp, bool invert_mate_p);

extern void
SAM_print_nomapping (Filestring_T fp, char *abbrev, Shortread_T queryseq, Stage3end_T mate, char *acc1, char *acc2,
		     Univ_IIT_T chromosome_iit, Resulttype_T resulttype, bool first_read_p,
		     int pathnum, int npaths, bool artificial_mate_p, int npaths_mate,
		     Chrpos_T mate_chrpos, int quality_shift, char *sam_read_group_id, bool invertp, bool invert_mate_p);

extern void
SAM_print (Filestring_T fp, Filestring_T fp_failedinput, char *abbrev,
	   Stage3end_T this, Stage3end_T mate, char *acc1, char *acc2, int pathnum, int npaths,
	   int absmq_score, int first_absmq, int second_absmq, int mapq_score, Univ_IIT_T chromosome_iit, Shortread_T queryseq,
	   Shortread_T queryseq2, int pairedlength, Chrpos_T chrpos, Chrpos_T mate_chrpos,
	   int clipdir, int hardclip5_low, int hardclip5_high, int hardclip3_low, int hardclip3_high,
	   Resulttype_T resulttype, bool first_read_p, bool artificial_mate_p, int npaths_mate, int quality_shift,
	   char *sam_read_group_id, bool invertp, bool invert_mate_p, bool merge_samechr_p);

extern void
SAM_print_paired (Filestring_T fp, Filestring_T fp_failedinput_1, Filestring_T fp_failedinput_2,
		  Result_T result, Resulttype_T resulttype, Univ_IIT_T chromosome_iit,
		  Shortread_T queryseq1, Shortread_T queryseq2, bool invert_first_p, bool invert_second_p,
		  bool nofailsp, bool failsonlyp, bool merge_samechr_p,
		  int quality_shift, char *sam_read_group_id);
#endif

#endif

