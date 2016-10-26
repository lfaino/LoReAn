/* $Id: shortread.h 157572 2015-01-28 00:05:22Z twu $ */
#ifndef SHORTREAD_INCLUDED
#define SHORTREAD_INCLUDED
#ifdef HAVE_CONFIG_H
#include <config.h>		/* For HAVE_ZLIB, HAVE_BZLIB, USE_MPI_FILE_INPUT  */
#endif

#include <stdio.h>
#include "bool.h"
#include "filestring.h"

#if defined(USE_MPI) && defined(USE_MPI_FILE_INPUT)
#include <mpi.h>
#endif

#ifdef HAVE_ZLIB
#include <zlib.h>
#endif

#ifdef HAVE_BZLIB
#include "bzip2.h"
#endif


#define T Shortread_T
typedef struct T *T;

extern void
Shortread_setup (int acc_fieldi_start_in, int acc_fieldi_end_in,
		 bool force_singled_end_p_in, bool filter_chastity_p_in,
		 bool allow_paired_end_mismatch_p_in, bool fastq_format_p_in,
		 int barcode_length_in, bool invert_first_p_in, bool invert_second_p_in);

extern char *
Shortread_accession (T this);
extern char *
Shortread_header (T this);
extern bool
Shortread_filterp (T this);
extern bool
Shortread_invertedp (T this);

#if 0
extern unsigned long long **
Shortread_input_divide (bool *fastq_format_p, char **files, int nfiles, int naliquots);
#endif

extern int
Shortread_input_init (int *nchars, FILE *fp);

#ifdef HAVE_ZLIB
extern int
Shortread_input_init_gzip (gzFile fp);
#endif

#ifdef HAVE_BZLIB
extern int
Shortread_input_init_bzip2 (Bzip2_T fp);
#endif

extern char *
Shortread_fullpointer (T this);
extern char *
Shortread_trimpointer (T this);

extern char *
Shortread_fullpointer_uc (T this);
extern char *
Shortread_contents_uc (T this);

extern int
Shortread_barcode_length (T this);
extern char *
Shortread_barcode (T this);

extern int
Shortread_choplength (T this);

extern char *
Shortread_quality_string (T this);

extern int
Shortread_fulllength (T this);

extern void
Shortread_free (T *old);
extern T
Shortread_new (char *acc, char *restofheader, bool filterp,
	       char *short_sequence, char *long_sequence, int sequence_length,
	       char *short_quality, char *long_quality, int quality_length,
	       int barcode_length, bool invertp, bool copy_acc_p, bool skipp);

extern bool
Shortread_chop_primers (T queryseq1, T queryseq2);
extern bool
Shortread_find_primers (T queryseq1, T queryseq2);
extern int
Shortread_max_overlap (T queryseq1, T queryseq2);
extern int
Shortread_find_overlap (T queryseq1, T queryseq2);

extern T
Shortread_new (char *acc, char *restofheader, bool filterp,
	       char *short_sequence, char *long_sequence, int sequence_length,
	       char *short_quality, char *long_quality, int quality_length,
	       int barcode_length, bool invertp, bool copy_acc_p, bool skipp);

extern T
Shortread_read_fastq_text (int *nextchar, int *nchars1, int *nchars2, T *queryseq2,
			   FILE **input1, FILE **input2,
			   char ***files, int *nfiles, bool skipp);

#ifdef HAVE_ZLIB
extern T
Shortread_read_fastq_gzip (int *nextchar, T *queryseq2,
#ifdef USE_MPI
			   Filestring_T filestring1, Filestring_T filestring2,
#endif
			   gzFile *input1, gzFile *input2,
			   char ***files, int *nfiles, bool skipp);
#endif

#ifdef HAVE_BZLIB
extern T
Shortread_read_fastq_bzip2 (int *nextchar, T *queryseq2,
#ifdef USE_MPI
			    Filestring_T filestring1, Filestring_T filestring2,
#endif
			    Bzip2_T *input1, Bzip2_T *input2,
			    char ***files, int *nfiles, bool skipp);
#endif

extern T
Shortread_read_fasta_text (int *nextchar, int *nchars1, int *nchars2, T *queryseq2,
			   FILE **input1, FILE **input2,
			   char ***files, int *nfiles, bool skipp);
#ifdef HAVE_ZLIB
extern T
Shortread_read_fasta_gzip (int *nextchar, T *queryseq2,
#ifdef USE_MPI
			   Filestring_T filestring1, Filestring_T filestring2,
#endif
			   gzFile *input1, gzFile *input2,
			   char ***files, int *nfiles, bool skipp);
#endif

#ifdef HAVE_BZLIB
extern T
Shortread_read_fasta_bzip2 (int *nextchar, T *queryseq2,
#ifdef USE_MPI
			    Filestring_T filestring1, Filestring_T filestring2,
#endif
			    Bzip2_T *input1, Bzip2_T *input2,
			    char ***files, int *nfiles, bool skipp);
#endif


extern T
Shortread_read (int *nextchar, int *nchars1, int *nchars2, T *queryseq2,
#ifdef USE_MPI
		Filestring_T filestring1, Filestring_T filestring2,
#endif
		FILE **input1, FILE **input2,
#ifdef HAVE_ZLIB
		gzFile *gzipped1, gzFile *gzipped2,
#endif
#ifdef HAVE_BZLIB
		Bzip2_T *bzipped1, Bzip2_T *bzipped2,
#endif
		char ***files, int *nfiles, bool skipp);

#ifdef USE_MPI
extern T
Shortread_read_filecontents (int *nextchar, T *queryseq2,
			     char **filecontents1, char **filecontents2,
#ifdef USE_MPI_FILE_INPUT
			     MPI_File *input1, MPI_File *input2, MPI_Comm workers_comm,
#else
			     FILE **input1, FILE **input2,
#endif
			     char ***files, int *nfiles, bool skipp);
#endif

extern void
Shortread_print_header (Filestring_T fp, T queryseq1, T queryseq2);

extern void
Shortread_stderr_query_singleend_fasta (T queryseq, T headerseq);
extern void
Shortread_stderr_query_pairedend_fasta (T queryseq1, T queryseq2,
					bool invert_first_p, bool invert_second_p);

extern void
Shortread_print_query_singleend_fastq (Filestring_T fp, T queryseq, T headerseq);
extern void
Shortread_print_query_pairedend_fastq (Filestring_T fp1, Filestring_T fp2, T queryseq1, T queryseq2,
				      bool invert_first_p, bool invert_second_p);

extern void
Shortread_print_query_singleend (Filestring_T fp, T queryseq, T headerseq);
extern void
Shortread_print_query_pairedend (Filestring_T fp1, Filestring_T fp2, T queryseq1, T queryseq2);

extern void
Shortread_print_oneline (Filestring_T fp, T this);
extern void
Shortread_print_oneline_revcomp (Filestring_T fp, T this);

extern void
Shortread_print_chopped_sam (Filestring_T fp, T this, int hardclip_low, int hardclip_high);
extern void
Shortread_print_chopped_revcomp_sam (Filestring_T fp, T this, int hardclip_low, int hardclip_high);
extern void
Shortread_print_chopped_end (Filestring_T fp, T this, int hardclip_low, int hardclip_high);
extern void
Shortread_print_chopped_end_revcomp (Filestring_T fp, T this, int hardclip_low, int hardclip_high);
extern void
Shortread_print_chopped_end_quality (Filestring_T fp, T this, int hardclip_low, int hardclip_high);
extern void
Shortread_print_chopped_end_quality_reverse (Filestring_T fp, T this, int hardclip_low, int hardclip_high);


extern void
Shortread_print_barcode (Filestring_T fp, T this);
extern void
Shortread_print_chop (Filestring_T fp, T this, bool invertp);
extern void
Shortread_print_chop_symbols (Filestring_T fp, T this);
extern void
Shortread_print_quality (Filestring_T fp, T this, int hardclip_low, int hardclip_high,
			int shift, bool show_chopped_p);
extern void
Shortread_print_quality_revcomp (Filestring_T fp, T this, int hardclip_low, int hardclip_high,
				int shift, bool show_chopped_p);
extern void
Shortread_print_oneline_uc (Filestring_T fp, T this);
extern void
Shortread_print_oneline_revcomp_uc (Filestring_T fp, T this);

#undef T
#endif
