/* $Id: sequence.h 170023 2015-07-17 16:47:21Z twu $ */
#ifndef SEQUENCE_INCLUDED
#define SEQUENCE_INCLUDED
#ifdef HAVE_CONFIG_H
#include <config.h>		/* For HAVE_ZLIB, HAVE_BZLIB */
#endif

#include <stdio.h>
#include "bool.h"
#include "filestring.h"

#ifdef HAVE_ZLIB
#include <zlib.h>
#endif

#ifdef HAVE_BZLIB
#include "bzip2.h"
#endif


#ifdef PMAP
#define MAXSEQLEN 300000
#else
#define MAXSEQLEN 1000000
#endif
#define HALFLEN MAXSEQLEN/2

#define T Sequence_T
typedef struct T *T;

extern bool
Sequence_firstp (T this);

extern char *
Sequence_accession (T this);

extern int
Sequence_input_init (FILE *fp);

#ifdef HAVE_ZLIB
extern int
Sequence_input_init_gzip (gzFile fp);
#endif

#ifdef HAVE_BZLIB
extern int
Sequence_input_init_bzip2 (Bzip2_T fp);
#endif

extern char *
Sequence_fullpointer (T this);
extern char *
Sequence_trimpointer (T this);
extern char *
Sequence_quality_string (T this);

extern int
Sequence_ntlength (T this);
extern int
Sequence_fulllength (T this);
extern char *
Sequence_subseq_pointer (T this, int querystart);
extern int
Sequence_subseq_length (T this, int querystart);
extern int
Sequence_trimlength (T this);
extern int
Sequence_fulllength_given (T this);
extern void
Sequence_trim (T this, int trim_start, int trim_end);
extern int
Sequence_trim_start (T this);
extern int
Sequence_trim_end (T this);
extern int
Sequence_subseq_offset (T this);
extern int
Sequence_skiplength (T this);

extern void
Sequence_free (T *old);
extern T
Sequence_genomic_new (char *contents, int length, bool copyp);
extern T
Sequence_read (int *nextchar, FILE *input);
extern T
Sequence_read_multifile (int *nextchar, FILE **input, char ***files, int *nfiles);
extern T
Sequence_read_unlimited (int *nextchar, FILE *input);


#ifdef PMAP
extern char
Sequence_codon_char (char aa, int codonpos);
extern T
Sequence_convert_to_nucleotides (T this);
#endif
extern int
Sequence_count_bad (T this, int pos, int max, int direction);

extern T
Sequence_subsequence (T this, int start, int end);
extern T
Sequence_revcomp (T this);
extern T
Sequence_uppercase (T this);
extern T
Sequence_alias (T this);


extern void
Sequence_print_digest (Filestring_T fp, T this);
extern void
Sequence_print_header (Filestring_T fp, T this, bool checksump);

extern void
Sequence_print (Filestring_T fp, T this, bool uppercasep, int wraplength, bool trimmedp);

extern void
Sequence_stdout (T this, bool uppercasep, int wraplength, bool trimmedp);
extern void
Sequence_stdout_alt (T ref, T alt, T snp, bool uppercasep, int wraplength);
extern void
Sequence_stdout_two (T ref, T alt, bool uppercasep, int wraplength);

extern void
Sequence_stdout_raw (T this);
extern void
Sequence_stdout_stream_chars (T this);
extern void
Sequence_stdout_stream_ints (T this);

extern T
Sequence_substring (T usersegment, unsigned int left, unsigned int length, 
		    bool revcomp);

#undef T
#endif
