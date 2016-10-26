/* $Id: samread.h 154089 2014-11-25 21:03:16Z twu $ */
#ifndef SAMREAD_INCLUDED
#define SAMREAD_INCLUDED
#include <stdio.h>
#include "samflags.h"
#include "genomicpos.h"
#include "intlist.h"
#include "uintlist.h"
#include "iit-read-univ.h"


extern char *
Samread_get_acc_fromfile (int *acclength, FILE *fp, int linelength);

extern char *
Samread_parse_acc_and_softclip_fromfile (int *acclength, unsigned int *flag, SAM_split_output_type *split_output,
					 char **hiti, Univcoord_T *genomicpos, int *initial_softclip, bool *query_lowp,
					 FILE *fp, Univ_IIT_T chromosome_iit, Univcoord_T *chroffsets, int linelength);

extern int
Samread_parse_linelen_fromfile (FILE *fp);

extern Univcoord_T
Samread_parse_genomicpos_fromfile (FILE *fp, unsigned int *flag, SAM_split_output_type *split_output,
				   Univ_IIT_T chromosome_iit, Univcoord_T *chroffsets, int linelength);

extern void
Samread_parse_read_fromfile (FILE *fp, unsigned int *flag, int *readlength, char **read, int linelength);

extern Univcoord_T
Samread_parse_mate_genomicpos_fromfile (FILE *fp, Univ_IIT_T chromosome_iit, Univcoord_T *chroffsets, int linelength);

extern char *
Samread_parse_aux_fromfile (FILE *fp, char *auxfield, int linelength);

extern void
Samread_print_as_duplicate_fromfile (FILE *fp, int linelength);

#endif

