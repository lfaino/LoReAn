/* $Id: samheader.h 157230 2015-01-22 18:49:34Z twu $ */
#ifndef SAMHEADER_INCLUDED
#define SAMHEADER_INCLUDED

#ifdef USE_MPI
#include <mpi.h>
#endif

#include <stdio.h>
#include "bool.h"
#include "samflags.h"
#include "filestring.h"

#ifdef USE_MPI
extern MPI_File
#else
extern FILE *
#endif
SAM_header_open_file (SAM_split_output_type split_output, char *split_output_root, bool appendp);

extern Filestring_T
SAM_header_change_HD_tosorted (FILE *input, int headerlen);

#ifdef USE_MPI
extern void
SAM_header_print_HD (MPI_File fp, int nworkers, bool orderedp);
#else
extern void
SAM_header_print_HD (FILE *fp, int nworkers, bool orderedp);
#endif

#ifdef USE_MPI
extern void
SAM_header_print_PG (MPI_File fp, int argc, char **argv, int optind);
#else
extern void
SAM_header_print_PG (FILE *fp, int argc, char **argv, int optind);
#endif

extern int
SAM_header_length (int *lastchar, FILE *fp);

#endif

