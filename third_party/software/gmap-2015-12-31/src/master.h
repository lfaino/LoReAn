/* $Id: master.h 162088 2015-03-26 18:29:04Z twu $ */
#ifndef MASTER_INCLUDED
#define MASTER_INCLUDED
#ifdef HAVE_CONFIG_H
#include <config.h>		/* For HAVE_ZLIB, HAVE_BZLIB, USE_MPI_FILE_INPUT */
#endif

#ifdef USE_MPI
#include <mpi.h>
#include "mpidebug.h"
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


#define T Master_T
typedef struct T *T;


extern int
Master_ntotal (T this);

extern T
Master_new (int n_slave_ranks, int nextchar, int nchars1, int nchars2,
	    FILE *input_parser, FILE *input2_parser,
#ifdef HAVE_ZLIB
	    gzFile gzipped, gzFile gzipped2,
#endif
#ifdef HAVE_BZLIB
	    Bzip2_T bzipped, Bzip2_T bzipped2,
#endif
	    char **files, int nfiles, int nspaces, int part_modulus, int part_interval);

extern void
Master_free (T *old);

extern void *
Master_write_stdout (void *data);

extern void *
Master_parser (void *data);

extern void
Master_self_interface (T this, int *nextchar_start, int *nextchar,
		       int *offset_start_1, int *offset_start_2,
		       int *offset_end_1, int *offset_end_2,
		       Filestring_T *filestring1, Filestring_T *filestring2,
		       bool *donep);

extern void *
Master_mpi_interface (void *data);

#undef T
#endif

