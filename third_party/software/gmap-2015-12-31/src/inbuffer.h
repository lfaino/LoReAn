/* $Id: inbuffer.h 159549 2015-02-25 22:23:14Z twu $ */
#ifndef INBUFFER_INCLUDED
#define INBUFFER_INCLUDED
#ifdef HAVE_CONFIG_H
#include <config.h>		/* For HAVE_ZLIB, HAVE_BZLIB, USE_MPI_FILE_INPUT */
#endif

#ifdef USE_MPI
#include <mpi.h>
#include "mpidebug.h"
#include "master.h"
#endif

#include <stdio.h>
#include "bool.h"
#include "outbuffer.h"
#include "request.h"

#ifdef GSNAP
#include "shortread.h"
#else
#include "sequence.h"
#endif

#ifdef HAVE_ZLIB
#include <zlib.h>
#endif

#ifdef HAVE_BZLIB
#include "bzip2.h"
#endif


#define T Inbuffer_T
typedef struct T *T;

extern void
Inbuffer_setup (bool filter_if_both_p_in, 
#if defined(USE_MPI) && defined(USE_MPI_FILE_INPUT)
		MPI_Comm workers_comm_in,
#endif
#ifndef GSNAP
		bool user_pairalign_p_in, Sequence_T global_usersegment_in,
#endif
		int part_modulus_in, int part_interval_in);

#ifndef GSNAP
extern T
Inbuffer_cmdline (char *contents, int length);
#endif

extern T
Inbuffer_new (int nextchar,
#ifdef USE_MPI
	      int myid,
#endif
#if defined(USE_MPI) && defined(USE_MPI_FILE_INPUT)
	      MPI_File input,
#else
	      FILE *input,
#endif
#ifdef GSNAP
#if defined(USE_MPI) && defined(USE_MPI_FILE_INPUT)
	      MPI_File input2,
#else
	      FILE *input2,
#endif
#ifdef HAVE_ZLIB
	      gzFile gzipped, gzFile gzipped2,
#endif
#ifdef HAVE_BZLIB
	      Bzip2_T bzipped, Bzip2_T bzipped2,
#endif
#endif
	      char **files, int nfiles, unsigned int nspaces);

#ifdef USE_MPI
extern void
Inbuffer_set_master (T this, Master_T master);
#endif

extern void
Inbuffer_set_outbuffer (T this, Outbuffer_T outbuffer);

extern void
Inbuffer_free (T *old);

#ifdef GSNAP
extern Shortread_T
Inbuffer_read (Shortread_T *queryseq2, T this, bool skipp);
#else
Sequence_T
Inbuffer_read (Sequence_T *pairalign_segment, T this, bool skipp);
#endif

#ifndef USE_MPI
extern unsigned int
Inbuffer_fill_init (T this);
#endif


extern Request_T
#ifdef GSNAP
Inbuffer_get_request (T this);
#else
Inbuffer_get_request (Sequence_T *pairalign_segment, T this);
#endif

#ifndef GSNAP
extern Request_T
Inbuffer_first_request (T this);
#endif

#ifdef USE_MPI
extern int
Inbuffer_master_process (int n_worker_ranks, int nextchar, int nchars1, int nchars2,
			 FILE *input, FILE *input2,
#ifdef HAVE_ZLIB
			 gzFile gzipped, gzFile gzipped2,
#endif
#ifdef HAVE_BZLIB
			 Bzip2_T bzipped, Bzip2_T bzipped2,
#endif
			 char **files, int nfiles, int nspaces, int part_modulus, int part_interval);
#endif

#undef T
#endif

