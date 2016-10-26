/* $Id: outbuffer.h 157571 2015-01-28 00:04:37Z twu $ */
#ifndef OUTBUFFER_INCLUDED
#define OUTBUFFER_INCLUDED

#include "types.h"
#include "bool.h"
#include "genomicpos.h"
#include "sequence.h"
#include "iit-read-univ.h"
#include "iit-read.h"
#include "samflags.h"
#include "filestring.h"

#include "request.h"
#include "mem.h"		/* To get MEMUSAGE */

#ifdef GSNAP
#include "resulthr.h"

#else
#include "stage3.h"		/* Has Printtype_T */
#include "genome.h"

#endif


#define T Outbuffer_T
typedef struct T *T;

extern void
Outbuffer_setup (int argc_in, char **argv_in, int optind_in,
		 Univ_IIT_T chromosome_iit_in, bool any_circular_p_in,
		 int nworkers_in, bool orderedp_in, bool quiet_if_excessive_p_in,
#ifdef GSNAP
		 bool output_sam_p_in, 
#else
		 Printtype_T printtype_in, Sequence_T usersegment_in,
#endif
		 bool sam_headers_p_in, char *sam_read_group_id_in, char *sam_read_group_name_in,
		 char *sam_read_group_library_in, char *sam_read_group_platform_in,
		 bool appendp_in, char *output_file_in, char *split_output_root_in, char *failedinput_root_in);

extern void
Outbuffer_cleanup ();

extern T
Outbuffer_new (unsigned int output_buffer_size, unsigned int nread);

extern void
Outbuffer_close_files ();

extern void
Outbuffer_free (T *old);

extern unsigned int
Outbuffer_nread (T this);

extern unsigned int
Outbuffer_nbeyond (T this);

extern void
Outbuffer_add_nread (T this, unsigned int nread);

#ifdef GSNAP
extern void
Outbuffer_put_filestrings (T this, Filestring_T fp, Filestring_T fp_failedinput_1, Filestring_T fp_failedinput_2);

extern void
Outbuffer_print_filestrings (Filestring_T fp, Filestring_T fp_failedinput_1, Filestring_T fp_failedinput_2);
#else
extern void
Outbuffer_put_filestrings (T this, Filestring_T fp, Filestring_T fp_failedinput);

extern void
Outbuffer_print_filestrings (Filestring_T fp, Filestring_T fp_failedinput_1);
#endif

extern void *
Outbuffer_thread_anyorder (void *data);

extern void *
Outbuffer_thread_ordered (void *data);

#ifdef USE_MPI
extern void
Outbuffer_mpi_process (T this, int n_worker_procs, int part_modulus, int part_interval);
#endif


#undef T
#endif

