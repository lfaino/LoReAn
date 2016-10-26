/* $Id: sarray-write.h 165971 2015-05-20 00:20:26Z twu $ */
#ifndef SARRAY_WRITE_INCLUDED
#define SARRAY_WRITE_INCLUDED
#include "types.h"
#include "genome.h"

/* If conversion is NULL, then no conversion is performed */
extern void
Sarray_write_array (char *sarrayfile, Genome_T genomecomp, UINT4 genomelength);
extern void
Sarray_write_array_from_genome (char *sarrayfile, unsigned char *gbuffer, UINT4 genomelength);

extern void
Sarray_write_index_separate (char *indexiptrsfile, char *indexicompfile, char *indexjptrsfile,char *indexjcompfile,
			     char *sarrayfile, Genome_T genomecomp, UINT4 genomelength, bool compressp,
			     char chartable[]);

extern void
Sarray_write_index_interleaved (char *indexptrsfile, char *indexcompfile,
				char *sarrayfile, Genome_T genomecomp, UINT4 genomelength, bool compressp,
				char chartable[]);

extern void
Sarray_write_csa (char **csaptrfiles, char **csacompfiles, char *sasampleqfile, char *sasamplesfile, char *saindex0file,
		  char *sarrayfile, char *rankfile, Genome_T genomecomp, UINT4 genomelength, char chartable[]);

extern UINT4 *
Sarray_compute_lcp (char *rankfile, char *permuted_sarray_file, char *sarrayfile, UINT4 n);
extern UINT4 *
Sarray_compute_lcp_from_genome (UINT4 *SA, unsigned char *gbuffer, UINT4 n);

extern unsigned char *
Sarray_discriminating_chars (UINT4 *nbytes, char *sarrayfile, Genome_T genome,
			     unsigned char *lcp_bytes, UINT4 *lcp_guide, UINT4 *lcp_exceptions, int guide_interval,
			     UINT4 n, char chartable[]);

extern UINT4 *
Sarray_compute_child (unsigned char *lcp_bytes, UINT4 *lcp_guide, UINT4 *lcp_exceptions, UINT4 n);


extern void
Sarray_array_uncompress (Genome_T genomecomp, char *sarrayfile, char *plcpptrsfile, char *plcpcompfile,
			 UINT4 genomelength, UINT4 start, UINT4 end);

extern void
Sarray_child_uncompress (Genome_T genomecomp, unsigned char *lcpchilddc, UINT4 *lcp_guide, UINT4 *lcp_exceptions,
			 int n_lcp_exceptions, UINT4 *child_guide, UINT4 *child_exceptions, int n_child_exceptions,
			 UINT4 *SA, UINT4 genomelength, UINT4 start, UINT4 end);

extern void
Sarray_child_test (char *sarrayfile,
		   char *childbpfile, char *childs_pagesfile, char *childs_ptrsfile, char *childs_compfile,
		   char *childr_ptrsfile, char *childr_compfile, char *childx_ptrsfile, char *childx_compfile,
		   char *pioneerbpfile, char *pior_ptrsfile, char *pior_compfile,
		   char *piom_ptrsfile, char *piom_compfile, UINT4 genomelength, int sampling_interval);

#endif

