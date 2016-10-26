/* $Id: mapq.h 133760 2014-04-20 05:16:56Z twu $ */
#ifndef MAPQ_INCLUDED
#define MAPQ_INCLUDED

#include "types.h"
#include "compress.h"
#include "genomicpos.h"


#define MAX_QUALITY_SCORE_INPUT 96	/* Was 40 */
#define MAX_QUALITY_SCORE 40

extern void
MAPQ_init (int quality_score_adj_in);
extern int
MAPQ_max_quality_score (char *quality_string, int querylength);
extern float
MAPQ_loglik_exact (char *quality_string, int querystart, int queryend);
extern float
MAPQ_loglik (Compress_T query_compress, Univcoord_T left, int querystart, int queryend,
	     int querylength, char *quality_string, bool plusp, int genestrand, bool first_read_p);

#endif

