/* $Id: parserange.h 157221 2015-01-22 18:38:57Z twu $ */
#ifndef PARSERANGE_INCLUDED
#define PARSERANGE_INCLUDED

#include "bool.h"
#include "genomicpos.h"
#include "types.h"
#include "iit-read-univ.h"

extern bool
Parserange_iscoordp (Univcoord_T *result, char *string);
extern bool
Parserange_islengthp (Chrpos_T *result, char *string);
extern bool
Parserange_israngep (Univcoord_T *left, Chrpos_T *length, bool *revcomp, char *string);

extern bool
Parserange_query (char **divstring, Univcoord_T *coordstart, Univcoord_T *coordend, bool *revcomp,
		  char *query, char *filename);


/* genomicstart is 0-based, chrstart is 1-based */
extern bool
Parserange_universal (char **div, bool *revcomp,
		      Univcoord_T *genomicstart, Chrpos_T *genomiclength,
		      Chrpos_T *chrstart, Chrpos_T *chrend,
		      Univcoord_T *chroffset, Chrpos_T *chrlength,
		      char *query, char *genomesubdir, char *fileroot);

extern bool
Parserange_universal_iit (char **div, bool *revcomp,
			  Univcoord_T *genomicstart, Chrpos_T *genomiclength,
			  Chrpos_T *chrstart, Chrpos_T *chrend,
			  Univcoord_T *chroffset, Chrpos_T *chrlength,
			  char *query, Univ_IIT_T chromosome_iit, Univ_IIT_T contig_iit);

extern bool
Parserange_simple (char **div, bool *revcomp, Chrpos_T *chrstart, Chrpos_T *chrend,
		   char *query);

#endif

