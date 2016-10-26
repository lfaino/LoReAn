/* $Id: bytecoding.h 170515 2015-07-23 23:03:24Z twu $ */
#ifndef BYTECODING_INCLUDED
#define BYTECODING_INCLUDED

#include "bool.h"
#include "types.h"

extern void
Bytecoding_write (char *bytesfile, char *excfile, char *guidefile, UINT4 *values,
		  UINT4 genomelength, int guide_interval);
extern unsigned char *
Bytecoding_write_exceptions_only (char *excfile, char *guidefile, UINT4 *values,
				  UINT4 genomelength, int guide_interval);

extern void
Bytecoding_write_lcpchilddc (char *bytesfile, char *excfile, char *guidefile, UINT4 *child,
			     unsigned char *discrim_chars, unsigned char *lcpbytes,
			     UINT4 genomelength, int guide_interval);
#if 0
extern void
Bytecoding_write_lcpchilddcn (char *bytesfile, char *excfile, char *guidefile, UINT4 *child,
			      unsigned char *discrim_chars, unsigned char *lcpbytes,
			      UINT8 *predictive_nextp, UINT4 genomelength, int guide_interval);
#endif

extern UINT4
Bytecoding_read (UINT4 key, unsigned char *bytes, UINT4 *exceptions, int nexceptions);
extern UINT4
Bytecoding_read_wguide (UINT4 key, unsigned char *bytes, UINT4 *guide, UINT4 *exceptions,
			int guide_interval);
extern UINT4
Bytecoding_lcpchilddc_lcp (UINT4 key, unsigned char *bytes, UINT4 *exceptions, int nexceptions);
extern char
Bytecoding_lcpchilddc_dc (char *c1, UINT4 key, unsigned char *bytes);
extern UINT4
Bytecoding_lcpchilddc_child_up (UINT4 key, unsigned char *bytes, UINT4 *guide, UINT4 *exceptions, int guide_interval);
extern UINT4
Bytecoding_lcpchilddc_child_next (UINT4 key, unsigned char *bytes, UINT4 *guide, UINT4 *exceptions, int guide_interval);
extern UINT4
Bytecoding_lcpchilddc_lcp_next (UINT4 *child_next, UINT4 key,
				unsigned char *lcpchilddc, UINT4 *child_guide,
				UINT4 *child_exceptions, int child_guide_interval,
				UINT4 *lcp_exceptions, int n_lcp_exceptions);

#if 0
extern UINT4
Bytecoding_lcpchilddcn_child_up (bool *nextp, UINT4 key, unsigned char *bytes, UINT4 *guide, UINT4 *exceptions, int guide_interval);
extern UINT4
Bytecoding_lcpchilddcn_child_next (bool *nextp, UINT4 key, unsigned char *bytes, UINT4 *guide, UINT4 *exceptions, int guide_interval);
#endif


#endif


