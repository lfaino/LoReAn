#ifndef BITPACK64_READ_INCLUDED
#define BITPACK64_READ_INCLUDED
#ifdef HAVE_CONFIG_H
#include <config.h>		/* For HAVE_64_BIT */
#endif

#include "types.h"

/* For reading differential-coded bitstreams */

extern UINT4
Bitpack64_read_one (Storedoligomer_T oligo,
		    UINT4 *bitpackptrs, UINT4 *bitpackcomp);

extern UINT8
Bitpack64_read_one_huge (Storedoligomer_T oligo, UINT4 *bitpackpages,
			 UINT4 *bitpackptrs, UINT4 *bitpackcomp);

#ifndef PMAP
extern void
Bitpack64_block_offsets (UINT4 *offsets, Storedoligomer_T oligo,
			 UINT4 *bitpackptrs, UINT4 *bitpackcomp);
#endif


#ifndef PMAP
#if defined(HAVE_64_BIT) && (defined(UTILITYP) || defined(LARGE_GENOMES))
extern void
Bitpack64_block_offsets_huge (UINT8 *offsets, Storedoligomer_T oligo,
			      UINT4 *bitpackpages, UINT4 *bitpackptrs, UINT4 *bitpackcomp);
#endif
#endif

#endif
