#ifndef BITPACK64_ACCESS_INCLUDED
#define BITPACK64_ACCESS_INCLUDED
#include "types.h"

/* For reading direct-coded bitstreams */
extern UINT4
Bitpack64_access (UINT4 position, UINT4 *ptrs, UINT4 *comp);

#endif
