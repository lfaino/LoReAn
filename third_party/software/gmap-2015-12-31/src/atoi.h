/* $Id: atoi.h 157222 2015-01-22 18:40:00Z twu $ */
#ifndef ATOI_INCLUDED
#define ATOI_INCLUDED

#include "indexdb.h"		/* For Storedoligomer_T */

extern Storedoligomer_T
Atoi_reduce_ag (Storedoligomer_T oligo);
extern Storedoligomer_T
Atoi_reduce_tc (Storedoligomer_T oligo);

#endif

