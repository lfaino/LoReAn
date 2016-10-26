/* $Id: cmet.h 157222 2015-01-22 18:40:00Z twu $ */
#ifndef CMET_INCLUDED
#define CMET_INCLUDED

#include "indexdb.h"		/* For Storedoligomer_T */

extern Storedoligomer_T
Cmet_reduce_ct (Storedoligomer_T oligo);
extern Storedoligomer_T
Cmet_reduce_ga (Storedoligomer_T oligo);

#endif

