/* $Id: maxent_hr.h 157221 2015-01-22 18:38:57Z twu $ */
#ifndef MAXENT_HR_INCLUDED
#define MAXENT_HR_INCLUDED

#include "genomicpos.h"
#include "types.h"

extern void
Maxent_hr_setup (Genomecomp_T *ref_blocks_in, Genomecomp_T *snp_blocks_in);

extern double
Maxent_hr_donor_prob (Univcoord_T splice_pos, Univcoord_T chroffset);

extern double
Maxent_hr_acceptor_prob (Univcoord_T splice_pos, Univcoord_T chroffset);

extern double
Maxent_hr_antidonor_prob (Univcoord_T splice_pos, Univcoord_T chroffset);

extern double
Maxent_hr_antiacceptor_prob (Univcoord_T splice_pos, Univcoord_T chroffset);

#endif

