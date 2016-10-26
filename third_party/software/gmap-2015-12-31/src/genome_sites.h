/* $Id: genome_sites.h 157221 2015-01-22 18:38:57Z twu $ */
#ifndef GENOME_SITES_INCLUDED
#define GENOME_SITES_INCLUDED

#include "bool.h"
#include "types.h"
#include "genomicpos.h"


extern void
Genome_sites_setup (Genomecomp_T *ref_blocks_in, Genomecomp_T *snp_blocks_in);

extern int
Genome_donor_positions (int *site_positions, int *site_knowni, int *knownpos, int *knowni,
			Univcoord_T left, int pos5, int pos3);

extern int
Genome_acceptor_positions (int *site_positions, int *site_knowni, int *knownpos, int *knowni,
			   Univcoord_T left, int pos5, int pos3);

extern int
Genome_antidonor_positions (int *site_positions, int *site_knowni, int *knownpos, int *knowni,
			    Univcoord_T left, int pos5, int pos3);

extern int
Genome_antiacceptor_positions (int *site_positions, int *site_knowni, int *knownpos, int *knowni,
			       Univcoord_T left, int pos5, int pos3);

extern bool
Genome_sense_canonicalp (Univcoord_T donor_rightbound, Univcoord_T donor_leftbound,
			 Univcoord_T acceptor_rightbound, Univcoord_T acceptor_leftbound,
			 Univcoord_T chroffset);

extern bool
Genome_antisense_canonicalp (Univcoord_T donor_rightbound, Univcoord_T donor_leftbound,
			     Univcoord_T acceptor_rightbound, Univcoord_T acceptor_leftbound,
			     Univcoord_T chroffset);


extern Univcoord_T
Genome_prev_donor_position (Univcoord_T pos, Univcoord_T prevpos, Univcoord_T chroffset);
extern Univcoord_T
Genome_prev_acceptor_position (Univcoord_T pos, Univcoord_T prevpos, Univcoord_T chroffset);
extern Univcoord_T
Genome_prev_antidonor_position (Univcoord_T pos, Univcoord_T prevpos, Univcoord_T chroffset);
extern Univcoord_T
Genome_prev_antiacceptor_position (Univcoord_T pos, Univcoord_T prevpos, Univcoord_T chroffset);

#endif

