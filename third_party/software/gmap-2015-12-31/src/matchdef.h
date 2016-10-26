/* $Id: matchdef.h 157221 2015-01-22 18:38:57Z twu $ */
#ifndef MATCHDEF_INCLUDED
#define MATCHDEF_INCLUDED

#include "bool.h"
#include "chrnum.h"
#include "genomicpos.h"

#define T Match_T
struct T {
  Univcoord_T position;
  Chrnum_T chrnum;
  Chrpos_T chrpos;
  double weight;		/* equal to 1/nentries */
  bool has_weight_p;
  int querypos;
  int npairings;		/* number of matchpairs made with
				   other matches */
  bool forwardp;
  bool fivep;
};

#undef T
#endif

