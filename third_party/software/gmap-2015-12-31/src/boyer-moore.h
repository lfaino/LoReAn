/* $Id: boyer-moore.h 157221 2015-01-22 18:38:57Z twu $ */
#ifndef BOYER_MOORE_INCLUDED
#define BOYER_MOORE_INCLUDED

#include "intlist.h"
#include "genomicpos.h"

extern Intlist_T
BoyerMoore (char *query, int querylen, char *text, int textlen);
extern Intlist_T
BoyerMoore_nt (char *query, int querylen, int textoffset, int textlen,
	       Univcoord_T chroffset, Univcoord_T chrhigh, bool watsonp);
extern void
BoyerMoore_bad_char_shift (int *bad_char_shift, char *query, int querylen);
extern int
BoyerMoore_maxprefix (char *query, int querylen, char *text, int textlen,
		      int *bad_char_shift);

#endif

