/* $Id: iit-write.h 157221 2015-01-22 18:38:57Z twu $ */
#ifndef IIT_WRITE_INCLUDED
#define IIT_WRITE_INCLUDED

#include "bool.h"
#include "list.h"
#include "uintlist.h"
#include "table.h"
#include "iitdef.h"

#define T IIT_T

extern void
IIT_output_direct (char *iitfile, T this, int version);
extern void
IIT_write (char *iitfile, List_T divlist, List_T typelist, List_T fieldlist, Table_T intervaltable,
	   Table_T valuetable, Table_T labeltable, Table_T annottable, Sorttype_T divsort, int version,
	   bool label_pointers_8p, bool annot_pointers_8p);
extern T
IIT_create (List_T divlist, List_T typelist, List_T fieldlist, Table_T intervaltable,
	    Table_T labeltable, Table_T datatable, Sorttype_T divsort, int version);
extern T
IIT_new (List_T intervallist);
extern void
IIT_backfill_sequence (T this, int index, int offset, char *Buffer);

#undef T
#endif

