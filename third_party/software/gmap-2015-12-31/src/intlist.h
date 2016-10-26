/* $Id: intlist.h 166641 2015-05-29 21:13:04Z twu $ */
#ifndef INTLIST_INCLUDED
#define INTLIST_INCLUDED

#include "bool.h"

#define T Intlist_T
typedef struct T *T;

extern T 
Intlist_push (T list, int x);
extern T 
Intlist_push_in (T list, int x);
extern T
Intlist_insert_second (T list, int x);
extern T 
Intlist_pop (T list, int *x);
extern void
Intlist_delete (T prev, T this);
extern int 
Intlist_head (T list);
extern T 
Intlist_next (T list);
extern void 
Intlist_head_set (T list, int x);
extern void 
Intlist_free (T *list);
extern void 
Intlist_free_in (T *list);
extern T 
Intlist_reverse (T list);
extern int 
Intlist_length (T list);
extern int
Intlist_max (T list);
extern int
Intlist_min (T list);
extern bool
Intlist_vary (T list);
extern bool
Intlist_exists_p (T list, int x);
extern int *
Intlist_to_array (int *n, T list);
extern void
Intlist_fill_array (int *array, T list);
extern void
Intlist_fill_array_and_free (int *array, T *list);
extern int *
Intlist_to_array_out (int *n, T list);
extern char *
Intlist_to_char_array (int *n, T list);
extern T
Intlist_from_array (int *array, int n);
extern T 
Intlist_copy (T list);
extern T 
Intlist_append (T list, T tail);
extern int 
Intlist_last_value (T this);
extern int 
Intlist_index (T this, int index);
extern T
Intlist_from_string (char *string);
extern char *
Intlist_to_string (T this);
extern int *
Intlist_array_ascending_by_key (int *n, T this, T key);
extern void
Intlist_array_dual_ascending_by_key (int *sorted, int *keyarray, int n, T this, T key);
extern T
Intlist_list_ascending_by_key (T this, T key);
extern T
Intlist_list_descending_by_key (T this, T key);
extern T
Intlist_sort_ascending (T this);
extern bool
Intlist_equal (T x, T y);

#undef T
#endif
