/* $Id: list.h 166641 2015-05-29 21:13:04Z twu $ */
#ifndef LIST_INCLUDED
#define LIST_INCLUDED

#define T List_T
typedef struct T *T;

extern T List_push (T list, void *x);
extern T List_push_keep (T list, void *x);
extern T List_push_out (T list, void *x);
extern T List_pop (T list, void **x);
extern T List_pop_out (T list, void **x);
extern void *List_head (T list);
extern T List_next (T list);
extern void List_head_set (T list, void *x);
extern void List_tail_set (T this, T rest);
extern void List_free (T *list);
extern void List_free_keep (T *list);
extern void List_free_out (T *list);
extern T List_reverse (T list);
extern int List_length (T list);
extern T
List_truncate (T list, int n);
extern void **List_to_array (T list, void *end);
extern void List_fill_array (void **array, T list);
extern void List_fill_array_and_free (void **array, T *list);
extern T List_fill_array_with_handle (struct T *new, void **array, int nelts);
extern void **List_to_array_out (T list, void *end);
extern void **List_to_array_n (int *n, T list);
extern T List_copy (T list);
extern void
List_dump (T list);
extern T List_append (T list, T tail);
extern void *
List_last_value (T this);
extern T
List_last_item (T this);
extern void *
List_index (T this, int index);
extern T
List_insert (T this, void *x);
extern void
List_reinsert (T *listptr, T cell);
extern T
List_transfer_one (T dest, T *source);
extern T
List_push_existing (T dest, T source);
extern T
List_from_string (char *string);

#undef T
#endif
