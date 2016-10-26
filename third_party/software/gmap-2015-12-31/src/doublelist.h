/* $Id: doublelist.h 166641 2015-05-29 21:13:04Z twu $ */
#ifndef DOUBLELIST_INCLUDED
#define DOUBLELIST_INCLUDED

#define T Doublelist_T
typedef struct T *T;

extern T Doublelist_push (T list, double index);
extern T Doublelist_pop (T list, double *index);
extern double Doublelist_head (T list);
extern T Doublelist_next (T list);
extern void Doublelist_free (T *list);
extern T Doublelist_reverse (T list);
extern int Doublelist_length (T list);
extern double *
Doublelist_to_array (int *n, T list);
extern double *
Doublelist_to_array_out (int *n, T list);
extern void
Doublelist_fill_array (double *array, T list);
extern T Doublelist_from_string (char *string);
extern T
Doublelist_from_array (double *array, int n);
extern double
Doublelist_max (T this);
extern double
Doublelist_min (T this);
extern void
Doublelist_print (T this);

#undef T
#endif
