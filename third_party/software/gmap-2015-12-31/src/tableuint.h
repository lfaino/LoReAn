/* $Id: tableuint.h 99737 2013-06-27 19:33:03Z twu $ */
#ifndef TABLEUINT_INCLUDED
#define TABLEUINT_INCLUDED

#define T Tableuint_T
typedef struct T *T;

extern T
Tableuint_new (int hint,
	      int (*cmp)(const void *x, const void *y),
	      unsigned int hash(const void *key));
extern void 
Tableuint_free (T *table);
extern int   
Tableuint_length (T table);
extern unsigned int
Tableuint_put (T table, const void *key, unsigned int value);
extern unsigned int
Tableuint_get (T table, const void *key);
extern unsigned int
Tableuint_remove (T table, const void *key);
extern void   
Tableuint_map (T table,
	      void (*apply)(const void *key, unsigned int *value, void *cl),
	      void *cl);
extern void **
Tableuint_keys (T table, void *end);
extern void **
Tableuint_keys_by_timeindex (T table, void *end);
extern unsigned int *
Tableuint_values (T table, int end);

#undef T
#endif
