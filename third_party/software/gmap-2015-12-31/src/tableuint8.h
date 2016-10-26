/* $Id: tableuint8.h 157221 2015-01-22 18:38:57Z twu $ */
#ifndef TABLEUINT8_INCLUDED
#define TABLEUINT8_INCLUDED

#include "types.h"

#define T Tableuint8_T
typedef struct T *T;

extern T
Tableuint8_new (int hint,
		int (*cmp)(const void *x, const void *y),
		unsigned int hash(const void *key));
extern void 
Tableuint8_free (T *table);
extern int   
Tableuint8_length (T table);
extern UINT8
Tableuint8_put (T table, const void *key, UINT8 value);
extern UINT8
Tableuint8_get (T table, const void *key);
extern UINT8
Tableuint8_remove (T table, const void *key);
extern void   
Tableuint8_map (T table,
	      void (*apply)(const void *key, UINT8 *value, void *cl),
	      void *cl);
extern void **
Tableuint8_keys (T table, void *end);
extern void **
Tableuint8_keys_by_timeindex (T table, void *end);
extern UINT8 *
Tableuint8_values (T table, int end);

#undef T
#endif
