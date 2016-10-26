static char rcsid[] = "$Id: tableuint8.c 153955 2014-11-24 17:54:45Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "tableuint8.h"

#include <stdio.h>
#include <limits.h>
#include <stddef.h>
#include <stdlib.h>		/* For qsort */
#include "mem.h"
#include "assert.h"
#include "types.h"

#define T Tableuint8_T
struct T {
  int size;
  int (*cmp)(const void *x, const void *y);
  unsigned int (*hash)(const void *key);
  int length;
  unsigned int timestamp;
  struct binding {
    struct binding *link;
    const void *key;
    UINT8 value;
    unsigned int timeindex;
  } **buckets;
};


static int 
cmpatom (const void *x, const void *y) {
  return x != y;
}

static unsigned int
hashatom (const void *key) {
  return (unsigned long)key>>2;
}

T 
Tableuint8_new (int hint,
	       int (*cmp)(const void *x, const void *y),
	       unsigned int hash(const void *key)) {
  T table;
  int i;
  static int primes[] = { 509, 509, 1021, 2053, 4093,
			  8191, 16381, 32771, 65521, INT_MAX };

  assert(hint >= 0);
  for (i = 1; primes[i] < hint; i++) {
  }
  table = (T) MALLOC(sizeof(*table) +
		     primes[i-1]*sizeof(table->buckets[0]));
  table->size = primes[i-1];
  table->cmp  = cmp  ?  cmp : cmpatom;
  table->hash = hash ? hash : hashatom;
  table->buckets = (struct binding **)(table + 1);
  for (i = 0; i < table->size; i++) {
    table->buckets[i] = NULL;
  }
  table->length = 0;
  table->timestamp = 0;
  return table;
}

UINT8
Tableuint8_get (T table, const void *key) {
  int i;
  struct binding *p;

  assert(table);
  /* assert(key); -- Doesn't hold for atomic 0 */
  i = (*table->hash)(key)%table->size;
  /* fprintf(stderr,"Doing Tableuint8_get on %p at bucket %d\n",(char *) key, i); */
  for (p = table->buckets[i]; p; p = p->link) {
    /* fprintf(stderr,"  Comparing keys %p and %p\n",key,p->key); */
    if ((*table->cmp)(key, p->key) == 0) {
      /* fprintf(stderr,"Success: keys are identical\n"); */
      break;
    }
  }

#if 0
  if (p == NULL) {
    fprintf(stderr,"p is NULL\n");
  } else {
    fprintf(stderr,"Found p with value %llu\n",p->value);
  }
#endif

  return p ? p->value : 0;
}

UINT8
Tableuint8_put (T table, const void *key, UINT8 value) {
  int i;
  struct binding *p;
  UINT8 prev;

  /* fprintf(stderr,"Doing Tableuint8_put of key = %p, value %llu\n",key,value); */

  assert(table);
  /* assert(key); -- Doesn't hold for atomic 0 */
  i = (*table->hash)(key)%table->size;
  for (p = table->buckets[i]; p; p = p->link) {
    if ((*table->cmp)(key, p->key) == 0) {
      break;
    }
  }
  if (p == NULL) {
    /* fprintf(stderr,"New entry\n"); */
    NEW(p);
    p->key = key;
    p->link = table->buckets[i];
    table->buckets[i] = p;
    table->length++;
    prev = 0;
  } else {
    /* fprintf(stderr,"Existing entry\n"); */
    prev = p->value;
  }
  p->value = value;
  p->timeindex = table->timestamp;
  table->timestamp++;
  return prev;
}

int 
Tableuint8_length (T table) {
  assert(table);
  return table->length;
}

void 
Tableuint8_map (T table,
	   void (*apply)(const void *key, UINT8 *value, void *cl),
	   void *cl) {
  int i;
  unsigned int stamp;
  struct binding *p;

  assert(table);
  assert(apply);
  stamp = table->timestamp;
  for (i = 0; i < table->size; i++)
    for (p = table->buckets[i]; p; p = p->link) {
      apply(p->key, &p->value, cl);
      assert(table->timestamp == stamp);
    }
}

UINT8
Tableuint8_remove (T table, const void *key) {
  int i;
  struct binding **pp;

  assert(table);
  /* assert(key); -- Doesn't hold for atomic 0 */
  table->timestamp++;
  i = (*table->hash)(key)%table->size;
  for (pp = &table->buckets[i]; *pp; pp = &(*pp)->link) {
    if ((*table->cmp)(key, (*pp)->key) == 0) {
      struct binding *p = *pp;
      UINT8 value = p->value;
      *pp = p->link;
      FREE(p);
      table->length--;
      return value;
    }
  }
  return 0;
}

void **
Tableuint8_keys (T table, void *end) {
  void **keyarray;
  int i, j = 0;
  struct binding *p;

  assert(table);
  keyarray = (void **) CALLOC(table->length+1,sizeof(void *));
  for (i = 0; i < table->size; i++) {
    for (p = table->buckets[i]; p; p = p->link) {
      keyarray[j++] = (void *) p->key;
    }
  }
  keyarray[j] = end;
  return keyarray;
}


static int
timeindex_cmp (const void *x, const void *y) {
  struct binding *a = * (struct binding **) x;
  struct binding *b = * (struct binding **) y;

  if (a->timeindex < b->timeindex) {
    return -1;
  } else if (a->timeindex > b->timeindex) {
    return +1;
  } else {
    return 0;
  }
}


void **
Tableuint8_keys_by_timeindex (T table, void *end) {
  void **keyarray;
  int i, j = 0;
  struct binding **buckets, *p;

  assert(table);
  buckets = (struct binding **) CALLOC(table->length,sizeof(struct binding *));
  for (i = 0; i < table->size; i++) {
    for (p = table->buckets[i]; p; p = p->link) {
      buckets[j++] = p;
    }
  }
  qsort(buckets,table->length,sizeof(struct binding *),timeindex_cmp);

  keyarray = (void **) CALLOC(table->length+1,sizeof(void *));
  for (j = 0; j < table->length; j++) {
    p = buckets[j];
    keyarray[j] = (void *) p->key;
  }
  keyarray[j] = end;

  return keyarray;
}

UINT8 *
Tableuint8_values (T table, int end) {
  UINT8 *valuearray;
  int i, j = 0;
  struct binding *p;

  assert(table);
  valuearray = (UINT8 *) CALLOC(table->length+1,sizeof(UINT8));
  for (i = 0; i < table->size; i++) {
    for (p = table->buckets[i]; p; p = p->link) {
      valuearray[j++] = p->value;
    }
  }
  valuearray[j] = end;
  return valuearray;
}

void 
Tableuint8_free (T *table) {
  assert(table && *table);
  if ((*table)->length > 0) {
    int i;
    struct binding *p, *q;
    for (i = 0; i < (*table)->size; i++) {
      for (p = (*table)->buckets[i]; p; p = q) {
	q = p->link;
	FREE(p);
      }
    }
  }
  FREE(*table);
}
