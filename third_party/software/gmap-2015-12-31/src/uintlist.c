static char rcsid[] = "$Id: uintlist.c 166641 2015-05-29 21:13:04Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "uintlist.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mem.h"

#define T Uintlist_T
struct T {
  UINT4 first;
  T rest;
};

T
Uintlist_push (T list, UINT4 x) {
  T new = (T) MALLOC(sizeof(*new));
  
  new->first = x;
  new->rest = list;
  return new;
}

T
Uintlist_pop (T list, UINT4 *x) {
  T head;

  if (list) {
    head = list->rest;
    *x = list->first;
    FREE(list);
    return head;
  } else {
    return list;
  }
}
  
UINT4
Uintlist_head (T list) {
  return list->first;
}

T
Uintlist_next (T list) {
  if (list) {
    return list->rest;
  } else {
    return NULL;
  }
}

void
Uintlist_head_set (T this, UINT4 x) {
  this->first = x;
  return;
}

void
Uintlist_free (T *list) {
  T prev;

  while ((prev = *list) != NULL) {
    *list = prev->rest;
    FREE(prev);
  }

  return;
}

T
Uintlist_reverse (T list) {
  T head = NULL, next;

  for ( ; list; list = next) {
    next = list->rest;
    list->rest = head;
    head = list;
  }
  return head;
}

int
Uintlist_length (T list) {
  int n;
  
  for (n = 0; list; list = list->rest) {
    n++;
  }
  return n;
}

UINT4
Uintlist_max (T list) {
  UINT4 m = 0;

  while (list) {
    if (list->first > m) {
      m = list->first;
    }
    list = list->rest;
  }

  return m;
}

UINT4
Uintlist_min (T list) {
  UINT4 m;

  if (list == NULL) {
    return 0;

  } else {
    m = list->first;
    list = list->rest;
    while (list) {
      if (list->first < m) {
	m = list->first;
      }
      list = list->rest;
    }

    return m;
  }
}

UINT4 *
Uintlist_to_array (int *n, T list) {
  UINT4 *array;
  int i;

  *n = Uintlist_length(list);
  array = (UINT4 *) CALLOC(*n,sizeof(UINT4));
  for (i = 0; i < *n; i++) {
    array[i] = list->first;
    list = list->rest;
  }
  return array;
}

void
Uintlist_fill_array (UINT4 *array, T list) {
  int i = 0;

  while (list) {
    array[i++] = list->first;
    list = list->rest;
  }

  return;
}

void
Uintlist_fill_array_and_free (UINT4 *array, T *list) {
  T prev;
  int i = 0;

  while ((prev = *list) != NULL) {
    array[i++] = prev->first;
    *list = prev->rest;
    FREE(prev);
  }

  return;
}

UINT4 *
Uintlist_to_array_out (int *n, T list) {
  UINT4 *array;
  int i;

  *n = Uintlist_length(list);
  if (*n == 0) {
    return NULL;
  } else {
    array = (UINT4 *) CALLOC_OUT(*n,sizeof(UINT4));
    for (i = 0; i < *n; i++) {
      array[i] = list->first;
      list = list->rest;
    }
    return array;
  }
}

T
Uintlist_from_array (UINT4 *array, int n) {
  T list = NULL, p;

  while (--n >= 0) {
    p = (T) MALLOC(sizeof(*p));
    p->first = array[n];
    p->rest = list;
    list = p;
  }

  return list;
}

T
Uintlist_copy (T list) {
  T head, *p = &head;

  for ( ; list; list = list->rest) {
    *p = (T) MALLOC(sizeof(**p));
    (*p)->first = list->first;
    p = &(*p)->rest;
  }
  *p = NULL;
  return head;
}

T
Uintlist_append (T list, T tail) {
  T *p = &list;

  while (*p) {
    p = &(*p)->rest;
  }
  *p = tail;
  return list;
}

UINT4
Uintlist_last_value (T this) {
  T last = NULL, r;

  for (r = this; r != NULL; r = r->rest) {
    last = r;
  }
  return last->first;
}

UINT4
Uintlist_index (T this, int index) {
  while (index-- > 0) {
    this = this->rest;
  }
  return this->first;
}


bool
Uintlist_find (T this, UINT4 value) {
  T r;

  for (r = this; r != NULL; r = r->rest) {
    if (r->first == value) {
      return true;
    }
  }
  return false;
}

char *
Uintlist_to_string (T this) {
  char *string, Buffer[256];
  T p;
  int n, i, strlength;

  if ((n = Uintlist_length(this)) == 0) {
    string = (char *) CALLOC(1,sizeof(char));
    string[0] = '\0';
  } else {
    strlength = 0;
    for (i = 0, p = this; i < n-1; i++, p = Uintlist_next(p)) {
      sprintf(Buffer,"%u,",Uintlist_head(p));
      strlength += strlen(Buffer);
    }
    sprintf(Buffer,"%u",Uintlist_head(p));
    strlength += strlen(Buffer);

    string = (char *) CALLOC(strlength + 1,sizeof(char));
    string[0] = '\0';
    for (i = 0, p = this; i < n-1; i++, p = Uintlist_next(p)) {
      sprintf(Buffer,"%u,",Uintlist_head(p));
      strcat(string,Buffer);
    }
    sprintf(Buffer,"%u",Uintlist_head(p));
    strcat(string,Buffer);
  }

  return string;
}
