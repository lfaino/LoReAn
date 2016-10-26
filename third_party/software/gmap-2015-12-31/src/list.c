static char rcsid[] = "$Id: list.c 166641 2015-05-29 21:13:04Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "list.h"
#include "listdef.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mem.h"

#define T List_T

T
List_push (T list, void *x) {
  T new = (T) MALLOC(sizeof(*new));
  
  new->first = x;
  new->rest = list;
  return new;
}

T
List_push_keep (T list, void *x) {
  T new = (T) MALLOC_KEEP(sizeof(*new));
  
  new->first = x;
  new->rest = list;
  return new;
}

T
List_push_out (T list, void *x) {
  T new = (T) MALLOC_OUT(sizeof(*new));
  
  new->first = x;
  new->rest = list;
  return new;
}

T
List_pop (T list, void **x) {
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

T
List_pop_out (T list, void **x) {
  T head;

  if (list) {
    head = list->rest;
    *x = list->first;
    FREE_OUT(list);
    return head;
  } else {
    return list;
  }
}

void *
List_head (T list) {
  return list->first;
}

T
List_next (T list) {
  if (list) {
    return list->rest;
  } else {
    return NULL;
  }
}

void
List_head_set (T this, void *x) {
  this->first = x;
  return;
}

void
List_tail_set (T this, T rest) {
  this->rest = rest;
  return;
}

void
List_free (T *list) {
  T prev;

  while ((prev = *list) != NULL) {
    *list = prev->rest;
    FREE(prev);
  }

  return;
}

void
List_free_keep (T *list) {
  T prev;

  while ((prev = *list) != NULL) {
    *list = (*list)->rest;
    FREE_KEEP(prev);
  }
}

void
List_free_out (T *list) {
  T prev;

  while ((prev = *list) != NULL) {
    *list = (*list)->rest;
    FREE_OUT(prev);
  }
}

T
List_reverse (T list) {
  T head = NULL, next;

  for ( ; list; list = next) {
    next = list->rest;
    list->rest = head;
    head = list;
  }
  return head;
}

int
List_length (T list) {
  int n;
  
  for (n = 0; list; list = list->rest) {
    n++;
  }
  return n;
}

T
List_truncate (T list, int n) {
  T head = list;

  while (--n > 0) {
    list = list->rest;
  }
  if (list) {
    list->rest = (T) NULL;
  }
  return head;
}

void **
List_to_array (T list, void *end) {
  void **array;
  int i, n = List_length(list);

#if 0
  if (n == 0) {
    return (void *) NULL;
  } else {
#endif
    array = (void **) CALLOC((n+1),sizeof(*array));
    for (i = 0; i < n; i++) {
      array[i] = list->first;
      list = list->rest;
    }
    array[i] = end;
    return array;
#if 0
  }
#endif
}

void
List_fill_array (void **array, T list) {
  int i = 0;

  while (list) {
    array[i++] = list->first;
    list = list->rest;
  }
  return;
}

void
List_fill_array_and_free (void **array, T *list) {
  T prev;
  int i = 0;

  while ((prev = *list) != NULL) {
    array[i++] = prev->first;
    *list = prev->rest;
    FREE(prev);
  }

  return;
}

T
List_fill_array_with_handle (struct T *new, void **array, int nelts) {
  T list = NULL;
  int i;

  for (i = nelts; i > 0; i--) {
    new[i].first = array[i-1];
    new[i].rest = list;
    list = &(new[i]);
  }

  /* Add initial list element as a handle */
  new[0].first = (void *) NULL;
  new[0].rest = list;

  return &(new[0]);
}
    

void **
List_to_array_out (T list, void *end) {
  void **array;
  int i, n = List_length(list);

#if 0
  if (n == 0) {
    return (void *) NULL;
  } else {
#endif
    array = (void **) CALLOC_OUT((n+1),sizeof(*array));
    for (i = 0; i < n; i++) {
      array[i] = list->first;
      list = list->rest;
    }
    array[i] = end;
    return array;
#if 0
  }
#endif
}

void **
List_to_array_n (int *n, T list) {
  void **array;
  int i;

  *n = List_length(list);
  if (*n == 0) {
    return NULL;
  } else {
    array = (void **) CALLOC(*n,sizeof(*array));
    for (i = 0; i < *n; i++) {
      array[i] = list->first;
      list = list->rest;
    }
    return array;
  }
}

T
List_copy (T list) {
  T head, *p = &head;

  for ( ; list; list = list->rest) {
    *p = (T) MALLOC(sizeof(**p));
    (*p)->first = list->first;
    p = &(*p)->rest;
  }
  *p = NULL;
  return head;
}

void
List_dump (T list) {
  while (list) {
    printf("%p\n",list);
    list = list->rest;
  }
  return;
}

T
List_append (T list, T tail) {
  T *p = &list;

  while (*p) {
    p = &(*p)->rest;
  }
  *p = tail;
  return list;
}

void *
List_last_value (T this) {
  T last = NULL, r;

  for (r = this; r != NULL; r = r->rest) {
    last = r;
  }
  return last->first;
}

T
List_last_item (T this) {
  T last = NULL, r;

  for (r = this; r != NULL; r = r->rest) {
    last = r;
  }
  return last;
}

void *
List_index (T this, int index) {
  while (index-- > 0) {
    this = this->rest;
  }
  return this->first;
}


#if 0
/* Old definition, which doesn't make sense */
void
List_insert (T *listptr, void *x) {
  T new = (T) MALLOC(sizeof(*new));
  
  new->first = x;
  new->rest = *listptr;
  *listptr = new;

  return;
}
#else
T
List_insert (T p, void *x) {
  T new = (T) MALLOC(sizeof(*new));
  
  new->first = x;
  new->rest = p->rest;
  p->rest = new;
  return new;
}
#endif

void
List_reinsert (T *listptr, T cell) {
  cell->rest = *listptr;
  *listptr = cell;

  return;
}

T
List_transfer_one (T dest, T *source) {
  T next;

  next = (*source)->rest;
  (*source)->rest = dest;
  dest = *source;
  *source = next;
  return dest;
}

T
List_push_existing (T dest, T source) {
  source->rest = dest;
  return source;
}

T
List_from_string (char *string) {
  T this = NULL;
  char *p = string, *scout = string, *substring;
  int substringlen;

  while (*++scout != '\0') {
    if (*scout == ',') {
      substringlen = (scout-p)/sizeof(char);
      substring = (char *) CALLOC(substringlen+1,sizeof(char));
      strncpy(substring,p,substringlen);
      this = List_push(this,substring);
      scout++;
      p = scout;
    }
  }

  substringlen = (scout-p)/sizeof(char);
  substring = (char *) CALLOC(substringlen+1,sizeof(char));
  strncpy(substring,p,substringlen);
  this = List_push(this,substring);

  return List_reverse(this);
}
