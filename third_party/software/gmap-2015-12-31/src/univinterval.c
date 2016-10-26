static char rcsid[] = "$Id: univinterval.c 153955 2014-11-24 17:54:45Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "univinterval.h"
#include <stdio.h>
#include <stdlib.h>		/* For qsort */
#include "mem.h"


#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif


#define T Univinterval_T

T
Univinterval_new (Univcoord_T low, Univcoord_T high, int type) {
  T new = (T) MALLOC(sizeof(*new));

  if (low < high) {
    new->low = low;
    new->high = high;
    new->sign = +1;
  } else if (low > high) {
    new->low = high;
    new->high = low;
    new->sign = -1;
  } else {
    new->low = low;
    new->high = high;
    new->sign = 0;
  }
  new->type = type;
  return new;
}

T
Univinterval_copy (T old) {
  T new = (T) MALLOC(sizeof(*new));

  new->low = old->low;
  new->high = old->high;
  new->sign = old->sign;
  new->type = old->type;
  return new;
}
  
void
Univinterval_free (T *old) {
  if (*old) {
    FREE(*old);
  }
  return;
}

void
Univinterval_table_free (void **x) {
  T *old = (T *) x;
  FREE(*old);
  return;
}


void
Univinterval_print (T this) {
  printf("%llu %llu %d",(unsigned long long) this->low,(unsigned long long) this->high,this->type);
  return;
}

Univcoord_T
Univinterval_low (T this) {
  return this->low;
}

Univcoord_T
Univinterval_high (T this) {
  return this->high;
}

void
Univinterval_store_length (T this, Chrpos_T length) {
  this->high = this->low - 1 + length;
  return;
}

int
Univinterval_sign (T this) {
  return this->sign;
}

Chrpos_T
Univinterval_length (T this) {
  return this->high - this->low + 1;
}

int
Univinterval_type (T this) {
  return this->type;
}

/* Have to subtract 1 because intervals array is zero-based */
Univcoord_T
Univinterval_array_low (struct T *intervals, int index) {
  return intervals[index-1].low;
}

/* Have to subtract 1 because intervals array is zero-based */
Univcoord_T
Univinterval_array_high (struct T *intervals, int index) {
  return intervals[index-1].high;
}


/* Have to subtract 1 because intervals array is zero-based */
bool
Univinterval_is_contained (Univcoord_T x, struct T *intervals, int index) {
  Univcoord_T low = intervals[index-1].low;
  Univcoord_T high = intervals[index-1].high;

  if (low <= x && x <= high) {
    return true;
  } else {
    return false;
  }
}

/* Have to subtract 1 because intervals array is zero-based */
bool
Univinterval_overlap_p (Univcoord_T x, Univcoord_T y, struct T *intervals, int index) {
  Univcoord_T low = intervals[index-1].low;
  Univcoord_T high = intervals[index-1].high;

  if (x <= high && y >= low) {
    return true;
  } else {
    return false;
  }
}




/************************************************************************/
/* These sorting procedures are accessed only by iit-write.c            */
/************************************************************************/

static struct T *current_intervals;

/* Have to subtract 1 because intervals array is zero-based */
static int 
sigma_compare (const void *i, const void *j) {
  int x = * (int *) i;
  int y = * (int *) j;

  Univcoord_T a = current_intervals[x-1].low;
  Univcoord_T b = current_intervals[y-1].low;

  if (a < b) {
    return -1;
  } else if (a > b) {
    return 1;
  } else {
    return 0;
  }
}

/* Have to subtract 1 because intervals array is zero-based */
static int 
omega_compare (const void *i, const void *j) {
  int x = * (int *) i;
  int y = * (int *) j;

  Univcoord_T a = current_intervals[x-1].high;
  Univcoord_T b = current_intervals[y-1].high;

  if (a < b) {
    return -1;
  } else if (a > b) {
    return 1;
  } else {
    return 0;
  }
}



/* These routines sort table[i..j] in place.  Assume that
   current_intervals has been set. */
void
Univinterval_qsort_by_sigma (int *table, int i, int j, struct T *intervals) {
  current_intervals = intervals;
  qsort(&(table[i]), j - i + 1, sizeof(int), sigma_compare);
  return;
}

void
Univinterval_qsort_by_omega (int *table, int i, int j, struct T *intervals) {
  current_intervals = intervals;
  qsort(&(table[i]), j - i + 1, sizeof(int), omega_compare);
  return;
}


int
Univinterval_cmp (const void *a, const void *b) {
  T x = * (T *) a;
  T y = * (T *) b;

  debug(printf("Comparing %llu..%llu with %llu..%llu => ",
	       (unsigned long long) x->low,(unsigned long long) x->high,(unsigned long long) y->low,(unsigned long long) y->high));
  if (x->low < y->low) {
    debug(printf("-1\n"));
    return -1;
  } else if (x->low > y->low) {
    debug(printf("+1\n"));
    return +1;
  } else if (x->high < y->high) {
    debug(printf("-1\n"));
    return -1;
  } else if (x->high > y->high) {
    debug(printf("+1\n"));
    return +1;
  } else if (x->type < y->type) {
    debug(printf("-1\n"));
    return -1;
  } else if (x->type > y->type) {
    debug(printf("-1\n"));
    return +1;
  } else {
    debug(printf("0\n"));
    return 0;
  }
}


int
Univinterval_cmp_low (const void *a, const void *b) {
  T x = * (T *) a;
  T y = * (T *) b;

  debug(printf("Comparing %llu..%llu with %llu..%llu => ",
	       (unsigned long long) x->low,(unsigned long long) x->high,(unsigned long long) y->low,(unsigned long long) y->high));
  if (x->low < y->low) {
    debug(printf("-1\n"));
    return -1;
  } else if (x->low > y->low) {
    debug(printf("+1\n"));
    return +1;
  } else if (x->type < y->type) {
    debug(printf("-1\n"));
    return -1;
  } else if (x->type > y->type) {
    debug(printf("-1\n"));
    return +1;
  } else {
    debug(printf("0\n"));
    return 0;
  }
}


int
Univinterval_cmp_high (const void *a, const void *b) {
  T x = * (T *) a;
  T y = * (T *) b;

  debug(printf("Comparing %llu..%llu with %llu..%llu => ",
	       (unsigned long long) x->low,(unsigned long long) x->high,(unsigned long long) y->low,(unsigned long long) y->high));
  if (x->high < y->high) {
    debug(printf("-1\n"));
    return -1;
  } else if (x->high > y->high) {
    debug(printf("+1\n"));
    return +1;
#if 0
    /* Not needed by Splicetrie_retrieve_via_introns */
  } else if (x->type < y->type) {
    debug(printf("-1\n"));
    return -1;
  } else if (x->type > y->type) {
    debug(printf("-1\n"));
    return +1;
#endif
  } else {
    debug(printf("0\n"));
    return 0;
  }
}


#if 0
int
Univinterval_table_cmp (const void *a, const void *b) {
  T x = (T) a;
  T y = (T) b;

  debug(printf("Comparing %u..%u with %u..%u => ",x->low,x->high,y->low,y->high));
  if (x->low < y->low) {
    debug(printf("-1\n"));
    return -1;
  } else if (x->low > y->low) {
    debug(printf("+1\n"));
    return +1;
  } else if (x->high < y->high) {
    debug(printf("-1\n"));
    return -1;
  } else if (x->high > y->high) {
    debug(printf("+1\n"));
    return +1;
  } else if (x->type < y->type) {
    debug(printf("-1\n"));
    return -1;
  } else if (x->type > y->type) {
    debug(printf("+1\n"));
    return +1;
  } else {
    debug(printf("0\n"));
    return 0;
  }
}
#endif

#if 0
unsigned int
Univinterval_table_hash (const void *a) {
  T x = (T) a;
  return (unsigned int) (x->low + x->high);
}
#endif

bool
Univinterval_equal (T x, T y) {
  debug(printf("Comparing %u..%u with %u..%u => ",x->low,x->high,y->low,y->high));
  if (x->low != y->low) {
    return false;
  } else if (x->high != y->high) {
    return false;
  } else if (x->type != y->type) {
    return false;
  } else {
    return true;
  }
}
