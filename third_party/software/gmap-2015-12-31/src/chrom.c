static char rcsid[] = "$Id: chrom.c 138522 2014-06-09 17:08:44Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "chrom.h"
#include <stdio.h>
#include <stdlib.h>		/* For atoi */
#include <string.h>
#include "mem.h"
#include "interval.h"


#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif


/* This is the order for chrom sort */
typedef enum {PURE_NUMERIC, SEX_CHROMOSOME, MITOCHONDRIAL, NUMERIC_ALPHA, PURE_ALPHA} Chromtype_T;

#ifdef DEBUG
static char *
Chromtype_string (Chromtype_T chromtype) {
  switch (chromtype) {
  case PURE_NUMERIC: return "numeric";
  case SEX_CHROMOSOME: return "sex";
  case MITOCHONDRIAL: return "mitochondrial";
  case NUMERIC_ALPHA: return "numeric_alpha";
  case PURE_ALPHA: return "alpha";
  default: abort();
  }
}
#endif


#define T Chrom_T
struct T {
  Univcoord_T order;            /* order used for sorting (can be universal coordinate) */
  bool numericp;
  char *string;			/* The original string */
  unsigned int num;		/* The initial numeric part; valid only if numericp == true */
  char *alpha;			/* The alphabetic part, possibly after the number; valid only if numericp == true */
  Chromtype_T chromtype;
  bool circularp;
};

void
Chrom_free (T *old) {
  if ((*old)->numericp == true) {
    FREE((*old)->alpha);
  }
  FREE((*old)->string);
  FREE(*old);
  return;
}


char *
Chrom_string (T this) {
  return this->string;
}

bool
Chrom_circularp (T this) {
  return this->circularp;
}


/* Largest number for an unsigned int is 4294967295, which is 10
   digits.  However, the organism with the most chromosomes is the
   Indian fern, with 1260.  Hence, more than 4 digits would suggest
   a non-chromosomal string.  Also, if first digit is '0', treat as
   a string. */

T
Chrom_from_string (char *string, char *mitochondrial_string, Univcoord_T order, bool circularp) {
  T new = (T) MALLOC(sizeof(*new));
  int ndigits = 0;
  char *p;
  bool mitochondrial_p = false, sex_p = false;

  debug(printf("string is %s\n",string));

  new->order = order;
  new->circularp = circularp;

  new->string = (char *) CALLOC(strlen(string)+1,sizeof(char));
  strcpy(new->string,string);

  if (mitochondrial_string != NULL && !strcmp(string,mitochondrial_string)) {
    mitochondrial_p = true;
  }

  if (!strncmp(string,"chr",3)) {
    /* Ignore leading chr for sorting purposes */
    string += 3;
    debug(printf("  => chop chr to yield %s\n",string));
  }

  if (!strcmp(string,"X")) {
    sex_p = true;
  } else if (!strcmp(string,"Y")) {
    sex_p = true;
  } else if (!strcmp(string,"M")) {
    mitochondrial_p = true;
  } else if (!strcmp(string,"MT")) {
    mitochondrial_p = true;
  } else if (mitochondrial_string != NULL && !strcmp(string,mitochondrial_string)) {
    mitochondrial_p = true;
  }

  p = string;
  while (p != '\0' && *p >= '0' && *p <= '9') {
    ndigits++;
    p++;
  }
  debug(printf("  => saw %d digits\n",ndigits));


  if (ndigits > 0 && ndigits <= 4 && string[0] != '0') {
    new->numericp = true;
    new->num = atoi(string);
    new->alpha = (char *) CALLOC(strlen(p)+1,sizeof(char));
    strcpy(new->alpha,p);

    if (mitochondrial_p == true) {
      new->chromtype = MITOCHONDRIAL;
    } else if (strlen(new->alpha) == 0) {
      new->chromtype = PURE_NUMERIC;
    } else {
      new->chromtype = NUMERIC_ALPHA;
    }

    debug(printf("  => numeric with value %d and then alpha %s, type %s\n",
		 new->num,new->alpha,Chromtype_string(new->chromtype)));

  } else {
    new->numericp = false;
    new->num = 0;
    new->alpha = (char *) NULL;

    if (mitochondrial_p == true) {
      new->chromtype = MITOCHONDRIAL;
    } else if (sex_p == true) {
      new->chromtype = SEX_CHROMOSOME;
    } else {
      new->chromtype = PURE_ALPHA;
    }

    debug(printf("  => alphabetical %s, type %s\n",
		 new->string,Chromtype_string(new->chromtype)));
  }

  return new;
}


int
Chrom_cmp_alpha (T a, T b) {

  return strcmp(a->string,b->string);
}


int
Chrom_cmp_numeric_alpha (T a, T b) {

  if (a->numericp == true && b->numericp == false) {
    /* 1 and X */
    return -1;
  } else if (a->numericp == false && b->numericp == true) {
    /* X and 1 */
    return +1;
  } else if (a->numericp == true && b->numericp == true) {
    if (a->num < b->num) {
      /* 1 and 2U */
      return -1;
    } else if (a->num > b->num) {
      /* 2U and 1 */
      return +1;
    } else {
      return strcmp(a->alpha,b->alpha);
    }
  } else {
    return strcmp(a->string,b->string);
  }
}


int
Chrom_cmp_chrom (T a, T b) {

  debug(printf("Comparing %s and %s => ",a->string,b->string));

  if (a->chromtype < b->chromtype) {
    debug(printf("chromtype %d < %d => -1\n",a->chromtype,b->chromtype));
    return -1;
  } else if (b->chromtype < a->chromtype) {
    debug(printf("chromtype %d > %d => +1\n",a->chromtype,b->chromtype));
    return +1;
  } else if (a->numericp == true && b->numericp == true) {
    if (a->num < b->num) {
      /* 1 and 2U */
      debug(printf("numeric %d < %d => -1\n",a->num,b->num));
      return -1;
    } else if (a->num > b->num) {
      /* 2U and 1 */
      debug(printf("numeric %d > %d => +1\n",a->num,b->num));
      return +1;
    } else {
      debug(printf("numeric_alpha %s cmp %s => %d\n",
		   a->alpha,b->alpha,strcmp(a->alpha,b->alpha)));
      return strcmp(a->alpha,b->alpha);
    }
  } else {
      debug(printf("alpha %s cmp %s => %d\n",
		   a->string,b->string,strcmp(a->string,b->string)));
    return strcmp(a->string,b->string);
  }
}


int
Chrom_cmp_order (T a, T b) {

  if (a->order < b->order) {
    return -1;
  } else if (b->order < a->order) {
    return +1;
  } else {
    return 0;
  }
}


/* For use by qsorting an array */
int
Chrom_compare_order (const void *x, const void *y) {
  T a = * (T *) x;
  T b = * (T *) y;

  if (a->order < b->order) {
    return -1;
  } else if (b->order < a->order) {
    return +1;
  } else {
    return 0;
  }
}

/* For use by qsorting an array */
int
Chrom_compare_alpha (const void *x, const void *y) {
  T a = * (T *) x;
  T b = * (T *) y;

  return strcmp(a->string,b->string);
}

/* For use by qsorting an array */
int
Chrom_compare_numeric_alpha (const void *x, const void *y) {
  T a = * (T *) x;
  T b = * (T *) y;

  return Chrom_cmp_numeric_alpha(a,b);
}

/* For use by qsorting an array */
int
Chrom_compare_chrom (const void *x, const void *y) {
  T a = * (T *) x;
  T b = * (T *) y;

  return Chrom_cmp_chrom(a,b);
}


/* For use in table comparisons */
int
Chrom_compare_table (const void *x, const void *y) {
  T a = (T) x;
  T b = (T) y;

  return Chrom_cmp_chrom(a,b);
}

/* This is the X31 hash function */
static unsigned int
string_hash (char *a) {
  unsigned int h = 0U;
  char *p;
  
  for (p = a; *p != '\0'; p++) {
    h = (h << 5) - h + *p;
  }
  return h;
}

unsigned int
Chrom_hash_table (const void *key) {
  T this = (T) key;

  return string_hash(this->string);
}
