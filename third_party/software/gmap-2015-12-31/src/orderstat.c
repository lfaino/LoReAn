static char rcsid[] = "$Id: orderstat.c 40271 2011-05-28 02:29:18Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#ifndef HAVE_MEMCPY
# define memcpy(d,s,n) bcopy((s),(d),(n))
#endif

#include "orderstat.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>		/* For memcpy */


#if 0
static void
dump_vector (double *vector, int n) {
  int i;

  for (i = 0; i < n; i++) {
    printf("%d ",(int) vector[i]);
  }
  printf("\n");
  return;
}
#endif


static double
quickselect_double_aux (double *set, int n, int k) {
  double x, elt;
  int lowcount = 0, eqcount = 1, j = 0, i;

  /*
  printf("n = %d, k = %d: ",n,k);
  dump_vector(set,n);
  */

  x = set[0];
  for (i = 1; i < n; i++) {
    if ((elt = set[i]) < x) {
      lowcount++;
#if 0
    } else if (elt == x) {
      eqcount++;
#else
    } else if (elt > x) {
      /* Do nothing */
    } else {
      eqcount++;
#endif
    }
  }

  if (k <= lowcount) {
    for (i = 1; i < n; i++) {
      if (set[i] < x) {
	set[j++] = set[i];
      }
    }
    return quickselect_double_aux(set,j,k);

  } else if (k > lowcount+eqcount) {
    for (i = 1; i < n; i++) {
      if (set[i] > x) {
	set[j++] = set[i];
      }
    }
    return quickselect_double_aux(set,j,k-lowcount-eqcount);

  } else {
    return x;
  }
}

static int
quickselect_int_aux (int *set, int n, int k) {
  int x, elt;
  int lowcount = 0, eqcount = 1, j = 0, i;

  /*
  printf("n = %d, k = %d: ",n,k);
  dump_vector(set,n);
  */

  x = set[0];
  for (i = 1; i < n; i++) {
    if ((elt = set[i]) < x) {
      lowcount++;
    } else if (elt == x) {
      eqcount++;
    }
  }

  if (k <= lowcount) {
    for (i = 1; i < n; i++) {
      if (set[i] < x) {
	set[j++] = set[i];
      }
    }
    return quickselect_int_aux(set,j,k);

  } else if (k > lowcount+eqcount) {
    for (i = 1; i < n; i++) {
      if (set[i] > x) {
	set[j++] = set[i];
      }
    }
    return quickselect_int_aux(set,j,k-lowcount-eqcount);

  } else {
    return x;
  }
}

static long int
quickselect_long_int_aux (long int *set, int n, int k) {
  long int x, elt;
  int lowcount = 0, eqcount = 1, j = 0, i;

  /*
  printf("n = %d, k = %d: ",n,k);
  dump_vector(set,n);
  */

  x = set[0];
  for (i = 1; i < n; i++) {
    if ((elt = set[i]) < x) {
      lowcount++;
    } else if (elt == x) {
      eqcount++;
    }
  }

  if (k <= lowcount) {
    for (i = 1; i < n; i++) {
      if (set[i] < x) {
	set[j++] = set[i];
      }
    }
    return quickselect_long_int_aux(set,j,k);

  } else if (k > lowcount+eqcount) {
    for (i = 1; i < n; i++) {
      if (set[i] > x) {
	set[j++] = set[i];
      }
    }
    return quickselect_long_int_aux(set,j,k-lowcount-eqcount);

  } else {
    return x;
  }
}

static double
quickselect_double (double *vector, int n, int k) {
  double result;
  double *temp;

  temp = (double *) calloc(n,sizeof(double));
  memcpy((void *) temp,vector,n*sizeof(double));

  result = quickselect_double_aux(temp,n,k);

  free(temp);
  return result;
}

static int
quickselect_int (int *vector, int n, int k) {
  int result;
  int *temp;

  temp = (int *) calloc(n,sizeof(int));
  memcpy((void *) temp,vector,n*sizeof(int));

  result = quickselect_int_aux(temp,n,k);

  free(temp);
  return result;
}

static long int
quickselect_long_int (long int *vector, int n, int k) {
  long int result;
  long int *temp;

  temp = (long int *) calloc(n,sizeof(long int));
  memcpy((void *) temp,vector,n*sizeof(long int));

  result = quickselect_long_int_aux(temp,n,k);

  free(temp);
  return result;
}


double
Orderstat_double_pct (double *set, int length, double pct) {
  int cutoff;

  cutoff = (int) (pct*length+1);
  if (cutoff > length) {
    cutoff = length;
  }
  return quickselect_double(set,length,cutoff);
}

double
Orderstat_double_pct_inplace (double *set, int length, double pct) {
  int cutoff;

  cutoff = (int) (pct*length+1);
  if (cutoff > length) {
    cutoff = length;
  }
  return quickselect_double_aux(set,length,cutoff);
}

int
Orderstat_int_pct (int *set, int length, double pct) {
  int cutoff;

  cutoff = (int) (pct*length+1);
  if (cutoff > length) {
    cutoff = length;
  }
  return quickselect_int(set,length,cutoff);
}
 
long int
Orderstat_long_int_pct (long int *set, int length, double pct) {
  int cutoff;

  cutoff = (int) (pct*length+1);
  if (cutoff > length) {
    cutoff = length;
  }
  return quickselect_long_int(set,length,cutoff);
}

int
Orderstat_int_pct_inplace (int *set, int length, double pct) {
  int cutoff;

  cutoff = (int) (pct*length+1);
  if (cutoff > length) {
    cutoff = length;
  }
  return quickselect_int_aux(set,length,cutoff);
}

  

#if 0
int
compare (const void *x, const void *y) {
  double a = * (double *) x;
  double b = * (double *) y;

  if (a < b) {
    return -1;
  } else if (a > b) {
    return 1;
  } else {
    return 0;
  }
}

int
main (int argc, char *argv[]) {
  int i, n;
  double elt, *vector, median, first, second, last;

  n = atoi(argv[1]);
  vector = (double *) calloc(n,sizeof(double));

  for (i = 0; i < n; i++) {
    elt = (double) rand();
    vector[i] = elt;
  }
  qsort(vector,n,sizeof(double),compare);
  dump_vector(vector,n);

  printf("Median true           = %d\n",(int) vector[n/2]);

  median = quickselect(vector,n,n/2+1); 
  printf("Median by quickselect = %d\n",(int) median);

  first = quickselect(vector,n,1);

  last = quickselect(vector,n,n);
  printf("First by quickselect = %d\n",(int) first);
  printf("Last by quickselect = %d\n",(int) last);

  return 0;
}

#endif



/*

static double
sample_select (double *set, int n, int k) {
  double x;
  int m, j, i, count = 0;

  m = rint(exp(2.0*log(n*log(n))/3.0));
  j = (k*m)/n;		
  x = sample_select(set,m,j);

  for (i = 0; i < n
  
*/  
