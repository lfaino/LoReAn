static char rcsid[] = "$Id: diagnostic.c 40271 2011-05-28 02:29:18Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "diagnostic.h"
#include <stdio.h>
#include <stdlib.h>
#include "mem.h"

#define T Diagnostic_T

void
Diagnostic_free (T *old) {
  FREE(*old);
  return;
}

T
Diagnostic_new () {
  T new = (T) MALLOC(sizeof(*new));

#ifndef PMAP
  new->query_oligodepth = 0.0;
#endif
  new->query_badoligos = 0;
  new->query_trimoligos = 0;
  new->query_repoligos = 0;
  new->query_trimoligos = 0;
  new->query_trim_start = 0;
  new->query_trim_end = 0;

  new->stage1_runtime = 0.0;
  new->firstpair_found_5 = 0;
  new->firstpair_found_3 = 0;
  new->stutter_searched_5 = 0;
  new->stutter_searched_3 = 0;
  new->stutter_nmatchpairs = 0;
  new->stutter_matches_5 = 0;
  new->stutter_matches_3 = 0;
  new->sampling_rounds = 0;
  new->sampling_nskip = 0;
  new->ngregions = 0;

  return new;
}

void
Diagnostic_print (T this) {
  if (this != NULL) {
#ifndef PMAP
    printf("Query check oligodepth: %f\n",this->query_oligodepth);
#endif
    printf("Query check badoligos: %d/%d\n",this->query_badoligos,this->query_trimoligos);
    printf("Query check repoligos: %d/%d\n",this->query_repoligos,this->query_trimoligos);
    printf("Query check trim bounds: %d..%d\n",this->query_trim_start,this->query_trim_end);

    printf("Stage 1 runtime: %.3f sec\n",this->stage1_runtime);
    printf("Stage 1 firstpair_found_5: %d\n",this->firstpair_found_5);
    printf("Stage 1 firstpair_found_3: %d\n",this->firstpair_found_3);
    printf("Stage 1 stutter_searched_5: %d\n",this->stutter_searched_5);
    printf("Stage 1 stutter_searched_3: %d\n",this->stutter_searched_3);
    printf("Stage 1 stutter_matchpairs: %d\n",this->stutter_nmatchpairs);
    printf("Stage 1 stutter_matches_5: %d\n",this->stutter_matches_5);
    printf("Stage 1 stutter_matches_3: %d\n",this->stutter_matches_3);
    printf("Stage 1 sampling rounds: %d\n",this->sampling_rounds);
    printf("Stage 1 sampling nskip: %d\n",this->sampling_nskip);
    printf("Stage 1 number gregions: %d\n",this->ngregions);
  }
    
  return;
}

