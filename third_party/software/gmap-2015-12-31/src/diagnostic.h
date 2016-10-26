/* $Id: diagnostic.h 40271 2011-05-28 02:29:18Z twu $ */
#ifndef DIAGNOSTIC_INCLUDED
#define DIAGNOSTIC_INCLUDED

#define T Diagnostic_T
typedef struct T *T;
struct T {
#ifndef PMAP
  double query_oligodepth;
#endif
  int query_badoligos;
  int query_repoligos;
  int query_trimoligos;
  int query_trim_start;
  int query_trim_end;

  double stage1_runtime;
  int firstpair_found_5;
  int firstpair_found_3;
  int stutter_searched_5;
  int stutter_searched_3;
  int stutter_nmatchpairs;
  int stutter_matches_5;
  int stutter_matches_3;
  int dangling_5;
  int dangling_3;
  int dangling_nmatchpairs;

  int sampling_rounds;
  int sampling_nskip;
  int ngregions;
};

extern void
Diagnostic_free (T *old);

extern T
Diagnostic_new ();

extern void
Diagnostic_print (T this);


#undef T
#endif

