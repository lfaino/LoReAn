/* $Id: intlistdef.h 40271 2011-05-28 02:29:18Z twu $ */
#ifndef INTLISTDEF_INCLUDED
#define INTLISTDEF_INCLUDED

#define T Intlist_T
struct T {
  int first;
  struct T *rest;
};

#undef T
#endif
