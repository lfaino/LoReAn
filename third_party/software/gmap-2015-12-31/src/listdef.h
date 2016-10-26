/* $Id: listdef.h 40271 2011-05-28 02:29:18Z twu $ */
#ifndef LISTDEF_INCLUDED
#define LISTDEF_INCLUDED

#define T List_T
struct T {
  void *first;
  struct T *rest;
};

#undef T
#endif

