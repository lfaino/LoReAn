/* $Id: assert.h 40271 2011-05-28 02:29:18Z twu $ */
#ifndef ASSERT_INCLUDED
#define ASSERT_INCLUDED

#include "except.h"

#undef assert

/* #define CHECK_ASSERTIONS 1 */
#ifdef CHECK_ASSERTIONS
extern void assert (int e);
#define assert(e) ((void) ((e) || (RAISE(Assert_Failed),0)))
#else
#define assert(e) ((void) 0)
#endif

#endif
