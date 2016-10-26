/* $Id: popcount.h 157225 2015-01-22 18:47:23Z twu $ */
#ifndef POPCOUNT_INCLUDED
#define POPCOUNT_INCLUDED
#ifdef HAVE_CONFIG_H
#include <config.h>		/* For HAVE_BUILTIN_CTZ, HAVE_BUILTIN_POPCOUNT, HAVE_BUILTIN_CLZ */
#endif

#ifndef HAVE_BUILTIN_CTZ
extern const int mod_37_bit_position[];
#endif

#ifndef HAVE_BUILTIN_POPCOUNT
extern const int count_bits[];
#endif

#ifndef HAVE_BUILTIN_CLZ
extern const int clz_table[];
#endif

#endif

