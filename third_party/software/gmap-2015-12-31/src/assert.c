static char rcsid[] = "$Id: assert.c 40271 2011-05-28 02:29:18Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "assert.h"

const Except_T Assert_Failed = { "Assertion Failed" };

/*
void
(assert) (int e) {
  assert(e);
}
*/

