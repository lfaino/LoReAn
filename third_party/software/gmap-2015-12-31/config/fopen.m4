

define([_ACX_FUNC_FOPEN], [
AC_CACHE_CHECK([whether fopen accepts "$1" mode],
               [cvx_cv_func_fopen_$1],
               [AC_TRY_RUN([
#include <stdio.h>
int
main () {
  FILE *fp = fopen("conftest.bin","w$1");
  fprintf(fp,"\n");
  fclose(fp);
  return 0;
}],
                  [cvx_cv_func_fopen_$1=yes],
                  [cvx_cv_func_fopen_$1=no],
                  [cvx_cv_func_fopen_$1=no])])
if test x$cvx_cv_func_fopen_$1 = xyes; then
  AC_DEFINE([$2], 1,
            [Define this if we can use the "$1" mode for fopen safely.])
fi[]dnl
])


AC_DEFUN([ACX_FUNC_FOPEN_BINARY], [_ACX_FUNC_FOPEN(b, USE_FOPEN_BINARY)])
AC_DEFUN([ACX_FUNC_FOPEN_TEXT], [_ACX_FUNC_FOPEN(t, USE_FOPEN_TEXT)])

