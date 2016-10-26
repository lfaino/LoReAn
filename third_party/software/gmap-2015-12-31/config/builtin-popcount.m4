
AC_DEFUN([ACX_BUILTIN_POPCOUNT], [
AC_REQUIRE([AC_CANONICAL_HOST])
AC_LANG_SAVE
AC_LANG(C)

CFLAGS_ORIG=$CFLAGS
CFLAGS="$CFLAGS -mpopcnt"

AC_MSG_CHECKING(whether -mpopcnt compiler flag works)
AC_RUN_IFELSE(
  [AC_LANG_PROGRAM([[#include <stdio.h>
#include <stdlib.h>]],
                   [[unsigned int x = rand();
printf("%08X ",x);
printf("clz=%d ",__builtin_clz(x));
printf("ctz=%d ",__builtin_ctz(x));
printf("popcount=%d ",__builtin_popcount(x));
]])],
  [AC_MSG_RESULT(yes)
   acx_mpopcnt_ok=yes],
  [AC_MSG_RESULT(no)])

if test "x$acx_mpopcnt_ok" != xyes; then
  CFLAGS=$CFLAGS_ORIG
fi

# Test for __builtin functions with or without the -mpopcnt compiler flag
AC_MSG_CHECKING(for __builtin_popcount)
AC_RUN_IFELSE(
  [AC_LANG_PROGRAM([[]],
                   [[return (__builtin_popcount(0xffffffffu) == 32) ? 0 : 9;]])],
  [AC_MSG_RESULT(yes)
   AC_DEFINE([HAVE_BUILTIN_POPCOUNT],[1],[Define to 1 if __builtin_popcount works.])],
  [AC_MSG_RESULT(no)])

AC_MSG_CHECKING(for __builtin_clz)
AC_RUN_IFELSE(
  [AC_LANG_PROGRAM([[]],
                   [[return (__builtin_clz(0x1u) == 31) ? 0 : 9;]])],
  [AC_MSG_RESULT(yes)
   AC_DEFINE([HAVE_BUILTIN_CLZ],[1],[Define to 1 if __builtin_clz works.])],
  [AC_MSG_RESULT(no)])

AC_MSG_CHECKING(for __builtin_ctz)
AC_RUN_IFELSE(
  [AC_LANG_PROGRAM([[]],
                   [[return (__builtin_ctz(0x80000000u) == 31) ? 0 : 9;]])],
  [AC_MSG_RESULT(yes)
   AC_DEFINE([HAVE_BUILTIN_CTZ],[1],[Define to 1 if __builtin_ctz works.])],
  [AC_MSG_RESULT(no)])

CFLAGS=$CFLAGS_ORIG

AC_LANG_RESTORE
])


