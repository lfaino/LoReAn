
AC_DEFUN([ACX_MADVISE_FLAGS], [
AC_LANG_SAVE
AC_LANG(C)

AC_MSG_CHECKING(for MADV_DONTNEED in madvise)
AC_COMPILE_IFELSE(
  [AC_LANG_PROGRAM([[#include <sys/types.h>
#include <sys/mman.h>]],
                   [[int flags = MADV_DONTNEED;]])],
  [AC_MSG_RESULT(yes)
   AC_DEFINE([HAVE_MADVISE_MADV_DONTNEED],[1],[Define to 1 if MADV_DONTNEED available for madvise.])],
  [AC_MSG_RESULT(no)])

AC_MSG_CHECKING(for MADV_WILLNEED in madvise)
AC_COMPILE_IFELSE(
  [AC_LANG_PROGRAM([[#include <sys/types.h>
#include <sys/mman.h>]],
                   [[int flags = MADV_WILLNEED;]])],
  [AC_MSG_RESULT(yes)
   AC_DEFINE([HAVE_MADVISE_MADV_WILLNEED],[1],[Define to 1 if MADV_WILLNEED available for madvise.])],
  [AC_MSG_RESULT(no)])

AC_MSG_CHECKING(for MADV_RANDOM in madvise)
AC_COMPILE_IFELSE(
  [AC_LANG_PROGRAM([[#include <sys/types.h>
#include <sys/mman.h>]],
                   [[int flags = MADV_RANDOM;]])],
  [AC_MSG_RESULT(yes)
   AC_DEFINE([HAVE_MADVISE_MADV_RANDOM],[1],[Define to 1 if MADV_RANDOM available for madvise])],
  [AC_MSG_RESULT(no)])

AC_MSG_CHECKING(for MADV_SEQUENTIAL in madvise)
AC_COMPILE_IFELSE(
  [AC_LANG_PROGRAM([[#include <sys/types.h>
#include <sys/mman.h>]],
                   [[int flags = MADV_SEQUENTIAL;]])],
  [AC_MSG_RESULT(yes)
   AC_DEFINE([HAVE_MADVISE_MADV_SEQUENTIAL],[1],[Define to 1 if MADV_SEQUENTIAL available for madvise])],
  [AC_MSG_RESULT(no)])

AC_LANG_RESTORE
])


