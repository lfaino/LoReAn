
AC_DEFUN([ACX_MMAP_FLAGS], [
AC_LANG_SAVE
AC_LANG(C)

AC_MSG_CHECKING(for MAP_FILE in mmap)
AC_COMPILE_IFELSE(
  [AC_LANG_PROGRAM([[#include <sys/types.h>
#include <sys/mman.h>]],
                   [[int flags = MAP_FILE;]])],
  [AC_MSG_RESULT(yes)
   AC_DEFINE([HAVE_MMAP_MAP_FILE],[1],[Define to 1 if MAP_FILE available for mmap.])],
  [AC_MSG_RESULT(no)])

AC_MSG_CHECKING(for MAP_VARIABLE in mmap)
AC_COMPILE_IFELSE(
  [AC_LANG_PROGRAM([[#include <sys/types.h>
#include <sys/mman.h>]],
                   [[int flags = MAP_VARIABLE;]])],
  [AC_MSG_RESULT(yes)
   AC_DEFINE([HAVE_MMAP_MAP_VARIABLE],[1],[Define to 1 if MAP_VARIABLE available for mmap.])],
  [AC_MSG_RESULT(no)])

AC_MSG_CHECKING(for MAP_SHARED in mmap)
AC_COMPILE_IFELSE(
  [AC_LANG_PROGRAM([[#include <sys/types.h>
#include <sys/mman.h>]],
                   [[int flags = MAP_SHARED;]])],
  [AC_MSG_RESULT(yes)
   AC_DEFINE([HAVE_MMAP_MAP_SHARED],[1],[Define to 1 if MAP_SHARED available for mmap])],
  [AC_MSG_RESULT(no)])

AC_MSG_CHECKING(for MAP_PRIVATE in mmap)
AC_COMPILE_IFELSE(
  [AC_LANG_PROGRAM([[#include <sys/types.h>
#include <sys/mman.h>]],
                   [[int flags = MAP_PRIVATE;]])],
  [AC_MSG_RESULT(yes)
   AC_DEFINE([HAVE_MMAP_MAP_PRIVATE],[1],[Define to 1 if MAP_PRIVATE available for mmap])],
  [AC_MSG_RESULT(no)])

AC_MSG_CHECKING(for MAP_FAILED in mmap)
AC_COMPILE_IFELSE(
  [AC_LANG_PROGRAM([[#include <sys/types.h>
#include <sys/mman.h>]],
                   [[void *p = MAP_FAILED;]])],
  [AC_MSG_RESULT(yes)],
  [AC_MSG_RESULT(no)
   AC_DEFINE([MAP_FAILED],[((void *)-1)],[Define MAP_FAILED here if not available otherwise.])])

AC_LANG_RESTORE
])


