AC_DEFUN([ACX_STRUCT_STAT64], [
AC_LANG_SAVE
AC_LANG(C)

AC_MSG_CHECKING([for struct stat64])
AC_COMPILE_IFELSE(
  [AC_LANG_PROGRAM([[#include <sys/stat.h>]],
	           [[struct stat64 st; stat64("/tmp/foo",&st);]])],
  [AC_MSG_RESULT(yes)
   AC_DEFINE([HAVE_STRUCT_STAT64],[1],[have struct stat64])],
  [AC_MSG_RESULT(no)])

AC_LANG_RESTORE
])

