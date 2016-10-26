
AC_DEFUN([ACX_PAGESIZE], [
AC_LANG_SAVE
AC_LANG(C)

AC_MSG_CHECKING(for pagesize via sysconf)
AC_COMPILE_IFELSE(
  [AC_LANG_PROGRAM([[#include <unistd.h>]],
                   [[sysconf(_SC_PAGESIZE);]])],
  [AC_MSG_RESULT(yes)
   AC_DEFINE([PAGESIZE_VIA_SYSCONF],[1],[pagesize is available via sysconf])],
  [AC_MSG_RESULT(no)])

AC_MSG_CHECKING(for pagesize via sysctl)
AC_COMPILE_IFELSE(
  [AC_LANG_PROGRAM([[#include <sys/types.h>
#include <sys/sysctl.h>]],
                   [[int mib[2], pagesize; size_t len = sizeof(pagesize); mib[0] = CTL_HW; mib[1] = HW_PAGESIZE; sysctl(mib,2,&pagesize,&len,NULL,0);]])],
  [AC_MSG_RESULT(yes)
   AC_DEFINE([PAGESIZE_VIA_SYSCTL],[1],[pagesize is available via sysctl])],
  [AC_MSG_RESULT(no)])

AC_LANG_RESTORE
])

