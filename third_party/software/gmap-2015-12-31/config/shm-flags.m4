
AC_DEFUN([ACX_SHM_FLAGS], [
AC_LANG_SAVE
AC_LANG(C)

AC_MSG_CHECKING(for SHM_NORESERVE in shmget)
AC_COMPILE_IFELSE(
  [AC_LANG_PROGRAM([[#include <sys/ipc.h>
#include <sys/shm.h>]],
                   [[int flags = SHM_NORESERVE;]])],
  [AC_MSG_RESULT(yes)
   AC_DEFINE([HAVE_SHM_NORESERVE],[1],[Define to 1 if SHM_NORESERVE available for shmget.])],
  [AC_MSG_RESULT(no)])

AC_LANG_RESTORE
])


