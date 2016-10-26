
AC_DEFUN([ACX_SSE2_SHIFT_DEFECT], [
  AC_REQUIRE([AC_CANONICAL_HOST])
  AC_LANG_SAVE
  AC_LANG_C

  AC_MSG_CHECKING(compiler is defective and requires an immediate in sse2 shift commands)
  AC_COMPILE_IFELSE(
  [AC_LANG_PROGRAM([[#include <stdio.h>
#include <stdlib.h>
#include <emmintrin.h>]],
                   [[int nshift = rand() % 32;
__m128i shifted;
shifted = _mm_slli_epi32(_mm_set1_epi32(1),nshift);
]])],
  [ax_cv_sse2_shift_defect=no],
  [ax_cv_sse2_shift_defect=yes])

  AC_MSG_RESULT($ax_cv_sse2_shift_defect)
  if test "$ax_cv_sse2_shift_defect" = yes; then
    AC_DEFINE(DEFECTIVE_SSE2_COMPILER,1,[Define to 1 if your compiler is defective and requires an immediate in sse2 shift commands.])
  fi

AC_LANG_RESTORE
])dnl

