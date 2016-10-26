# ===========================================================================
#          http://www.gnu.org/software/autoconf-archive/ax_ext.html
# ===========================================================================
#
# SYNOPSIS
#
#   AX_EXT
#
# DESCRIPTION
#
#   Find supported SIMD extensions by requesting cpuid. When an SIMD
#   extension is found, the -m"simdextensionname" is added to SIMD_CFLAGS if
#   compiler supports it. For example, if "sse2" is available, then "-msse2"
#   is added to SIMD_CFLAGS.
#
#   This macro calls:
#
#     AC_SUBST(SIMD_CFLAGS)
#
#   And defines:
#
#     HAVE_MMX / HAVE_SSE / HAVE_SSE2 / HAVE_SSE3 / HAVE_SSSE3 / HAVE_SSE4_1 / HAVE_SSE4_2 / HAVE_AVX
#
# LICENSE
#
#   Copyright (c) 2007 Christophe Tournayre <turn3r@users.sourceforge.net>
#   Copyright (c) 2013 Michael Petch <mpetch@capp-sysware.com>
#
#   Copying and distribution of this file, with or without modification, are
#   permitted in any medium without royalty provided the copyright notice
#   and this notice are preserved. This file is offered as-is, without any
#   warranty.

#serial 13

AC_DEFUN([AX_EXT],
[
  AC_REQUIRE([AC_CANONICAL_HOST])

  CFLAGS_ORIG=$CFLAGS

  case $host_cpu in
    powerpc*)
      AC_CACHE_CHECK([whether altivec is enabled and supported], [ax_cv_cpu_have_altivec_ext],
          [
            if test `/usr/sbin/sysctl -a 2>/dev/null| grep -c hw.optional.altivec` != 0; then
                if test `/usr/sbin/sysctl -n hw.optional.altivec` = 1; then
                  ax_cv_cpu_have_altivec_ext=yes
                fi
            fi
          ])

          if test "$ax_cv_cpu_have_altivec_ext" = yes; then
            AC_DEFINE(HAVE_ALTIVEC,1,[Define to 1 if you support Altivec instructions])
            AX_CHECK_COMPILE_FLAG(-faltivec, SIMD_CFLAGS="$SIMD_CFLAGS -faltivec", [])
          fi
    ;;


    i[[3456]]86*|x86_64*|amd64*)

      AC_REQUIRE([AX_GCC_X86_CPUID])
      AC_REQUIRE([AX_GCC_X86_AVX_XGETBV])

      AX_GCC_X86_CPUID(0x00000001)
      ecx=0
      edx=0
      if test "$ax_cv_gcc_x86_cpuid_0x00000001" != "unknown";
      then
        ecx=`echo $ax_cv_gcc_x86_cpuid_0x00000001 | cut -d ":" -f 3`
        edx=`echo $ax_cv_gcc_x86_cpuid_0x00000001 | cut -d ":" -f 4`
      fi

      AC_CACHE_CHECK([whether mmx is enabled and supported by CPU], [ax_cv_cpu_have_mmx_ext],
      [
        ax_cv_cpu_have_mmx_ext=no
        if test "$((0x$edx>>23&0x01))" = 1; then
          ax_cv_cpu_have_mmx_ext=yes
          ax_cv_cpu_features="$ax_cv_cpu_features mmx"

        fi
      ])

      AC_CACHE_CHECK([whether sse is enabled and supported by CPU], [ax_cv_cpu_have_sse_ext],
      [
        ax_cv_cpu_have_sse_ext=no
        if test "$((0x$edx>>25&0x01))" = 1; then
          ax_cv_cpu_have_sse_ext=yes
          ax_cv_cpu_features="$ax_cv_cpu_features sse"
        fi
      ])

      AC_CACHE_CHECK([whether sse2 is enabled and supported by CPU], [ax_cv_cpu_have_sse2_ext],
      [
        ax_cv_cpu_have_sse2_ext=no
        if test "$((0x$edx>>26&0x01))" = 1; then
          ax_cv_cpu_features="$ax_cv_cpu_features sse2"
          if test "$ax_cv_want_sse2_ext" = yes; then
            ax_cv_cpu_have_sse2_ext=yes
          fi
        fi
      ])

      AC_CACHE_CHECK([whether sse3 is enabled and supported by CPU], [ax_cv_cpu_have_sse3_ext],
      [
        ax_cv_cpu_have_sse3_ext=no
        if test "$((0x$ecx&0x01))" = 1; then
          ax_cv_cpu_have_sse3_ext=yes
          ax_cv_cpu_features="$ax_cv_cpu_features sse3"
        fi
      ])

      AC_CACHE_CHECK([whether ssse3 is enabled and supported by CPU], [ax_cv_cpu_have_ssse3_ext],
      [
        ax_cv_cpu_have_ssse3_ext=no
        if test "$((0x$ecx>>9&0x01))" = 1; then
          ax_cv_cpu_features="$ax_cv_cpu_features ssse3"
          if test "$ax_cv_want_ssse3_ext" = yes; then
            ax_cv_cpu_have_ssse3_ext=yes
          fi
        fi
      ])

      AC_CACHE_CHECK([whether sse4.1 is enabled and supported by CPU], [ax_cv_cpu_have_sse41_ext],
      [
        ax_cv_cpu_have_sse41_ext=no
        if test "$((0x$ecx>>19&0x01))" = 1; then
          ax_cv_cpu_features="$ax_cv_cpu_features sse4.1"
          if test "$ax_cv_want_sse41_ext" = yes; then
            ax_cv_cpu_have_sse41_ext=yes
          fi
        fi
      ])

      AC_CACHE_CHECK([whether sse4.2 is enabled and supported by CPU], [ax_cv_cpu_have_sse42_ext],
      [
        ax_cv_cpu_have_sse42_ext=no
        if test "$((0x$ecx>>20&0x01))" = 1; then
          ax_cv_cpu_features="$ax_cv_cpu_features sse4.2"
          if test "$ax_cv_want_sse42_ext" = yes; then
            ax_cv_cpu_have_sse42_ext=yes
          fi
        fi
      ])

      AC_CACHE_CHECK([whether avx is enabled and supported by CPU], [ax_cv_cpu_have_avx_ext],
      [
        ax_cv_cpu_have_avx_ext=no
        if test "$((0x$ecx>>28&0x01))" = 1; then
          ax_cv_cpu_features="$ax_cv_cpu_features avx"
          if test "$ax_cv_want_avx_ext" = yes; then
            ax_cv_cpu_have_avx_ext=yes
          fi
        fi
      ])

      if test x"$ax_cv_cpu_have_avx_ext" = x"yes"; then
        AX_GCC_X86_AVX_XGETBV(0x00000000)

        xgetbv_eax="0"
        if test x"$ax_cv_gcc_x86_avx_xgetbv_0x00000000" != x"unknown"; then
          xgetbv_eax=`echo $ax_cv_gcc_x86_avx_xgetbv_0x00000000 | cut -d ":" -f 1`
        fi

        AC_CACHE_CHECK([whether avx is supported by operating system], [ax_cv_have_avx_os_ext],
        [
          ax_cv_have_avx_os_ext=no

          if test "$((0x$ecx>>27&0x01))" = 1; then
            if test "$((0x$xgetbv_eax&0x6))" = 6; then
              ax_cv_have_avx_os_ext=yes
            fi
          fi
        ])
        if test x"$ax_cv_have_avx_os_ext" = x"no"; then
          AC_MSG_WARN([Your processor supports AVX, but your operating system doesn't])
        fi
      fi


      AC_CACHE_CHECK([whether popcnt is enabled and supported by CPU], [ax_cv_cpu_have_popcnt_ext],
      [
        ax_cv_cpu_have_popcnt_ext=no
        if test "$((0x$ecx>>23&0x01))" = 1; then
          ax_cv_cpu_have_popcnt_ext=yes
        fi
      ])


      AX_GCC_X86_CPUID(0x80000001)
      ecx=`echo $ax_cv_gcc_x86_cpuid_0x80000001 | cut -d ":" -f 3`
      edx=`echo $ax_cv_gcc_x86_cpuid_0x80000001 | cut -d ":" -f 4`

      AC_CACHE_CHECK([whether lzcnt is enabled and supported by CPU], [ax_cv_cpu_have_lzcnt_ext],
      [
        ax_cv_cpu_have_lzcnt_ext=no
        if test "$((0x$ecx>>5&0x01))" = 1; then
          ax_cv_cpu_have_lzcnt_ext=yes
        fi
      ])


      AX_GCC_X86_CPUID(0x00000007)
      ebx=`echo $ax_cv_gcc_x86_cpuid_0x00000007 | cut -d ":" -f 2`
      ecx=`echo $ax_cv_gcc_x86_cpuid_0x00000007 | cut -d ":" -f 3`
      edx=`echo $ax_cv_gcc_x86_cpuid_0x00000007 | cut -d ":" -f 4`

      AC_CACHE_CHECK([whether avx2 is enabled and supported by CPU], [ax_cv_cpu_have_avx2_ext],
      [
        ax_cv_cpu_have_avx2_ext=no
        if test "$((0x$ebx>>5&0x01))" = 1; then
          ax_cv_cpu_features="$ax_cv_cpu_features avx2"
          if test "$ax_cv_want_avx2_ext" = yes; then
            ax_cv_cpu_have_avx2_ext=yes
          fi
        fi
      ])

      AC_CACHE_CHECK([whether bmi1 is enabled and supported by CPU], [ax_cv_cpu_have_bmi1_ext],
      [
        ax_cv_cpu_have_bmi1_ext=no
        if test "$((0x$ebx>>3&0x01))" = 1; then
          ax_cv_cpu_have_bmi1_ext=yes
        fi
      ])

      AC_CACHE_CHECK([whether bmi2 is enabled and supported by CPU], [ax_cv_cpu_have_bmi2_ext],
      [
        ax_cv_cpu_have_bmi2_ext=no
        if test "$((0x$ebx>>8&0x01))" = 1; then
          ax_cv_cpu_have_bmi2_ext=yes
        fi
      ])

      if test "$ax_cv_cpu_have_mmx_ext" = yes; then
        AX_CHECK_COMPILE_FLAG(-mmmx, [ax_cv_compile_mmx_ext=yes], [ax_cv_ext_compile_problem=yes])
        if test x"$ax_cv_compile_mmx_ext" != x"yes"; then
          AC_MSG_WARN([Your CPU supports MMX instructions but not your compiler.  Can you try another compiler or update yours?])
        else
          AC_MSG_CHECKING(for mmintrin.h header file)
          CFLAGS=-mmmx
          AC_LINK_IFELSE([AC_LANG_PROGRAM([#include <mmintrin.h>])],
                         [ax_cv_link_mmintrin_h=yes],
  		         [ax_cv_ext_linker_problem=yes])
          if test x"$ax_cv_link_mmintrin_h" != x"yes"; then
            AC_MSG_RESULT([no])
            AC_MSG_WARN([Your compiler supports MMX instructions but not your linker.])
          else
            AC_MSG_RESULT([yes])
            SIMD_CFLAGS="$SIMD_CFLAGS -mmmx"
            AC_DEFINE(HAVE_MMX,1,[Define to 1 if you support MMX instructions])
          fi            
        fi
      fi


      if test "$ax_cv_cpu_have_sse_ext" = yes; then
        AX_CHECK_COMPILE_FLAG(-msse, [ax_cv_compile_sse_ext=yes], [ax_cv_ext_compile_problem=yes])
        if test x"$ax_cv_compile_sse_ext" != x"yes"; then
          AC_MSG_WARN([Your CPU supports SSE instructions but not your compiler.  Can you try another compiler or update yours?])
        else
          AC_MSG_CHECKING(for xmmintrin.h header file)
          CFLAGS=-msse
          AC_LINK_IFELSE([AC_LANG_PROGRAM([#include <xmmintrin.h>])],
                         [ax_cv_link_xmmintrin_h=yes],
		         [ax_cv_ext_linker_problem=yes])
          if test x"$ax_cv_link_xmmintrin_h" != x"yes"; then
            AC_MSG_RESULT([no])
            AC_MSG_WARN([Your compiler supports SSE instructions but not your linker.])
          else
            AC_MSG_RESULT([yes])
            SIMD_CFLAGS="$SIMD_CFLAGS -msse"
            AC_DEFINE(HAVE_SSE,1,[Define to 1 if you support SSE (Streaming SIMD Extensions) instructions])
          fi            
        fi
      fi


      if test "$ax_cv_cpu_have_sse2_ext" = yes; then
	AX_CHECK_COMPILE_FLAG(-msse2, [ax_cv_compile_sse2_ext=yes], [ax_cv_ext_compile_problem=yes])
	if test x"$ax_cv_compile_sse2_ext" != x"yes"; then
	  AC_MSG_WARN([Your CPU supports SSE2 instructions but not your compiler.  Can you try another compiler or update yours?])
	else
          AC_MSG_CHECKING(for emmintrin.h header file)
          CFLAGS=-msse2
          AC_LINK_IFELSE([AC_LANG_PROGRAM([#include <emmintrin.h>])],
                         [ax_cv_link_emmintrin_h=yes],
		         [ax_cv_ext_linker_problem=yes])
          if test x"$ax_cv_link_emmintrin_h" != x"yes"; then
            AC_MSG_RESULT([no])
            AC_MSG_WARN([Your compiler supports SSE2 instructions but not your linker.])
          else
            AC_MSG_RESULT([yes])
  	    SIMD_CFLAGS="$SIMD_CFLAGS -msse2"
  	    AC_DEFINE(HAVE_SSE2,1,[Define to 1 if you support SSE2 (Streaming SIMD Extensions 2) instructions])
          fi            
	fi
      fi


      if test "$ax_cv_cpu_have_sse3_ext" = yes; then
        AX_CHECK_COMPILE_FLAG(-msse3, [ax_cv_compile_sse3_ext=yes], [ax_cv_ext_compile_problem=yes])
        if test x"$ax_cv_compile_sse3_ext" != x"yes"; then
          AC_MSG_WARN([Your CPU supports SSE3 instructions but not your compiler.  Can you try another compiler or update yours?])
        else
          AC_MSG_CHECKING(for pmmintrin.h header file)
          CFLAGS=-msse3
          AC_LINK_IFELSE([AC_LANG_PROGRAM([#include <pmmintrin.h>])],
                         [ax_cv_link_pmmintrin_h=yes],
		         [ax_cv_ext_linker_problem=yes])
          if test x"$ax_cv_link_pmmintrin_h" != x"yes"; then
            AC_MSG_RESULT([no])
            AC_MSG_WARN([Your compiler supports SSE3 instructions but not your linker.])
          else
            AC_MSG_RESULT([yes])
            SIMD_CFLAGS="$SIMD_CFLAGS -msse3"
            AC_DEFINE(HAVE_SSE3,1,[Define to 1 if you support SSE3 (Streaming SIMD Extensions 3) instructions])
          fi            
        fi
      fi


      if test "$ax_cv_cpu_have_ssse3_ext" = yes; then
        AX_CHECK_COMPILE_FLAG(-mssse3, [ax_cv_compile_ssse3_ext=yes], [ax_cv_ext_compile_problem=yes])
        if test x"$ax_cv_compile_ssse3_ext" != x"yes"; then
          AC_MSG_WARN([Your CPU supports SSSE3 instructions but not your compiler.  Can you try another compiler or update yours?])
        else
          AC_MSG_CHECKING(for tmmintrin.h header file)
          CFLAGS=-mssse3
          AC_LINK_IFELSE([AC_LANG_PROGRAM([#include <tmmintrin.h>])],
                         [ax_cv_link_tmmintrin_h=yes],
		         [ax_cv_ext_linker_problem=yes])
          if test x"$ax_cv_link_tmmintrin_h" != x"yes"; then
            AC_MSG_RESULT([no])
            AC_MSG_WARN([Your compiler supports SSSE3 instructions but not your linker.])
          else
            AC_MSG_RESULT([yes])
            SIMD_CFLAGS="$SIMD_CFLAGS -mssse3"
            AC_DEFINE(HAVE_SSSE3,1,[Define to 1 if you support SSSE3 (Supplemental Streaming SIMD Extensions 3) instructions])
          fi            
        fi
      fi


      if test "$ax_cv_cpu_have_sse41_ext" = yes; then
	AX_CHECK_COMPILE_FLAG(-msse4.1, [ax_cv_compile_sse41_ext=yes], [ax_cv_ext_compile_problem=yes])
	if test x"$ax_cv_compile_sse41_ext" != x"yes"; then
	  AC_MSG_WARN([Your CPU supports SSE4.1 instructions but not your compiler.  Can you try another compiler or update yours?])
	else
          AC_MSG_CHECKING(for smmintrin.h header file)
          CFLAGS=-msse4.1
          AC_LINK_IFELSE([AC_LANG_PROGRAM([#include <smmintrin.h>])],
                         [ax_cv_link_smmintrin_h=yes],
		         [ax_cv_ext_linker_problem=yes])
          if test x"$ax_cv_link_smmintrin_h" != x"yes"; then
            AC_MSG_RESULT([no])
            AC_MSG_WARN([Your compiler supports SSE4.1 instructions but not your linker.])
          else
            AC_MSG_RESULT([yes])
  	    SIMD_CFLAGS="$SIMD_CFLAGS -msse4.1"
  	    AC_DEFINE(HAVE_SSE4_1,1,[Define to 1 if you support SSE4.1 (Streaming SIMD Extensions 4.1) instructions])
          fi            
	fi
      fi


      if test "$ax_cv_cpu_have_sse42_ext" = yes; then
        AX_CHECK_COMPILE_FLAG(-msse4.2, [ax_cv_compile_sse42_ext=yes], [ax_cv_ext_compile_problem=yes])
        if test x"$ax_cv_compile_sse42_ext" != x"yes"; then
          AC_MSG_WARN([Your CPU supports SSE4.2 instructions but not your compiler.  Can you try another compiler or update yours?])
        else
          AC_MSG_CHECKING(for nmmintrin.h header file)
          CFLAGS=-msse4.2
          AC_LINK_IFELSE([AC_LANG_PROGRAM([#include <nmmintrin.h>])],
                         [ax_cv_link_nmmintrin_h=yes],
		         [ax_cv_ext_linker_problem=yes])
          if test x"$ax_cv_link_nmmintrin_h" != x"yes"; then
            AC_MSG_RESULT([no])
            AC_MSG_WARN([Your compiler supports SSE4.2 instructions but not your linker.])
          else
            AC_MSG_RESULT([yes])
            SIMD_CFLAGS="$SIMD_CFLAGS -msse4.2"
            AC_DEFINE(HAVE_SSE4_2,1,[Define to 1 if you support SSE4.2 (Streaming SIMD Extensions 4.2) instructions])

            if test "$ax_cv_cpu_have_popcnt_ext" = yes; then
              AC_RUN_IFELSE(
                [AC_LANG_PROGRAM([[#include <nmmintrin.h>]],
                                 [[return (_mm_popcnt_u32(0xffffffffu) == 32) ? 0 : 9;]])],
			         [ax_cv_run_mm_popcnt_ext=yes])
              if test x"$ax_cv_run_mm_popcnt_ext" = x"yes"; then
                AC_DEFINE(HAVE_MM_POPCNT,1,[Define to 1 if you support Intel intrinsic _mm_popcnt_u32/64 instructions])
              fi
            fi
          fi            
        fi
      fi


      if test "$ax_cv_cpu_have_avx_ext" = yes; then
        AX_CHECK_COMPILE_FLAG(-mavx, [ax_cv_compile_avx_ext=yes], [ax_cv_ext_compile_problem=yes])
        if test x"$ax_cv_compile_avx_ext" != x"yes"; then
          AC_MSG_WARN([Your CPU supports AVX instructions but not your compiler.  Can you try another compiler or update yours?])
        else
          AC_MSG_CHECKING(for immintrin.h header file)
          CFLAGS=-mavx
          AC_LINK_IFELSE([AC_LANG_PROGRAM([#include <immintrin.h>])],
                         [ax_cv_link_immintrin_h=yes],
		         [ax_cv_ext_linker_problem=yes])
          if test x"$ax_cv_link_immintrin_h" != x"yes"; then
            AC_MSG_RESULT([no])
            AC_MSG_WARN([Your compiler supports AVX instructions but not your linker.])
          else
            AC_MSG_RESULT([yes])
            SIMD_CFLAGS="$SIMD_CFLAGS -mavx"
            AC_DEFINE(HAVE_AVX,1,[Define to 1 if you support AVX (Advanced Vector Extensions) instructions])

            if test "$ax_cv_cpu_have_popcnt_ext" = yes; then
              AC_RUN_IFELSE(
                [AC_LANG_PROGRAM([[#include <immintrin.h>]],
                                 [[return (_popcnt32(0xffffffffu) == 32) ? 0 : 9;]])],
			         [ax_cv_run_popcnt_ext=yes])
              if test x"$ax_cv_run_popcnt_ext" = x"yes"; then
                AC_DEFINE(HAVE_POPCNT,1,[Define to 1 if you support Intel intrinsic _popcnt instruction])
              fi
            fi

            if test "$ax_cv_cpu_have_lzcnt_ext" = yes; then
              AC_RUN_IFELSE(
                [AC_LANG_PROGRAM([[#include <immintrin.h>]],
                                 [[return (_lzcnt_u32(0xffffffffu) == 32) ? 0 : 9;]])],
			         [ax_cv_run_lzcnt_ext=yes])
              if test x"$ax_cv_run_lzcnt_ext" = x"yes"; then
                AC_DEFINE(HAVE_LZCNT,1,[Define to 1 if you support Intel intrinsic _lzcnt instruction])
              fi
            fi

            if test "$ax_cv_cpu_have_bmi1_ext" = yes; then
              CFLAGS=-mbmi
              AC_RUN_IFELSE(
                [AC_LANG_PROGRAM([[#include <immintrin.h>]],
                                 [[return (_tzcnt_u32(0xffffffffu) == 32) ? 0 : 9;]])],
			         [ax_cv_run_tzcnt_ext=yes])
              if test x"$ax_cv_run_tzcnt_ext" = x"yes"; then
                SIMD_CFLAGS="$SIMD_CFLAGS -mbmi"
                AC_DEFINE(HAVE_TZCNT,1,[Define to 1 if you support Intel intrinsic _tzcnt instruction])
              fi
            fi

            if test "$ax_cv_cpu_have_bmi2_ext" = yes; then
              SIMD_CFLAGS="$SIMD_CFLAGS -mbmi2"
              AC_DEFINE(HAVE_BMI2,1,[Define to 1 if you support BMI2 (Bit Manipulation Instruction set 2)])
            fi

          fi            
        fi
      fi


      if test "$ax_cv_cpu_have_avx2_ext" = yes; then
        AX_CHECK_COMPILE_FLAG(-mavx2, [ax_cv_compile_avx2_ext=yes], [ax_cv_ext_compile_problem=yes])
        if test x"$ax_cv_compile_avx2_ext" != x"yes"; then
          AC_MSG_WARN([Your CPU supports AVX2 instructions but not your compiler.  Can you try another compiler or update yours?])
        else
          AC_MSG_CHECKING(for immintrin.h header file)
          CFLAGS=-mavx2
          AC_LINK_IFELSE([AC_LANG_PROGRAM([#include <immintrin.h>])],
                         [ax_cv_link_immintrin_h=yes],
		         [ax_cv_ext_linker_problem=yes])
          if test x"$ax_cv_link_immintrin_h" != x"yes"; then
            AC_MSG_RESULT([no])
            AC_MSG_WARN([Your compiler supports AVX2 instructions but not your linker.])
          else
            AC_MSG_RESULT([yes])
            SIMD_CFLAGS="$SIMD_CFLAGS -mavx2"
            AC_DEFINE(HAVE_AVX2,1,[Define to 1 if you support AVX2 (Advanced Vector Extensions 2) instructions])

          fi            
        fi
      fi

  ;;
  esac

  CFLAGS=$CFLAGS_ORIG

  AC_SUBST(SIMD_CFLAGS)
])
