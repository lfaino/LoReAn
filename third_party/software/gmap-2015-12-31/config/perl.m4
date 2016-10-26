
AC_DEFUN([ACX_PATH_PERL],[
  AC_MSG_CHECKING([for a working version of perl])
  if test x"$PERL" != x; then
    # 'PERL' is specified by user
    if ( $PERL -e 'use IO::File; use Getopt::Std;' ) > /dev/null 2>&1; then
      found=yes
    else
      $as_unset PERL
      found=no
    fi
  else
    found=no
  fi

  if test $found = no; then
    candidate_perl_names='perl perl5'
    as_save_IFS=$IFS; IFS=$PATH_SEPARATOR
    for as_dir in $PATH; do
      IFS=$as_save_IFS
      test -z "$as_dir" && as_dir=.
      for perl in $candidate_perl_names; do
        if ( $as_dir/$perl -e 'use IO::File; use Getopt::Std;' ) > /dev/null 2>&1; then
          ac_cv_path_PERL=$as_dir/$perl
          found=yes
          break 2
        fi
      done
    done
    PERL=$ac_cv_path_PERL
  fi

  if test -n "$PERL"; then
    AC_MSG_RESULT($PERL)
  else
    AC_MSG_RESULT([not found])
    perl_warning=yes
  fi

  AC_SUBST(PERL)
])

