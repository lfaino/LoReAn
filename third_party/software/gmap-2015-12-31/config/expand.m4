
AC_DEFUN([ACX_EXPAND], [
  EXP_VAR=[$1]
  FROM_VAR=[$2]

  dnl first expand prefix and exec_prefix
  prefix_save=$prefix
  exec_prefix_save=$exec_prefix

  dnl set defaults if necesary
  if test "x$prefix" = xNONE; then
    prefix="$ac_default_prefix"
  fi
  if test "x$exec_prefix" = xNONE; then
    exec_prefix=$prefix
  fi

  full_var="$FROM_VAR"
  while true; do
    new_full_var="`eval echo $full_var`"
    if test "x$new_full_var" = "x$full_var"; then break; fi
    full_var=$new_full_var
  done

  full_var=$new_full_var
  AC_SUBST([$1],"$full_var")

  prefix=$prefix_save
  exec_prefix=$exec_prefix_save
])dnl ACX_EXPAND

