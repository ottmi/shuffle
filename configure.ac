AC_INIT
AC_CONFIG_SRCDIR([src/shuffle.cpp])
AC_LANG(C++)

if test -z $CXXFLAGS; then
    CXXFLAGS='-O3'
fi

# Checks for programs.
AC_PROG_CXX

AC_ARG_WITH([libmpfr],
            [AS_HELP_STRING([--with-libmpfr],
              [use GNU MPFR Library @<:@default=yes@:>@])],
            [],
            [with_libmpfr=yes])

AC_ARG_WITH([openmp],
            [AS_HELP_STRING([--with-openmp],
              [Compile OpenMP parallel code @<:@default=yes@:>@])],
            [],
            [with_openmp=yes])


AC_CHECK_HEADERS([cmath], , AC_MSG_ERROR([This package needs libm]))
AC_CHECK_LIB(m, log, , AC_MSG_ERROR([This package needs libm]))

MPFR=0
if test "$with_libmpfr" != "no"; then
  if test "$with_libmpfr" != "yes"; then
    MPFR_INC=`find $with_libmpfr -name "mpfr.h"|head -n1`
    if test "$MPFR_INC" != ""; then
      MPFR_INC=`dirname $MPFR_INC`
    fi
    MPFR_LIB=`find $with_libmpfr -name "libmpfr.so"|head -n1`
    if test "$MPFR_LIB" != ""; then
      MPFR_LIB=`dirname $MPFR_LIB`
    fi

    CPPFLAGS="-I$MPFR_INC $CPPFLAGS"
    CXXFLAGS="-I$MPFR_INC $CXXFLAGS"
    LDFLAGS="-L$MPFR_LIB $LDFLAGS"
  fi

  AC_CHECK_HEADERS([mpfr.h], MPFR=$(($MPFR+1)), MPFR=0)
  AC_CHECK_LIB(mpfr, mpfr_log, MPFR=$(($MPFR+1)), MPFR=0)
fi

if test $MPFR -ge 2; then
  CXXFLAGS="$CXXFLAGS -D_MPFR"
  LIBS="$LIBS -lmpfr"
fi


if test "$with_openmp" != "no"; then
  m4_include([ax_openmp.m4])
  AX_OPENMP
  CXXFLAGS="$CXXFLAGS $OPENMP_CXXFLAGS"
  LDFLAGS="$LDFLAGS $OPENMP_CXXFLAGS"
fi

AC_OUTPUT(Makefile)

echo
echo "================================================================================"
echo "Summary:"
if test $MPFR -ge 2; then
  echo "  Using libmpfr: yes"
else
  echo "  Using libmpfr: no"
fi
if test "$OPENMP_CXXFLAGS" != ""; then
  echo "  Using OpenMP: yes"
else
  echo "  Using OpenMP: no"
fi
echo "================================================================================"
