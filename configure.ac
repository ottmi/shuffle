AC_INIT
AC_CONFIG_SRCDIR([src/shuffle.cpp])
AC_LANG(C++)

# Checks for programs.
AC_PROG_CXX

AC_ARG_WITH([libgmp],
            [AS_HELP_STRING([--with-libgmp],
              [use GNU Multiple Precision Library @<:@default=check@:>@])],
            [],
            [with_libgmp=check])

AC_ARG_WITH([openmp],
            [AS_HELP_STRING([--with-openmp],
              [Compile OpenMP parallel code @<:@default=check@:>@])],
            [],
            [with_openmp=check])


AC_CHECK_HEADERS([cmath], , AC_MSG_ERROR([This package needs libm]))
AC_CHECK_LIB(m, log, , AC_MSG_ERROR([This package needs libm]))

GMP=0
if test "$with_libgmp" != "no"; then
  AC_CHECK_HEADERS([gmp.h], GMP=1, GMP=0)
  AC_CHECK_LIB(gmp, __gmpz_init, GMP=$(($GMP+1)), GMP=0)

  AC_CHECK_HEADERS([mpfr.h], GMP=$(($GMP+1)), GMP=0)
  AC_CHECK_LIB(mpfr, mpfr_log, GMP=$(($GMP+1)), GMP=0)
fi

if test $GMP -ge 4; then
  CXXFLAGS="$CXXFLAGS -D_GMP"
  LIBS="$LIBS -lgmp -lmpfr"
fi


if test "$with_openmp" != "no"; then
  m4_include([ax_openmp.m4])
  AX_OPENMP
  CXXFLAGS="$CXXFLAGS $OPENMP_CXXFLAGS"
fi

AC_OUTPUT(Makefile)

echo
echo "================================================================================"
echo "Summary:"
echo -n "  Using libgmp: "
if test $GMP -ge 4; then
  echo "yes"
else
  echo "no"
fi
echo -n "  Using OpenMP: "
if test "$OPENMP_CXXFLAGS" != ""; then
  echo "yes"
else
  echo "no"
fi
echo "================================================================================"
