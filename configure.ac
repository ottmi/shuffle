AC_INIT
AC_CONFIG_SRCDIR([src/shuffle.cpp])
AC_LANG(C++)

# Checks for programs.
AC_PROG_CXX


AC_CHECK_HEADERS([cmath], , AC_MSG_ERROR([This package needs libm]))
AC_CHECK_LIB(m, log, , AC_MSG_ERROR([This package needs libm]))

AC_CHECK_HEADERS([gmp.h], GMP=1, GMP=0)
AC_CHECK_LIB(gmp, __gmpz_init, GMP=$(($GMP+1)), GMP=0)

AC_CHECK_HEADERS([gmpxx.h], GMP=$(($GMP+1)), GMP=0)
AC_CHECK_LIB(gmpxx, _init, GMP=$(($GMP+1)), GMP=0)

AC_CHECK_HEADERS([mpfr.h], GMP=$(($GMP+1)), GMP=0)
AC_CHECK_LIB(mpfr, mpfr_log, GMP=$(($GMP+1)), GMP=0)

if test $GMP -ge 6; then
  CXXFLAGS="$CXXFLAGS -D_GMP"
  LIBS="$LIBS -lgmp -lgmpxx -lmpfr"
fi


m4_include([ax_openmp.m4])
AX_OPENMP
CXXFLAGS="$CXXFLAGS $OPENMP_CXXFLAGS"


AC_OUTPUT(Makefile)
