AC_INIT
AC_CONFIG_SRCDIR([src/shuffle.cpp])
AC_LANG(C++)

#if test -z $CXXFLAGS; then
#    CXXFLAGS='-O3'
#fi

AC_ARG_ENABLE([debug],
              [AS_HELP_STRING([--enable-debug],[Enable debug mode: creates reproducible randomizations for POC scores])])
AC_ARG_WITH([openmp],
            [AS_HELP_STRING([--with-openmp],
            [Compile OpenMP parallel code @<:@default=auto@:>@])],
            [],
            [with_openmp=auto])
AC_ARG_WITH([mpi],
            [AS_HELP_STRING([--with-mpi],
            [Compile MPI parallel code @<:@default=no@:>@])],
            [],
            [with_mpi=no])
AC_ARG_WITH([libmpfr],
            [AS_HELP_STRING([--with-libmpfr],
            [use GNU MPFR Library @<:@default=auto@:>@])],
            [],
            [with_libmpfr=auto])

if test "$with_mpi" = "yes" && test "$with_openmp" = "yes"; then
  AC_MSG_FAILURE([Cannot compile for OpenMP and MPI at the same time])
fi

m4_include([ax_prog_cxx_mpi.m4])
AX_PROG_CXX_MPI([test "$with_mpi" != "no"],[use_mpi=yes],[use_mpi=no])
if test "$use_mpi" = "no"; then
  if test "$with_mpi" = "yes"; then
    AC_MSG_FAILURE([MPI requested, but failed to enable it])
  fi
else
  CXXFLAGS="$CXXFLAGS -D_MPI"
fi

use_openmp=no
if test "$with_openmp" != "no" && test "$use_mpi" = "no"; then
  m4_include([ax_openmp.m4])
  AX_OPENMP([use_openmp=yes],)
  if test "$use_openmp" = "no"; then  
    if test "$with_openmp" = "yes"; then
      AC_MSG_FAILURE([OpenMP requested, but failed to enable it])
    fi
  else
    CXXFLAGS="$CXXFLAGS $OPENMP_CXXFLAGS"
    LDFLAGS="$LDFLAGS $OPENMP_CXXFLAGS"
  fi
fi

if test "$enable_debug" = yes; then
  CXXFLAGS="$CXXFLAGS -D_DEBUG"
fi

AC_CHECK_HEADERS([cmath], , AC_MSG_ERROR([This package needs libm]))
AC_CHECK_LIB(m, log, , AC_MSG_ERROR([This package needs libm]))

use_mpfr=no
if test "$with_libmpfr" != "no"; then
  if test "$with_libmpfr" != "yes" && test "$with_libmpfr" != "auto"; then
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

  AC_CHECK_HEADERS([mpfr.h], use_mpfr=yes, )
  if test "$use_mpfr" = "yes"; then
    AC_CHECK_LIB(mpfr, mpfr_log, , use_mpfr=no)
  fi

  if test "$use_mpfr" = "yes"; then
    CXXFLAGS="$CXXFLAGS -D_MPFR"
  else
    if test "$with_libmpfr" = "yes"; then
      AC_MSG_FAILURE([libmpfr requested, but failed to enable it])
    fi
  fi
fi

AC_OUTPUT(Makefile)


echo
echo "================================================================================"
echo "Summary:"

if test "$enable_debug" = "yes"; then
  echo "  Debug mode: yes"
else
  echo "  Debug mode: no"
fi

if test "$use_mpfr" = "yes"; then
  echo "  Using libmpfr: yes"
else
  echo "  Using libmpfr: no"
fi

if test "$use_openmp" = "yes"; then
  echo "  Using OpenMP: yes"
else
  echo "  Using OpenMP: no"
fi

if test "$use_mpi" = "yes"; then
  echo "  Using MPI: yes"
else
  echo "  Using MPI: no"
fi

echo "================================================================================"
