AC_INIT([spokes], [0.1], [andrew@acooke.org])
AC_CONFIG_MACRO_DIR([m4])
AM_INIT_AUTOMAKE([foreign -Wall -Werror])
AM_PROG_AR()
LT_INIT
AC_PROG_CC_C99
dnl following required by check
AM_PROG_CC_C_O
PKG_CHECK_MODULES([CHECK], [check >= 0.9.0])
AC_CHECK_LIB([lu], [ludebug], [], [AC_MSG_ERROR([No liblu found])])
AC_CHECK_LIB([cairo], [cairo_create], [], [AC_MSG_ERROR([No libcairo found])])
dnl AC_CHECK_LIB([cblas], [cblas_dgemm], [LIBS="$LIBS -lopenblas"], [AC_MSG_ERROR([No libcblas found])], [-lopenblas])
dnl AC_CHECK_LIB([lapacke], [LAPACKE_dposv], [LIBS="$LIBS -llapacke -llapack -lopenblas -lm -lgfortran"], [AC_MSG_ERROR(["No liblapacke found"])], [-llapack -lopenblas -lm -lgfortran])
AC_CHECK_LIB([m], [cos], [], [AC_MSG_ERROR([No libm found])])
AC_CHECK_LIB([gslcblas], [cblas_dgemm], [], [AC_MSG_ERROR([No libgslcblas found])])
AC_CHECK_LIB([gsl], [gsl_blas_dgemm], [], [AC_MSG_ERROR([No libgsl found])])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([Makefile src/Makefile tests/Makefile])
AC_OUTPUT
