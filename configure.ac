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
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([Makefile src/Makefile tests/Makefile])
AC_OUTPUT
