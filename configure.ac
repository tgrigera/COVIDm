#
# configure.ac -- autoconf script for COVIDm
#
# COVIDm is copyright (c) 2020 by the authors (see AUTHORS)
#
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT(COVIDm, 0.3, tgrigera@iflysib.unlp.edu.ar)
AC_CONFIG_SRCDIR([qdrandom.cc])
AC_CONFIG_HEADERS([config.h])
AC_LANG(C++)
AC_CONFIG_MACRO_DIR([m4])
AM_INIT_AUTOMAKE

# Checks for programs.
AC_PROG_MAKE_SET
AC_PROG_CPP
AC_PROG_CC
AC_PROG_CXX
AC_PROG_LN_S
AC_PROG_INSTALL

# Checks for typedefs, structures, and compiler characteristics.
AC_CHECK_HEADER_STDBOOL
AC_C_INLINE

# Checks for library functions.
AC_CHECK_FUNCS([strerror])

# Checks for header files.
AC_CHECK_HEADERS([stdlib.h string.h])
AC_CHECK_LIB([gslcblas],[cblas_dgemm])
AC_CHECK_HEADER([gsl/gsl_rng.h],,exit)
AC_CHECK_HEADER([lemon/list_graph.h],,exit)

# Checks for libraries.
AC_CHECK_LIB([gsl], [gsl_rng_alloc],,exit)
dnl AX_CXX_CHECK_LIB([lemon], [lemon::ListDigraph])
dnl AS_IF([test -z $HAVE_LEMON],AC_MSG_FAILURE([cannot build without liblemon]))


# Generate Makefile
AC_CONFIG_FILES([Makefile graph/Makefile])
AC_OUTPUT
