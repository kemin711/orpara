#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([orpara], [0.1.0], [kmzhou4@yaho.com])
AM_INIT_AUTOMAKE([-Wall -Werror])
AM_PROG_AR
AC_CONFIG_SRCDIR([src/interval.cpp])
AC_CONFIG_HEADERS([config.h])
LT_PREREQ([2.4])
LT_INIT([dlopen])
AC_CONFIG_MACRO_DIR([m4])

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC
AC_PROG_INSTALL
AC_LANG([C], [C++])
# Checks for libraries.
AC_LANG_PUSH([C++])
AX_CXX_COMPILE_STDCXX_11([noext], [mandatory])
AC_LANG_POP([C++])

save_LDFLAGS=$LDFLAGS
save_LIBS=$LIBS
LDFLAGS="-L/usr/lib/x86_64-linux-gnu $LDFLAGS"
LIBS="$LIBS -lpthread"
AX_PTHREAD([AC_DEFINE([HAVE_PTHREAD], [1], [Define if you have POSIX threads libraries and headers])
	AC_MSG_NOTICE([pthread_libs: $PTHREAD_LIBS $PTHREAD_LDFLAGS])]
	,[AC_MSG_ERROR([could not find pthread])])
AC_SUBST([PTHREAD_LIBS])
LDFLAGS=$save_LDFLAGS
LIBS=$save_LIBS
# could not figure out why PTHREAD_xxx variables are set to empty string!

AC_ARG_VAR([GTESTINC], [gtest PATH for gtest include dir])
AC_ARG_VAR([GTESTLIB], [gtest PATH for gtest lib dir])
if (test -n "$GTESTINC" && test -z "$GTESTLIB") || \
	(test -z "$GTESTINC" && test -n "$GTESTLIB"); then 
   AC_MSG_ERROR([please specify both lib and include dir for gtest])
fi
if test "$prefix" = "NONE"; then
	prefix=$ac_default_prefix
fi
AC_LANG_PUSH([C++])
save_LIBS=$LIBS
save_LDFLAGS=$LDFLAGS
save_CXXFLAGS=$CXXFLAGS
LDFLAGS="-L$prefix/lib $LDFLAGS"
CXXFLAGS="-I$prefix/include $CXXFLAGS"
LIBS="$LIBS -lgtest -lgtest_main -lgmock -lgmock_main -lpthread"
if test -n "$GTESTINC"; then
	LDFLAGS="-L$GTESTLIB $LDFLAGS"
	CXXFLAGS="-I$GTESTINC $CXXFLAGS"
fi
CXXFLAGS="$PTHREAD_CXXFLAGS $CXXFLAGS"
LDFLAGS="$PTHREAD_LDFLAGS $LDFLAGS"
AC_LINK_IFELSE(
	[AC_LANG_PROGRAM([#include <gmock/gmock.h>], [testing::Cardinality dummy])],
	[gtest_LIBS=$LIBS] [gtest_CXXFLAGS=$CXXFLAGS] [gtest_LDFLAGS=$LDFLAGS] [HAVE_GTEST=1] [HAVE_GMOCK=1],
	[AC_MSG_WARN([libgtest and libgmock not found, will not be able to make check.])])
AC_SUBST([gtest_CFLAGS])
AC_SUBST([gtest_CXXFLAGS])
AC_SUBST([gtest_LIBS])
LIBS=$save_LIBS
LDFLAGS=$save_LDFLAGS
CXXFLAGS=$save_CXXFLAGS
AC_LANG_POP([C++])

# Checks for typedefs, structures, and compiler characteristics.
AC_CHECK_HEADER_STDBOOL
AC_TYPE_SIZE_T

# check headers

# Checks for library functions.
AC_CHECK_FUNCS([pow sqrt strchr strstr])
AC_CHECK_PROGS([DOXYGEN], [doxygen])
if test -z "$DOXYGEN"; then
	AC_MSG_WARN([Doxygen not found - continue without Doxygen support])
else
        AC_MSG_NOTICE([Doxygen installed check dot and pdflatex ...])
	AC_CHECK_PROGS([DOT], [dot])
	if test -z "$DOT"; then
		AC_MSG_ERROR([Doxygen needs dot part of graphviz, please install dot first])
	fi
	AC_CHECK_PROGS([PDFLATEX], [pdflatex])
	if test -z "$PDFLATEX"; then
		AC_MSG_ERROR([Doxygen needs pdflatex program, it is part of TeX http://www.tug.org/texlive/acquire-netinstall.html])
	fi
fi
AM_CONDITIONAL([HAVE_DOXYGEN], [test -n "$DOXYGEN"])
AM_COND_IF([HAVE_DOXYGEN], [AC_CONFIG_FILES([docs/Doxyfile])])
AM_COND_IF([HAVE_DOXYGEN], [AC_CONFIG_FILES([docs/Makefile])])

AC_CONFIG_FILES([Makefile
                 src/Makefile
                 src/seqaln/Makefile
                 src/test/Makefile])
AC_OUTPUT
