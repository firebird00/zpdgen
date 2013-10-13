#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.61)

AC_INIT([ieps_gpdf], [0.1], [ogurcan@lpp.polytechnique.fr])

AM_INIT_AUTOMAKE([-Wall -Werror foreign])

# Checks for programs.
AC_PROG_CC
AC_PROG_F77
AM_PROG_AR
AC_PROG_RANLIB
AM_PROG_CC_C_O
AC_CONFIG_MACRO_DIR([m4])

AC_ARG_ENABLE(debug, AC_HELP_STRING([--enable-debug],[enables debugging]), ok=$enableval, ok=no)
if test "$ok" = "yes"; then
	AC_DEFINE(ZPDGEN_DEBUG,1,[define to enable debugging code.])
	CFLAGS=" -g -O0 -pg"
	FFLAGS=" -g -O0 -pg"
else
	CFLAGS="-O0"
	FFLAGS="-O0"
fi

AC_ARG_ENABLE(python, AC_HELP_STRING([--enable-python],[enables python wrappers]), ok=$enableval, ok=no)
if test "$ok" = "yes"; then
	AC_DEFINE(ZPDGEN_PYTHON,1,[define to enable python wrappers.])
	AM_CONDITIONAL(ZPDGEN_PYTHON_ON,true)
	ZPDGEN_PYTHON_DEF=true
else
	ZPDGEN_PYTHON_DEF=false
	AM_CONDITIONAL(ZPDGEN_PYTHON_ON,false)
fi

if $ZPDGEN_PYTHON_DEF ; then
   AM_CHECK_PYMOD([numpy],[ndarray],,
   [
	AC_MSG_FAILURE([could not find the Python module "numpy"])
   ])
   
   AC_PATH_PROG(F2PY, f2py, false)
   if test x"$F2PY" = xfalse ; then
      AC_MSG_ERROR([F2PY not found; see README on f2py])
   fi
   AC_SUBST(F2PY)
fi

# Checks for libraries.

AC_CHECK_LIB([pthread],[main],[],[echo "Error! You need libpthread."
                        echo "see http://xmlsoft.org/"
                        exit -1
                        ])

AC_CHECK_LIB([m],[cos],[],[echo "Error! You need libm."
                        echo  "see www.gnu.org/software/libc/"
                        exit -1
                        ])

AC_CHECK_LIB([hdf5],[main],[],[echo "Error! You need hdf5."
                        echo "see http://hdf.ncsa.uiuc.edu/HDF5/"
                        exit -1
                        ])

AC_CHECK_LIB([gfortran],[main],[],[echo "Error! You need gfortran runtime library."
                        echo  "see http://gcc.gnu.org/fortran/"
                        exit -1
                        ])

# check for MPI stuff
#  if test -z "$CC" && test -n "$MPICC"; then
#     CC="$MPICC"
#  fi
#  AC_PROG_CC
#  AX_MPI
#  CC="$MPICC"
#  LIBS="$MPILIBS $LIBS"

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([stdlib.h])
#AC_CHECK_HEADERS([omp.h],[CFLAGS+=" -fopenmp"],[])

# Checks for typedefs, structures, and compiler characteristics.

# Checks for library functions.
AC_FUNC_MALLOC

#AC_CONFIG_FILES([Makefile])
AC_OUTPUT(Makefile amos/Makefile quadp/Makefile)
