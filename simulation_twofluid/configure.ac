#Process this with autoreconf -vfi

AC_INIT([runsim],[0.1],[kchiou@sas.upenn.edu]) #initialize
AC_CONFIG_SRCDIR([threshold_wave_simulation.c]) #main program
AC_PROG_CC #compiler check

AM_INIT_AUTOMAKE

# Check for libraries like math and gsl
AC_CHECK_LIB([gsl],[gslblas_dgemm])
AC_CHECK_LIB([gsl],[gsl_blas_dgemm])
AC_CHECK_LIB([m],[cos])

# Check for header files
AC_HEADER_STDC
AC_CHECK_HEADERS([stdlib.h]);

# Check for typedefs, structures, compiler characteristics via const
AC_C_CONST

# Check for other library (base) functions
AC_FUNC_MALLOC
AC_FUNC_REALLOC

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
