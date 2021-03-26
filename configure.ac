# Process this file with autoconf to produce a configure script.

AC_INIT(opium, ac1.2, ejwalt@wm.edu)

AC_CONFIG_SRCDIR([src/uniPPlib/xmcomplex.h])

AC_ARG_WITH(f2c,
[  --with-f2c=PFX             Prefix where fort77 and f2c executables,
                             includes, and libraries are found],
f2c_prefix="$withval")
AC_ARG_WITH(f2c-libraries,
[  --with-f2c-libraries=DIR   Directory where libf2c is located],
f2c_libraries="$withval")
AC_ARG_WITH(f2c-includes,
[  --with-f2c-includes=DIR    Directory where f2c.h is located],
f2c_includes="$withval")

dnl dont need a common header....
dnl AC_CONFIG_HEADER([config.h])

AC_CANONICAL_HOST

if test "${f2c_prefix+set}" = set; then
  PATH=$f2c_prefix/bin:$PATH; export PATH
  F77=fort77
fi

dnl Checks for programs.
AC_PROG_F77(ifort ifc efc gfortran pgf77 pgf90 fort fort77 pathf90)
AC_F77_LIBRARY_LDFLAGS

if test x"$F77" = x; then
  AC_MSG_ERROR([
*** a suitable fortran compiler is required, see the README.fortran
])
fi

AC_PROG_CC(ccc pgcc icc ecc pathcc gcc cc)

if test x"$CC" = x; then
  AC_MSG_ERROR([
*** You need to install a working C compiler, check with your vendor
])
fi

AC_PROG_CPP

AC_SUBST(LD,[$CC])
AC_SUBST(FC,[$F77])
AC_SUBST(FPP,[$FPP])
AC_SUBST(FPPFLAGS,[$FPPFLAGS])
AC_SUBST(F2C_H,[$F2C_H])
AC_SUBST(NOOPT,[$NOOPT])


FOO_MONKEY="$PATH:/lib:/usr/lib"
AC_PATH_PROG([FPP], [cpp], [no], [$FOO_MONKEY])

FPPFLAGS="$FPPFLAGS"
AC_MSG_CHECKING([if $FPP needs -traditional])
need_traditional=no
cat <<_ACEOF > conftest.f
c well do'i need -traditional?????
_ACEOF

if $FPP $FPPFLAGS conftest.f 1>/dev/null 2>/dev/null; then
  AC_MSG_RESULT(no)
else
  FPPFLAGS="-traditional $FPPFLAGS"
  if $FPP $FPPFLAGS conftest.f 1>/dev/null 2>/dev/null; then
    AC_MSG_RESULT(yes)
  else
    AC_MSG_ERROR([*** Could not figure out how to preprocess fortran files ***])
  fi
fi
rm -f conftest.f

AC_MSG_CHECKING([if $FPP is broken with unnecessary pragma statements])
need_pragma_removed=no
cat <<_ACEOF > conftest.c
int main();
_ACEOF

if $FPP $FPPFLAGS conftest.c 2>/dev/null | grep '#pragma' >/dev/null 2>/dev/null; then
  cat <<_ACEOF > fpp.fixed
#!/bin/sh

$FPP \$* | grep -v '^#pragma'
_ACEOF
  chmod +x fpp.fixed
  FPP="`pwd`/fpp.fixed"
  FPPFLAGS=''
  if $FPP $FPPFLAGS conftest.c >/dev/null 2>/dev/null; then
    AC_MSG_RESULT(yes)
  else
    AC_MSG_ERROR([*** Could not figure out how to preprocess fortran files ***])
  fi
else
  AC_MSG_RESULT(no)
fi
rm -f conftest.c


AC_PROG_MAKE_SET
AC_PROG_RANLIB

dnl dunno why I cant just do this directly, but...
FOO_MONKEY="$PATH:/usr/ccs/bin"
AC_PATH_PROG([AR], [ar], [no], [$FOO_MONKEY])

if test x"$AR" = xno; then
  AC_MSG_ERROR([
*** ar could not be found on your system, make sure you have a development package installed
])
fi

AC_CHECK_LIB([m], [cos],[LDFLAGS="$LDFLAGS -lm"])

FCT=`echo "$F77" |sed -e s,^.*/,,`
CCT=`echo "$F77" |sed -e s,^.*/,,`
dnl echo $FCT $CCT
dnl echo $F77 $CC
dnl echo "$FFLAGS"



case $ac_cv_host in

dnl Sun's f77 and C compiler
  *sun*)
    if test "$FCT" = "fort77"; then
      :
    else
      if test x"$CFLAGS" = "x-g"; then
	CFLAGS="-fast -dalign"
      fi
      if test x"$FFLAGS" = "x-g"; then
	FFLAGS="-fast -dalign -xtypemap=real:64,double:64"
      fi
      LD="$F77 -fast -dalign"
      NOOPT="-O0 -xtypemap=real:64,double:64"
    fi
    ;;

dnl SGI's native C and fortran compilers
  *sgi*)
    if test "$FCT" = "fort77"; then
      :
    else
      LD=$F77
      if test x"$FFLAGS" = "x-g"; then
	FFLAGS="-O2 -r8"
      fi
      if test x"$CFLAGS" = "x-g"; then
	CFLAGS="-O3"
      fi
      NOOPT="-O0"
    fi
    ;;

dnl alpha linux with Compaq (HP)'s ccc and fort compilers
  *alpha*linux*)
    if test "$CCT" = "ccc"; then
      if test x"$CFLAGS" = "x-g"; then
	CFLAGS="-tune host"
      fi
    fi
    if test "$FCT" = "fort"; then
      if test x"$FFLAGS" = "x-g"; then
	FFLAGS="-tune host -r8"
      fi
      if echo $LDFLAGS | grep Futil >/dev/null; then
	:
      else
	LDFLAGS="$LDFLAGS -lUfor -lFutil -lfor"
      fi
      NOOPT="-O0"
    fi
    ;;

dnl 32 bit PC running linux
  *86-pc-linux-*)
    NOOPT="-O0"	
    if test $FCT = pgf77 || test $FCT = pgf90; then
      dnl do this to get the library path correct
      AC_PATH_PROG([XPGF77],[pgf77],[huh],[$PATH])
      PGF77_LIB_DIR=`echo "$XPGF77/lib" | sed s,/bin/pgf77,,`
      LDFLAGS="-L$PGF77_LIB_DIR $LDFLAGS"
      test x"$FFLAGS" = x-g && FFLAGS=" -fast -r8"
    fi    
    
    if test $FCT = ifc || test $FCT = ifort; then
      LD=$F77
      LDFLAGS="$LDFLAGS -Vaxlib"
      test x"$FFLAGS" = x-g && FFLAGS="-O3 -save -r8 -fp-model strict"
      LDFLAGS="$LDFLAGS -L$i -static-libcxa"

      if test $FCT = ifort; then
        LDFLAGS="$LDFLAGS -nofor_main"
      fi
    fi

    if test x"$CCT" = xpgcc; then
      test x"$CFLAGS" = x-g && CFLAGS="-O2"
    fi

    if test x"$CCT" = xicc; then
      test x"$CFLAGS" = x-g && CFLAGS="-O2"
    fi

    if test $FCT = gfortran ; then
    echo ""
    echo ""
    echo "*******************************************************************"
    echo "*******************************************************************"
    echo "Using gfortran for OPIUM is experimental... please be careful      "
    echo "*******************************************************************"
    echo "*******************************************************************"
    echo ""
    echo ""

    LD=$F77
    LDFLAGS="-fdefault-real-8 -fdefault-double-8 -fno-automatic $LDFLAGS"
    FFLAGS="$FFLAGS -fdefault-real-8 -fdefault-double-8 -fno-automatic"
    NOOPT="-O0"		     
    fi
    if test $FCT = gcc ; then
      test x"$CFLAGS" = x-g && CFLAGS="-O2"
    fi
    ;;

  *-apple-*)
    NOOPT="-O0"		     

    if test $FCT = ifc || test $FCT = ifort; then
      LD=$FCT
      LDFLAGS="$LDFLAGS -Vaxlib"
      test x"$FFLAGS" = x-g && FFLAGS="-O3 -save -r8 -fp-model strict"
      LDFLAGS="$LDFLAGS -L$i -static-libcxa"

      if test $FCT = ifort; then
        LDFLAGS="$LDFLAGS -nofor_main"
      fi

    fi

    if test x"$CCT" = xicc; then
      test x"$CFLAGS" = x-g && CFLAGS="-O2"
    fi

    if test $FCT = gfortran ; then
    echo ""
    echo ""
    echo "*******************************************************************"
    echo "*******************************************************************"
    echo "Using gfortran for OPIUM is experimental... please be careful      "
    echo "*******************************************************************"
    echo "*******************************************************************"
    echo ""
    echo ""
    LD=$F77
    LDFLAGS="-lm -fdefault-real-8 -fdefault-double-8 -fno-automatic $LDFLAGS"
    FFLAGS="-fdefault-real-8 -fdefault-double-8 -fno-automatic $FFLAGS"
    fi
    if test $FCT = gcc ; then
      test x"$CFLAGS" = x-g && CFLAGS="-O2"
    fi
    ;;


dnl 64bit PC running linux 
  *x86_64*linux-*)
    NOOPT="-O0"	
    if test $FCT = pgf77 || test $FCT = pgf90; then
      dnl do this to get the library path correct
      AC_PATH_PROG([XPGF77],[pgf77],[huh],[$PATH])
      PGF77_LIB_DIR=`echo "$XPGF77/lib" | sed s,/bin/pgf77,,`
      LDFLAGS="-L$PGF77_LIB_DIR $LDFLAGS"
      test x"$FFLAGS" = x-g && FFLAGS=" -fast -r8"
    fi    

    if test $FCT = ifc || test $FCT = ifort; then
      LD=$FCT
      LDFLAGS="$LDFLAGS -Vaxlib"
      test x"$FFLAGS" = x-g && FFLAGS="-O3 -save -r8 -fp-model strict"
      LDFLAGS="$LDFLAGS -L$i -static-libcxa"

      if test $FCT = ifort; then
        LDFLAGS="$LDFLAGS -nofor_main"
      fi

    fi

    if test $FCT = pathf77; then
       FFLAGS="-fast -r8"
    fi

    if test x"$CCT" = xpgcc; then
      test x"$CFLAGS" = x-g && CFLAGS="-O2"
    fi

    if test x"$CCT" = xicc; then
      test x"$CFLAGS" = x-g && CFLAGS="-O2"
    fi

    if test $FCT = gfortran ; then
    echo ""
    echo ""
    echo "*******************************************************************"
    echo "*******************************************************************"
    echo "Using gfortran for OPIUM is experimental... please be extra careful"
    echo "*******************************************************************"
    echo "*******************************************************************"
    echo ""
    echo ""
    LD=$F77
    LDFLAGS="-lm -fdefault-real-8 -fno-automatic $LDFLAGS"
    test x"$FFLAGS" = x-g && FFLAGS="-g -O2 -fdefault-real-8 -fno-automatic"
    fi
    if test $FCT = gcc ; then
      test x"$CFLAGS" = x-g && CFLAGS="-O2"
    fi
    ;;

esac


AC_PATH_PROG([XFORT77],[fort77],[huh],[$PATH])

dnl looks like we are stuck with fort77 and f2c
if test $FCT = fort77; then
  NOOPT="-O0"
  FFLAGS="-r8 $FFLAGS"
  if test "${f2c_libraries+set}" = set; then
    LDFLAGS="-L$f2c_libraries $LDFLAGS"
  elif test "${f2c_prefix+set}" = set; then
    LDFLAGS="-L$f2c_prefix/lib $LDFLAGS"
  else
    FORT77_LIB_DIR=`echo "$XFORT77/lib" | sed s,/bin/fort77,,`
    LDFLAGS="-L$FORT77_LIB_DIR $LDFLAGS"
  fi
  CHECK_FAILED=no
  AC_CHECK_LIB([f2c],[do_fio],[LDFLAGS="$LDFLAGS -lf2c"],[CHECK_FAILED=yes])
  if test $CHECK_FAILED = yes; then
    AC_MSG_ERROR([
*** could not find libf2c, run ./configure --help to find the flags for
*** specifying where the f2c library is located
])
  fi

  xLDFLAGS_SAVE="$LDFLAGS"
  xSAVE_CFLAGS="$CFLAGS"
  xSAVE_CPPFLAGS="$CPPFLAGS"
  if test "${f2c_includes+set}" = set; then
    CPPFLAGS="-I$f2c_includes"
  elif test "${f2c_prefix+set}" = set; then
    CPPFLAGS="-I$f2c_prefix/include"
  else
    LDFLAGS=
    FORT77_H_DIR=`echo $XFORT77/include | sed s,/bin/fort77,,`
    CPPFLAGS="-I$FORT77_H_DIR"
  fi
  CHECK_FAILED=no
  AC_CHECK_HEADER([f2c.h],:,CHECK_FAILED=yes)
  if test $CHECK_FAILED = yes; then
    AC_MSG_ERROR([
*** could not find f2c.h, run --configure --help to find the flags for
*** specifying where the f2c include is located
])
  fi

  F2C_H="$CPPFLAGS"
  CPPFLAGS="$xSAVE_CPPFLAGS"
  CFLAGS="$xSAVE_CFLAGS"
  LDFLAGS="$xLDFLAGS_SAVE"

dnl  test x"$FFLAGS" = x-g && FFLAGS="-O3 -r8"
dnl  test x"$FFLAGS" = "x-g -O2" && FFLAGS="-O3 -r8"
dnl Boy! This is a hack EJW
fi

dnl Checks for header files.
xLDFLAGS_SAVE="$LDFLAGS"
LDFLAGS=
AC_HEADER_STDC
AC_CHECK_HEADERS([stdlib.h string.h unistd.h])
LDFLAGS="$xLDFLAGS_SAVE"

dnl Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST

dnl Checks for library functions.
dnl AC_FUNC_MALLOC
dnl AC_FUNC_STRFTIME
dnl AC_CHECK_FUNCS([memset pow sqrt strchr])

CFLAGS="-I../include -I../FlexiLib -I../uniPPlib $CFLAGS"
CDATE=\'`date`\'
CHOST=`hostname`
CSYS=\'`uname -srp`\'
CFLAGS="$CFLAGS -DCSYS=\\\"${CSYS}\\\" -DCHOST=\\\"${CHOST}\\\" -DCDATE=\\\"${CDATE}\\\" "
AC_CONFIG_FILES([makevars])

AC_OUTPUT

echo ""
echo ""
echo "Building with the following options:"
echo ""
echo "Fortran compiler     :       $F77"
echo "Fortan flags         :       $FFLAGS"
echo "Fortran preprocessor :       $FPP"
#echo "Fortran preproc flags:       $FPPFLAGS"
echo "C compiler           :       $CC"
echo "C flags              :       $CFLAGS"
echo "Linker               :       $LD"
echo "Linker flags         :       $LDFLAGS $FLIBS"
echo ""
echo ""

# vim: ai sw=2
