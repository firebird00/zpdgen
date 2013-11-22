lib_LIBRARIES = libgpdf.a
libgpdf_a_SOURCES = wofzmw.f wofzh.f wofzwh16.f wofzwh8.f wofzwh4.f cpolyev.f inmzpd.f epszpd.f inmkur.f utils.f cbj0.f pldisp_utils.c pldisp_specfun.c
libgpdf_a_CFLAGS = -DH5_USE_18_API -fPIC -I./
libgpdf_a_FFLAGS = -I./ -fPIC -fimplicit-none
libgpdf_a_INCLUDES=$(all_includes)
EXTRA_DIST=gpdf.py
include_HEADERS=gpdf.h
noinst_HEADERS=doxy.cfg
SUBDIRS = quadp
ACLOCAL_AMFLAGS = -I m4

noinst_PROGRAMS = timetest test

timetest_SOURCES = timetest.f
timetest_LDADD = -L./ -L./quadp -lgpdf -lquadp -lgfortran
timetest_CFLAGS = -I./
timetest_FFLAGS = -I./ -fimplicit-none
timetest_DEPENDENCIES = ./libgpdf.a

test_SOURCES = test.c
test_LDADD = -L./ -L./quadp -lgpdf -lquadp -lgfortran
test_CFLAGS = -I./
test_FFLAGS = -I./ -fimplicit-none
test_DEPENDENCIES = ./libgpdf.a

if ZPDGEN_PYTHON_ON
bin_SCRIPTS=inmzpd.pyf
endif
CLEANFILES=$(bin_SCRIPTS) *.so
inmzpd.pyf: inmzpd.f Makefile.am
	$(F2PY) --overwrite-signature -m inmzpd -h inmzpd.pyf only: inmweid : inmzpd.f
	sed 's/:: res/,intent(out) :: res/' <inmzpd.pyf | sed 's/:: za/, intent(in) :: za/'| sed 's/:: b/, intent(in) :: b/'| sed 's/:: zb/, intent(in) :: zb/'| sed 's/:: n/, intent(in) :: n/'| sed 's/:: m/, intent(in) :: m/'| sed 's/:: z1/, intent(in) :: z1/'| sed 's/:: z2/, intent(in) :: z2/'>tmpout
	mv tmpout inmzpd.pyf
	$(F2PY) -c inmzpd.pyf $(timetest_LDADD)
