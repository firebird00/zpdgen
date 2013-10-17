lib_LIBRARIES = libgpdf.a
libgpdf_a_SOURCES = gpdf_inmzpd.c wofzwh.f rpolyev.f inmzpd.f
libgpdf_a_CFLAGS = -fPIC -I./
libgpdf_a_FFLAGS = -I./ -fPIC -fimplicit-none
libgpdf_a_INCLUDES=$(all_includes)
EXTRA_DIST=gpdf.py
include_HEADERS=gpdf.h
noinst_HEADERS=doxy.cfg
SUBDIRS = amos quadp
ACLOCAL_AMFLAGS = -I m4

noinst_PROGRAMS = ftest ctest

ftest_SOURCES = ftest.f
ftest_LDADD = -L./ -L./amos/ -L./quadp -lgpdf -lamos -lquadp -lgfortran
ftest_CFLAGS = -I./
ftest_FFLAGS = -I./ -fimplicit-none
ftest_DEPENDENCIES = ./libgpdf.a ./amos/libamos.a

ctest_SOURCES = ctest.c
ctest_LDADD = -L./ -L./amos/ -L./quadp -lgpdf -lamos -lquadp -lgfortran
ctest_CFLAGS = -I./
ctest_FFLAGS = -I./ -fimplicit-none
ctest_DEPENDENCIES = ./libgpdf.a ./amos/libamos.a

if ZPDGEN_PYTHON_ON
bin_SCRIPTS=inmzpd.pyf
endif
CLEANFILES=$(bin_SCRIPTS) *.so
inmzpd.pyf: inmzpd.f Makefile.am
	$(F2PY) --overwrite-signature -m inmzpd -h inmzpd.pyf only : inmweid : inmzpd.f
	sed 's/:: res/,intent(out) :: res/' <inmzpd.pyf | sed 's/:: za/, intent(in) :: za/'| sed 's/:: b/, intent(in) :: b/'| sed 's/:: zb/, intent(in) :: zb/'| sed 's/:: n/, intent(in) :: n/'| sed 's/:: m/, intent(in) :: m/'| sed 's/:: z1/, intent(in) :: z1/'| sed 's/:: z2/, intent(in) :: z2/'>tmpout
	mv tmpout inmzpd.pyf
	$(F2PY) -c inmzpd.pyf $(ftest_LDADD)
