lib_LIBRARIES = libgpdf.a
libgpdf_a_SOURCES = pldisp_utils.c pldisp_specfun.c wofzh.f wofzwh.f rpolyev.f inmzpd.f epszpd.f inmkur.f
libgpdf_a_CFLAGS = -DH5_USE_18_API -fPIC -I./ 
libgpdf_a_FFLAGS = -I./ -fPIC
libgpdf_a_INCLUDES=$(all_includes)
EXTRA_DIST=gpdf.py
noinst_HEADERS=gpdf.h doxy.cfg
SUBDIRS = amos quadp
ACLOCAL_AMFLAGS = -I m4

bin_PROGRAMS = test

test_SOURCES = test.c
test_LDADD = -L./ -L./amos/ -L./quadp -lgpdf -lamos -lquadp -lgfortran
test_CFLAGS = -I./
test_DEPENDENCIES = ./libgpdf.a ./amos/libamos.a
if ZPDGEN_PYTHON_ON
bin_SCRIPTS=inmzpd.pyf
endif
CLEANFILES=$(bin_SCRIPTS) *.so
inmzpd.pyf: inmzpd.f
	$(F2PY) --overwrite-signature -m inmzpd -h inmzpd.pyf only: inmweid : inmzpd.f
	sed 's/:: res/,intent(out) :: res/' <inmzpd.pyf | sed 's/:: numza/, intent(in) :: numza/' | sed 's/:: za/, intent(in) :: za/'| sed 's/:: b/, intent(in) :: b/'| sed 's/:: zb/, intent(in) :: zb/'| sed 's/:: n/, intent(in) :: n/'| sed 's/:: m/, intent(in) :: m/'>tmpout
	mv tmpout inmzpd.pyf
	$(F2PY) -c inmzpd.pyf $(test_LDADD)
