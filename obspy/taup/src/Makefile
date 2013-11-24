CC = gcc
FC = gfortran
PYTHON = /usr/bin/python

CFLAGS = -pthread -fno-strict-aliasing -DNDEBUG -g -fwrapv -O2 -Wall -Wstrict-prototypes -fPIC
FFLAGS = -Wall -ffixed-form -fno-second-underscore -fPIC -O3 -funroll-loops

INCLUDES = -I/usr/share/pyshared/numpy/core/include -I/usr/include/python2.7 


_libtau.so: _libtaumodule.c emdlv.f libtau.f ttimes_subrout.f _libtau.pyf
    # copy f2py dependenceis from source tree fortranobject.c fortranobject.h 
	cp /usr/share/pyshared/numpy/f2py/src/fortranobject.h .
	cp /usr/share/pyshared/numpy/f2py/src/fortranobject.c .

	# play around to generate _libtaumodule.c, by defaul, this includes to
	# many subroutins, the skip only flags only work with -c options...
	# http://cens.ioc.ee/pipermail/f2py-users/2006-April/001250.html
	#f2py -c -m _libtau _libtau.pyf emdlv.f libtau.f ttimes_subrout.f

	# building wrappers: libtaumodule and fortran objects
	$(CC) $(CFLAGS) $(INCLUDES) -c fortranobject.c -o fortranobject.o
	$(CC) $(CFLAGS) $(INCLUDES) -c _libtaumodule.c -o _libtaumodule.o
	
	# building emdv.f and libtau.f
	$(FC) $(FFLAGS) -c -c emdlv.f -o emdlv.o
	$(FC) $(FFLAGS) -c -c libtau.f -o libtau.o
	$(FC) $(FFLAGS) -c -c ttimes_subrout.f -o ttimes_subrout.o
	
	# linking shared library
	$(FC) -Wall -Wall -shared _libtaumodule.o fortranobject.o emdlv.o libtau.o ttimes_subrout.o -lgfortran -o _libtau.so
	# test that it worked
	$(PYTHON) -c "import _libtau; print _libtau.__doc__ + '\n----\n\nOK'"
