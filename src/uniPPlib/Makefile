# This is the Makefile for 'uniPPlib' distributed with OPIUM
#

SHELL=/bin/sh

include ../../makevars

objects = uniPP_read.o uniPP_write.o uniPP_writefhi.o uniPP_free.o uniPP_alpha.o uniPP_ffloc.o uniPP_ffnlcc.o uniPP_kbden.o uniPP_kbfnl.o

all_objects = $(objects)

this_dir = uniPPlib

all: $(all_objects)
	-rm -f objects
	for i in $(all_objects); do echo $(this_dir)/$$i >> objects; done

uniPP_read.o: uniPPlib.h
uniPP_write.o: uniPPlib.h
uniPP_free.o: uniPPlib.h 
uniPP_alpha.o: uniPPlib.h
uniPP_ffloc.o: uniPPlib.h
uniPP_kbden.o: uniPPlib.h
uniPP_kbfnl.o: uniPPlib.h

clean: 
	rm -f $(objects) libuniPP.a core* *~ objects

distclean: clean

