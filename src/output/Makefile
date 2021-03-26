# Makefile for somthin'

SHELL=/bin/sh

include ../../makevars

fobjects = klby.o radina.o writepwf.o writerecpot.o

cobjects = cwritepwf.o do_fhi.o do_ncpp.o do_pwf.o do_qeupf.o do_spinor.o nwrite.o do_recpot.o do_siesta.o do_cpmd.o do_casino.o do_champ.o do_teter.o do_qso.o

all_objects = $(fobjects) $(cobjects)

this_dir = output

.F.o:
	$(FC) $(FFLAGS) $(OPTIONS)  -I../include $*.F -c 


all:  $(all_objects)
	-rm -f objects
	for i in $(all_objects); do echo $(this_dir)/$$i >> objects; done

clean:
	rm -f $(fobjects) $(cobjects) liboutput.a core* *~ *.f.f objects

distclean: clean


