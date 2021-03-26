# Makefile for somthin'

SHELL=/bin/sh

include ../../makevars

fobjects = klby.o radina.o writepwf.o writerecpot.o

cobjects = cwritepwf.o do_fhi.o do_ncpp.o do_pwf.o do_upf.o nwrite.o do_recpot.o do_siesta.o do_cpmd.o do_casino.o do_champ.o

all_objects = $(fobjects) $(cobjects)

this_dir = output

.F.o:
	if [ $(OSTYPE) = darwin ] ; then \
         $(FPP) $(FPPFLAGS) -I../include $< > $<.f ;\
	 $(FC) $(FFLAGS) $(OPTIONS) -c -I../include $<.f -o $@ ; \
	 rm -f $<.f ; \
        else \
         $(FC) $(FFLAGS) $(OPTIONS)  -I../include $*.F -c ; \
        fi

all:  $(all_objects)
	-rm -f objects
	for i in $(all_objects); do echo $(this_dir)/$$i >> objects; done

clean:
	rm -f $(fobjects) $(cobjects) liboutput.a core* *~ *.f.f objects

distclean: clean


