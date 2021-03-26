# Makefile for src/atom

SHELL=/bin/sh

include ../../makevars

fobjects = atm.o cubint.o difrel.o dsolv2.o etotal.o excorr.o ghostnl.o \
	hrtree.o interp.o logder.o nllogd.o pbe.o getpcc.o pratt.o radin.o \
	schsl.o scpot.o splift.o spliq.o velect.o denkcomp.o interp2.o

cobjects = do_ae.o do_nl.o do_tc.o do_fc.o

all_objects = $(fobjects) $(cobjects)

this_dir = atom

.F.o:
	if [ $(OSTYPE) = darwin ] ; then \
         $(FPP) $(FPPFLAGS) -I../include $< > $<.f ;\
	 $(FC) $(FFLAGS) $(OPTIONS) -c -I../include $<.f -o $@ ; \
	 rm -f $<.f ; \
        else \
         $(FC) $(FFLAGS) $(OPTIONS)  -I../include $*.F -c ; \
        fi

all: $(all_objects)
	-rm -f objects
	for i in $(all_objects); do echo $(this_dir)/$$i >> objects; done

clean:
	rm -f $(fobjects) $(cobjects) libatom.a core* *~ *.f.f objects

distclean: clean


