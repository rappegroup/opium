# Makefile for somthin'

SHELL=/bin/sh

include ../../makevars

fobjects = ai.o besfn.o biq.o conmax.o cq.o descreen.o \
	dij.o ei.o finteg.o fitwv.o gaussq.o ghost.o ginteg.o \
	guess.o optinvert.o optim.o optsolve.o conmaxdriver.o\
	xroots.o kpot.o kerker.o lambda.o newmin.o fnset2.o tmsub.o tmpot.o \
        kcomp.o hfinvdesc.o hfsmooth.o

cobjects = do_ps.o do_ke.o

all_objects = $(fobjects) $(cobjects)

this_dir = pseudo

.F.o:
	$(FC) $(FFLAGS) $(OPTIONS)  -I../include $*.F -c


all:  $(all_objects)
	-rm -f objects
	for i in $(all_objects); do echo $(this_dir)/$$i >> objects; done

clean:
	rm -f $(fobjects) $(cobjects) libpseudo.a core* *~ *.f.f objects

distclean: clean


