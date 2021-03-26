# Makefile for src/atom

SHELL=/bin/sh

include ../../makevars

fobjects = cubint.o excorr.o pbe.o ghostnl.o schslrel.o \
	hrtree.o average.o logder.o nllogd.o getpcc.o pratt.o radin.o \
	schsl.o dftsolve.o denkcomp.o interp2.o corpz.o corvwn.o hfcoef.o \
	potl.o xch.o ykf.o zkf.o dfsolve.o soldir.o \
        setuphf.o intdir.o norm.o dentfa.o \
	cwig3j.o potrdf.o ortdat.o cofcon.o dsordf.o akeato.o aprdev.o \
	yzkrdf.o yzkteg.o lagdat.o fdrirk.o etotal.o hfmethd.o hfde.o \
	hfsolve.o hfenergy.o potrdf2.o hfcoef_rel.o hfint.o hfuncs.o

cobjects = do_ae.o do_nl.o do_tc.o do_fc.o writeae.o startae.o relorbae.o nrelorbae.o relorbnl.o nrelorbnl.o

all_objects = $(fobjects) $(cobjects)

this_dir = atom

.F.o:
	$(FC) $(FFLAGS) $(OPTIONS)  -I../include $*.F -c 

all: $(all_objects)
	-rm -f objects
	for i in $(all_objects); do echo $(this_dir)/$$i >> objects; done

clean:
	rm -f $(fobjects) $(cobjects) libatom.a core* *~ *.f.f objects

distclean: clean


