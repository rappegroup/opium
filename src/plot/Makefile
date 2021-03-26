# Makefile for somthin'

SHELL=/bin/sh

include ../../makevars

cobjects =  do_plot.o do_pccplot.o do_logplot.o do_wplot.o do_vplot.o do_qplot.o do_keplot.o
fobjects = btrans.o

all_objects = $(fobjects) $(cobjects)

this_dir = plot

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
	rm -f $(cobjects) $(fobjects)  *~ core* objects

distclean: clean

