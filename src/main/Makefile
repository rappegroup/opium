# Makefile for somthin'

SHELL=/bin/sh

include ../../makevars

objects = opium.o nlm.o parameter.o

all_objects = $(objects)

this_dir = main

all: $(all_objects)
	-rm -f objects
	for i in $(all_objects); do echo $(this_dir)/$$i >> objects; done

clean:
	rm -f $(objects) libmain.a core* *~ *.f.f objects

distclean: clean


