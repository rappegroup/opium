# This is the Makefile for 'FlexiLib' distributed with OPIUM
#

SHELL=/bin/sh

include ../../makevars

this_dir = FlexiLib

objects = flexi.o

all_objects = $(objects)

all: $(all_objects)
	-rm -f objects
	for i in $(all_objects); do echo $(this_dir)/$$i >> objects; done

flexi.o: flexi.h

clean: 
	rm -f $(objects) libflexi.a core* *~ objects

distclean: clean

