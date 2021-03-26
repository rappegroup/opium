SHELL=/bin/sh

include ../makevars

subdirs=main atom pseudo output plot FlexiLib uniPPlib lapack

.PHONEY=all-subdirs

SHELL=/bin/sh

all: all-subdirs opium

opium: 
	$(LD) -o ../opium `for i in $(subdirs); do cat $$i/objects; done` $(LDFLAGS)

all-subdirs:
	for i in $(subdirs); do cd $$i && $(MAKE) || exit 1; cd ..; done

clean:
	for i in $(subdirs); do cd $$i && $(MAKE) clean; cd ..; done

distclean:
	for i in $(subdirs); do \
		test -f $$i/Makefile && cd $$i && $(MAKE) distclean; cd ..; done
	-rm -rf autom4te.cache
	-rm -f opium core* config.status config.log *~

