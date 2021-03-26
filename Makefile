SHELL=/bin/sh

all:
	cd src && $(MAKE)  

clean:
	cd src && $(MAKE) clean 

distclean:
	test -f src/Makefile && cd src && $(MAKE) distclean
	rm -fr makevars
	rm -rf au*.cache fpp.fixed
	rm -f opium core* config.status config.log *~
	rm -f cu.*
	rm -f fe.*
	rm -f h.*
	rm -f o.*
	rm -f ti.*
	rm -f c.*
	rm -f pt.*
	rm -f rh.*
	rm -f ca.*
	rm -f *.agr
	rm -f *.par
	rm -f b.*
	rm -f au.*
	rm -f ne.*
	rm -f li.*
	rm -f s.*
	rm -f y.*

