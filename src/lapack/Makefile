SHELL=/bin/sh

include ../../makevars

objects = dbdsqr.o dcopy.o ddot.o dgebd2.o dgebrd.o dgelq2.o dgelqf.o dgemm.o dgemv.o dgeqr2.o dgeqrf.o dger.o dgesvd.o dlabrd.o dlacpy.o dlamch.o dlange.o dlapy2.o dlarf.o dlarfb.o dlarfg.o dlarft.o dlartg.o dlas2.o dlascl.o dlaset.o dlasq1.o dlasq2.o dlasq3.o dlasq4.o dlasq5.o dlasq6.o dlasr.o dlasrt.o dlassq.o dlasv2.o dlauu2.o dlauum.o dnrm2.o dorg2r.o dorgbr.o dorgl2.o dorglq.o dorgqr.o dorm2r.o dormbr.o dorml2.o dormlq.o dormqr.o dpotf2.o dpotrf.o dpotri.o drot.o dscal.o dswap.o dsyrk.o dtrmm.o dtrmv.o dtrsm.o dtrti2.o dtrtri.o ieeeck.o ilaenv.o xerbla.o lsame.o dgesv.o dgetf2.o dgetrf.o dgetrs.o dlaswp.o idamax.o

all_objects = $(objects)

this_dir = lapack

all: $(all_objects)
	-rm -f objects
	for i in $(all_objects); do echo $(this_dir)/$$i >> objects; done

clean: 
	rm -f $(objects) opium_lapack.a core* *~ objects

distclean: clean

dlamch.o: dlamch.f 
	${FC} -c ${NOOPT} $<
