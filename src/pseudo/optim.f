      subroutine optim(ixc)
      
c*************************************************************************
c generates optimized pseudopotentials [Rappe, et al. PRB 41, 1227 (1990)]
c*************************************************************************

      implicit double precision(a-h,o-z)
      
#include "PARAMOPT"

c -------------------------------------------------------------------------
c     External (shared between C and Fortran) common blocks
c -------------------------------------------------------------------------
      common /grid/ h,r1,z,r(npdm),np
      common /atomic/ xion,rnl(npdm,n0),nlm(n0),wnl(n0),en(n0),norb
      common /rcall/ rcall(nvale0)
      common /filenames/ file_log
      common /lparam/ qcl(10),nbl(10)
      common /np/ ncores,nvales
      common /ibound/ ibd(n0)
c -------------------------------------------------------------------------

c -------------------------------------------------------------------------
c     Internal (Fortran only) common blocks
c -------------------------------------------------------------------------
      common /transum/ transumry
      common /angm/ ll
      common /cuts/ qc,rc
      common /nnn/ nnn
      common /numfn/ numfn
      common /atom3/ rsvale(npdm)
c -------------------------------------------------------------------------

      character*80 file_log
      
      open(unit=7,file=file_log,form='formatted',access='append')
      
      write (7,*) 'Optimized Pseudopotential Generation'
      transumry=0
      
      do i=1,np
         rsvale(i)=0.0
      enddo
      
      do nnn = 1,nvales

         write(7,*) 
         write(7,*) '=================='

         ic1=nlm(nnn)/100
         ll=(nlm(nnn)-ic1*100)/10
         rc = rcall(nnn)
         qc = qcl(nnn)
         numfn = nbl(nnn)

         if (numfn.gt.numfn0) then
            write (7,*) "Too many basis functions requested"
            write (7,*) "Execution terminates."
            stop
         endif
         write(7,9000) nlm(nnn) 
         write(7,9001) en(nnn)
         write(7,9002) qc
         write(7,9003) numfn
         
         call fitwv
         if (ibd(nnn).gt.0) then

            call xroots
         
            call ai
            call gaussq
            call biq
            call cq
            call dij
            call ei
            call finteg
            call ginteg
            
            call optsolve
            
         else
            call kpot
         endif

      enddo

c      write(7,*) 
c      write(7,*) '======Potentials created======'
c      write(7,9011) transumry*1000.0,transumry*13.6058*1000.0
      
      call descreen(ixc)
      call ghost


 9000 format(1x,'Pseudizing state : |',i3,'>')
 9001 format(1x,'eigenvalue       :',f10.6)
 9002 format(1x,'qc               :',f10.6)
 9003 format(1x,'# bessel fxns    :',i3)
 9011 format(1x,'Total Convergence Error: ',f10.6,1x,
     $       'mRy',3x,f10.6,1x,'meV') 
      
      
      close(unit=7)
      
      return      
      end
