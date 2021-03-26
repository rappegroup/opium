      subroutine pseudo(ixc)
      
c     *************************************************************************
c     generate optimized pseudopotentials [Rappe, et al. PRB 41, 1227 (1990)]
c     *************************************************************************

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
      
      write (7,*) 'Optimized Pseudopotential Generation Program'
      transumry=0
      
      do i=1,np
         rsvale(i)=0.0
      enddo
      
      do nnn = 1,nvales
         
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
         write(7,9001) numfn
         write(7,9002) rc
         write(7,9003) qc


            call fitwv
         if (en(nnn).lt.0.0) then

            call xroots
         
            call ai
            call gaussq
            call biq
            call cq
            call dij
            call ei
            call finteg
            call ginteg
            
            call scr
            
         else
            call kpot
         endif

      enddo

      write(7,*) 
      write(7,*) '======Potentials created======'
      write(7,9011) transumry*1000.0,transumry*13.6058*1000.0
      
      call descreen(ixc)
c      call ghost


 9000 format(1x,'--------------------Pseudizing state: ',
     $       '|',i3,'>',3x,'--------------------')
 9001 format(1x,'# basis functions        : ',i4)
 9002 format(1x,'rc                       : ',f8.4)
 9003 format(1x,'qc                       : ',f8.4)

 9011 format(1x,'Total Convergence Error: ',f10.6,1x,
     $       'mRy',3x,f10.6,1x,'meV') 
      
      
      close(unit=7)
      
      return      
      end
