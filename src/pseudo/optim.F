c
c Copyright (c) 1998-2012 The OPIUM Group
c
c This program is free software; you can redistribute it and/or modify
c it under the terms of the GNU General Public License as published by
c the Free Software Foundation; either version 2 of the License, or
c (at your option) any later version.
c
c This program is distributed in the hope that it will be useful,
c but WITHOUT ANY WARRANTY; without even the implied warranty of
c MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
c GNU General Public License for more details.
c
c You should have received a copy of the GNU General Public License
c along with this program; if not, write to the Free Software
c Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
c
c
      subroutine optim(ixc,irel,exccut)
      
c*************************************************************************
c generates optimized pseudopotentials [Rappe, et al. PRB 41, 1227 (1990)]
c*************************************************************************

      implicit double precision(a-h,o-z)
      
#include "fortdim.h"

c -------------------------------------------------------------------------
c     External (shared between C and Fortran) common blocks
c -------------------------------------------------------------------------

      common /totpot/ rvcore(npdm,n0),rvps(npdm,n0),rvcoul(npdm)

      common /consts/ etol,vtol,maxit,isoft
      common /grid/ h,r1,z,r(npdm),np

      common /filenames/ file_log

      common /rpcc/ rpcc,rpccz

      common /rscore/ rscore(npdm),rdd(npdm),rddd(npdm),rscoretot(npdm),
     $     rsval(npdm)

      common /aorb/ ncore,nval,norb,nlm(n0),
     $     no(n0),lo(n0),nmax(n0),maxim
      common /adat/ wnl(n0),en(n0),so(n0),xion
      common /aval/ rcall(n0),rvap(n0),rnorm(n0),ibd(n0),etot

      common /wfn/ rnl(npdm,n0)

      common /optparam/ qcl(n0),nbl(n0)
c -------------------------------------------------------------------------

c -------------------------------------------------------------------------
c     Internal (Fortran only) common blocks
c -------------------------------------------------------------------------
      common /transum/ transumry
      common /angm/ ll
      common /cuts/ qc,rc
      common /nnn/ nnn
      common /numfn/ numfn
c -------------------------------------------------------------------------

      character*80 file_log
      
      open(unit=7,file=file_log,form='formatted',access='append')
      
      write (7,*) 'Optimized Pseudopotential Generation ...'
      transumry=0
      
      do i=1,np
         rsval(i)=0.0
      enddo

      zeff=xion
      do nnn=1,nval
         zeff=zeff+wnl(nnn)
      enddo

      do nnn = 1,nval
         write(7,*) 
         write(7,*) '=================='
         write(7,*) 'zeff=',zeff

         ic1=nlm(nnn)/100
         ll=(nlm(nnn)-ic1*100)/10
         rc = rcall(nnn)
         qc = qcl(nnn)
         numfn = nbl(nnn)
         lo(nnn)=ll
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
            
            call optsolve(ixc)
            
         else
            write(7,*) "Using Kerker method for unbound state"
            call kpot(ixc)
         endif

      enddo

      if (ixc.ge.0) then
         call descreen(ixc,exccut)
         call ghost
      else          
c     do this after in do_ps
         call hfinvdesc(ixc,irel,zeff)
      endif

 9000 format(1x,'Pseudizing state : |',i3,'>')
 9001 format(1x,'eigenvalue       :',f10.6)
 9002 format(1x,'qc               :',f10.6)
 9003 format(1x,'# bessel fxns    :',i3)
 9011 format(1x,'Total Convergence Error: ',f10.6,1x,
     $       'mRy',3x,f10.6,1x,'meV') 
      

      close(unit=7)
      
      return      
      end
