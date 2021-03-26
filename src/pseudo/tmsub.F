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
      subroutine tmsub(ixc,irel,exccut)
      
c     *************************************************************************
c     generate TM style pseudopotentials
c     *************************************************************************

      implicit double precision(a-h,o-z)
      
#include "fortdim.h"

      common /grid/ h,r1,z,r(npdm),np      
      common /lall/ lall(n0)

      common /rscore/ rscore(npdm),rdd(npdm),rddd(npdm),rscoretot(npdm),
     $     rsval(npdm)
      common /aorb/ ncore,nval,norb,nlm(n0),
     $     no(n0),lo(n0),nmax(n0),maxim
      common /adat/ wnl(n0),en(n0),so(n0),xion
      common /aval/ rcall(n0),rvap(n0),rnorm(n0),ibd(n0),etot
      common /wfn/ rnl(npdm,n0)

      common /angm/ ll
      common /cuts/ qc,rc
      common /nnn/ nnn
      common /numfn/ numfn

      common /totpot/ rvcore(npdm,n0),rvps(npdm,n0),rvcoul(npdm)
      common /nonloc/ nonloc,ifnl,numloc
      common /nlall/ nlall(n0)
      
      common /filenames/ file_log      
      character*80 file_log
      
      open(unit=7,file=file_log,form='formatted',access='append')
      
      write (7,*) 'Troullier Martin  pseudopotential ...'

      do i=1,np
         rsval(i)=0.0
      enddo

      zeff=xion
      do nnn=1,nval
         zeff=zeff+wnl(nnn)
      enddo

      do nnn = 1,nval

         ic1=nlm(nnn)/100
         lall(nnn)=(nlm(nnn)-ic1*100)/10
         
         ll = lall(nnn)
         rc = rcall(nnn)

         write(7,9000) nlm(nnn) 
         write(7,9002) rc

c     Here, rnl is r * R_nl(r)         
c     1 - modfiy fitwv for 3rd and 4th deriv and 1st and 2nd deriv of V_scr
         call fitwv

c     Here, rnl is R_nl(r)

         call flush(7)
         call tmpot(ixc)
         call flush(7)
      enddo      


      if (ixc.ge.0) then
         call descreen(ixc,exccut)
         call ghost
      else          
         call hfinvdesc(ixc,irel,zeff)
      endif


 9000 format(1x,'--------------------Pseudizing state: ',
     $       '|',i3,'>',3x,'--------------------')
 9002 format(1x,'rc                       : ',f8.4)
      
      close(unit=7)
      
      return      
      end
