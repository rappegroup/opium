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
      subroutine fnset2(ipr,rc,numfn)
      implicit double precision(a-h,o-z)

#include "fortdim.h"

      common /rke/ rkmat(numfn0,numfn0),rkvec(numfn0)
      common /rn/ rnmat(numfn0,numfn0),rnvec(numfn0)
      common /rlm/ rlmat(3,numfn0),rlvec(numfn0)
      common /rconst/ rkcon,rncon,rlcon(3)
      common /re/ rkin,rnorm
      common /emat/ emat(3,3),rfl(numfn0)
      common /rf/ afor(3),rlam(3)
      common /rb/ rb(numfn0,numfn0)
      common /ncon/ ncon

      dimension rc(numfn0)

      rkin=0
      rnorm=0
      do i=1,numfn
         do j=1,numfn
            rkin=rkin+rc(i)*rkmat(i,j)*rc(j) 
            rnorm=rnorm+rc(i)*rnmat(i,j)*rc(j)
         enddo
         rkin=rkin+rkvec(i)*rc(i)
      enddo

      do i=1,ncon
         rlcon(i)=0.0
         do j=1,numfn
            rlcon(i)=rlcon(i)+rlmat(i,j)*rc(j)
         enddo
      enddo

      if (ipr.ne.0) then
         write(7,9004)
         write(7,9005) (rc(i),i=1,numfn)
         write(7,9004)
         write(7,9000) rkin+rkcon
         write(7,9001) rnorm+rncon
         write(7,9002) rlcon(1)+rlvec(1)
         write(7,9003) rlcon(2)+rlvec(2)
         if (ncon.eq.3) write(7,9007) rlcon(3)+rlvec(3)
         write(7,9004)
      endif
 
 9000 format(1x,'Resid KE (Ry)     :    ',2f16.10)
 9001 format(1x,'Norm error        :    ',e10.3)
 9002 format(1x,'Continuity error  :    ',e10.3)
 9003 format(1x,'Curvature error   :    ',e10.3)
 9007 format(1x,'Sum       error   :    ',2f10.6)
 9005 format(1x,'Coeff:    ',20f10.6)
 9006 format(1x,'Coeff SUM         :    ',f10.6)
 9004 format(1x,'-------------------------------------',
     $     '--------------------------')
      
      return
      end
