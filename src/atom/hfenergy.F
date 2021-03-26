c--------------------------------------------------------------------
c The following code is adapted from:
c
c A GENERAL HARTREE-FOCK PROGRAM
c
c     by C. Froese Fischer
c        Vanderbilt University
c        Nashville, TN 37235 USA
c
c     April, 1987
c     Published in Computer Physics Communications, 43 (1987) 355
c
c     .... Modified for for arbitrarily filled f-shells
c     by Gediminas GAIGALAS 
c        Institute of Theoretical Physics and Astronomy
c        A. Gostauto str. 12
c        Vilnius, 2600, Lithuania
c       
c     January, 1996
c
c     Pulished in Computer Physics Communications, 98 (1996) 255.
c
c
c
c
c Obtained from http://atoms.vuse.vanderbilt.edu/Elements/CompMeth/hf.f
c-----------------------------------------------------------------------


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
*
*     ------------------------------------------------------------------
*                       e n e r g y
*     ------------------------------------------------------------------
*     
*     determines the position of the electron in the electron list
*     
      subroutine hfenergy(ipsp,etotal,ehart,exx,ekn)
      implicit double precision(a-h,o-z)
#include "fortdim.h"
      common /wfn/ rnl(npdm,n0)
      common /grid/ h,r1,z,r(npdm),np
      common /totpot/ rvcore(npdm,n0),rvps(npdm,n0),rvcoul(npdm)
      common /aorb/ ncore,nval,norb,nlm(n0),
     $     no(n0),lo(n0),nmax(n0),maxim
      common /adat/ wnl(n0),en(n0),so(n0),xion
      common /edi/ e(n0,n0)
      common /ahf/ s(n0),az(n0),acc(n0),dpm(n0)

      dimension rm1(n0),ek(n0),yk(npdm)

      ekn=0.0
      etotal = 0.d0
      exx = 0.d0
      ehart = 0.d0
      eloc = 0.d0

      do  i = 1,norb
         ek(i) = -0.5*hl(i,i)
         if (ipsp.eq.0) then
            rm1(i) = z*quadr(i,i,-1)
         else
            do ii=1,np
               yk(ii)=rvcore(ii,i)/r(ii)
            enddo
            rm1(i) = quads2(yk,i,i,0)
         endif
         ekn = ekn+wnl(i)*(ek(i) + rm1(i))
      enddo

      do i = 1,norb
         etotal = etotal + wnl(i)*(ek(i))

         do j = 1,i
            do k = 0,2*min0(lo(i),lo(j)),2
               c = aaa(i,j,k)*wnl(i)
               if (i .eq. j) c = c/2.d0
               if (abs(c).ne.0.d0) etotal = etotal + c*fk(i,j,k)
               if (abs(c).ne.0.d0) ehart = ehart + c*fk(i,j,k)               
            enddo
         enddo

         do j = 1,i-1
            do k = abs(lo(i)-lo(j)),lo(i)+lo(j),2
               c = bbb(i,j,k)*wnl(i)
               if (abs(c).ne.0.d0) etotal=etotal+c*gk(i,j,k)
               if (abs(c).ne.0.d0) exx=exx+c*gk(i,j,k)
            enddo
         enddo
      enddo

      end
      
