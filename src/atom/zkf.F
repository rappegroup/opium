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
      subroutine zkf(yk,i,j,k)
     
c     *************************************************************************
c     get the zk function
c     *************************************************************************
     
      implicit double precision (a-h,o-z)

#include "fortdim.h"
      common /grid/ h,r1,z,r(npdm),np
      common /aorb/ ncore,nval,norb,nlm(n0),
     $     no(n0),lo(n0),nmax(n0),maxim
      common /adat/ wnl(n0),en(n0),so(n0),xion
      common /aval/ rcall(n0),rvap(n0),rnorm(n0),ibd(n0),etot
      common /wfn/ rnl(npdm,n0)

      dimension yk(npdm)

      do in=1,np
         yk(in)=0.d0
      enddo

      den=lo(i)+lo(j)+3+k
      fact = (1.d0/(lo(i)+1) + 1.d0/(lo(j)+1))/(den + 1.d0)
      a = dexp(-h)**k
      a2 = a*a
      h90 = h/90.d0
      h3=h/3.d0
      a3 = a2*a*h90
      ai = h90/a
      an = 114.d0*a*h90
      a34 = 34.d0*h90
      f1 = r(1)**2*rnl(1,i)*rnl(1,j)
      f2 = r(2)**2*rnl(2,i)*rnl(2,j)
      f3 = r(3)**2*rnl(3,i)*rnl(3,j)
      f4 = r(4)**2*rnl(4,i)*rnl(4,j)
      yk(1) = f1*(1.d0 + z*r(1)*fact)/den
      yk(2) = f2*(1.d0 + z*r(2)*fact)/den
      yk(3) = yk(1)*a2 + h3*(f3 + 4.d0*a*f2 + a2*f1)
      mx = (min(nmax(i),nmax(j))/2)*2
      do m = 5,mx
         f5 = (r(m)**2*rnl(m,i))*rnl(m,j)
         yk(m-1) = yk(m-3)*a2 + ( an*f3 + a34*(f4+a2*f2)-f5*ai-f1*a3)
         f1 = f2
         f2 = f3
         f3 = f4
         f4 = f5
      enddo

      m1 = mx - 1
      if (iabs(i-j)  +  iabs(k) .eq. 0) then
*     
*     *****  for y0(i,i) set the limit to 1 and remove oscillations
*     *****  introduced by the use of simpson's rule


         m2 = m1 - 1
         c1 = 1.d0 - yk(m1)
         c2 = 1.d0 - yk(m2)
         do m = 1,m1,2
            yk(m) = yk(m) + c1
            yk(m+1) = yk(m+1) + c2
         enddo 
      endif

      do m = m1+1,np
         yk(m) = a*yk(m-1)
      enddo

      return
      end

      
