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
      subroutine ykf(yk,i,j,k)
     
c     *************************************************************************
c     get the yk function
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

      call zkf(yk,i,j,k)

      a =  dexp(-h)**(k+1)
      ccc = 2*k+1
      a2 = a*a
      h3=h/3.d0
      h90 = ccc*h3/30.d0
      a3 = a2*a*h90
      ai = h90/a
      an = 114.d0*a*h90
      a34 = 34.d0*h90
      mx = (min(nmax(i),nmax(j))/2)*2
      f1 = yk(mx)*dexp(-h)**k
      f2 = yk(mx)
      f3 = yk(mx-1)
      f4 = yk(mx-2)

      do  m = mx-2,2,-1
         f5 = yk(m-1)
         yk(m) = yk(m+2)*a2 + ( an*f3 + a34*(f4+a2*f2)-f5*ai-f1*a3)
         f1 = f2
         f2 = f3
         f3 = f4
         f4 = f5
      enddo
      yk(1) = yk(3)*a2+ccc*h3*(f4 + 4.d0*a*f3 + a2*f2)

      return
      end

