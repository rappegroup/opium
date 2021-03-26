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
      subroutine finteg
      implicit double precision (a-h,o-z)
c Here we integrate  phi(r)*delsq(phi(r)) from rc to infinity.
c We do this in one interval (relatively low accuracy) since it is
c a costly calculation and is relatively unimportant to the calc.
      
#include "fortdim.h"
c -------------------------------------------------------------------------
c     External (shared between C and Fortran) common blocks
c -------------------------------------------------------------------------
      common /grid/ h,r1,z,r(npdm),np
      common /aorb/ ncore,nval,norb,nlm(n0),
     $     no(n0),lo(n0),nmax(n0),maxim
      common /adat/ wnl(n0),en(n0),so(n0),xion
      common /aval/ rcall(n0),rvap(n0),rnorm(n0),ibd(n0),etot
      common /wfn/ rnl(npdm,n0)
c -------------------------------------------------------------------------

c -------------------------------------------------------------------------
c     Internal (Fortran only) common blocks
c -------------------------------------------------------------------------
      common /gauss/ xquad(nquad0),wquad(nquad0)
      common /nnn/ nnn
      common /cuts/ qc,rc
      common /angm/ ll
      common /f/ fint
c -------------------------------------------------------------------------
      dimension xf(nquad0),f(npdm)

      nord = 10
      do i = 1,np
         f(i) = rnl(i,nnn)
      enddo

      nint = 1
      width = (r(nmax(nnn)) - rc)/float(nint)
      zsum = 0.0

      do i = 1,nint
         xlend = float(i-1) * width + rc
         do j = 1,nquad0
            xf(j) = (xquad(j) + 1.0)*width/2.0 + xlend

            vv = val(f,r,np,xf(j),nord)
            v2 = val2(f,r,np,xf(j),nord)
            v3 = val3(f,r,np,xf(j),nord)

            xl = float(ll * (ll + 1))
            dl = v3+2.0*v2/xf(j)-xl*vv/xf(j)**2
            zsum=zsum+vv*dl*xf(j)**2*wquad(j)
         enddo
      enddo

      zsum = zsum * width/2.0
      fint = zsum

      return
      end
