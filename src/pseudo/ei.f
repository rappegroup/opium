c
c Copyright (c) 1998-2004 The OPIUM Group
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
      subroutine ei
      implicit double precision (a-h,o-z)

c This code calculates the integral from 0 to qc of
c q**4 * bi(q) * c(q) using gaussian quadrature.

#include "PARAMOPT"
c -------------------------------------------------------------------------
c     Internal (Fortran only) common blocks
c -------------------------------------------------------------------------
      common/gauss/ xquad(nquad0),wquad(nquad0)
      common/quads/ qq(nquad0)
      common/cuts/ qc,rc
      common /b/ b(numfn0,nquad0)
      common /c/ c(nquad0)
      common /e/ e(numfn0)
      common /numfn/ numfn
c -------------------------------------------------------------------------

      do i = 1,numfn
         sum = 0.0
         do j = 1,nquad0
            sum = sum + b(i,j) * c(j) * qq(j)**4 * wquad(j)
         enddo
         e(i) = sum * qc/2.0
      enddo

      return
      end
