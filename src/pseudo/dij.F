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
      subroutine dij
c This subroutine calculates the definite integral from 0 to qc of
c q**4 * b(i,...) * b(j,...) dq by gaussian quadrature.
c That is why b(i,...) is actually a 2-d array stored for each 
c b(i,...) at all the gaussian quadrature points.

      implicit double precision(a-h,o-z)
      
#include "fortdim.h"

c -------------------------------------------------------------------------
c     Internal (Fortran only) common blocks
c -------------------------------------------------------------------------
      common/gauss/ xquad(nquad0),wquad(nquad0)
      common/quads/ qq(nquad0)
      common/b/ b(numfn0,nquad0)
      common/cuts/ qc,rc
      common/d/ d(numfn0,numfn0)
      common /numfn/ numfn
c -------------------------------------------------------------------------

      do i = 1,numfn
         do j = i,numfn
            sum = 0.0
            do k = 1,nquad0
               sum = sum + wquad(k) * b(i,k) * b(j,k) * qq(k)**4
            enddo
            d(i,j) = sum * qc/2.0
            d(j,i) = d(i,j)
         enddo
      enddo

      return
      end

