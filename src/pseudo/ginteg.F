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
      subroutine ginteg
      implicit double precision (a-h,o-z)
c This code calculates the integral from 0 to qc of
c q**4 * c(q)**2 using gaussian quadrature.

#include "fortdim.h"

c -------------------------------------------------------------------------
c     Internal (Fortran only) common blocks
c -------------------------------------------------------------------------
      common/gauss/ xquad(nquad0),wquad(nquad0)
      common/quads/ qq(nquad0)
      common/cuts/ qc,rc
      common /c/ c(nquad0)
      common /g/ gint
c -------------------------------------------------------------------------

      sum = 0.0

      do j = 1,nquad0
         sum = sum + c(j)**2 * qq(j)**4 * wquad(j)
      enddo

      gint = sum*qc/2.0

      return
      end
