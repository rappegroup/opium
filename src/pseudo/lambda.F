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
      subroutine lambda(rc,numfn)
      implicit double precision(a-h,o-z)

#include "fortdim.h"

      common /rlm/ rlmat(3,numfn0),rlvec(numfn0)

      dimension rc(numfn0)

      rc(numfn-1)=0.0
      rc(numfn)=0.0

      alpha = rlmat(2,numfn)*rlmat(1,numfn-1)
     $     - rlmat(2,numfn-1)*rlmat(1,numfn)
      beta  = (-rlmat(2,numfn)*rlvec(1)+rlmat(1,numfn)*rlvec(2))

      do i=1,numfn-2
         rc(numfn-1) = rc(numfn-1) + rc(i)*(rlmat(1,numfn)*rlmat(2,i)
     $        - rlmat(2,numfn)*rlmat(1,i))
      enddo

      rc(numfn-1)=(rc(numfn-1)+beta)/alpha

      do i=1,numfn-1
         rc(numfn) = rc(numfn) + rlmat(1,i)*rc(i)
      enddo
      rc(numfn)=(-rlvec(1)-rc(numfn))/rlmat(1,numfn)

      end

               
