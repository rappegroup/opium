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
      subroutine ai
      
c     *************************************************************************
c     This subroutine computes the definite integral from 0 to rc of
c     r**2 * (jl(qr))**2 dr analytically.
c     Now if q<0 integrand is r**2 * (il(abs(q)r))**2 , where
c     il(x) = i**-l * jl(ix)
c     *************************************************************************

      implicit double precision(a-h,o-z)
      
#include "fortdim.h"

      common /a/ a(numfn0)
      common /cuts/ qc,rc
      common /roots/ xroot(numfn0)
      common /angm/ ll
      common /numfn/ numfn
      common /bs/ bs(numfn0)
      common /bd/ bd(numfn0)
      
      do i = 1,numfn
         x = xroot(i) * rc
         xs = 1.0
         if (xroot(i).lt.0.0) xs = -1.0
         bs(i) = besfn(x,ll)
         bd(i) = besder(x,ll)
         t1 = 1.0 - xs * dfloat(ll * (ll + 1))/x/x
         a(i) = (t1 * bs(i)**2 + xs * bd(i) * bs(i)/abs(x) + 
     $        xs * bd(i)**2) * rc**3/2.0
      enddo
      
      return
      end
      
