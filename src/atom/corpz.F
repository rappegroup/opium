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
      subroutine corpz(rs,ec,uc)
      implicit double precision (a-h,o-z)
      
#include "fortdim.h"
      parameter(thrd=1.d0/3.d0,thrd2=2.d0*thrd)
      parameter(pi32=29.608813203268075856503472999628d0)
      parameter(pi=3.1415926535897932384626433832795d0)
      parameter(alpha=1.91915829267751300662482032624669d0)

      b11 = 1.0529
      b22 = 0.3334
      aa = 0.0622
      bb =-0.0960
      cc = 0.0040
      dd =-0.0232
      ga =-0.2846
      t11 = 7./6.
      t22 = (4./3.)*b22      

      if (rs.lt.1) goto 506

c     Low density correlation
      
      root = sqrt(rs) * b11
      xd = 1. + root + b22 * rs
      ec = ga/xd
      uc = (1.+t11*root+t22*rs)*ec/xd
      goto 509
      
c     High density correlation
      
 506  dl = log(rs)
      t33 = cc*rs
      t44 = dd*rs
      t55 = t33*dl
      ec = aa*dl +bb + t55 + t44
      uc = ec - (aa+ t55 + t44 + t33)*thrd
      
 509  continue

      ec=ec*0.5d0
      uc=uc*0.5d0

      return
      end
