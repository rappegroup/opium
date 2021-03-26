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
      subroutine corvwn(rs,ec,uc,rho)
      implicit double precision (a-h,o-z)
      
#include "fortdim.h"
      parameter(thrd=1.d0/3.d0,thrd2=2.d0*thrd)
      parameter(pi32=29.608813203268075856503472999628d0)
      parameter(pi=3.1415926535897932384626433832795d0)
      parameter(alpha=1.91915829267751300662482032624669d0)
      
      va=0.0310907
      vb=3.72744
      vc=12.9352
      x0=-0.10498

      xrs=sqrt(rs)

      x2=xrs*xrs
      fx=x2+vb*xrs+vc

      q=4.d0*vc-vb*vb
      q=sqrt(q)

      xb2=2.d0*xrs+vb

      atn=q/xb2
      atn=atan(atn)

      al1=x2/fx
      al1=log(al1)

      al2=(xrs-x0)**2/fx
      al2=log(al2)

      coef=vb*x0/(x0*x0+vb*x0+vc)

      ct1=2.d0*vb/q
      ct2=ct1+4.d0*x0/q

      ec=va*(al1+ct1*atn-coef*(al2+ct2*atn))

      dal=-xb2/fx
      dal1=dal+2.d0/xrs
      dal2=dal+2.d0/(xrs-x0)
      datn=-2.d0*q/(xb2*xb2+q*q)
      dec=va*(dal1+ct1*datn-coef*(dal2+ct2*datn))

      ec=ec
      dec=-dec*xrs/(6.d0*rho)
      uc=ec+rho*dec
     
      ec=ec

      return
      end


