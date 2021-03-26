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
      subroutine hrtree(np,h,r,q,rvh)
      implicit real*8  (a-h,o-z)

#include "fortdim.h"      

c compute hartree term of scf potential from charge dens. q
c this version uses corrected trapezoid for indef. integrals
c this version improves implementation of trap. rule corrections, using
c asymptotic expansions and explicit cancellations.
c np=number of grid points; h=spacing of uniform y grid where r(y)=const
c exp(y) (makes r a logarithmic grid); r = r(y); q = charge dens on same
c grid; and finally, rvh = r(y) * vhartree(y).
c integrations use asymptotic forms near r=0.

      dimension q(npdm),r(npdm),rvh(npdm)
      iout = 8
      h2 = h+h
      t1 = -h*(1.+h)/2880.
      t2 =  h*(2.+h)/ 360.
      t3 = -h*(67.5-h*h)/720.
      t4 =  h*(2.-h)/ 360.
      t5 = -h*(1.-h)/2880.
      s1 = q(1)*r(1)
      s2 = q(2)*r(2)
      s3 = q(3)*r(3)
      s4 = q(4)*r(4)
      rvh(1) = -s1/(6.*h)
      rvh(2) = -s2/(6.*h)
c a and b accumulate corrected trap. from origin (analytic there)outward
      a = s2  *(0.5 + 1./(3.*h) +  h*(0.25-0.0375*h*h)        )
      b = q(2)*(0.5 + 1./(2.*h) +  h*( 15.-         h*h)/90.  )
c
      m = np - 2
      do 10 j=3,m
            s5 = q(j+2)*r(j+2)
            rvh(j) =   a + t1*s1+t2*s2+t3*s3+t4*s4+t5*s5 - r(j) * b
            a = a + s3
            b = b + q(j)
            s1 = s2
            s2 = s3
            s3 = s4
            s4 = s5
   10 continue
      izero = 1
      r1 = r(np-1)
      r2 = r( np )
      if(q(np-1).le.1.e-25.or.q(np).le.1.e-25) izero = 0
      corr = 0.
      decay = 0.
      if(izero.ne.0) decay = log(q(np-1)/q(np))/(r2-r1)
      d1 = decay*r1
      corr = - h*s3*(3. - 0.05*(1.-d1*(6.-3.*d1))*h*h)/36.
      rvh(np-1) = a - r1*b + corr
      a = a + s3
      b = b + q(np-1)
      d2 = decay*r2
      corr = - h*s4*(3. - 0.05*(1.-d2*(6.-3.*d2))*h*h)/36.
      rvh( np ) = a - r2*b + corr
      bcorr = h*s4*decay*(3.-0.05*(1.-d2*(3.-d2))  *h*h)/36.
      b = b + 0.5*q(np) + bcorr
      do 20 i=1,np
   20 rvh(i) = h2 * (rvh(i) + b * r(i))
c following lines were used for testing integrals.
c     acorr = corr + r2 * bcorr
c     a = ( a + 0.5*s4 + acorr ) * h
c     decay = (0.5*decay)**2
c     write(iout,1000) a,decay,q(np)
c1000 format(1x,'hrtreed finds ',f15.11,' electrons from charge',
c    $' density integral.',/,1x,'(.5*decay)**2 =',e20.7,
c    $  2x,'r**2*rho(r) at max=',e12.3)
      return
      end
