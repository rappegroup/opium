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
      subroutine logder(en,z,l,isoft,rv,g,r,np,h,rlog,dl,nnt1,nr,p,
     $     dlwf,ider,de,qp,dld)
     
      implicit double precision (a-h,o-z)
      
#include "fortdim.h"

c     -------------------------------------------------------------------------
c     Internal (Fortran only) common blocks
c     -------------------------------------------------------------------------
      common /pofy/ pofy
c     -------------------------------------------------------------------------
      
      dimension rv(npdm),r(npdm),dl(n0),p(npdm),dld(npdm),g(npdm)
      
c     *************************************************************************
c     local variables
c     *************************************************************************

      dimension w(4),wi(4)
      
c     ensure nr is multiple of 2 integer and nr >= 4

      nr = 2*(nr/2)
      if(nr.lt.4) nr = 4
      
      nrh = nr/2

c     ensure that rlog is within range of grid
      if(rlog.ge.r(np-nrh-1)) then
        write(7,1000) rlog,r(np)
 1000   format(' ****rlog outside range:rlog,r(np)=',2e20.7)
        return
      endif
   
      ylog = 1. + log(rlog/r(1))/h
      ilog = int(ylog+0.000001)
      imax = ilog + nrh + 1
      imin = ilog - nrh
      imp1 = imin + 1
      i1 = ilog - 1
      i2 = ilog
      i3 = ilog + 1
      i4 = ilog + 2
      yz = float(ilog+ilog+1)*0.5
      dy = ylog - yz
      x2 = dy
      x3 = dy**2
      x4 = dy**3
      w(1) = (- 3.      + 2.*x2+12.*x3- 8.*x4)/48.
      w(2) = ( 27.      -54.*x2-12.*x3+24.*x4)/48.
      w(3) = ( 27.      +54.*x2-12.*x3-24.*x4)/48.
      w(4) = (- 3.      - 2.*x2+12.*x3+ 8.*x4)/48.
c
      lp = l + 1
      e = en
      ie = 3
      if(ider.eq.0) goto 95
      x1 = dy + 0.5
      x2 = (dy*dy-0.25)*0.5
      x3 = (dy*dy*dy+0.125)/3.
      x4 = (dy*dy*dy*dy-0.0625)*0.25
      wi(1) = (- 3.*x1+ 2.*x2+12.*x3- 8.*x4)/48.
      wi(2) = ( 27.*x1-54.*x2-12.*x3+24.*x4)/48.
      wi(3) = ( 27.*x1+54.*x2-12.*x3-24.*x4)/48.
      wi(4) = (- 3.*x1- 2.*x2+12.*x3+ 8.*x4)/48.
      ie = 1
   94 e = en + float(ie-3)*de
   95 continue
   
c     power series near r = 0 (remove z/r for soft core vps)
c     **********************************************************************
c     beginning of changes made by njr
c       addition of inhomogeneous term g(i)
c       local potential was taken as x*r**-1 + y + w*r
c       non-local potential was taken as g1*r**(l+2.5)+g2*r**(l+3.5)+
c       g3*r**(l+4.5)
c     **********************************************************************
      ztemp = z
      if(isoft.eq.1) ztemp = 0.
      z2 = ztemp+ztemp
c     z2 = -x
c     a = 1 (assumption of wavefunction leading term coefficient for
c            inhomogeneous solution only)
      a0 = ( rv(1) + z2)/r(1)
c     a0 = y + w * r(1)
      a1 = ((rv(2) + z2)/r(2) - a0)/(r(2)-r(1))
c     a1 = w
      b0 = g(1)/(r(1)**(l + 2.5))
      b1 = g(2)/(r(2)**(l + 2.5))
      g3 = (b1 - b0)/(r(2) - r(1))
      g2 = b0 - g3*r(1)
      g1 = 0
      at = e - a0 + a1 * r(1)
c     at = e - y
      x3 = float(l+l+1)
      d1 = 1.0
c     d1 = a/a
      d2 = (-ztemp*d1+g1)/float(lp)
c     d2 = b/a
      d3 = (-(z2*d2+at*d1)+g2)/float(4*l+6)
c     d3 = c/a
      d4 = (-(z2*d3+at*d2)-a1*d1+g3)/float(6*l+12)
c     d4 = d/a
      p(1) = d1 + r(1)*(d2+r(1)*(d3+r(1)*d4))
      p(2) = d1 + r(2)*(d2+r(2)*(d3+r(2)*d4))
      j = 4
  100 ap = float(j)
      d5 = -(z2*d4 + at*d3 - a1*d2)/(ap*(ap+x3))
      ncross = -1
      do i=1,2
        d2 = d5 * r(i)**j
        if(abs(d2/p(i)).gt.1.e-13) ncross = 1
        p(i) = p(i) + d2
      enddo
      if(ncross.ge.0) then
        d2 = d3
        d3 = d4
        d4 = d5
        j = j + 1
        if(j.le.100) goto 100
        write(7,1010) e,l,rlog
 1010   format('0****logder fails to converge at origin: e=',e20.7,2x,
     $         'l=',i3,2x,'rlog=',e20.7 )
        return
      endif
      val = float(l) + 0.5
      p(1) = p(1)*r(1)**val
      p(2) = p(2)*r(2)**val
c     outward integration by noumerov method
      x2 = val**2
      h2 = (h**2)/12.
      d1 = ((rv(1)-e*r(1))*r(1)+x2)*h2*p(1)+h2*g(1)
      d2 = ((rv(2)-e*r(2))*r(2)+x2)*h2*p(2)+h2*g(2)
      do i=3,imax
        d3 = ((rv(i)-e*r(i))*r(i)+x2)*h2
        p(i) = (p(i-1)+p(i-1)-p(i-2)+10.*d2+d1+h2*g(i))/(1.-d3)
        d1 = d2
        d2 = d3 * p(i)+h2*g(i)
      enddo
      d1 = ((rv(imin)-e*r(imin))*r(imin)+x2)*h2*p(imin)+h2*g(imin)
      d2 = ((rv(imp1)-e*r(imp1))*r(imp1)+x2)*h2*p(imp1)+h2*g(imp1)
      do j=1,nr
        i = imin + j
        ip = i + 1
        d3 = ((rv(ip)-e*r(ip))*r(ip)+x2)*h2*p(ip)+h2*g(ip)
        dl(j) = (0.5*(p(ip)-p(i-1))-(d3-d1))/h
        d1 = d2
        d2 = d3
      enddo
c     **********************************************************************
c     end of changes made by njr
c     **********************************************************************

c     **********************************************************************
c     cubic interpolation of both p(y) and dpdy on uniform y grid
c     **********************************************************************
      j = 0
      pofy = 0.
      dpdy = 0.
      do i=i1,i4
        j = j + 1
        pofy  = pofy + w(j)*p(i)
        dpdy  = dpdy + w(j)*dl(nrh+j-2)
      enddo
      if(abs(pofy).lt.1.d-15) pofy = sign(1.d-15,pofy)
      dlp  = dpdy/pofy
      dlwf = (0.5 + dlp)/rlog
      if(ie.eq.3) then
        do j=1,nr
          i = imin + j
          dl(j) = (0.5 + dl(j)/p(i))/r(i)
        enddo
      endif
      if(ider.eq.0) return
      do i=1,i4
        p(i) = p(i)*p(i)*r(i)
      enddo
      q = float(lp+lp)
      call radin(r,p,0,ilog,h,q)
      q = q + (wi(1)*r(i1)*p(i1) + wi(2)*r(i2)*p(i2) + wi(3)*r(i3)*p(i3)
     $      +  wi(4)*r(i4)*p(i4) )*h
      q = q / (pofy*pofy*rlog)
      dld(ie) = - q
      ie = ie + 1
      if(ie.le.5) goto 94
      
c     **********************************************************************
c     energy derivative evaluated by 5-pt symmetric interpolation
c     **********************************************************************
      qp = (dld(1)-dld(5)+8.*(dld(4)-dld(2)))/(12.*de)

      return
      end
