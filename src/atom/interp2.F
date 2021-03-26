      subroutine interp2(nnn,f,f2,f3)
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

c     Interpolate onto the 0 origin grid      

      implicit double precision (a-h,o-z)
      
#include "fortdim.h"

      common /grid/ h,r1,z,r(npdm),np
      common /rgrid/ a,b,rrg(nrmax),rab(nrmax)
      common /nrgrid/ nr

      dimension f(npdm),f2(npdm),fp(npdm),f3(npdm)

      yp1 = 2e33
      ypn = 2e33

      call spline(r,f,np,yp1,ypn,fp)
      
      do j=1,nnn
         call splint(r,f,fp,np,f3(j),yy)
         f2(j)=yy
      enddo

 9595 format(4f20.10)
      return
      end

      SUBROUTINE spline(x,y,n,yp1,ypn,y2)
      INTEGER n,NMAX
      REAL*8 yp1,ypn,x(n),y(n),y2(n)
      PARAMETER (NMAX=5000)

c     Given arrays x(1:n) and y(1:n) containing a tabulated function,
c     i.e., yi = f(xi), with x1 < x2 < .. . < xN, and given values yp1
c     and ypn for the first derivative of the interpolating function at
c     points 1 and n, respectively, this routine returns an array
c     y2(1:n) of length n which contains the second derivatives of the
c     interpolating function at the tabulated points xi. If yp1 and/or
c     ypn are equal to 1 × 1030 or larger, the routine is signaled to
c     set the corresponding boundary condition for a natural spline,
c     with zero second derivative on that boundary.
c     Parameter: NMAX is the largest anticipated value of n.

      INTEGER i,k
      REAL*8 p,qn,sig,un,u(NMAX)
      if (yp1.gt..99e30) then 
         y2(1)=0. 
         u(1)=0.
      else 
         y2(1)=-0.5
         u(1)=(3./(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      endif
      do i=2,n-1 
         sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
         p=sig*y2(i-1)+2.
         y2(i)=(sig-1.)/p
         u(i)=(6.*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))
     $        /(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
      enddo 
      if (ypn.gt..99e30) then 
         qn=0.
         un=0.
      else 
         qn=0.5
         un=(3./(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      endif
      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.)
      do k=n-1,1,-1 
         y2(k)=y2(k)*y2(k+1)+u(k) 
      enddo 
      return
      END
      
      SUBROUTINE splint(xa,ya,y2a,n,x,y)
      INTEGER n
      REAL*8 x,y,xa(n),y2a(n),ya(n)
      INTEGER k,khi,klo
      REAL*8 a,b,h
      klo=1 
      khi=n
 1    if (khi-klo.gt.1) then
         k=(khi+klo)/2
         if (xa(k).gt.x)then
            khi=k
         else
            klo=k
         endif
         goto 1
      endif 
      h=xa(khi)-xa(klo)
      if (h.eq.0.) stop "bad input"
      a=(xa(khi)-x)/h 
      b=(x-xa(klo))/h
      y=a*ya(klo)+b*ya(khi)+
     $     ((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.
      return
      END


     
