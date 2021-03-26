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
      subroutine optinvert(psi,r,h,ll,en0,indrc,np,maxim,
     $     g,coef,numfn,xion,v,ixc)
      implicit double precision(a-h,o-z)
      
#include "fortdim.h"

      common /totpot/ rvcore(npdm,n0),rvps(npdm,n0),rvcoul(npdm)
      common /nnn/ nnn

      dimension psi(npdm),r(npdm),v(npdm),d0(npdm),pv(npdm),p(npdm)
      dimension g(numfn0),coef(numfn0)

cEJW Most of the work in this routine is to make the PS and AE 
cEJW pots exactly the same.  Just inverting the wavefunction analytically
cEJW leads to a mall (< 1mRy) discontinuity at rc, therefore, the 
cEJW reference NL calculation will take more than 1 iteration.
cEJW this is basically a compromise between making the potential
cEJW match exactly and making the eigenvalues and norms match close to
cEJW exactly.

c     good compromise
      ioff=100

c     This is for NO analytic section, 0 -> indrc is done numerically
c      if (ixc.lt.0) ioff=indrc

c     This is for NO numerical section,0 -> indrc is done anlytically 
c      ioff=0

cEJW Therefore, there are three regions:  analytic | numerical | AE
cEJW the boundary between the analytic and numerical section is indrc - ioff
cEJW the boundary between the AE and numerical section is indrc + ioff
cEJW So, if ioff=0, there is no numerical section, the boundary between
cEJW the analytic and AE sections is indrc, the AE and PS potentials will match
cEJW exactly at rc, but the eigenvalues and norms will be in error by 
cEJW 10e-4 Ry vs. 10e-8 for ioff around 100.


c     Section I - Numerical procedure
      bbb = (psi(2)/r(2)**ll - psi(1)/r(1)**ll) / (r(2)**2 - r(1)**2)
      aaa = psi(2)/r(2)**ll - bbb * r(2)**2
      vvv = (bbb/aaa) * (4 * ll + 6) + en0

      do i = 1,np
         pv(i) = 0.0
         d0(i) = 0.0
         p(i) = psi(i) * sqrt(r(i))
      enddo

      hsq12 = h * h/12.0

      sqlp = (float(ll) + 0.5)**2
      
      xmult = exp(h)

c We assume that by the far end of the interval, the potential
c is a constant + a 1/r term.

      i0=maxim
      do i = 1,i0
         d0(i) = p(i+2)-p(i+1)-p(i+1)+p(i)
      enddo

c Invert the Numerov difference equivalent of the Schrodinger equation.
c Make j loop go for more than 20 iterations if the potential does not
c give the proper eigenvalue out of scheq.

      do j = 1,20

         do i = 2,i0-1,2
            pv(i) = (d0(i-1) - pv(i-1) - pv(i+1))/10.0
         enddo

         do i = 3,i0-1,2
            pv(i) = (d0(i-1) - pv(i-1) - pv(i+1))/10.0
         enddo

         pv(1) = ((vvv - en0)*r(1)**2+sqlp)*hsq12*p(1)

         if (p(i0-2)*p(i0-1).ne.0.0) then
            tem1 = (pv(i0-2)/p(i0-2)/hsq12-sqlp)/r(i0-2)**2 
            tem2 = (pv(i0-1)/p(i0-1)/hsq12-sqlp)/r(i0-1)**2 
            tem3 = tem2*(xmult+1.0)/xmult - tem1/xmult
            pv(i0) = (r(i0)**2 * tem3 + sqlp) * hsq12 * p(i0)
         else
            tem3 = (-xion-xion)/r(i0) - en0
            pv(i0) = (r(i0)**2 * tem3 + sqlp) * hsq12 * p(i0)
         endif
      enddo

      do i = 1,np
         if (psi(i).gt.0.0) then
            v(i) = (pv(i)/p(i)/hsq12-sqlp)/r(i)/r(i) + en0             
         endif
      enddo

c     section II - analytical procedure
      do i=1,indrc-ioff
         v(i)=0.0
         do j=1,numfn
            gr=g(j)*r(i)
            v(i)=v(i) + coef(j)*g(j)*abs(g(j))*besfn(gr,ll)  
         enddo
         v(i)=en0-v(i)/psi(i)       
      enddo

c     section III - AE -> oo

      if (ixc.ge.0) then

         do i=indrc+ioff,np
            v(i)=(rvcore(i,nnn)+rvcoul(i))/r(i)
         enddo
         
         do i = 1,np
            rvps(i,nnn) = v(i) * r(i)
            if (i.gt.maxim) rvps(i,nnn) = -xion-xion
         enddo
         
      else

         do i=indrc+ioff,np
            v(i)=rvps(i,nnn)/r(i)
         enddo
         
         do i = 1,np
            rvps(i,nnn) = v(i) * r(i)
         enddo

      endif

      return
      end                                                               
