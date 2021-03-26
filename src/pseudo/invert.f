c ****************************************************************************
c Copyright (C) 2002 Andrew M Rappe group, University of Pennsylvania
c This file is distributed under the terms of the GNU General Public
c License as described in the file 'License' in the current directory.
c ****************************************************************************

      subroutine potinv(psi,r,h,ll,en0,i0,np,xion,v)
      implicit real(a-h,o-z)
      
#include "PARAMOPT"

c This program takes a function and finds the potential which gives 
c rise to that wavefunction.
      dimension psi(npdm),r(npdm),v(npdm),d0(npdm),pv(npdm),p(npdm)
camr Part 1 of bug fix
      bbb = (psi(2)/r(2)**ll - psi(1)/r(1)**ll) / (r(2)**2 - r(1)**2)
      aaa = psi(2)/r(2)**ll - bbb * r(2)**2
      vvv = (bbb/aaa) * (4 * ll + 6) + en0
camr Part 1 of bug fix
      do 5 i = 1,np
         pv(i) = 0.0
         d0(i) = 0.0
         p(i) = psi(i) * sqrt(r(i))
 5    continue
      hsq12 = h * h/12.0
      sqlp = (float(ll) + 0.5)**2
      xmult = exp(h)
c We assume that at r=0, potential is linear.
c We assume that by the far end of the interval, the potential
c is a constant + a 1/r term.
      do 14 i = 1,i0
         d0(i) = p(i+2)-p(i+1)-p(i+1)+p(i)
 14   continue
c Invert the Numerov difference equivalent of the Schrodinger equation.
c Make j loop go for more than 20 iterations if the potential does not
c give the proper eigenvalue out of scheq.
      do 16 j = 1,20
         do 17 i = 2,i0-1,2
            pv(i) = (d0(i-1) - pv(i-1) - pv(i+1))/10.0
 17      continue
         do 18 i = 3,i0-1,2
            pv(i) = (d0(i-1) - pv(i-1) - pv(i+1))/10.0
 18      continue
camr Part 2 of bug fix
c$$$         tem2 = (pv(2)/p(2)/hsq12-sqlp)/r(2)**2 
c$$$         tem3 = (pv(3)/p(3)/hsq12-sqlp)/r(3)**2 
c$$$         tem1 = tem2*(xmult+1.0)/xmult - tem3/xmult
c$$$         pv(1) = (r(1)**2 * tem1 + sqlp) * hsq12 * p(1)
         pv(1) = ((vvv - en0)*r(1)**2+sqlp)*hsq12*p(1)
camr Part 2 of bug fix
         if (p(i0-2)*p(i0-1).ne.0.0) then
            tem1 = (pv(i0-2)/p(i0-2)/hsq12-sqlp)/r(i0-2)**2 
            tem2 = (pv(i0-1)/p(i0-1)/hsq12-sqlp)/r(i0-1)**2 
            tem3 = tem2*(xmult+1.0)/xmult - tem1/xmult
            pv(i0) = (r(i0)**2 * tem3 + sqlp) * hsq12 * p(i0)
         else
            tem3 = (-xion-xion)/r(i0) - en0
            pv(i0) = (r(i0)**2 * tem3 + sqlp) * hsq12 * p(i0)
         endif
 16   continue
      do 19 i = 1,np
         if (psi(i).eq.0.0) then
            v(i) = (-xion-xion)/r(i)
            goto 19
         endif
         v(i) = (pv(i)/p(i)/hsq12-sqlp)/r(i)/r(i) + en0 
camr Part 3 of bug fix
camr This method seems to have problems for small r and large ll.
         if (r(i).lt.0.01) v(i) = vvv
camr Part 3 of bug fix
 19   continue
      return
      end                                                               
