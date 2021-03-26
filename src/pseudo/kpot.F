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
      subroutine kpot(ixc)
      implicit double precision (a-h,o-z)
      
#include "fortdim.h"
      common /grid/ h,r1,z,r(npdm),np

      common /totpot/ rvcore(npdm,n0),rvps(npdm,n0),rvcoul(npdm)
      common /rscore/ rscore(npdm),rdd(npdm),rddd(npdm),rscoretot(npdm),
     $     rsval(npdm)
      common /aorb/ ncore,nval,norb,nlm(n0),
     $     no(n0),lo(n0),nmax(n0),maxim
      common /adat/ wnl(n0),en(n0),so(n0),xion
      common /aval/ rcall(n0),rvap(n0),rnorm(n0),ibd(n0),etot
      common /wfn/ rnl(npdm,n0)

      common /wavrc/ wavrc, slope, curvrc, indrc
      common /angm/ ll
      common /nnn/ nnn
      common /numfn/ numfn
      
      dimension rpsi(npdm),rtemp(npdm)

      pi=acos(-1.0)

      astep = 0.2
      idir = 1
      iter = 1
      itmax = 10000
      pow=ll

      do i = 1,np
         rpsi(i) = rnl(i,nnn)*r(i)
      enddo

      do i=1,indrc
         rtemp(i)=(rnl(i,nnn)*r(i))**2
      enddo

      rrnorm=2*ll+2
      call radin(r,rtemp,0,indrc,h,rrnorm)
      headnorm=rrnorm

      do i=indrc+1,np
         rtemp(i)=(rnl(i,nnn)*r(i))**2
      enddo

c     Section I : Fit rpsi to Kerker form
      
      delta=0.0
      rc=r(indrc)
      rc2=rc*rc
      rc3=rc2*rc
      rc4=rc3*rc
      rr=rc**(ll)

      curvrc2 = curvrc - 2.0/rc*slope + float(ll*(ll+1))/rc/rc*wavrc

      slope=(wavrc+rc*slope)/(wavrc*rc)
      wavrc=wavrc*rc

      wavrc=abs(wavrc)
      fkr=log(wavrc/rr)
      fkdr=slope/wavrc - ll/rc
      c1=wavrc*ll*(ll-1)/rc2
      c2=fkdr*2*ll*wavrc/rc
      c3=wavrc*fkdr*fkdr
      fkddr=(curvrc2-c1-c2-c3)/wavrc
      lp=ll+1
      vrc=rvps(indrc,nnn)/rc

      ei=en(nnn)

      alpha = ( 3*log(wavrc/rc**lp) - 2*(rc*slope-lp)
     1     + (rc2*vrc+lp*lp-rc2*(ei+slope*slope))/2 ) / rc4
      beta  = (-8*log(wavrc/rc**lp) + 5*(rc*slope-lp)
     1     - (rc2*vrc+lp*lp-rc2*(ei+slope*slope))   ) / rc3
      gamma = ( 6*log(wavrc/rc**lp) - 3*(rc*slope-lp)
     1     + (rc2*vrc+lp*lp-rc2*(ei+slope*slope))/2 ) / rc2

 911  continue

      do i = 1,indrc
         rpsi(i) = r(i)**lp * exp(delta+gamma*r(i)**2 
     $        + beta*r(i)**3 + alpha*r(i)**4)
      enddo

      rrnorm = float(ll+ll+2)

      do i=1,indrc
         rtemp(i)=(rpsi(i))**2
      enddo

      if (ibd(nnn).eq.0) then
         call radin(r,rtemp,0,indrc,h,rrnorm)
      else
         call radin(r,rtemp,0,np,h,rrnorm)
         headnorm=1.0
      endif

      fdnew=(headnorm-rrnorm)
      if (abs(fdnew).gt.1e-10) then
         if (iter.eq.1) then
            ddd=0.5
         else
            ddd= -fdnew * ddd / (fdnew-fdold)
         endif
         alpha = alpha - 3*ddd/rc4
         beta  = beta  + 8*ddd/rc3
         gamma = gamma - 6*ddd/rc2
         delta = delta + ddd
         fdold = fdnew
         iter=iter+1
         if (iter.eq.itmax) then
            write (7,*) "Couldn't find new rpsi"
            stop
         endif
         goto 911
      endif

      do i=1,indrc
         xlamda=(4*alpha*r(i)+3*beta)*r(i)+2*gamma
         rvps(i,nnn) = (ei + xlamda * (2 * lp + xlamda * r(i)**2)
     $        + (12 * alpha * r(i) + 6 * beta) * r(i) + 2 * gamma)*r(i)
       enddo

      do i = 1,np
         rnl(i,nnn) = rpsi(i) 
         if (ibd(nnn).ne.0) rsval(i)=rsval(i) + wnl(nnn)*rpsi(i)**2
         if (ixc.ge.0) then
            if (i.gt.indrc+50) rvps(i,nnn) = -z-z+rvcoul(i)
            if (i.gt.maxim) rvps(i,nnn) = -xion-xion
         endif
      enddo

      return
      end
