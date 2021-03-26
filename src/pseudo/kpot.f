      subroutine kpot
      implicit double precision (a-h,o-z)
      
#include "PARAMOPT"
      common /grid/ h,r1,z,r(npdm),np
      common /atomic/ xion,rnl(npdm,n0),nlm(n0),wnl(n0),en(n0),norb
      common /angm/ ll
      common /nnn/ nnn
      common /ibound/ ibd(n0)
      common /numfn/ numfn
      common /wavrc/ wavrc, slope, curvrc, indrc
      common /totpot/ rvcore(npdm,n0),rvps(npdm,n0),rvcoul(npdm)
      common /atom3/ rsvale(npdm)
      common /nmax/ nmax(nvale0),maxim
      
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
        rsvale(i)=rsvale(i) + wnl(nnn)*rpsi(i)**2
        if (i.gt.indrc+50) rvps(i,nnn) = -z-z+rvcoul(i)
        if (i.gt.maxim) rvps(i,nnn) = -xion-xion
      enddo

      return
      end
