      subroutine difrel(iorb,v,ar,br,
     +     znuc,cdd,cdu,cdc,
     +     viod,viou,vid,viu,vod,vou,
     +     etot,ev,ek,ep)
      implicit double precision(a-h,o-z)
c     
c     integrate the relativistic Dirac equation:
c     find eigenvalue ev and the major and minor
c     components of the wavefunction, ar and br
c
c     Revision 1.1  89/10/26  19:53:21  sverre
c     Initial revision
c     

#include "param.h"

      common /rcrelatm/ rcrelatm(30)
      common /rcrel/ rcrel(30),relnorm(30)

      common /reli/ norb,ncore,nval,no(norbp),lo(norbp)
      common /reld/ zcore,zval,zo(norbp),so(norbp)
      common /rgrid/ aa,bb,r(nrmax),rab(nrmax)
      common /nrgrid/ nr

      dimension v(nrmax),ar(nrmax),br(nrmax)
      dimension 
     +     cdd(nrmax),cdu(nrmax),cdc(nrmax),
     +     viod(lmax,nrmax),viou(lmax,nrmax),vid(nrmax),viu(nrmax),
     $     vod(nrmax),vou(nrmax),
     +     etot(10),ev(norbp),ek(norbp),ep(norbp)

c     eigenvalue tolerance and max number of iterations

      tol = 1.D-10
      itmax = 50
c     
c     fine-structure constant etc
c     
      ai = 2*137.04D0
      az = znuc/(2*ai)
      ka = lo(iorb)+1
      if (so(iorb) .lt. 0.1D0 .and. lo(iorb) .ne. 0) ka=-lo(iorb)
c     
c     integration coefficients
c     
      abc1 = 1901.D0/720.D0
      abc2 = -1387.D0/360.D0
      abc3 = 109.D0/30.D0
      abc4 = -637.D0/360.D0
      abc5 = 251.D0/720.D0
      amc0 = 251.D0/720.D0
      amc1 = 323.D0/360.D0
      amc2 = -11.D0/30.D0
      amc3 = 53.D0/360.D0
      amc4 = -19.D0/720.D0
c     
c     determine effective charge and vzero for
c     startup of outward integration
c     
c     ar = r**s * (1  + a1 r + a2 r**2 + ... )
c     br = r**s * (b0 + b1 r + b2 r**2 + ... )
c     
c     s = sqrt (ka**2 - az**2)    b0 = - az / (s + ka)
c     
c     an = (az (v0 - e) a(n-1) - (s + n + ka) (v0 - e - ai**2) b(n-1))
c     .    / (n ai (2 s + n))
c     
c     bn = ((v0 - e) a(n-1) - 2 znuc an ) / ( ai (s + n + ka))
c     
      s = sqrt(ka*ka-az*az)
      if (ka .gt. 0) b0 = -az/(s+ka)
      if (ka .le. 0) b0 = (s-ka)/az
      if (so(iorb) .lt. 0.1D0) vzero = vid(2)
      if (so(iorb) .gt. 0.1D0) vzero = viu(2)
c     
c     these are used to bracket eigenvalue
c     
      emax = +1.D+20
      emin = -1.D+20
c     
c     max step size for eigenvalue changes
c     
      devmax = -ev(iorb) / 5
      if (devmax .lt. 0.3D0) devmax = 0.3D0
c     
c     begin iteration loop
c     
      do 190 i=1,itmax
c     
c     find closest point inside rwell - nr,
c     practical infinity ninf, and
c     classical turning point nctp
c     
         nr = nr
         ninf = nr
         nctp = nr
         do 100 jj=2,nr
            j = nr-jj+2
            ar(j) = 0.D0
            br(j) = 0.D0
            idone = 1
c            if (r(j) .gt. rwell) then
c               nr = j - 1
c               idone = 0
c            end if
            if (r(j)*r(j)*(v(j)-ev(iorb)) .gt. log(tol)**2) then
               ninf = j
               idone = 0
            end if
            if (v(j) .gt. ev(iorb)) then
               nctp = j
               idone = 0
            end if
            if (idone .eq. 1) goto 110
  100    continue
c     
c     three possibilities (nr is normally equal to nr)
c     
c     nctp < ninf < nr  -- normal case, exponetial inward startup
c     nctp < nr < ninf  -- bounded case, linear inward startup
c     nr < nctp         -- bounded case, no inward integration
c     
c     reset ninf and nctp to allow at least two inward startup points
c     


  110    if (ninf .gt. nr) ninf = nr
         if (nctp .gt. nr - 1) nctp = nr - 1
c     
c     outward integration from 1 to nctp -- startup
c     
         a1 = (az*(vzero-ev(iorb))-(s+1+ka)
     +        *(vzero-ev(iorb)-ai**2)*b0) / (ai*(2*s+1))
         b1 = ((vzero-ev(iorb))-2*znuc*a1) / (ai*(s+1+ka))
         a2 = (az*(vzero-ev(iorb))*a1-(s+2+ka)
     +        *(vzero-ev(iorb)-ai**2)*b1) / (2*ai*(2*s+2))
         b2 = ((vzero-ev(iorb))*a1-2*znuc*a2) / (ai*(s+2+ka))
         ar(1) = 0.D0
         br(1) = 0.D0
         do 120 j=2,5
            ar(j) = r(j)**s * (1 +(a1+a2*r(j))*r(j))
            br(j) = r(j)**s * (b0+(b1+b2*r(j))*r(j))
  120    continue
         fa5 = 0.D0
         fb5 = 0.D0
         fa4 = rab(2)*(+ka*ar(2)/r(2)+(ev(iorb)-v(2)+ai*ai)*br(2)/ai)
         fb4 = rab(2)*(-ka*br(2)/r(2)-(ev(iorb)-v(2))*ar(2)/ai)
         fa3 = rab(3)*(+ka*ar(3)/r(3)+(ev(iorb)-v(3)+ai*ai)*br(3)/ai)
         fb3 = rab(3)*(-ka*br(3)/r(3)-(ev(iorb)-v(3))*ar(3)/ai)
         fa2 = rab(4)*(+ka*ar(4)/r(4)+(ev(iorb)-v(4)+ai*ai)*br(4)/ai)
         fb2 = rab(4)*(-ka*br(4)/r(4)-(ev(iorb)-v(4))*ar(4)/ai)
         fa1 = rab(5)*(+ka*ar(5)/r(5)+(ev(iorb)-v(5)+ai*ai)*br(5)/ai)
         fb1 = rab(5)*(-ka*br(5)/r(5)-(ev(iorb)-v(5))*ar(5)/ai)
c     
c     outward integration loop
c     
         nodes = 0
         do 130 j=6,nctp
c     
c     predictor (Adams-Bashforth)
c     
            arp = ar(j-1) + abc1*fa1+abc2*fa2+abc3*fa3+abc4*fa4+abc5*fa5
            brp = br(j-1) + abc1*fb1+abc2*fb2+abc3*fb3+abc4*fb4+abc5*fb5
            fa0 = rab(j)*(+ka*arp/r(j)+(ev(iorb)-v(j)+ai*ai)*brp/ai)
            fb0 = rab(j)*(-ka*brp/r(j)-(ev(iorb)-v(j))*arp/ai)
c     
c     corrector (Adams-Moulton)
c     
            arc = ar(j-1) + amc0*fa0+amc1*fa1+amc2*fa2+amc3*fa3+amc4*fa4
            brc = br(j-1) + amc0*fb0+amc1*fb1+amc2*fb2+amc3*fb3+amc4*fb4
            fa5 = fa4
            fb5 = fb4
            fa4 = fa3
            fb4 = fb3
            fa3 = fa2
            fb3 = fb2
            fa2 = fa1
            fb2 = fb1
            fa1 = rab(j)*(+ka*arc/r(j)+(ev(iorb)-v(j)+ai*ai)*brc/ai)
            fb1 = rab(j)*(-ka*brc/r(j)-(ev(iorb)-v(j))*arc/ai)
            ar(j) = arc + amc0*(fa1-fa0)
            br(j) = brc + amc0*(fb1-fb0)
            fa1 = rab(j)*(+ka*ar(j)/r(j)+(ev(iorb)-v(j)+ai*ai)*br(j)/ai)
            fb1 = rab(j)*(-ka*br(j)/r(j)-(ev(iorb)-v(j))*ar(j)/ai)
c     
c     count nodes
c     
            if (ar(j)*ar(j-1) .le. 0) nodes = nodes + 1
  130    continue
c     
c     end outward integration
c     
c     if incorrect number of nodes modify energy stepwise
c     
         if (nodes .gt. no(iorb)-lo(iorb)-1) then
c     
c     too many nodes -- decrease ev
c     
            if (ev(iorb) .lt. emax) emax = ev(iorb)
            if (devmax .gt. 0.D0) devmax = -devmax / 2
            ev(iorb) = ev(iorb) + devmax
            goto 190
         else if (nodes .lt. no(iorb)-lo(iorb)-1) then
c     
c     too few nodes -- increase ev
c     
            if (ev(iorb) .gt. emin) emin = ev(iorb)
            if (devmax .lt. 0.D0) devmax = -devmax / 2
            ev(iorb) = ev(iorb) + devmax
            goto 190
         end if
c     
c     correct number of nodes
c     
         arout = ar(nctp)
         arpout = fa1
c     
c     inward integration from ninf to nctp -- startup
c     
         if (ninf .eq. nr) then
            ar0 = 0.25D0/ai
            ar1 = -1.D0 + 0.25D0*ka/(ai*r(nr))
            ar2 = 0.25D0*ka*(ka-1)/(ai*r(nr)*r(nr))
         end if
         istart = nr - nctp + 1
         if (istart .gt. 5) istart = 5
         do 140 jj=1,istart
            j = ninf-jj+1
            if (ninf .lt. nr) then
               alf = v(j) - ev(iorb)
               if (alf .lt. 0.D0) alf = 0.D0
               alf = sqrt(alf)
               ar(j) = exp(-alf*r(j))
               arp = -alf * ar(j)
            else
c               dr = r(j) - rwell
               dr=r(j) - r(nr)
               ar(j) = ar0 + (ar1 + 0.5D0 * ar2 * dr) * dr
               arp = ar1 + ar2 * dr
            end if
            br(j) = ai*(-arp+ka*ar(j)/r(j))/(v(j)-ev(iorb)-ai*ai)
  140    continue
         fa5 = rab(ninf)*(+ka*ar(ninf)/r(ninf)
     +        +(ev(iorb)-v(ninf)+ai*ai)*br(ninf)/ai)
         fb5 = rab(ninf)*(-ka*br(ninf)/r(ninf)
     +        -(ev(iorb)-v(ninf))*ar(ninf)/ai)
         fa4 = rab(ninf-1)*(+ka*ar(ninf-1)/r(ninf-1)
     +        +(ev(iorb)-v(ninf-1)+ai*ai)*br(ninf-1)/ai)
         fb4 = rab(ninf-1)*(-ka*br(ninf-1)/r(ninf-1)
     +        -(ev(iorb)-v(ninf-1))*ar(ninf-1)/ai)
         fa3 = rab(ninf-2)*(+ka*ar(ninf-2)/r(ninf-2)
     +        +(ev(iorb)-v(ninf-2)+ai*ai)*br(ninf-2)/ai)
         fb3 = rab(ninf-2)*(-ka*br(ninf-2)/r(ninf-2)
     +        -(ev(iorb)-v(ninf-2))*ar(ninf-2)/ai)
         fa2 = rab(ninf-3)*(+ka*ar(ninf-3)/r(ninf-3)
     +        +(ev(iorb)-v(ninf-3)+ai*ai)*br(ninf-3)/ai)
         fb2 = rab(ninf-3)*(-ka*br(ninf-3)/r(ninf-3)
     +        -(ev(iorb)-v(ninf-3))*ar(ninf-3)/ai)
         fa1 = rab(ninf-4)*(+ka*ar(ninf-4)/r(ninf-4)
     +        +(ev(iorb)-v(ninf-4)+ai*ai)*br(ninf-4)/ai)
         fb1 = rab(ninf-4)*(-ka*br(ninf-4)/r(ninf-4)
     +        -(ev(iorb)-v(ninf-4))*ar(ninf-4)/ai)
c     
c     integration loop
c     
         istop = ninf - nctp
         do 150 jj=5,istop
            j = ninf - jj
c     
c     predictor (Adams-Bashforth)
c     
            arp = ar(j+1)
     +           - (abc1*fa1+abc2*fa2+abc3*fa3+abc4*fa4+abc5*fa5)
            brp = br(j+1)
     +           - (abc1*fb1+abc2*fb2+abc3*fb3+abc4*fb4+abc5*fb5)
            fa0 = rab(j)*(+ka*arp/r(j)+(ev(iorb)-v(j)+ai*ai)*brp/ai)
            fb0 = rab(j)*(-ka*brp/r(j)-(ev(iorb)-v(j))*arp/ai)
c     
c     corrector (Adams-Moulton)
c     
            arc = ar(j+1)
     +           - (amc0*fa0+amc1*fa1+amc2*fa2+amc3*fa3+amc4*fa4)
            brc = br(j+1)
     +           - (amc0*fb0+amc1*fb1+amc2*fb2+amc3*fb3+amc4*fb4)
            fa5 = fa4
            fb5 = fb4
            fa4 = fa3
            fb4 = fb3
            fa3 = fa2
            fb3 = fb2
            fa2 = fa1
            fb2 = fb1
            fa1 = rab(j)*(+ka*arc/r(j)+(ev(iorb)-v(j)+ai*ai)*brc/ai)
            fb1 = rab(j)*(-ka*brc/r(j)-(ev(iorb)-v(j))*arc/ai)
c           ar(j) = arc
c           br(j) = brc
            ar(j) = arc + amc0*(fa1-fa0)
            br(j) = brc + amc0*(fb1-fb0)
            fa1 = rab(j)*(+ka*ar(j)/r(j)+(ev(iorb)-v(j)+ai*ai)*br(j)/ai)
            fb1 = rab(j)*(-ka*br(j)/r(j)-(ev(iorb)-v(j))*ar(j)/ai)
  150    continue
         arin = ar(nctp)
         arpin = rab(nctp)*(+ka*ar(nctp)/r(nctp)
     +        +(ev(iorb)-v(nctp)+ai*ai)*br(nctp)/ai)
c     
c     end inward integration
c     
c     rescale ar and br outside nctp to match
c     ar(nctp) from outward integration
c     
         factor = arout/arin
         do 160 j=nctp,ninf
            ar(j) = factor * ar(j)
            br(j) = factor * br(j)
  160    continue
         arpin = factor * arpin
c     
c     find normalization
c     
         factor = 0.D0
         ll = 4
         do 170 j=2,ninf
            factor = factor + ll*(ar(j)*ar(j)+br(j)*br(j))*rab(j)
            ll = 6 - ll
  170    continue
         factor = factor / 3
c     
c     modify eigenvalue ev
c     
         dev = arout * (arpout-arpin) / (factor * rab(nctp))
c
c     resort to bisection if dev too large
c
         if (abs(dev) .gt. abs(devmax)) then
            if (devmax*dev .lt. 0.D0) devmax = -devmax / 2
            dev = devmax
         end if
         evold = ev(iorb)
         ev(iorb) = ev(iorb) + dev
         if (ev(iorb) .gt. emax) ev(iorb) = (evold + emax) / 2
         if (ev(iorb) .lt. emin) ev(iorb) = (evold + emin) / 2
         if (abs(dev) .lt. tol*(1+abs(ev(iorb)))) goto 220
  190 continue
c     
c     eigenpar not converged in itmax iterations
c     
c     if missing -- find normalization
c     
      if (nodes .ne. no(iorb)-lo(iorb)-1) then
         factor = 0.D0
         ll = 4

         do j=2,ninf
            factor = factor + ll*(ar(j)*ar(j)+br(j)*br(j))*rab(j)
            ll = 6 - ll
         enddo

      end if
      factor = factor / 3
c     
c     error message
c     
      write(nout,210) iorb,ev(iorb),nodes,dev
  210 format(' orb #',i3,' did not converge',/,
     +     ' ev =',e18.10,' nodes =',i2,' dev =',e18.10)
c     
c     normalize wavefunction
c     

  220 factor = 1 / sqrt(factor)

      do j=1,ninf
         ar(j) = factor*ar(j)
         br(j) = factor*br(j)
      enddo

      return
      end
