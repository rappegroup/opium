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
      subroutine schslrel(iorb,v,ar,br,
     +     znuc,vid,ev,ninf,ibound)
      implicit double precision(a-h,o-z)
c     
c     integrate the relativistic Dirac equation:
c     find eigenvalue ev and the major and minor
c     components of the wavefunction, ar and br
c
c     Orignal author info:
c     **Revision 1.1  89/10/26  19:53:21  sverre**
c     

#include "fortdim.h"

      common /aorb/ ncore,nval,norb,nlm(n0),
     $     no(n0),lo(n0),nmax(n0),maxim
      common /adat/ wnl(n0),en(n0),so(n0),xion
      common /aval/ rcall(n0),rvap(n0),rnorm(n0),ibd(n0),etot
      common /grid/ h,r1,z,r(npdm),np
      
      common /ipos/ ipos(n0),itermcount

      dimension v(nrmax),ar(nrmax),br(nrmax)
      dimension 
     +     vid(nrmax),
     +     etot(10),ev(n0),f(npdm),rab(npdm)

c     eigenvalue tolerance and max number of iterations

      tol = 1.D-12
      itmax = 50
c     
c     fine-structure constant etc
c     
      ai = 2*137.04D0
      az = znuc/(2*ai)
      ka = lo(iorb)+1
      if (so(iorb) .lt. 0.1D0 
     $     .and. lo(iorb) .ne. 0) ka=-lo(iorb)
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
c      if (so(iorb) .lt. 0.1D0) vzero = vid(1)/r(1)
c      if (so(iorb) .gt. 0.1D0) vzero = viu(1)/r(1)
      vzero = vid(1)/r(1)

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
      do k=1,np
         rab(k)=r(k)*h
      enddo

      do 190 i=1,itmax

c         write(7,4949) i,iorb,ev(iorb),
c     $        lo(iorb),so(iorb),wnl(iorb)
 4949    format(2i8,f10.6,i5,2f10.6)
c     
c     find closest point inside rwell - np,
c     practical infinity ninf, and
c     classical turning point nctp
c     
         np = np
         ninf = np
         nctp = np
         do 100 jj=1,np
            j = np-jj+1
            ar(j) = 0.D0
            br(j) = 0.D0
            idone = 1
c            if (r(j) .gt. rwell) then
c               np = j - 1
c               idone = 0
c            end if
c            if (r(j)*r(j)*(v(j)-ev(iorb)) .gt. 750) then
            if (r(j)*r(j)*(v(j)/r(j)-ev(iorb)) .gt. log(tol)**2) then
               ninf = j
               idone = 0
            end if
            if (v(j)/r(j) .gt. ev(iorb)) then
               nctp = j
               idone = 0
            end if
            if (idone .eq. 1) goto 110
  100    continue
c     
c     three possibilities (np is normally equal to np)
c     
c     nctp < ninf < np  -- normal case, exponetial inward startup
c     nctp < np < ninf  -- bounded case, linear inward startup
c     np < nctp         -- bounded case, no inward integration
c     
c     reset ninf and nctp to allow at least two inward startup points
c     

  110    if (ninf .gt. np) ninf = np
         if (nctp .gt. np - 1) nctp = np - 1
         
         if (ibound.eq.0) nctp=np-200

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
         do 120 j=1,4
            ar(j) = r(j)**s * (1 +(a1+a2*r(j))*r(j))
            br(j) = r(j)**s * (b0+(b1+b2*r(j))*r(j))
c            write(7,*) j,ar(j),br(j)
  120    continue
         fa5 = 0.D0
         fb5 = 0.D0
         fa4 = rab(1)*(+ka*ar(1)/r(1)+(ev(iorb)
     $        -v(1)/r(1)+ai*ai)*br(1)/ai)
         fb4 = rab(1)*(-ka*br(1)/r(1)-(ev(iorb)
     $        -v(1)/r(1))*ar(1)/ai)
         fa3 = rab(2)*(+ka*ar(2)/r(2)+(ev(iorb)
     $        -v(2)/r(2)+ai*ai)*br(2)/ai)
         fb3 = rab(2)*(-ka*br(2)/r(2)-(ev(iorb)-v(2)/r(2))*ar(2)/ai)
         fa2 = rab(3)*(+ka*ar(3)/r(3)+(ev(iorb)
     $        -v(3)/r(3)+ai*ai)*br(3)/ai)
         fb2 = rab(3)*(-ka*br(3)/r(3)-(ev(iorb)-v(3)/r(3))*ar(3)/ai)
         fa1 = rab(4)*(+ka*ar(4)/r(4)+(ev(iorb)
     $        -v(4)/r(4)+ai*ai)*br(4)/ai)
         fb1 = rab(4)*(-ka*br(4)/r(4)-(ev(iorb)-v(4)/r(4))*ar(4)/ai)
c     
c     outward integration loop
c     
         nodes = 0
         do 130 j=4,nctp
c     
c     predictor (Adams-Bashforth)
c     
            arp = ar(j-1) + abc1*fa1+abc2*fa2+abc3*fa3+abc4*fa4+abc5*fa5
            brp = br(j-1) + abc1*fb1+abc2*fb2+abc3*fb3+abc4*fb4+abc5*fb5
            fa0 = rab(j)*(+ka*arp/r(j)+(ev(iorb)
     $           -v(j)/r(j)+ai*ai)*brp/ai)
            fb0 = rab(j)*(-ka*brp/r(j)-(ev(iorb)-v(j)/r(j))*arp/ai)
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
            fa1 = rab(j)*(+ka*arc/r(j)+(ev(iorb)
     $           -v(j)/r(j)+ai*ai)*brc/ai)
            fb1 = rab(j)*(-ka*brc/r(j)-(ev(iorb)-v(j)/r(j))*arc/ai)
            ar(j) = arc + amc0*(fa1-fa0)
            br(j) = brc + amc0*(fb1-fb0)
            fa1 = rab(j)*(+ka*ar(j)/r(j)+(ev(iorb)
     $           -v(j)/r(j)+ai*ai)*br(j)/ai)
            fb1 = rab(j)*(-ka*br(j)/r(j)-(ev(iorb)-v(j)/r(j))*ar(j)/ai)
c     
c     count nodes
c     
            if (ar(j)*ar(j-1) .le. 0) nodes = nodes + 1
  130    continue
c     
c     end outward integration

         if (ibound.eq.0) return
c     
c     if incorrect number of nodes modify energy stepwise
c     
         if (nodes .gt. no(iorb)-lo(iorb)-1) then
c     
c     too many nodes -- decrease ev
c     
            if (ev(iorb) .lt. emax) emax = ev(iorb)
            if (devmax .gt. 0.D0) devmax = -devmax / 2
c            write(7,*) 'EV ',iorb,ev(iorb),devmax
            ev(iorb) = ev(iorb) + devmax
            goto 190
         else if (nodes .lt. no(iorb)-lo(iorb)-1) then
c     
c     too few nodes -- increase ev
c     
            if (ev(iorb) .gt. emin) emin = ev(iorb)
            if (devmax .lt. 0.D0) devmax = -devmax / 2
c            write(7,*) 'EV ',iorb,ev(iorb),devmax
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
         if (ninf .eq. np) then
            ar0 = 0.25D0/ai
            ar1 = -1.D0 + 0.25D0*ka/(ai*r(np))
            ar2 = 0.25D0*ka*(ka-1)/(ai*r(np)*r(np))
         end if
         istart = np - nctp + 1
         if (istart .gt. 4) istart = 4
         do 140 jj=1,istart
            j = ninf-jj+1
            if (ninf .lt. np) then
               alf = v(j)/r(j) - ev(iorb)
               if (alf .lt. 0.D0) alf = 0.D0
               alf = sqrt(alf)
               ar(j) = exp(-alf*r(j))
               arp = -alf * ar(j)
            else
c               dr = r(j) - rwell
               dr=r(j) - r(np)
               ar(j) = ar0 + (ar1 + 0.5D0 * ar2 * dr) * dr
               arp = ar1 + ar2 * dr
            end if
            br(j) = ai*(-arp+ka*ar(j)/r(j))/(v(j)/r(j)-ev(iorb)-ai*ai)
  140    continue
         fa5 = rab(ninf)*(+ka*ar(ninf)/r(ninf)
     +        +(ev(iorb)-v(ninf)/r(ninf)+ai*ai)*br(ninf)/ai)
         fb5 = rab(ninf)*(-ka*br(ninf)/r(ninf)
     +        -(ev(iorb)-v(ninf)/r(ninf))*ar(ninf)/ai)
         fa4 = rab(ninf-1)*(+ka*ar(ninf-1)/r(ninf-1)
     +        +(ev(iorb)-v(ninf-1)/r(ninf-1)+ai*ai)*br(ninf-1)/ai)
         fb4 = rab(ninf-1)*(-ka*br(ninf-1)/r(ninf-1)
     +        -(ev(iorb)-v(ninf-1)/r(ninf-1))*ar(ninf-1)/ai)
         fa3 = rab(ninf-2)*(+ka*ar(ninf-2)/r(ninf-2)
     +        +(ev(iorb)-v(ninf-2)/r(ninf-2)+ai*ai)*br(ninf-2)/ai)
         fb3 = rab(ninf-2)*(-ka*br(ninf-2)/r(ninf-2)
     +        -(ev(iorb)-v(ninf-2)/r(ninf-2))*ar(ninf-2)/ai)
         fa2 = rab(ninf-3)*(+ka*ar(ninf-3)/r(ninf-3)
     +        +(ev(iorb)-v(ninf-3)/r(ninf-3)+ai*ai)*br(ninf-3)/ai)
         fb2 = rab(ninf-3)*(-ka*br(ninf-3)/r(ninf-3)
     +        -(ev(iorb)-v(ninf-3)/r(ninf-3))*ar(ninf-3)/ai)
         fa1 = rab(ninf-4)*(+ka*ar(ninf-4)/r(ninf-4)
     +        +(ev(iorb)-v(ninf-4)/r(ninf-4)+ai*ai)*br(ninf-4)/ai)
         fb1 = rab(ninf-4)*(-ka*br(ninf-4)/r(ninf-4)
     +        -(ev(iorb)-v(ninf-4)/r(ninf-4))*ar(ninf-4)/ai)
c     
c     integration loop
c     
         istop = ninf - nctp
         do 150 jj=4,istop
            j = ninf - jj
c     
c     predictor (Adams-Bashforth)
c     
            arp = ar(j+1)
     +           - (abc1*fa1+abc2*fa2+abc3*fa3+abc4*fa4+abc5*fa5)
            brp = br(j+1)
     +           - (abc1*fb1+abc2*fb2+abc3*fb3+abc4*fb4+abc5*fb5)
            fa0 = rab(j)*(+ka*arp/r(j)+(ev(iorb)
     $           -v(j)/r(j)+ai*ai)*brp/ai)
            fb0 = rab(j)*(-ka*brp/r(j)-(ev(iorb)-v(j)/r(j))*arp/ai)
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
            fa1 = rab(j)*(+ka*arc/r(j)+(ev(iorb)
     $           -v(j)/r(j)+ai*ai)*brc/ai)
            fb1 = rab(j)*(-ka*brc/r(j)-(ev(iorb)-v(j)/r(j))*arc/ai)
c           ar(j) = arc
c           br(j) = brc
            ar(j) = arc + amc0*(fa1-fa0)
            br(j) = brc + amc0*(fb1-fb0)
            fa1 = rab(j)*(+ka*ar(j)/r(j)+(ev(iorb)
     $           -v(j)/r(j)+ai*ai)*br(j)/ai)
            fb1 = rab(j)*(-ka*br(j)/r(j)-(ev(iorb)-v(j)/r(j))*ar(j)/ai)
  150    continue
         arin = ar(nctp)
         arpin = rab(nctp)*(+ka*ar(nctp)/r(nctp)
     +        +(ev(iorb)-v(nctp)/r(nctp)+ai*ai)*br(nctp)/ai)
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
c         factor = 0.D0
c         ll = 4
c         do 170 j=2,ninf
c            factor = factor + ll*(ar(j)*ar(j)+br(j)*br(j))*rab(j)
c            ll = 6 - ll
c  170    continue
c         factor = factor / 3

         do j=1,ninf
            f(j)=ar(j)**2+br(j)**2
         enddo
         pow=2
         call radin(r,f,0,ninf,h,pow)       
         factor=pow
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
         do j=1,ninf
            f(j)=ar(j)**2+br(j)**2
         enddo
         pow=2
         call radin(r,f,0,ninf,h,pow)       
         factor=pow
      endif
c      write(7,*) iorb,pow
c      write(7,*) 'factor 2', iorb,factor
c         ll = 4
c         do j=2,ninf
c            factor = factor + ll*(ar(j)*ar(j)+br(j)*br(j))*rab(j)
c            ll = 6 - ll
c         enddo

c      end if
c      factor = factor / 3
c     
c     error message
c     
      write(nout,210) iorb,ev(iorb),nodes,dev
c         ipos(iorb)=ipos(iorb)+1
c         if (nodes.gt.12) then
c            write(7,*) 'Doing reset of this eigenvalue',ev(iorb)
c            ev(iorb)=-z*z/(no(iorb)**2)
c            write(7,*) 'New value', en(iorb)
c         endif

  210 format(' orb #',i3,' did not converge',/,
     +     ' ev =',e18.10,' nodes =',i8,' dev =',e18.10)
c     
c     normalize wavefunction
c     

  220 factor = 1 / sqrt(factor)

      do j=1,ninf
         ar(j) = factor*ar(j)
         br(j) = factor*br(j)
      enddo

c      rewind(400+iorb)
c      rewind(450+iorb)
c      do j=1,np
c         read(400+iorb,*) rjunk,ar(j)
c         read(450+iorb,*) rjunk,br(j)
c      enddo

      do j=1,ninf
         f(j)=ar(j)**2+br(j)**2
      enddo
      pow=2
      call radin(r,f,0,ninf,h,pow)       
      factor=pow
c      write(7,*) iorb,pow

      return
      end
