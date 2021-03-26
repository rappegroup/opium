c
c Copyright (c) 1998-2004 The OPIUM Group
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
      subroutine dsolv2(iter,iconv,
     +     znuc,cdd,cdu,cdc,
     +     viod,viou,vid,viu,vod,vou,
     +     etot,ev,ek,ep)
      implicit double precision(a-h,o-z)
c     
c     find the (non) relativistic wave function using
c     difnrl to intgrate the Scroedinger equation or
c     difrel to intgrate the Dirac equation
c     the energy levels from the previous iteration are used
c     as initial guesses and must therefore be reasonable accurate
c
c     Revision 1.1  89/10/26  19:53:26  sverre
c     Initial revision
c     

#include "param.h"

      common /reli/ norb,ncore,nval,no(norbp),lo(norbp)
      common /reld/ zcore,zval,zo(norbp),so(norbp)
      common /rgrid/ aa,bb,r(nrmax),rab(nrmax)
      common /nrgrid/ nr

      common /atmwave/ nnr,numorb,rr(nrmax),aar(nrmax,30),bbr(nrmax,30),
     $                 ccdc(nrmax),eev(norbp),llo(norbp),sso(norbp)

      dimension cdd(nrmax),cdu(nrmax),cdc(nrmax),
     $ viod(lmax,nrmax),viou(lmax,nrmax),vid(nrmax),viu(nrmax),
     $ vod(nrmax),vou(nrmax),
     $ etot(10),ev(norbp),ek(norbp),ep(norbp)

      dimension v(nrmax),ar(nrmax),br(nrmax)

c     initialize arrays for charge density
c     
      do 100 i=1,nr
         cdd(i) = 0.D0
         cdu(i) = 0.D0
         cdc(i) = 0.0
  100 continue
c     
c     start loop over orbitals
c     note that spin zero is treated as down
c     
      do 140 i=1,norb
         if (no(i) .le. 0) goto 140
         if (zo(i) .eq. 0.D0 .and. iconv .eq. 0) goto 140
c     
c     set up potential
c     
         lp  = lo(i)+1
         llp = lo(i)*lp
         do 110 j=2,nr
            if (so(i) .lt. 0.1D0) v(j) = viod(lp,j)/r(j) + vid(j)
            if (so(i) .gt. 0.1D0) v(j) = viou(lp,j)/r(j) + viu(j)
  110    continue

         call difrel(i,v,ar,br,
     +        znuc,cdd,cdu,cdc,
     +        viod,viou,vid,viu,vod,vou,
     +        etot,ev,ek,ep)

c - gjt: store the current orbital wavefunctions
        if (i-ncore.gt.0) then
          do j=1,nr
            aar(j,i-ncore) = ar(j)
            bbr(j,i-ncore) = br(j)
          enddo
        endif

c     add to the charge density
        do 120 j=1,nr
           denr = zo(i) * ar(j) * ar(j)
           denr = denr + zo(i) * br(j) * br(j)
           if (so(i) .lt. 0.1D0) cdd(j) = cdd(j) + denr
           if (so(i) .gt. 0.1D0) cdu(j) = cdu(j) + denr
           if (i .le. ncore) cdc(j)=cdc(j)+denr
 120    continue
        if (iconv .ne. 1) goto 140
c     
c     find orbital kinetic and potential energy
c     potential part includes only interaction with the nucleus
c     
        ek(i) = 0.D0
        ekt = 0.D0
        ep(i) = 0.D0
        ll = 2
        if (2*(nr/2) .eq. nr) ll=4
        do 130 jj=2,nr
           j = nr-jj+2
           ar2 = ar(j)*ar(j)
           br2 = br(j)*br(j)
           denj = ar2
           denj = denj + br2
           ek(i) = ek(i) + ll * (ev(i) - v(j)) * denj * rab(j)
           ekt   = ekt   + ll * (br2 + ar2*llp/r(j)**2) * rab(j)
           if (so(i) .lt. 0.1D0) ep(i) = ep(i)
     +          + ll * denj*viod(lp,j)*rab(j)/r(j)
           if (so(i) .gt. 0.1D0) ep(i) = ep(i)
     +          + ll * denj*viou(lp,j)*rab(j)/r(j)
           ll = 6 - ll
 130    continue
        ek(i) = ek(i) / 3
        ekt   = ekt   / 3
        ep(i) = ep(i) / 3

 140  continue
      
      return
      end
