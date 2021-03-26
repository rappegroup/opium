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
      subroutine atm(znuc,ixc)

      implicit double precision(a-h,o-z)
#include "param.h"

      dimension cdd(nrmax),cdu(nrmax),cdc(nrmax),viod(lmax,nrmax),
     $ viou(lmax,nrmax),vid(nrmax),viu(nrmax),vod(nrmax),vou(nrmax),
     $ etot(10),ev(norbp),ek(norbp),ep(norbp)

      common /reli/ norb,ncore,nval,no(norbp),lo(norbp)
      common /reld/ zcore,zval,zo(norbp),so(norbp)
      common /rgrid/ a,b,r(nrmax),rab(nrmax)
      common /nrgrid/ nr
      common /rcrel/ rcrel(30),relnorm(30)
      common /rcrelatm/ rcrelatm(30)
      common /results/ rnorm(30), eetot, lghost(30)
      common /rcall/ rcall(16)

      common /atmwave/ nnr,numorb,rr(nrmax),aar(nrmax,30),bbr(nrmax,30),
     $     ccdc(nrmax),eev(norbp),llo(norbp),sso(norbp),
     $     nno(norbp)

      common /filenames/ file_log
      character*80 file_log

      open(unit=7,file=file_log,form='formatted',access='append')

      tol = 1.D-8
      dvold = 1.D10
      zion = znuc - zcore - zval
      zel = zval+zcore

      if (ixc.eq.2) write(nout,332)
      if (ixc.eq.1) write(nout,335)
      if (ixc.eq.0) write(nout,334)
      write(nout,340) znuc,ncore,nval,zel,zion
      write(nout,360)
      call flush(nout) 

      xji = 0.0
      do i=1,norb
         xji = lo(i) + so(i)
         write(nout,370) i,no(i),lo(i),so(i),xji,zo(i)
      enddo
      write(nout,390) r(2),nr,r(nr),-log(znuc*a*b)/2.3025851,1.0/b
      call flush(nout) 

c     set up ionic potential
      do i=1,lmax
         do j=1,nrmax
            viod(i,j) = -2*znuc
            viou(i,j) = -2*znuc
         enddo
      enddo

c     set up initial electronic terms
      call velect(0,0,zel,ixc,cdd,cdu,cdc,
     +     viod,viou,vid,viu,vod,vou,etot,ev,ek,ep)

      do ii=1,nr
         i = nr - ii + 1
         vid(i) = vod(i)
         viu(i) = vou(i)
      enddo

c     find approximate energylevels
      do k=1,norb
         ev(k)=-1.0 * znuc*znuc/(no(k)**2)
      enddo

c     start iteration loop
      iconv = 0
      maxit = 200
      xmixo = 0.4D0
      do 170 iter=1,maxit
         if (iter .eq. maxit) iconv=1

c     compute orbitals
         call dsolv2(iter,iconv,
     +        znuc,cdd,cdu,cdc,
     +        viod,viou,vid,viu,vod,vou,
     +        etot,ev,ek,ep)

c     set up output electronic potential from charge density
         call velect(iter,iconv,zel,ixc,
     +        cdd,cdu,cdc,
     +        viod,viou,vid,viu,vod,vou,
     +        etot,ev,ek,ep)
         if (iconv .gt. 0) goto 190
         dvmax = 0.D0
         do 150 i=1,nr
            dv = (vod(i)-vid(i))/(1.D0+vod(i)+vou(i))
            if (abs(dv) .gt. dvmax) dvmax=abs(dv)
            dv = (vou(i)-viu(i))/(1.D0+vou(i)+vod(i))
            if (abs(dv) .gt. dvmax) dvmax=abs(dv)
  150    continue
         iconv = 1
         if (dvmax .gt. tol) iconv=0
         if (dvmax .ge. dvold) xmixo=0.8D0*xmixo
         if (xmixo .lt. 0.25D0) xmixo=0.25D0
         dvold = dvmax

c     mix input and output electronic potentials
         do i=1,nr
            vid(i)=vid(i)+xmixo*(vod(i)-vid(i))
            viu(i)=viu(i)+xmixo*(vou(i)-viu(i))
         enddo
 170  continue

      write(nout,180) dvmax,xmixo
  180 format(/,' potential not converged - dvmax =',e10.4,
     +     '  xmixo =',f5.3)
      stop


c     find total energy
 190  continue
      write(7,*) 'Converged in ', iter, ' iterations'

      call etotal(etot,ev,ek,ep)
c     update eetot in the results COMMON block  
      eetot = etot(10)   

c     set up necessary info for interp -- should be removed
      numorb = norb - ncore
      nnr = nr
      do j = 1,nr
         rr(j) = r(j)
         ccdc(j) = cdc(j)
      enddo

      do i=1,numorb
         eev(i)=ev(i+ncore)
         llo(i)=lo(i+ncore)
         sso(i)=so(i+ncore)
         nno(i)=no(i+ncore)
      enddo

      k=1
      lastlo=-20
      do i=1,numorb

         factor = 0.D0
         ll = 4
         rc=rcall(k)
         if (llo(i).eq.0.or.llo(i).eq.lastlo) k=k+1  
         lastlo=llo(i)

         do j=nr,0,-1
           if (abs(aar(j,i)).gt.1e-6) goto 913
         enddo
 913     continue
         if (aar(j,i).lt.0) then
           do j=1,nr
             aar(j,i)=-aar(j,i)
           enddo
         endif

         do j=nr,0,-1
           if (abs(bbr(j,i)).gt.1e-6) goto 916
         enddo
 916     continue

         if (bbr(j,i).lt.0) then
           do j=1,nr
             bbr(j,i)=-bbr(j,i)
           enddo
         endif

         do j=2,nr
            if (r(j).gt.rc) goto 911
            factor = factor + ll*(aar(j,i)*aar(j,i)+bbr(j,i)*bbr(j,i))
     $           *rab(j)
            ll = 6 - ll
         enddo
 911     continue
         factor = factor / 3
         relnorm(i)=1.0-factor
      enddo

      close(unit=7)

 370  format(1x,i2,2i5,2f6.1,f10.4)
 390  format(//,' radial grid parameters',//,
     +     ' r(1) = .0 , r(2) =',e8.2,' , ... , r(',i5,') =',f6.2,
     +     /,' aa =',f5.2,'  bb =',f6.2,/)
      
 360  format(' input data for orbitals',//,
     +     '  i    n    l    s     j     occ',/)
      
 340  format(' nuclear charge             =',f10.6,/,
     +     ' number of core orbitals    =',i3,/,
     +     ' number of valence orbitals =',i3,/,
     +     ' electronic charge          =',f10.6,/,
     +     ' ionic charge               =',f10.6,//)
      
 332  format(' XC functional is GGA (Perdew-Burke-Ernzerhof)')
 334  format(' XC functional is LDA (Perdew-Zunger)')
 335  format(' XC functional is LDA (Perdew-Wang)')

      end
