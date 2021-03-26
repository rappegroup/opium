c--------------------------------------------------------------------
c The following code is adapted from:
c
c A GENERAL HARTREE-FOCK PROGRAM
c
c     by C. Froese Fischer
c        Vanderbilt University
c        Nashville, TN 37235 USA
c
c     April, 1987
c     Published in Computer Physics Communications, 43 (1987) 355
c
c     .... Modified for for arbitrarily filled f-shells
c     by Gediminas GAIGALAS 
c        Institute of Theoretical Physics and Astronomy
c        A. Gostauto str. 12
c        Vilnius, 2600, Lithuania
c       
c     January, 1996
c
c     Pulished in Computer Physics Communications, 98 (1996) 255.
c
c
c
c
c Obtained from http://atoms.vuse.vanderbilt.edu/Elements/CompMeth/hf.f
c-----------------------------------------------------------------------


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
      subroutine hfint(m,ipsp,p,pp,zeff,ed,ntp,azd,deltae,
     $     first,eu,em,iprint)
      implicit double precision (a-h,o-z)
      
#include "fortdim.h"

c -------------------------------------------------------------------------
c     External (shared between C and Fortran) common blocks
c -------------------------------------------------------------------------
      common /aorb/ ncore,nval,norb,nlm(n0),
     $     no(n0),lo(n0),nmax(n0),maxim
      common /adat/ wnl(n0),en(n0),so(n0),xion
      common /aval/ rcall(n0),rvap(n0),rnorm(n0),ibd(n0),etot
      common /consts/ etol,vtol,maxit,isoft
      common /wfn/ rnl(npdm,n0)
      common /grid/ h,r1,z,r(npdm),np
      common /edi/ e(n0,n0)
      common /ahf/ s(n0),az(n0),acc(n0)

c -------------------------------------------------------------------------
c     Internal (Fortran only) common blocks                
c -------------------------------------------------------------------------
      common /totpot/ rvcore(npdm,n0),rvps(npdm,n0),rvcoul(npdm)
      common /iterm/ iterm
      common /elim/ emore,eless,elim
      common /ipos/ ipos(n0),itermcount
      common /itoo/ ietoolo,ietoohi
      common /iskip/ iskip
      common /ameth/ meth(n0),imeth
      common /ighost/ ighost,ighost2

c -------------------------------------------------------------------------
      dimension p(npdm), yk(npdm)
      dimension x(npdm)
      dimension xz(npdm),w(npdm),p1(npdm),p2(npdm)
      dimension yr(npdm),hq(npdm),xx(npdm)
      logical v2,first

      if (.not.first) goto 17

      if (ighost*ighost2.gt.1) goto 17
      call potl(yr,m)
      call xch(x,m,3)
      ighost2=ighost2+1

      cd=(lo(m)+1.d0/2.d0)**2
      mk=0
      
      do i=1,np
         yk(i)=-2.d0*(rvcore(i,m)-yr(i))*r(i)+cd
      enddo

      zinf=max(0.05d0,rvcore(nmax(m),m)-yr(nmax(m)))
      v=yr(1)/r(1)
      c=4.d0*lo(m)+6.d0

      ch=h*h/12.d0
      cn=(2.d0*z/no(m))**(lo(m)+1)

      xy=x(1)
      xp=x(2)
      ed=e(m,m)
      x1=x(1)
      x2=x(2)
      x3=x(3)
      x4=x(4)
      
      do j=1,np
         x5=x(j+2)
         x(j)=ch*(-x5+24.d0*(x2+x4)+194.d0*x3-x1)/20.d0
         x1=x2
         x2=x3
         x3=x4
         x4=x5
      enddo
      x(nmax(m)-1)=ch*(x4+10.d0*x3+x2)

      x1=ch*rnl(1,m)*(yk(1)+ed*r(1)*r(1))
      x2=ch*rnl(2,m)*(yk(2)+ed*r(2)*r(2))
      x3=ch*rnl(3,m)*(yk(3)+ed*r(3)*r(3))
      x4=ch*rnl(4,m)*(yk(4)+ed*r(4)*r(4))

      do j=3,np
         x5=ch*rnl(j+2,m)*(yk(j+2)+ed*r(j+2)*r(j+2))
         x(j)=x(j)-(x5-4.d0*(x2+x4)+6.d0*x3+x1)/20.d0
         x1=x2
         x2=x3
         x3=x4
         x4=x5
      enddo
      r1=lo(m)+2.5d0
      x(2)=r(2)**r1*(x(5)/r(5)**r1 - 3.d0*(x(4)/r(4)**r1-x(3)/r(3)**r1))

      if (imeth.ne.3) goto 80

      do jj=15,np
         j=np-jj
         if (yk(j).lt.0.d0) goto 63
      enddo
      write(7,*) "POTENTIAL TOO SMALL!"
      goto 80
 63   em=-yk(j)/r(j)**2
      goto 81
 80   continue
      
      if (ipsp.eq.0) then
         em=(zinf/(no(m)+0.5))**2
      endif
 81   fm=em


      fl=0.d0
      cden=0.d0
      np2=np-2
      do j=2,np2
         fnum=rnl(j+1,m)-2.d0*rnl(j,m)+rnl(j-1,m)
         fnum=fnum-ch*(yk(j+1)*rnl(j+1,m)+10.d0*yk(j)*rnl(j,m)
     $        +yk(j-1)*rnl(j-1,m))
     $        -x(j)
         ddum=r(j+1)*r(j+1)*rnl(j+1,m) + 10.d0*r(j)*r(j)*rnl(j,m)
     $        + r(j-1)*r(j-1)*rnl(j-1,m)
         fl=fl+rnl(j,m)*fnum
         cden=cden+rnl(j,m)*ddum
      enddo
      ed=fl/(cden*ch)

      if (ed.lt.em) then
         ed=em
         if (dabs(fm-e(m,m)).gt.1e-6.or.kk.eq.3) goto 19
         if (ipsp.eq.0) then
            zzz=z
         else
            zzz=zeff
         endif
         pnf=hnorm(no(m),lo(m),zzz-s(m))
         do k=1,np
            rnl(k,m)=pnf*hwf(no(m),lo(m),zzz-s(m),r(k))/sqrt(r(k))
         enddo
      endif

 19   continue
      if (10.d0*ed.gt.eu) then
         eu=10.d0*ed
         fu=eu
      endif
      azd=az(m)

 17   continue
      
      do j=1,np
         yr(j)=(yk(j)+ed*r(j)*r(j))*ch
         xz(j)=0.d0
      enddo

      ia=100
      il=np
 4    if (yr(ia).lt.0.d0) goto 3
      ia=ia+2
      if (ia.lt.il) goto 4
      ntp=max(100,nmax(m)-100)
      goto 101
 3    nk=(ia+il)/2
      if (yr(nk).lt.0.d0) goto 1
      il=nk
      goto 2
 1    ia=nk
 2    if (il-ia.gt.1) goto 3
      ntp=il-7
 101  continue

      b3=(v+v+ed-(z/no(m))**2)/c
      b4=z*(lo(m)+4.d0/3.d0)/((lo(m)+1.d0)*(lo(m)+2.d0))

      do k=1,2
         if (ipsp.eq.0) then
            hw=hwf(no(m),lo(m),z,r(k))/cn
         else
            hw=rnl(k,m)*sqrt(r(k))
         endif
         hq(k)=azd*(hw+r(k)**(lo(m)+3)*b3*(1.d0-r(k)*b4))/sqrt(r(k))
      enddo

      call nmrvs(ntp,delh,mh,hq,xz,yr,m,iprint)
      p1(1)=hq(1)+xy/c
      p1(2)=hq(2)+xp/c
      call nmrvs(ntp,del1,m1,p1,x,yr,m,iprint)

      nmax(m)=max(m1,mh)
      pnorm=0.d0
      do j=1,nmax(m)
         pnorm=pnorm+r(j)**2*hq(j)*p1(j)
      enddo

      y1=p1(ntp-1)
      y2=p1(ntp)
      y3=p1(ntp+1)
      delta=y2-y1+y2-y3+yr(ntp-1)*y1+10.d0*yr(ntp)*y2
     $     +yr(ntp+1)*y3+x(ntp)
      deltae=hq(ntp)*delta/(h*h*pnorm)
      pp=-del1/delh
      
      do j=1,np
         p1(j)=p1(j)+pp*hq(j)
      enddo

      if (imeth.eq.2) then
c         write(7,*) 'Hey in var section'

         x1 = rnl(1,m)*r(1)**2
         x2 = rnl(2,m)*r(2)**2
         p2(1) = x1/c
         p2(2) = x2/c
         do  j = 3,np
            x3 = rnl(j,m)*r(j)**2
            xx(j-1) = (10.d0*x2 + x1 + x3)*ch
            x1 = x2
 8          x2 = x3
         enddo
         call nmrvs(ntp,del2,m2,p2,xx,yr,m,iprint)
         aa = -del2/delh
         nmax(m) = max0(nmax(m),m2)
         do j = 1,np
            p2(j) = p2(j) + aa*hq(j)
         enddo
         a11 = quad(m,nmax(m),p2,p2)
         b11 = quad(m,nmax(m),p1,p2)
         c11 = quad(m,nmax(m),p1,p1) - 1.d0
         disc = b11*b11 - a11*c11
         if ( disc .lt. 0.d0 ) go to 70

         de1 = -(b11+dsqrt(disc))/a11
         de2 = c11/a11/de1
         if ( p1(3)+de1*p2(3) .lt. 0.d0) de1 = de2
         go to 71
 70      de1 = c11/a11
 71      do  j = 1,np
            p1(j) = p1(j) + de1*p2(j)
         enddo
         pp = pp + de1*aa
      endif

      do j=1,np
         p(j)=p1(j)
      enddo
         
c      write(7,9955) m,pp,del1,delh
 9955 format(1x,'pp',i5,3f20.15)
      
      return
      end


*     ------------------------------------------------------------------
*               n m r v s
*     ------------------------------------------------------------------
*
*       given two starting values, pde(1) and pde(2), values of pde(j),
*   j=3,4,...,nj+1 are obtained by outward integration of
*               y" = yr y + f
*   using the discretization  of  eq.  (6-27 )  with  the  difference
*   correction.  with pde(nj) given, the tail procedure is applied to
*   pde(j),j=nj+1,  nj+2,...,mm, where mm is determined automatically
*   and delta is the difference between  pde(nj+1)  for  outward  and
*   inward integration. (see eq 6-32, 6-33, and 6-37 for further
*   details.)
*
      subroutine nmrvs(nj,delta,mm,pdef,f,yr,morb,iprint)
      implicit double precision(a-h,o-z)
#include "fortdim.h"
      common /grid/ h,r1,z,r(npdm),np
      dimension pdef(npdm),f(npdm),a(npdm),d(npdm),yr(npdm)
      equivalence (g,g3)

      ch=h*h/12.d0
      tol=1e-7
      eh=dexp(-h)
*  *****  integrate outward to nj+1

      y1 = pdef(1)
      y2 = pdef(2)
      g1 = yr(1)
      g2 = yr(2)
      m = nj + 1
      do i = 3,m
         g3 = yr(i)
         y3 = (y2+y2-y1 + (10.d0*g2*y2 + g1*y1) + f(i-1)) / (1.d0 - g3)
         pdef(i) = y3
         y1 = y2
         y2 = y3
         g1 = g2
         g2 = g3
      enddo
      delta = y3

*  *****  apply the tail procedure
      k = 1
      pdef(m) = -(1.d0 - g1)*y1 + f(m)
      a(1) = 1.d0 - g
      d(1) = -(2.d0 + 10.d0*g)
22    ratio = a(k)/d(k)

      
      if (k .ge. (np)-1 .or. m .eq. np) then 
         go to 23
      endif
      k = k +1
      m = m+1
      g = yr(m)
      a(k) = 1.d0 - g
      d(k) = -(2.d0 + 10.d0*g) - a(k)*ratio
      pdef(m) = -pdef(m-1)*ratio + f(m)

      if (dabs(pdef(m))+dabs(pdef(m-1)) .gt. tol .or. k .lt. 9) go to 22
20    con =dsqrt(eh)*dexp(-dsqrt(dabs(g/ch-.25)/r(m)**2)*(r(m+1)-r(m)))
      pdef(m) = pdef(m)/(d(k) + con*(1.d0-  yr(m+1)))

      j = m+1
      do  i= j,np
         pdef(i) = 0.d0
      enddo
      do j = 2,k
         i = m-j+1
         ii = k-j+1
         pdef(i) = (pdef(i)-a(ii+1)*pdef(i+1))/d(ii)
      enddo

*  *****  set delta = difference of the two solutions at nj+1
*  *****         mm = number of points in the range of the solution

      delta = delta - pdef(i)
      mm = m

      return
23    if (iprint.gt.0) write (7,24) morb,r(np),pdef(np)+pdef(np-1),
     $     pdef(np-1)+pdef(np-2),tol
24    format(i3,6x,f20.10,3e10.3
     $     '  : function truncated in nmrvs',i6)
      go to 20
      end

