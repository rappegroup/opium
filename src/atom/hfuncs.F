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
*    1/n*sqrt(Z*(n-l-1)!*(n+l)!) the loop is amazing!!
      double precision function hnorm(n,l,zs)
      implicit double precision(a-h,o-z)
#include "fortdim.h"

      mt = l + l + 1
      a = n + l
      b = mt
      t = a
      d = b
      mt = mt - 1
      if (mt .eq. 0) go to 2
      do 1 i = 1,mt
      a = a - 1.d0
      b = b - 1.d0
      t = t*a
1     d = d*b
2     hnorm = dsqrt(zs*t)/( n*d)
      return
      end

*    MCHF_HF (Part2 of 2)
*     ------------------------------------------------------------------
*               H W F
*     ------------------------------------------------------------------
*
*       Returns the value of an unnormalized (nl) hydrogenic function
*   with nuclear charge ZZ and radius r.
*   hwf stands for hydrogenic wavefunction, equations are in page 4 of 
*   Fisher's book.
*   the sum part of the Laguerre polynomials P(nl;r) computed by the
*   loop  
*
      double precision function hwf(n,l,zs,r)
      implicit double precision(a-h,o-z)
#include "fortdim.h"

      kval = n-l-1
      pt = 1.d0
      a = 1.d0
      b = kval
      cc = n+ l
      xp = -2.d0*zs*r/n
*
*  *****  test if underflow may occur, if so set hwf = 0
*
      if ( xp .lt. -150.d0 ) go to 5
      if (kval) 1,2,3
3     do 4 i = 1,kval
      pt = 1.d0 + a/b*pt/cc*xp
      a = a + 1.d0
      b = b - 1.d0
4     cc = cc - 1.d0
2     hwf = pt*dexp(xp/2.d0)*(-xp)**(l+1)
      return
1     write(6,7) n,l,zs,r
7     format(51h forbidden combination of n and l in hwf subprogram/
     :    4h n = ,i4,6h   l = ,i4,6h   z = ,f6.1,6h   r = ,f8.4)
      stop
5     hwf = 0.d0
      return
      end


      double precision function quadr(i,j,k)
      implicit double precision(a-h,o-z)
#include"fortdim.h"
      common /grid/ h,r1,z,r(npdm),np
      common /aorb/ ncore,nval,norb,nlm(n0),
     $     no(n0),lo(n0),nmax(n0),maxim
      common /adat/ wnl(n0),en(n0),so(n0),xion
      common /aval/ rcall(n0),rvap(n0),rnorm(n0),ibd(n0),etot
      common /wfn/ rnl(npdm,n0)
      common /ahf/ s(n0),az(n0),acc(n0),dpm(n0)

      k12 = k + 2
      li = lo(i)
      lj = lo(j)
      h1=h/1.5
      den = li + lj + 1 + k12
      zr = z*r(4)
      r24=sqrt(r(4))
      bi = (rnl(4,i)/(az(i)*r24*r(4)**li) - 1.d0+zr/(li+1) )/zr**2
      bj = (rnl(4,j)/(az(j)*r24*r(4)**lj) - 1.d0+zr/(lj+1) )/zr**2
      alpha= (1.d0/(li + 1) + 1.d0/(lj + 1))/(den + 1.d0)
      zr = z*r(1)
      beta = (den+1.d0)*alpha**2 
     $     - 2.d0*(bi+bj+1.d0/((li+1)*(lj+1)))/(den+2.d0)
      d = rnl(1,i)*rnl(1,j)*r(1)**k12
     $     *(((beta*zr+alpha)*zr+1.d0)/(den*h1)+0.5d0)
      dd = 0.d0
      m = min0(nmax(i),nmax(j)) - 1
      do 1 jj = 2,m,2
      jp = jj + 1
      d = d +rnl(jp,i)*rnl(jp,j)*r(jp)**k12
      dd = dd + rnl(jj,i)*rnl(jj,j)*r(jj)**k12
1     continue
      quadr = h1*(d + 2.d0*dd)
      return
      end


*     ------------------------------------------------------------------
*                 q u a d s
*     ------------------------------------------------------------------
*
*                                       kk
*       evaluates the integral of  (1/r)   yk(r) p (r) p (r)  with
*                                                 i     j
*   respect to r.
*
      double precision function  quads2(yk,i,j,k)
      implicit double precision(a-h,o-z)
#include "fortdim.h"
      common /grid/ h,r1,z,r(npdm),np
      common /aorb/ ncore,nval,norb,nlm(n0),
     $     no(n0),lo(n0),nmax(n0),maxim
      common /adat/ wnl(n0),en(n0),so(n0),xion
      common /aval/ rcall(n0),rvap(n0),rnorm(n0),ibd(n0),etot
      common /wfn/ rnl(npdm,n0)
      common /ahf/ s(n0),az(n0),acc(n0),dpm(n0)
      dimension yk(npdm)

      h1=h/1.5
      den = lo(i) + lo(j) + 3
      kk = 2 - k
      cd = 1.d0 + z*r(1)*(den-1.d0)/
     $     ((den+1.d0)*((lo(i)+1)*(lo(j)+1)))
      d = yk(1)*rnl(1,i)*rnl(1,j)*r(1)**kk*( cd/(den*h1)+ 0.5)
      dd = 0.d0

      mx = min0(nmax(i),nmax(j)) - 1
      mx=np-1
      do m = 2,mx,2
         dd = dd + yk(m)*rnl(m,i)*rnl(m,j)*r(m)**kk
         d= d+  yk(m+1)*rnl(m+1,i)*rnl(m+1,j)*r(m+1)**kk
      enddo
      quads2 = h1*(d + 2.d0*dd)
      return
      end

*     ------------------------------------------------------------------
*                       E K I N
*     ------------------------------------------------------------------
*
*       Returns the value of the integral of
*
*         (2/r)P (Y P  + X )
*               j  i i    i
*
*   integrated with respect to r.
*
*
      double precision function ekin(i,ii)
      implicit double precision(a-h,o-z)
#include "fortdim.h"
      common /grid/ h,r1,z,r(npdm),np
      common /aorb/ ncore,nval,norb,nlm(n0),
     $     no(n0),lo(n0),nmax(n0),maxim
      common /adat/ wnl(n0),en(n0),so(n0),xion
      common /aval/ rcall(n0),rvap(n0),rnorm(n0),ibd(n0),etot
      common /wfn/ rnl(npdm,n0)
      common /ahf/ s(n0),az(n0),acc(n0),dpm(n0)
      dimension yk(npdm),yr(npdm),x(npdm)

      call xch(x,i,2)
      call potl(yr,i)
      do 1 j=1,np
         yk(j) = yr(j)
         yr(j) = rnl(j,ii)
1     continue
      ekin = 2.d0*quads2(yk,i,ii,1) + quad(ii,np,yr,x)
      return
      end

*
*     ------------------------------------------------------------------
*               q u a d
*     ------------------------------------------------------------------
*
*       evaluates the integral of f(r)g(r) with respect to r , where
*   f(r) and g(r) have the same asymptotic properties as p (r).   the
*                                                         i
*   composite simpson's rule is used.   the integrand is zero for r >
*   r  .
*    m
*
      double precision function quad(i,m,yr,x)
      implicit double precision (a-h,o-z)

#include"fortdim.h"

      common /grid/ h,r1,z,r(npdm),np
      common /aorb/ ncore,nval,norb,nlm(n0),
     $     no(n0),lo(n0),nmax(n0),maxim
      common /adat/ wnl(n0),en(n0),so(n0),xion
      common /aval/ rcall(n0),rvap(n0),rnorm(n0),ibd(n0),etot
      common /wfn/ rnl(npdm,n0)
      common /ahf/ s(n0),az(n0),acc(n0),dpm(n0)

      dimension yr(npdm),x(npdm)

      h1=h/1.5
      d = (1.d0 + 0.5*z*r(1))/(h1*(2*lo(i) + 3))
      quad = r(1)**2* yr(1)*x(1)*( d -0.5)
      quad2 = 0.d0

      do j = 2,m,2
         rrm=r(j-1)**2
         rr=r(j)**2

         quad = quad + rrm*yr(j-1)*x(j-1)
         quad2 = quad2 + rr*yr(j)*x(j)

      enddo

      quad = h1*(quad + 2.d0*quad2)

      return
      end


*     ------------------------------------------------------------------
*               g r a n g e
*     ------------------------------------------------------------------
*
*       controls the calculation of off-diagonal energy parameters.
*   it searches for all pairs (i,j) which are constrained through  an
*   orthogonality requirement.   eq. (7-10) is used to calculate the
*   parameter.
*
      subroutine grange
      implicit double precision(a-h,o-z)
#include "fortdim.h"

      common /aorb/ ncore,nval,norb,nlm(n0),
     $     no(n0),lo(n0),nmax(n0),maxim
      common /adat/ wnl(n0),en(n0),so(n0),xion
      common /aval/ rcall(n0),rvap(n0),rnorm(n0),ibd(n0),etot
      common /edi/ e(n0,n0)
      common /ahf/ s(n0),az(n0),acc(n0),dpm(n0)
      common /wfn/ rnl(npdm,n0)
      common /consts/ etol,vtol,maxit,isoft
      common /grid/ h,r1,z,r(npdm),np
      common /coeff/coef(100),ijptr(n0,n0)

*  *****  rotate pairs connected by orthogonality but not when one of
*         the orbitals is simultaneously orthogonal to a non-orthogonal
*         pair
*
      do i = 1,norb-1
         do j = i+1,norb
            if (dabs(e(i,j)) .gt. 1.d-10) then
               write(7,*) 'Rotation needed ', i,j,e(i,j)
               if (dabs(e(i,j)) .gt. 1.d-10) call rotate(i,j)
            endif
         enddo
      enddo

*  *****   compute off-diagonal energy parameters

      do i = 2,norb
         do j = 1,i-1
            if (dabs(e(i,j)) .gt. 1.d-10) then
               if ( j .lt. 1) then
                  e(i,j) = hl(i,j) - ekin(i,j)
                  e(j,i) = 0.d0
               else if (wnl(i) .eq. wnl(j)) then
                  c=hl(i,j)-(ekin(i,j)+ekin(j,i))/2.d0
                  e(i,j) = c
                  e(j,i) = c
               else
                  res = 0.d0
                  do  ii = 1,norb
                     if (ii.eq.i .or. ii.eq.j) then
                        do  k = 0,2*lo(i),2
                           if(ii.eq.i) then
                              c= aaa(i,i,k)-aaa(j,i,k)-bbb(j,i,k)
                              if (dabs(c).gt.1.d-10)
     :                             res = res + c*rk(i,i,i,j,k)
                           else if (ii.eq.j) then
                              c= aaa(j,j,k)-aaa(i,j,k)-bbb(i,j,k)
                              if (dabs(c).gt.1.d-10)
     :                             res = res - c*rk(j,j,j,i,k)
                           end if
                        enddo
                     else
                        do k = 0,2*min0(lo(i),lo(ii)),2
                           c = aaa(i,ii,k) - aaa(j,ii,k)
                           if (dabs(c).gt.1.d-10)
     :                          res = res + c*rk(i,ii,j,ii,k)
                           kk = abs(lo(i)-lo(ii)) + k
                           c = bbb(i,ii,kk) - bbb(j,ii,kk)
                           if(dabs(c) .gt. 1.d-10)
     :                          res = res + c*rk(i,ii,ii,j,kk)
                        enddo
                     end if
                  enddo
                  e(i,j) = 2.d0*wnl(j)*res/(wnl(i)-wnl(j))
                  e(j,i) = wnl(i)*e(i,j)/wnl(j)
               end if
            end if
            if (dabs(e(i,j)) .gt. 1.d-10) then
c               write(7,35) nlm(i),nlm(j),e(i,j),nlm(j),nlm(i),e(j,i)
            end if
         enddo
      enddo
 35   format(7x,2(3x,'e(',i3,'|',i3') =',f12.5))
      return
      end

      !this includes excorr terms in off-diag energy params
      subroutine grange1
      implicit double precision(a-h,o-z)
#include "fortdim.h"

      common /aorb/ ncore,nval,norb,nlm(n0),
     $     no(n0),lo(n0),nmax(n0),maxim
      common /adat/ wnl(n0),en(n0),so(n0),xion
      common /aval/ rcall(n0),rvap(n0),rnorm(n0),ibd(n0),etot
      common /edi/ e(n0,n0)
      common /ahf/ s(n0),az(n0),acc(n0),dpm(n0)
      common /wfn/ rnl(npdm,n0)
      common /consts/ etol,vtol,maxit,isoft
      common /grid/ h,r1,z,r(npdm),np
      common /coeff/coef(100),ijptr(n0,n0)

*  *****  rotate pairs connected by orthogonality but not when one of
*         the orbitals is simultaneously orthogonal to a non-orthogonal
*         pair
*
      do i = 1,norb-1
         do j = i+1,norb
            if (dabs(e(i,j)) .gt. 1.d-10) then
               write(7,*) 'Rotation needed ', i,j,e(i,j)
               !if (dabs(e(i,j)) .gt. 1.d-10) call rotate(i,j)
            endif
         enddo
      enddo

*  *****   compute off-diagonal energy parameters

      do i = 2,norb
         do j = 1,i-1
            if (dabs(e(i,j)) .gt. 1.d-10) then
               if ( j .lt. 1) then
                  e(i,j) = hl(i,j) - ekin(i,j)
                  e(j,i) = 0.d0
               else if (wnl(i) .eq. wnl(j)) then
                  c=hl(i,j)-(ekin(i,j)+ekin(j,i))/2.d0
                  e(i,j) = c
                  e(j,i) = c
               else
                  res = ekin(i,j)-ekin(j,i)
                  
                  e(i,j) = 1.d0*wnl(j)*res/(wnl(i)-wnl(j))
                  e(j,i) = wnl(i)*e(i,j)/wnl(j)
               end if
            end if
         enddo
      enddo

      !testing: zeroing off-diags
      !do i = 2,norb
      !   do j = 1,i-1
      !      e(i,j) = 0.0
      !      e(j,i) = 0.0
      !   enddo
      !enddo

      return
      end



      subroutine orthog(ilast)
      implicit double precision (a-h,o-z)
#include "fortdim.h"
      common /aorb/ ncore,nval,norb,nlm(n0),
     $     no(n0),lo(n0),nmax(n0),maxim
      common /adat/ wnl(n0),en(n0),so(n0),xion
      common /aval/ rcall(n0),rvap(n0),rnorm(n0),ibd(n0),etot
      common /edi/ e(n0,n0)
      common /ahf/ s(n0),az(n0),acc(n0),dpm(n0)
      common /wfn/ rnl(npdm,n0)
      common /consts/ etol,vtol,maxit,isoft
      common /grid/ h,r1,z,r(npdm),np
      common /coeff/coef(100),ijptr(n0,n0)
      logical change,ilast

      if (norb .eq. 1) return
      do i = 2,norb
         change = .false.
         azz = az(i)
         do j = 1,i-1
            if (e(i,j) .ne. 0.d0 ) then
*     
*     orthogonality condition applies

               c = quadr(i,j,0)
               if (dabs(c) .gt. 1.d-10.or.ilast) then
                  if (ilast) write(7,63) nlm(j),nlm(i),c
 63               format(6x,'{',i3,'|',i3,'}=',1pd8.1)
                  m = max0(m,nmax(j))
                  do jj = 1,m
                     rnl(jj,i) = rnl(jj,i) - c*rnl(jj,j)
                  enddo
                  azz = azz - c*az(j)
                  change = .true.
               end if
            end if
         enddo
         
         if (change) then
            pnn = dsqrt(quadr(i,i,0))
            if (rnl(1,i) .lt. 0.d0) pnn = - pnn
            do jj = 1,m
               rnl(jj,i) = rnl(jj,i)/pnn
            enddo
            azz = azz/pnn
            m = np
 67         if (dabs(rnl(m,i)) .lt. 1.d-15) then
               rnl(m,i) = 0.d0
               m = m-1
               go to 67
            end if
            nmax(i) = m
            az(i) = azz
         end if
      enddo
      
      return
      end



*
*     ------------------------------------------------------------------
*               r o t a t e
*     ------------------------------------------------------------------
*
*        this routine analyses the energy expression to determine the
*   stationary condition with respect to rotation of orbials i and j.
*   if the condition is zero, the off-diagonal energy parameters may
*   be set to zero;  otherwise the orbials are rotated so as to satisfy
*   the stationay condition to first order in the perturbation.
*
*
      subroutine rotate(i,j)
      implicit double precision (a-h,o-z)
#include "fortdim.h"
      common /aorb/ ncore,nval,norb,nlm(n0),
     $     no(n0),lo(n0),nmax(n0),maxim
      common /adat/ wnl(n0),en(n0),so(n0),xion
      common /aval/ rcall(n0),rvap(n0),rnorm(n0),ibd(n0),etot
      common /edi/ e(n0,n0)
      common /ahf/ s(n0),az(n0),acc(n0),dpm(n0)
      common /wfn/ rnl(npdm,n0)
      common /consts/ etol,vtol,maxit,isoft
      common /grid/ h,r1,z,r(npdm),np
      common /totpot/ rvcore(npdm,n0),rvps(npdm,n0),rvcoul(npdm)
      common /coeff/coef(100),ijptr(n0,n0)
      logical all

      all = .true.
      g = 0.0
      dg = 0.0
      qi = wnl(i)
      qj = wnl(j)
      if (qi .eq. 2.0*(2*lo(i)+1) .and. qj .eq. 2.0*(2*lo(j)+1)) then
         go to 44
      endif

      if (dabs(qi - qj) .lt. 1.d-14) then
         go to 16
      endif


      c = 0.50*(qi - qj)
      g = g -c*hl(i,j)
      dg = dg -c*(hl(i,i) - hl(j,j))

 16   do 13 kf = 0,2*lo(i),2
         c = qi*(aaa(i,i,kf) - aaa(i,j,kf) - bbb(i,j,kf))

         if (dabs(c) .lt. 1.d-8) then 
            go to 21
         endif

         g = g  + c*rk(i,i,i,j,kf)
         fkii = fk(i,i,kf)
         fkij = fk(i,j,kf)
         gkij = gk(i,j,kf)
         dg = dg +c*(fkii - fkij - 2.0*gkij)
 21      cj = qj*(aaa(j,j,kf) - aaa(j,i,kf) - bbb(j,i,kf))

         if (dabs(cj) .lt. 1.d-8) then
            go to 13
         endif

         fkjj = fk(j,j,kf)

         if (dabs(c) .ge. 1.d-8) then
            go to 22
         endif

         fkij = fk(i,j,kf)
         gkij = gk(i,j,kf)
 22      g = g - cj*rk(j,j,j,i,kf)
         dg = dg + cj*(fkjj -fkij - 2.0*gkij)


 13   continue

      do 14 mf = 1,norb
         if (mf .eq. i .or. mf.eq. j) then
            go to 14
         endif

         do 15 kf = 0,2*min0(lo(i),lo(mf)),2
            c = aaa(i,mf,kf)*qi - aaa(j,mf,kf)*qj

            if (dabs(c) .lt. 1.d-8) then
               go to 23
            endif

            g = g + c*rk(i,mf,j,mf,kf)
            dg = dg + c*(fk(i,mf,kf) - fk(j,mf,kf))
 23         kkf = iabs(lo(i)-lo(mf)) + kf
            c = bbb(i,mf,kkf)*qi - bbb(j,mf,kkf)*qj

            if  (dabs(c) .lt. 1.d-8) then 
               go to 15
            endif
            g = g + c*rk(i,j,mf,mf,kkf)
            dg = dg + c*(gk(i,mf,kkf) - gk(j,mf,kkf))
 15      continue
 14   continue

      if (dabs(qi-qj) + dabs(g) + dabs(dg).
     $     le. 2.d-8 ) then 
         go to 44
      endif

 101  if (dabs(g)+dabs(dg) .gt. 1.d-8 .or. dabs(e(i,j)) .gt. 2.d-5) then

         eps = g/dg
         eps = dsign(dmin1(dabs(eps),0.2d0),eps)
         dd = dsqrt(1.d0 + eps*eps)

         do jj = 1,np
            pi = (rnl(jj,i) + eps*rnl(jj,j))/dd
            rnl(jj,j) = (rnl(jj,j) - eps*rnl(jj,i))/dd
            rnl(jj,i) = pi
         enddo
      else
         eps = 0.d0
      end if
c      write (7,100) nlm(i),nlm(j),g,nlm(i),nlm(j),dg,eps
 100  format(10x,'c(',i3,'|',i3,') =',f12.5,3x,'v(',i3,'|',i3') ='
     $     ,f12.5,3x,'eps =',f9.6)
      return

 44   e(i,j) = 1.d-10
      e(j,i) = 1.d-10
      return
      end


      double precision function fk(i,j,k)
      implicit double precision(a-h,o-z)
#include "fortdim.h"

      dimension yk(npdm)

      call ykf(yk,i,i,k)
      fk = quads2(yk,j,j,1)
      return
      end

      double precision function gk(i,j,k)
      implicit double precision(a-h,o-z)
#include "fortdim.h"

      dimension yk(npdm)

      call ykf(yk,i,j,k)
      gk = quads2(yk,i,j,1)
      return
      end

      double precision function rk(i,j,ii,jj,k)
      implicit double precision(a-h,o-z)
#include "fortdim.h"

      dimension yk(npdm)

      call ykf(yk,i,ii,k)
      rk = quads2(yk,j,jj,1)
      return
      end


*
*     ------------------------------------------------------------------
*                 H L
*     ------------------------------------------------------------------
*
*       Returns the value of <i^L^j>, using a special formula to
*  preserve symmetry.
*
      double precision function hl(i,j)
      implicit double precision(a-h,o-z)
#include "fortdim.h"
      common /aorb/ ncore,nval,norb,nlm(n0),
     $     no(n0),lo(n0),nmax(n0),maxim
      common /adat/ wnl(n0),en(n0),so(n0),xion
      common /aval/ rcall(n0),rvap(n0),rnorm(n0),ibd(n0),etot
      common /edi/ e(n0,n0)
      common /ahf/ s(n0),az(n0),acc(n0),dpm(n0)
      common /wfn/ rnl(npdm,n0)
      common /consts/ etol,vtol,maxit,isoft
      common /grid/ h,r1,z,r(npdm),np
      common /totpot/ rvcore(npdm,n0),rvps(npdm,n0),rvcoul(npdm)
      common /coeff/coef(100),ijptr(n0,n0)


      h1=h/1.5
      IF (IABS(LO(I)-LO(J)) .NE. 0) STOP
 3    LI = LO(I)
      C = 2*LI + 1

      A1 = -2.D0/(C*(LI+1))
      A2 = A1/((C+2.D0)*(LI+1))
      A3 = A2/((LI+2)*(LI+1))

      ZR = Z*R(1)

      HL = H*C*RNL(1,I)*RNL(1,J)*(1.D0+ZR*(A1+ZR*(A2+ZR*A3)))

      MM = MIN0(NMAX(I)+3,NMAX(J)+3,NP-1)
      K = 2
      C = 4.D0/3.D0
      DI1 = RNL(K+1,I) - RNL(K-1,I)
      DI2 = RNL(K+1,I) - 2.D0*RNL(K,I) + RNL(K-1,I)
      DJ1 = RNL(K+1,J) - RNL(K-1,J)
      DJ2 = RNL(K+1,J) - 2.D0*RNL(K,J) + RNL(K-1,J)
      HL = HL + DI1*DJ1 + C*DI2*DJ2

      DO 1 K = 4,MM,2
      DI1 = RNL(K+1,I) - RNL(K-1,I)
      DI2 = RNL(K+1,I) - 2.D0*RNL(K,I) + RNL(K-1,I)
      DI4 = RNL(K+2,I) - 4.D0*(RNL(K+1,I)+RNL(K-1,I)) 
     $     + 6.D0*RNL(K,I) +RNL(K-2,I)
      DI3 = RNL(K+2,I) - RNL(K-2,I) - 2.D0*DI1
      DI5 = RNL(K+3,I)-RNL(K-3,I) - 4.D0*(RNL(K+2,I)-RNL(K-2,I))
     :   + 5.D0*(RNL(K+1,I)-RNL(K-1,I))
      DI6 = RNL(K+3,I)+RNL(K-3,I) - 6.D0*(RNL(K+2,I)+RNL(K-2,I))
     :   + 15.D0*(RNL(K+1,I)+RNL(K-1,I)) - 20.D0*RNL(K,I)
      DJ1 = RNL(K+1,J) - RNL(K-1,J)
      DJ2 = RNL(K+1,J) - 2.D0*RNL(K,J) + RNL(K-1,J)
      DJ4 = RNL(K+2,J) - 4.D0*(RNL(K+1,J)+RNL(K-1,J)) 
     $     + 6.D0*RNL(K,J) +RNL(K-2,J)
      DJ3 = RNL(K+2,J) - RNL(K-2,J) - 2.D0*DJ1
      DJ5 = RNL(K+3,J)-RNL(K-3,J) - 4.D0*(RNL(K+2,J)-RNL(K-2,J))
     :   + 5.D0*(RNL(K+1,J)-RNL(K-1,J))
      DJ6 = RNL(K+3,J)+RNL(K-3,J) - 6.D0*(RNL(K+2,J)+RNL(K-2,J))
     :   + 15.D0*(RNL(K+1,J)+RNL(K-1,J)) - 20.D0*RNL(K,J)
1     HL = HL + DI1*DJ1 + C*DI2*DJ2 + (DI3*DJ3 + DI2*DJ4+DI4*DJ2)/45.D0
     :  -(DI3*DJ5+DI5*DJ3)/252.D0 - (DI2*DJ6+DI6*DJ2-1.1*DI4*DJ4)/378.D0
      TZ = 2.D0*RVCORE(1,I)
      C = (LI + 0.5D0)**2
      HL2 = 0.5D0*(TZ*R(1) - C)*RNL(1,I)*RNL(1,J)
      DO 2 K = 2,MM,2
         TZ = 2.D0*RVCORE(K,I)
2     HL2 = HL2 + 2.D0*(TZ*R(K) - C)*RNL(K,I)*RNL(K,J)
     :   + (TZ*R(K+1) - C)*RNL(K+1,I)*RNL(K+1,J)
      HL = -HL/(2.D0*H) + HL2*H1

      return
      end
