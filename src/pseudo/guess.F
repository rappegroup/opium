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
      subroutine guess(xg)
      implicit double precision (a-h,o-z)
      external fact
      
#include "fortdim.h"

      common /grid/ h,r1,z,r(npdm),np

      common /aorb/ ncore,nval,norb,nlm(n0),
     $     no(n0),lo(n0),nmax(n0),maxim
      common /adat/ wnl(n0),en(n0),so(n0),xion
      common /wfn/ rnl(npdm,n0)

      common /gauss/ xquad(nquad0),wquad(nquad0)
      common /roots/xroot(numfn0)
      common /angm/ ll
      common /nnn/ nnn
      common /np/ ncores,nvales
      common /numfn/ numfn
      common /wavrc/ wavrc, slope, curvrc, indrc
      common /a/ a(numfn0)
      common /bir/ bir(numfn0,npdm)
      common /d/ d(numfn0,numfn0)
      common /e/ e(numfn0)
      common /f/ fint
      common /g/ gint

c     EJW newopt variables
      common /rke/ rkmat(numfn0,numfn0),rkvec(numfn0)
      common /rn/ rnmat(numfn0,numfn0),rnvec(numfn0)
      common /rlm/ rlmat(3,numfn0),rlvec(numfn0)
      common /rconst/ rkcon,rncon,rlcon(3)
      common /re/ rkin,rnorm
      common /emat/ emat(3,3),rfl(numfn0)
      common /rf/ afor(3),rlam(3)

      dimension p(npdm),xg(numfn0),xpt(nquad0)
      dimension ptail(npdm),rov(npdm)

      do 15 i = 1,np
         p(i) = (rnl(i,nnn) * r(i))**2
         ptail(i)=p(i)
 15   continue
      
      tov = float(ll+ll+2)
      call radin(r,p,0,np,h,tov)
c      write (7,*) 'test integral ',tov
      xxx = fact(ll+ll+2)
c      write (7,*) 'fact ',xxx,ll,ll+ll+2
      astep = 0.2
      idir = 1
      iter = 1
      itmax = 100
      p0 = 1.0
      p1 = r(indrc)
      p2 = p1 * p1
      p3 = p2 * p1
      g0 = wavrc
      g1 = slope
      g2 = curvrc - 2.0 * g1/p1 + float(ll * (ll+1)) * g0/p2
      h0 = log(g0/p1**ll)
      h1 = g1/g0 - float(ll)/p1
      h2 = g2/g0 - g1*g1/g0/g0 + float(ll)/p2
      a3 = 1.0
 1    continue
      a2 = (h2 - 6.0 * p1 * a3)/2.0
      a1 = h1 - 3.0 * p2 * a3 - 2.0 * p1 * a2
      a0 = h0 - p1 * a1 - p2 * a2 - p3 * a3
c For this choice, find the normalization from 0 to rc.
      do 2 i = 1,indrc
         p(i) = exp((a0 + a1 * r(i) + a2 * r(i) * r(i) + a3 * r(i)**3)
     $        *2.0) * r(i) ** (ll + ll + 2)
 2    continue
      tov = float(ll+ll+2)
      call radin(r,p,0,np,h,tov)
c      write (7,*) 'integral ',tov,a0,a1,a2,a3
      if (abs(tov-1.0).gt.1.0e-6) then
         if (tov.gt.1.0) idir2 = -1
         if (tov.lt.1.0) idir2 =  1
         if (idir2.ne.idir) astep = astep/2.0
         if (idir2.eq.idir) astep = astep * 1.2
         idir = idir2
         a3 = a3 - float(idir) * astep
         iter = iter + 1
         if (iter.eq.itmax) then
            write (7,*) 'Could not find starting guess in '
            write (7,*) 'guess.f Therefore we stop. '
         endif
         goto 1
      endif
      do 11 i = 1,np
         p(i) = sqrt(p(i))/r(i)
c         write (7,*) 'wavefcn ',i,r(i),p(i)
 11   continue
c Now p(i) contains a nodeless wavefunction with norm of 1.
c Expand from 0 to rc in the set of Bessel functions given.

      do 105 i = 1,nquad0
         xpt(i) = (xquad(i) + 1.0) * r(indrc)/2.0
 105  continue
      do 100 i = 1,numfn
         sum = 0.0
         do 101 j = 1,nquad0
            y = xpt(j) * xroot(i)
            x = xpt(j)
            fx = exp(a0 + a1 * x + a2 * x * x + a3 * x**3) * x**ll
            sum = sum + besfn(y,ll) * fx * x**2 * wquad(j)
 101     continue
         sum = sum * r(indrc)/2.0
         xg(i) = sum/a(i)
 100  continue
      sum = 0.0
      sum2 = 0.0
      do 110 i = 1,numfn
         sum = sum + xg(i) * bir(i,indrc)
         sum2 = sum2 + xg(i) * bir(i,indrc) * xroot(i) * abs(xroot(i))
 110  continue

      do 102 i = 1,indrc
         sum = 0.0
         do 103 j = 1,numfn
            sum = sum + xg(j) * bir(j,i)
 103     continue
 102  continue

      do i=indrc+1,np
         rov(i)=0.0
      enddo
      do i=1,numfn
         do j=1,numfn
            do n=1,indrc
               rov(n)=r(n)*r(n)*bir(i,n)*bir(j,n)
            enddo
            tov = float(ll+ll+2)
            call radin(r,rov,0,np,h,tov)
            rnmat(i,j)=tov
         enddo
      enddo

      do i=1,indrc
         ptail(i)=0.0
      enddo
      tov = float(ll+ll+2)
      call radin(r,ptail,0,np,h,tov)
      rncon=tov-1.0

      do i=1,numfn
         do j=1,numfn
            rkmat(i,j)=-d(i,j)
            if (i.eq.j) then
               rkmat(i,j)=rkmat(i,j)+xroot(i)*abs(xroot(i))*a(i)
            endif
         enddo
         rkvec(i)=-2.0*e(i)
      enddo

      rkcon=-fint-gint

      do i=1,numfn
         rnvec(i)=0.0
         rlmat(1,i)=bir(i,indrc)
         rlmat(2,i)=bir(i,indrc)*xroot(i)*abs(xroot(i))
      enddo

      rlvec(1)=-wavrc
      rlvec(2)=curvrc

      return
      end

      function fact(i)
      implicit double precision(a-h,o-z)
      fact = 1
      do 1 j = 1,i
         fact = fact * j
 1    continue
      return
      end
