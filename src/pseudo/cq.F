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
      subroutine cq
      
c     *************************************************************************
c     Here we integrate sqrt(2/pi) * phi(r) * jl(qr) from rc to infinity
c     Because the Bessel functions have many wiggles between rc and 
c     r(nmax(nnn)) (the point where phi(r) = 0), the integral must be done
c     carefully.  We monitor convergence as a function of the number of
c     regions in which quadratures are applied.
c     *************************************************************************

      implicit double precision (a-h,o-z)
      
#include "fortdim.h"

      common /grid/ h,r1,z,r(npdm),np
      common /aorb/ ncore,nval,norb,nlm(n0),
     $     no(n0),lo(n0),nmax(n0),maxim
      common /adat/ wnl(n0),en(n0),so(n0),xion
      common /wfn/ rnl(npdm,n0)

      common /gauss/ xquad(nquad0),wquad(nquad0)
      common /nnn/ nnn
      common /cuts/ qc,rc
      common /angm/ ll

      common /quads/ qq(nquad0)
      common /c/ c(nquad0)
      
c     internal common block      
      common /vv/xf(nquad0,nintm0,nfam0),f(npdm),vv(nquad0,nintm0,nfam0)


      pi = 3.14159265358979323
      xcon1 = sqrt(2.0/pi)
      do 4 i = 1,np
         f(i) = rnl(i,nnn)
 4    continue

      nfammx = 1
      call valstor(nfammx)
      do 99 k = 1,nquad0
         zsumo = 0.0
         nfam = 1
 3       continue
         nint = 2**(nfam-1)
         width = (r(nmax(nnn)) - rc)/float(nint)
         zsum = 0.0
         do 1 i = 1,nint
            do 2 j = 1,nquad0
               xf0 = xf(j,i,nfam)
               bb = besfn(qq(k)*xf0,ll)
               zsum = zsum + vv(j,i,nfam) * bb * xf0 * xf0 * wquad(j)
 2          continue
 1       continue
         zsum = zsum * width/2.0
          if (abs(zsum - zsumo).gt.tol2) then
            zsumo = zsum
            nfam = nfam + 1
            if (nfam.gt.nfammx) then
               nfammx = nfam
               call valstor(nfammx)
            endif
            goto 3
         endif
         c(k) = zsum * xcon1
 99   continue
      zsumo = 0.0
      nfam = 1
 7    continue
      nint = 2**(nfam-1)
      width = (r(nmax(nnn)) - rc)/float(nint)
      zsum = 0.0
      do 5 i = 1,nint
         do 6 j = 1,nquad0
            xf0 = xf(j,i,nfam)
            vv0 = vv(j,i,nfam)
            zsum = zsum + vv0 * vv0 * xf0 * xf0 * wquad(j)
 6       continue
 5    continue
      zsum = zsum * width/2.0
      if (abs(zsum - zsumo).gt.tol2) then
         zsumo = zsum
         nfam = nfam + 1
         if (nfam.gt.nfammx) then
            nfammx = nfam
            call valstor(nfammx)
         endif
         goto 7
      endif
      xnorm = 1.0 - zsum
      return
      end

      
c     #########################################################################


      subroutine valstor(nfammx)
      implicit double precision(a-h,o-z)
      
#include "fortdim.h"

      common /grid/ h,r1,z,r(npdm),np
      common /aorb/ ncore,nval,norb,nlm(n0),
     $     no(n0),lo(n0),nmax(n0),maxim
      common /adat/ wnl(n0),en(n0),so(n0),xion
      common /wfn/ rnl(npdm,n0)
      common/gauss/ xquad(nquad0),wquad(nquad0)
      common /nnn/ nnn
      common/cuts/ qc,rc

      common /vv/xf(nquad0,nintm0,nfam0),f(npdm),vv(nquad0,nintm0,nfam0)

      nord=10
      if (nfammx.gt.nfam0) then
         write (7,*) 'xf and vv not dimensioned large enough'
         write (7,*) 'nfam0 nfammx',nfam0,nfammx
         stop
      endif
      nint = 2**(nfammx-1)
      width = (r(nmax(nnn))-rc)/float(nint)
      do 10 i = 1,nint
         xlend = float(i-1) * width + rc
         do 11 j = 1,nquad0
            xf0 = (xquad(j) + 1.0)*width/2.0 + xlend
            xf(j,i,nfammx) = xf0
            vv(j,i,nfammx) = val(f,r,np,xf0,nord)
 11      continue
 10   continue
      
      return
      end
