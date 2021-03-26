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
      subroutine hfcoef(iflag,ipsp)
      implicit double precision (a-h,o-z)

#include "fortdim.h"
      common /grid/ h,r1,z,r(npdm),np
      common /aorb/ ncore,nval,norb,nlm(n0),
     $     no(n0),lo(n0),nmax(n0),maxim
      common /adat/ wnl(n0),en(n0),so(n0),xion
      common /aval/ rcall(n0),rvap(n0),rnorm(n0),ibd(n0),etot
      common /wfn/ rnl(npdm,n0)
      common /coeff/ coef(100),ijptr(n0,n0)
      common /eav/cca(10),ccb(35)


      ip=0
*  *****  average interactions for equivalent electrons
*
*  *****  p - p
*
      cca(1) = 2.d0/25.d0
*
*  *****  d - d
*
      cca(2) = 2.d0/63.d0
      cca(3) = 2.d0/63.d0
*
*  *****  f - f
*
      cca(4) =   4.d0/ 195.d0
      cca(5) =   2.d0/ 143.d0
      cca(6) = 100.d0/5577.d0
*
*  *****  g - g
*
      cca(7) =   20.d0/  1309.d0
      cca(8) =  162.d0/ 17017.d0
      cca(9) =   20.d0/  2431.d0
      cca(10) = 4410.d0/371943.d0
*
*
*  ***** average interactions for non-equivalent electrons
*
*  *****  s - ( s, p, d, f, g )
*
      ccb(1) = 1.d0/ 2.d0
      ccb(2) = 1.d0/ 6.d0
      ccb(3) = 1.d0/10.d0
      ccb(4) = 1.d0/14.d0
      ccb(5) = 1.d0/18.d0
*
*  *****  p - ( p, d, f, g )
*
      ccb(6) = 1.d0/  6.d0
      ccb(7) = 1.d0/ 15.d0
      ccb(8) = 1.d0/ 15.d0
      ccb(9) = 3.d0/ 70.d0
      ccb(10) = 3.d0/ 70.d0
      ccb(11) = 2.d0/ 63.d0
      ccb(12) = 2.d0/ 63.d0
      ccb(13) = 5.d0/198.d0
*
*  *****  d - ( d, f, g )
*
      ccb(14) =  1.d0/ 10.d0
      ccb(15) =  1.d0/ 35.d0
      ccb(16) =  1.d0/ 35.d0
      ccb(17) =  3.d0/ 70.d0
      ccb(18) =  2.d0/105.d0
      ccb(19) =  5.d0/231.d0
      ccb(20) =  1.d0/ 35.d0
      ccb(21) = 10.d0/693.d0
      ccb(22) =  5.d0/286.d0
*
*  *****  f - ( f, g )
*
      ccb(23) =  1.d0/  14.d0
      ccb(24) =  2.d0/ 105.d0
      ccb(25) =  1.d0/  77.d0
      ccb(26) = 50.d0/3003.d0
      ccb(27) =  2.d0/  63.d0
      ccb(28) =  1.d0/  77.d0
      ccb(29) = 10.d0/1001.d0
      ccb(30) = 35.d0/2574.d0
*
*  *****  g - ( g )
*
      ccb(31) =   1.d0/   18.d0
      ccb(32) =  10.d0/  693.d0
      ccb(33) =   9.d0/ 1001.d0
      ccb(34) =  10.d0/ 1287.d0
      ccb(35) = 245.d0/21879.d0

      do i=ncore+1,norb
         isumi = wnl(i)
         dsumi = wnl(i) - isumi
         do j = ncore+1,norb
            isumj = wnl(j)
            dsumj = wnl(j) - isumj

            if ( i .ne. j) then
               c = wnl(j)
               if (dsumi .ne. 0.d0 .and. dsumj .ne. 0.d0)
     :              c = (dsumi*(isumi+1)*isumj +
     :              dsumj*(isumj+1)*isumi)/wnl(i)
            else
               c = wnl(i) - 1.d0
               if (dsumi .ne. 0.d0)
     :              c = (isumi*(wnl(i)+dsumi-1))/wnl(i)
               if (wnl(i).lt.1e-8) c=0.d0
            end if

            ijptr(i-ncore,j-ncore) = ip
       
c     Hartree part

            do k = 0,2*min(lo(i),lo(j)),2
               ip = ip + 1
               if (ip .gt. (100))
     :              stop ' coefh array too small: max = (100)'
               coef(ip) = 0.d0
               if (k .eq. 0) then
                  coef(ip) = c
               else if (i .eq. j) then
                  coef(ip) = -c*ca(lo(i),k)
               end if
            enddo

c     Exchange part

            if (i .ne. j) then
               do k = abs(lo(i)-lo(j)),lo(i)+lo(j),2
                  ip = ip + 1
                  if (ip .gt. (100))
     :                 stop ' coefx array too small: max = (100)'
                  coef(ip) = -c*cb(lo(i),lo(j),k)
               enddo
            endif

      
         enddo
      enddo

      return
      end


      double precision function aaa(i,j,k)
      implicit double precision(a-h,o-z)
#include "fortdim.h"

      common /grid/ h,r1,z,r(npdm),np
      common /aorb/ ncore,nval,norb,nlm(n0),
     $     no(n0),lo(n0),nmax(n0),maxim
      common /adat/ wnl(n0),en(n0),so(n0),xion
      common /aval/ rcall(n0),rvap(n0),rnorm(n0),ibd(n0),etot
      common /coeff/coef(100),ijptr(n0,n0)

      if (i.gt.ncore .and. j.gt.ncore) then
         istart = ijptr(i-ncore,j-ncore) + 1
         aaa = coef(istart + k/2)
      else if (i.eq.j) then
         c = wnl(i) - 1.d0

         if (k.eq.0) then
            aaa = c
         else
            aaa = -c*ca(lo(i),k)
         end if
      else if (k.eq.0) then
         aaa= wnl(j)
      else
         aaa = 0.d0
      endif
c      if ((i.eq.2).and.(j.eq.2).and.(k.eq.2)) aaa=-0.4
c      if ((i.eq.3).and.(j.eq.3).and.(k.eq.2)) aaa=-0.4

      end


      double precision function bbb(i,j,k)
      implicit double precision(a-h,o-z)
#include "fortdim.h"
      common /grid/ h,r1,z,r(npdm),np
      common /aorb/ ncore,nval,norb,nlm(n0),
     $     no(n0),lo(n0),nmax(n0),maxim
      common /adat/ wnl(n0),en(n0),so(n0),xion
      common /aval/ rcall(n0),rvap(n0),rnorm(n0),ibd(n0),etot
      common /coeff/coef(100),ijptr(n0,n0)

      if (i.eq.j) then
         b = 0.d0
      else if (i.gt.ncore .and. j.gt.ncore) then
*
*   ..... ll is the number of direct terms
*         istart the beginning of the exchange terms
*
         ll = min(lo(i),lo(j)) + 1
         istart = ijptr(i-ncore,j-ncore) + 1 + ll
         kkk = (k - abs(lo(i)-lo(j)))/2
         bbb = coef(istart + kkk)
      else
         bbb = -wnl(j)*cb(lo(i),lo(j),k)
      endif
      end

      double precision function ca(l,k)
      implicit double precision(a-h,o-z)
      common /eav/cca(10),ccb(35)
      

      if (l .le. 4) then
         ca = cca((l*(l-1) + k)/2)
      else
         stop 'max l = 4'
      endif

      end


      double precision function cb(l,lp,k)
      implicit double precision(a-h,o-z)
      common /eav/cca(10),ccb(35)
      integer icbptr(0:4)
      data    icbptr/1,6,14,23,31/

      if (l .le. lp) then
         l1 = l
         l2 = lp
      else
         l1 = lp
         l2 = l
      end if

      if ( l2 .le. 4) then
         cb = ccb(icbptr(l1)+(k+l1-l2)/2+(l1+1)*(l2-l1))
      else
         stop 'max l = 4'
      endif

      end

          
