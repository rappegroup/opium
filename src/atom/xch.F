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
      subroutine xch(x,i,iopt)
     
c     *************************************************************************
c     get the hartree potential function
c     *************************************************************************
     
      implicit double precision (a-h,o-z)

#include "fortdim.h"
      common /grid/ h,r1,z,r(npdm),np
      common /aorb/ ncore,nval,norb,nlm(n0),
     $     no(n0),lo(n0),nmax(n0),maxim
      common /adat/ wnl(n0),en(n0),so(n0),xion
      common /aval/ rcall(n0),rvap(n0),rnorm(n0),ibd(n0),etot
      common /edi/ e(n0,n0)
      common /wfn/ rnl(npdm,n0)

      dimension yk(npdm),yr(npdm),x(npdm)

      dimension meth(10)

      do  in=1,np
         x(in) = 0.0
      enddo

      do  j=1,norb
         if (j .ne. i) then
            do  k=iabs(lo(i)-lo(j)),lo(i)+lo(j),2
               if (dabs(bbb(i,j,k)) .ge. 1.d-10) then
                  call ykf(yk,i,j,k)
c                  write(7,*) 'xc', i,j,k,bbb(i,j,k)
                  do in=1,np
                     x(in) = x(in)+ 2.d0*bbb(i,j,k)*yk(in)*rnl(in,j)
                  enddo
               endif
            enddo
            
         endif
      enddo

      if (iopt.eq.2) then
         do in=1,np
            x(in) = x(in)/r(in)
         enddo
      endif

      if (iopt.eq.3) then
         do in=1,np
            x(in) = x(in)*r(in)
         enddo
      endif
      
      if (iopt.eq.3) then
         do j = 1,norb
            c = e(i,j)
            if (dabs(c).gt.1.d-10.and.(j.ne.i)) then
               do in = 1,np
                  x(in) = x(in) + c*rnl(in,j)*r(in)*r(in)
               enddo
            endif
         enddo
      endif

*  *****  check if exchange is zero: if so, method 2 should be used.
      if (meth(i) .eq. 2) return
      if ( dabs(x(1)) + dabs(x(2)) + dabs(x(3)) .eq. 0.d0 ) meth(i) = 2

      return
      end
      
