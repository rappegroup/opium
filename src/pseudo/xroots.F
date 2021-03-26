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
      subroutine xroots
      
c     *************************************************************************
c     find the wavevectors at which the logarithmic derivative of the 
c     spherical Bessel function of proper angular momentum is the same as
c     the all-electron wavefunction.
c     Use Newton's method.  
c     *************************************************************************
      
      implicit double precision(a-h,o-z)
      
#include "fortdim.h"

      common /angm/ ll
      common /roots/ xroot(numfn0)
      common /wavrc/ wavrc, slope, curvrc, indrc
      common /cuts/ qc,rc
      common /numfn/ numfn
      
c     internal common block      
      common /xcon1/ xcon1

      irtgr = numfn * 10
      xlder = slope/wavrc
      xcon1 = xlder * rc
      pi = 3.14159265358979323d0
      j = 0
      
      if (xlder.gt.float(ll)/rc) then
c       ***********************************************************************
c       We must look for one additional root for xguess<0.
c       This corresponds to looking for il(x)=i**-l*jl(ix).
c       The exponentially growing solution exists for these cases.
c       ***********************************************************************
        xtrial = -10. * xlder - 0.1
        xnow = -tol2
        call findrt(xtrial,xnow)
        if (xtrial.eq.-9999.0d0) then
          write (7,*) 'Newton solver jumping boundaries.'
          stop
        endif
        j = j + 1
        xroot(j) = xtrial/rc
      endif
      
      xtop = dfloat(numfn + 1 + ll) * pi
      xstep = xtop/dfloat(irtgr)
      xnow = 0.0
      test = fun0(tol2)
      
      do i = 1,irtgr
        if (j.eq.numfn) goto 1
        xnow = xnow + xstep
        test2 = fun0(xnow)
        xtrial = xnow - xstep
        if (test*test2.lt.0.0) then
          call findrt(xtrial,xnow)
          if (xtrial.eq.-9999.0d0) then
            write (7,*) 'Newton solver jumping boundaries.'
            stop
          endif
          j = j + 1
          xroot(j) = xtrial/rc
        endif
        test = test2
      enddo
 1    continue      
      
      if (j.ne.numfn) then
        write (7,*) 'j = ',j,'numfn = ',numfn,'so we missed roots'
        stop
      endif
      
      return
      end
      
      
c     #########################################################################
      
      
      subroutine findrt(xlo,xhi)
      
c     *************************************************************************
c     perform the Newton's method search.  If a large jump is recommended,
c     we scale it down so that we converge to the right root.
c     *************************************************************************

      implicit double precision(a-h,o-z)
      
#include "fortdim.h"

      delmax = (xhi - xlo)/10.0
      xguess = (xlo + xhi)/2.0
      iter = 0
 2    x = fun0(xguess)
      y = fun1(xguess)
      iter = iter + 1
      if (iter.gt.500) then
         write (7,*) 'cannot find Bessel function wavevectors.'
         stop
      endif
      if (abs(x).gt.tol2) then
         del = x/y
         if (abs(del).gt.delmax) del = abs(del)/del*delmax
         xguess = xguess - del * abs(xguess)/xguess
         goto 2
      endif
      if ((xguess - xlo) * (xguess - xhi).lt.0.0d0) then
         xlo = xguess
      else
         xlo = -9999.0d0
      endif
      return
      end
      
      
C**********************************************************************C
      function fun0(x)
C**********************************************************************C
      implicit double precision(a-h,o-z)
      common/xcon1/ xcon1
      common/angm/ ll
      bl = besfn(x,ll)
      bd = besder(x,ll)
      fun0 = -bl * xcon1 + bd * abs(x)
      if (abs(xcon1).gt.1.0d0) then
         fun0 = -bl + bd * abs(x)/xcon1
      endif
      return
      end
      
      
C**********************************************************************C
      function fun1(x)
C**********************************************************************C
      implicit double precision(a-h,o-z)
      common/xcon1/ xcon1
      common/angm/ ll
      bd = besder(x,ll)
      bdd = besder2(x,ll)
      fun1 = -bd * xcon1 + bdd * abs(x) + bd
      if (abs(xcon1).gt.1.0d0) then
         fun1 = -bd + bdd * abs(x)/xcon1 + bd/xcon1
      endif
      return
      end
