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
      subroutine radin(r,f,m,np,h,p)
      implicit double precision (a-h,o-z)
#include "fortdim.h"      
c     *************************************************************************
c     RADIAL INTEGRAL ON LOGARITHMIC GRID- CORRECTED TRAPEZOIDAL RULE
c     THIS METHOD IS PARTICULARLY CONVENIENT AND EFFECTIVE FOR DEF INT--
c     USES ENDPT CORRECTIONS IN TERMS OF DIFFERENCES UP TO FOURTH ORDER.
c     R = R(Y) WHERE Y IS UNIFORM VARIABLE; F = INPUT FUNCTION ON Y GRID;
c     M = POWER OF R USED IF MOMENTS ARE DESIRED (OTHERWISE=0); NP=NUMBER
c     OF GRID POINTS USED IN INTEGRAL; H = UNIFORM SPACING OF Y VARIABLE;
c     P = POWER OF R IN ASYMPTOTIC FORM OF F(R) NEAR R=0: F(R)=A*R**P.
c     ON RETURN, P IS REDEFINED AS VALUE OF DEFINITE INTEGRAL.
c     ASYMPTOTIC FORM ALLOWS ANALYTIC INTEGRAL FROM R=0 TO R=R(1).
c                                                                   DOUG ALLAN
c     *************************************************************************


      
      dimension r(npdm), f(npdm)

c     check if there are enough grid points available      

      if(np.lt.10) then
        write(7,1000) np
 1000   format(1x,'****terminal error in radin:np.lt.10,=',i5)
        stop
      endif
      
c     since this is on logarithmic grid increment the power of r by one      
      n = m+1

c     the first gridpoint is used twice -> hold it in a variable      
      s1 = (r(1)**n)*f(1)
      
c     use asymptotic behavior to approximate the interval [0,r(1)]
      p = s1/(h*(n+p))
            
c     use trapezoidal rule for all inner points [r(6),r(np-5)]     
      i2 = np - 5
      do i=6,i2
        p = p + (r(i)**n)*f(i)
      enddo

c     determine the endpoint corrections to fourth order      
      endpt  = (23.75*(     s1       +f(np  )*(r(np  )**n))
     $        + 95.10*(f(2)*(r(2)**n)+f(np-1)*(r(np-1)**n))
     $        + 55.20*(f(3)*(r(3)**n)+f(np-2)*(r(np-2)**n))
     $        + 79.30*(f(4)*(r(4)**n)+f(np-3)*(r(np-3)**n))
     $        + 70.65*(f(5)*(r(5)**n)+f(np-4)*(r(np-4)**n)))/ 72.
      
c     add end point correction and multiply with width of log. grid
      p = (p + endpt) * h
      
      return
      end
