# 1 "cubint.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "cubint.F"
c
c Copyright (c) 1998-2010 The OPIUM Group
c
c This program is free software; you can redistribute it and/or modify
c it under the terms of the GNU General Public License as published by
c the Free Software Foundation; either version 2 of the License, or
c (at your option) any later version.
c
c This program is distributed in the hope that it will be useful,
c but WITHOUT ANY WARRANTY; without even the implied warranty of
c MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
c GNU General Public License for more details.
c
c You should have received a copy of the GNU General Public License
c along with this program; if not, write to the Free Software
c Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
c
c
      subroutine cubint (iorb,pder,rv,g,e,lang,wfld)
c---------------------------------------------------------------------
c This program will compute the logarithmic derivative for either
c homogeneous or inhomogeneous solutions to the integro-differential
c equation. This program uses a 4 point fit around the rphas, the
c direct grid distance. It interpolates both the derivative of the
c wavefunction and wavefunction.
c NJR 5/28/97
c---------------------------------------------------------------------
      implicit double precision (a-h,o-z)
# 1 "../include/fortdim.h" 1
# 1 "../include/cdim.h" 1
# 2 "../include/fortdim.h" 2
      parameter(npdm = (10000))
      parameter(n0 = (100))
      parameter(npl0 = (3001))
      parameter(nvale0=(100))

      parameter(nout=7)
      parameter(numfn0 = 50)
      parameter(nquad0 = 50)
      parameter(tol2 = 1.0e-9)
      parameter(nfam0 = 10)
      parameter(nintm0 = 2**(nfam0-1))
      parameter(maxflq = 2001)
      parameter(npspt0 = 4002)
      parameter(nrmax = 10000)
      parameter(lmax = 4)
      parameter(norbp = 100)
# 30 "cubint.F" 2

c -------------------------------------------------------------------------
c External (shared between C and Fortran) common blocks
c -------------------------------------------------------------------------
      common /grid/ h,r1,z,r(npdm),np
      common /nlpot2/ inl,indrc(n0),IB,ID,IG
      common /logarith/ rphas,elogmax,elogmin,dlwf(npl0,n0)
      common /aval/ rcall(n0),rvap(n0),rnorm(n0),ibd(n0),etot
c -------------------------------------------------------------------------

      dimension rv(npdm),dl(6),pder(npdm),w(5),g(npdm)
      dimension rvtemp(6),rtemp(6),pdertemp(6),gtemp(6)
c ---------------------------------------------------------------------
      irc = indrc(iorb)
      if (rphas.lt.rcall(iorb)) then
         write(7,*) 'rlog is less than rc for state ',iorb,'.',rphas,
     $ rcall(iorb)
         return
      endif
 5 if(rphas.lt.r(np-3)) goto 15
      write (7,1000) rphas,r(np)
 1000 format(' ****rphas outside of range:rphas,r(np)=',2E20.7)
      return
 15 continue
c---------------------------------------------------------------------
c This section will determine preliminary information regarding the
c exponential grid and determine the 4 nearest linear grid points
c that surround rphas.
c---------------------------------------------------------------------
c plog is the numerical representation of rphas on the linear grid.
c iplog is the nearest integer linear grid point.
c---------------------------------------------------------------------
      plog = 1. + log(rphas/r(1))/h
      iplog = int(plog+0.000001)
c---------------------------------------------------------------------
c determination of the 6 linear grid points (ilog) for the cubic fit.
c array is index from (1,6). (iplog-2) [point 1] and (iplog+3) [6]
c are needed only to seed derivative determination. [point before
c iplog (iplog-1) [2], iplog [3] and 2 points after iplog (iplog+1)
c [4] and (iplog+2) [5] are used for fit].
c---------------------------------------------------------------------
      do 25 j = 1,6
         rvtemp (j) = rv(iplog+j-3)
         gtemp (j) = g(iplog+j-3)
         pdertemp (j) = pder(iplog+j-3)
         rtemp (j) = r(iplog+j-3)
 25 continue
c---------------------------------------------------------------------
c determination of the distance on the linear grid between the actual
c point corresponding to rphas and the nearest integer grid point.
c---------------------------------------------------------------------
      diplog = real(iplog+0.5)
      dplog = plog - diplog
      x1 = 1
      x2 = dplog
      x3 = dplog**2
      x4 = dplog**3
c---------------------------------------------------------------------
c determination of cooefficients of each value of the function to be
c interpolated (pder). The form of the function is assumed according
c to a Lagrange formulation:
c
c f(j) = f_(-1)*(j-j0)*(j-j1)*(j-j2)
c ----------------------------------- +
c (j_(-1)-j0)*(j_(-1)-j1)*(j_(-1)-j2)
c
c f0*(j-j_(-1))*(j-j1)*(j-j2)
c ----------------------------------- +
c (j0-j_(-1))*(j0-j1)*(j0-j2)
c
c f1*(j-j_(-1))*(j-j0)*(j-j2)
c ----------------------------------- +
c (j1-j_(-1))*(j1-j0)*(j1-j2)
c
c f2*(j-j_(-1))*(j-j0)*(j-j1)
c ----------------------------------- +
c (j2-j_(-1))*(j2-j0)*(j2-j1)
c
c The points are taken to be at (-3/2)[iplog-1], (-1/2)[iplog],
c (1/2)[iplog+1], and (3/2)[iplog+2] so that actual linear grid
c value of rphas is near 0.
c---------------------------------------------------------------------
      w(2) = (- 3.*x1+ 2.*x2+12.*x3- 8.*x4)/48.
      w(3) = ( 27.*x1-54.*x2-12.*x3+24.*x4)/48.
      w(4) = ( 27.*x1+54.*x2-12.*x3-24.*x4)/48.
      w(5) = (- 3.*x1- 2.*x2+12.*x3+ 8.*x4)/48.
c---------------------------------------------------------------------
c determination of the derivative at 4 points.
c---------------------------------------------------------------------
      xlang = (real(lang) + 0.5)**2
      h2 = (h**2)/12.
      d1 = ((rvtemp(1)-e*rtemp(1))*rtemp(1)+xlang)*h2*pdertemp(1)+
     $ h2*gtemp(1)
      d2 = ((rvtemp(2)-e*rtemp(2))*rtemp(2)+xlang)*h2*pdertemp(2)+
     $ h2*gtemp(2)
      do 100 j = 2,5
         d3 = ((rvtemp(j+1)-e*rtemp(j+1))*rtemp(j+1)+xlang)*h2*
     $ pdertemp(j+1)+h2*gtemp(j+1)
         dl(j) = (0.5*(pdertemp(j+1)-pdertemp(j-1))-(d3-d1))/h
         d1 = d2
         d2 = d3
  100 continue
c---------------------------------------------------------------------
c cubic interpolation of pder and dl at rphas.
c---------------------------------------------------------------------
      prphas = 0.
      dprphas = 0.
      do 150 j = 2,5
         prphas = prphas + w(j)*pdertemp(j)
         dprphas = dprphas + w(j)*dl(j)
  150 continue
c---------------------------------------------------------------------
c dlp is the logarithmic derivative on the linear grid.
c wfld is the logarithmic derivative on the direct grid.
c---------------------------------------------------------------------
      dlp = dprphas/prphas
      wfld = (0.5 + dlp)/rphas
      end
