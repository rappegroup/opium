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
      subroutine conmaxdriver(xguess)
c     *************************************************************************
c     This code finds an optimized pseudowavefunction following the paper
c       of Rappe et al. (PRB 41,1227 (1990)).
c     Modified from Andrew Rappe's subroutine scr to use the CONMAX
c       minimization routines rather than imsl's DNCONF.
c     The logarithmic derivative matching has been removed on Andrew's
c       advice.
c                                                           06/22/98. N. Hill
c     *************************************************************************

      implicit double precision (a-h,o-z)
      
#include "fortdim.h"

      parameter(maxitn=1000)

      common /roots/ xroot(numfn0)
      common /angm/ ll
      common /numfn/ numfn
      common /wavrc/ wavrc, slope, curvrc, indrc
      common /nnn/ nnn
      common /opt/ meth
      common /grid/ h,r1,z,r(npdm),np
      common /bir/ bir(numfn0,npdm)
      common /frrv/ fr(npdm), rv(npdm)
      common /transum/ transumry

      common /totpot/ rvcore(npdm,n0),rvps(npdm,n0),rvcoul(npdm)
      common /rscore/ rscore(npdm),rdd(npdm),rddd(npdm),rscoretot(npdm),
     $     rsval(npdm)
      common /aorb/ ncore,nval,norb,nlm(n0),
     $     no(n0),lo(n0),nmax(n0),maxim
      common /adat/ wnl(n0),en(n0),so(n0),xion

      common /wfn/ rnl(npdm,n0)
      
      common /scrconmax/ enc,tolc,lim,nstp,isw1,isw2,isw3
      
c     internal common block 
      common /xlog/ ffunc,sumo,sumt,sumn
      
      
c     *************************************************************************
c     local variables
c     *************************************************************************

      dimension xguess(numfn0)
      dimension iwork(7*7 + 7*numfn0 +3)
      dimension work(2*numfn0**2+4*7*numfn0+11*7+27*numfn0+13)
      dimension error(10),fun(1),pttbl(1,1),icntyp(7),confun(7,numfn0+1)

c     *************************************************************************
c     Start the CONMAX section
c     For info. on the CONMAX routines, see the comments in conmax.f or look at
c       http://www.netlib.org/opt/conmax.f
c     numgr is the number of constraints, which will always be equal to 7 for
c       Andrew's method (CONMAX thinks in a very strange way - the main
c       condition is treated as a constraint, and equality constraints
c       are imposed by setting + and - the expression .LE. 0.)
c     *************************************************************************

      nparm = numfn
      numgr = 7
      itlim = maxitn
      ifun = 1
      iptb = 1
      indm = 1
      liwrk = (7*numgr + 7*numfn0 +3)
      lwrk = 2*numfn0**2 + 4*NUMGR*numfn0+ 11*NUMGR + 27*numfn0+ 13
      
      ioptn = 1
      if (isw1.gt.0) ioptn=ioptn+10
      if (isw2.gt.0) ioptn=ioptn+100
      if (isw3.gt.0) ioptn=ioptn+200
      work(1)  = enc
      iwork(1) = lim
      work(2)  = tolc
      iwork(2) = nstp

      call conmax(ioptn,nparm,numgr,itlim,fun,ifun,pttbl,iptb,
     $     indm,iwork,liwrk,work,lwrk,iter,xguess,error)

      call fnset(nparm,numgr,pttbl,iptb,indm,xguess,7,1,icntyp,confun)
     
      write(7,*)
      write(7,*) 'Error in constraints'
      do i=1,numgr
        write(7,9013) i,error(i)
      enddo
      write(7,9013) numgr+1,error(numgr+1)
      write(7,9013) numgr+2,error(numgr+2)
      write(7,9013) numgr+3,error(numgr+3)

 9013 format(1x,i3,e20.5)

c     *************************************************************************
c     end CONMAX section.
c     *************************************************************************

      return
      end

c     #########################################################################
      
      
      subroutine fnset(
     $          NPARM,NUMGR,PTTBL,IPTB,INDM,x,IPT,INDFN,ICNTYP,CONFUN)
     
c     *************************************************************************
c     Subroutine to calculate the primary and secondary minimization
c       constraints to send to CONMAX package.
c     Modified from Andrew Rappe's subroutine fcn2 in file scr.f
c                                                             06/22/98. N. Hill
c     *************************************************************************

      implicit double precision(a-h,o-z)
      
#include "fortdim.h"


      common /bir/ bir(numfn0,npdm)
      common /numfn/ numfn
      common /wavrc/ wavrc, slope, curvrc, indrc
      common /nnn/ nnn
      common /grid/ h,r1,z,r(npdm),np
      common /wfn/ rnl(npdm,n0)
      common /angm/ ll
      common /frrv/ fr(npdm), rv(npdm)
      common /roots/ xroot(numfn0)
      common /a/ a(numfn0)
      common /d/ d(numfn0,numfn0)
      common /e/ e(numfn0)
      common /f/ fint
      common /g/ gint
      common /aorb/ ncore,nval,norb,nlm(n0),
     $     no(n0),lo(n0),nmax(n0),maxim
      common /adat/ wnl(n0),en(n0),so(n0),xion
      common /aval/ rcall(n0),rvap(n0),rnorm(n0),ibd(n0),etot

c     internal common block      
      common /xlog/ ffunc,sumo,sumt,sumn

c     *************************************************************************
c     local variables
c     *************************************************************************

      dimension x(nparm), pttbl(iptb,indm)
      dimension icntyp(7), confun(7,numfn0+1)
      dimension fr2(npdm)
      
c     *************************************************************************
c     Find the function to be minimized, including kinetic, and
c       norm-conservation terms. (Note that the logarithmic derivative
c       terms have been removed since Andrew advised against using them!)
c     x(j) are the Bessel functions alpha_i. In contrast to the original
c       paper, the pseudo-wavefunction is not split into F and C parts.
c     *************************************************************************

c     *************************************************************************
c     First calculate the pseudo-wavefunction (sum) and its Laplacian (sum2)
c       at rc, qi. (xroot is qi - the root of the transcendental equation).
c     *************************************************************************

      sum = 0.0
      sum2 = 0.0
      do j = 1,nparm
        sum = sum + x(j) * bir(j,indrc)
        sum2 = sum2 + x(j) * bir(j,indrc) * xroot(j) * abs(xroot(j))
c        write(7,*) j, x(j), bir(j,indrc), xroot(j)
      enddo
 
c -- gjt: this is debugging stuff ---
c      write(7,*) 'sum, sum2', sum, sum2    
c      stop

c     *************************************************************************
c     The first two constraints are the continuity of the wavefunction and
c       its derivatives. This is achieved by requiring that the pseudo-
c       and all-electron terms are equal at rc. For the weird formatting
c       of the CONMAX routine, this is achieved by setting the difference
c       and minus the difference each .LE. 0.0
c     *************************************************************************

      g2 = sum - wavrc
      g3 = sum2 + curvrc
      if (ipt.eq.1) then
        confun(ipt,1)=g2
        icntyp(ipt)=-2
        return
      end if
      if (ipt.eq.2) then
        confun(ipt,1)=-g2
        icntyp(ipt)=-2
        return
      end if
      if (ipt.eq.3) then
        confun(ipt,1)=g3
        icntyp(ipt)=-2
        return
      end if
      if (ipt.eq.4) then
        confun(ipt,1)=-g3
        icntyp(ipt)=-2
        return
      end if

c -- gjt: this is debugging stuff ---
c      write(7,*) 'g2, g3', g2, g3

c     *************************************************************************
c     Now construct the trial pseudo-wavefunction fr over the whole grid.
c     *************************************************************************

c     Up to rc this is the sum over Bessel functions:
      
      f = 0.0
      do i = 1,indrc
        fr(i) = 0.0
        do j = 1,numfn
          fr(i) = fr(i) + x(j) * bir(j,i)
        enddo
        fr2(i) = fr(i) * fr(i) * r(i) * r(i)
      enddo
      
c     Beyond rc this is the all-electron wavefunction:

      do i = indrc+1,np
        fr(i) = rnl(i,nnn)
        fr2(i) = fr(i) * fr(i) * r(i) * r(i)
      enddo
      
      xn9 = float(ll+ll+2)
      call radin(r,fr2,0,np,h,xn9)
      
c -- gjt: this is debugging stuff ---
c      do i=1, 20
c        write(7,*) 'fr2(',i,')=',fr2(i)
c      enddo
c      write(7,*) 'xn9=',xn9
      
c     *************************************************************************
c     xn9 is the norm of fr - this creates an additional constraint
c       (xn9 must equal 1!)
c     *************************************************************************

      g1 = xn9 - 1.0
      sumn = g1**2
      if (ipt.eq.5) then
        confun(ipt,1)=g1
        icntyp(ipt)=-2
        return
      end if
      if (ipt.eq.6) then
        confun(ipt,1)=-g1
        icntyp(ipt)=-2
        return
      end if

c     *************************************************************************
c     Here we find the part kinetic energy.
c     *************************************************************************

      sumt = 0.0
      do i = 1,numfn
        sum = 0.0
        do j = 1,numfn
          sum = sum + x(j) * d(i,j)
        enddo
        
c       Here xroot*abs(xroot) again corrects for xroot(1)<0 exception.
        fact = x(i) * xroot(i)*abs(xroot(i)) * a(i) - 2.0 * e(i) - sum
     
        sumt = sumt + fact * x(i)
      enddo
 
      sumt = sumt - fint - gint
      sumo = f
      f = f + sumt 

      if (f.lt.0.0) then
        write (7,*) 'less than zero',f
c        write (7,*) fint,gint,(x(i),i=1,numfn),(xroot(i),i=1,numfn),
c     $              (e(i),i=1,numfn)
c        write(7,*) (a(i),i=1,numfn),((d(i,j),i=1,numfn),j=1,numfn)


        write(7,9595)
        write(7,9596) sum2,sum2-f,f
        write(7,9597) 
        write(7,*)
c        write(6,9595)
c        write(6,9596) sum2,sum2-f,f
c        write(6,9597) 
c        write(6,*)
c        stop
 9595   format(1x,
     $   '!!WARNING!! Resid became negative (Total,lowq,highq)')
 9596   format(1x,3f20.10)
 9597   format(1x,
     $       'Consider reducing qc or changing the # bess fxns',
     $       ' for this state...')

      endif
      
      ffunc = f

c -- gjt: this is debugging stuff ---
c      write(7,*) 'stuff=',sumt,wke,fint,gint,ffunc
      if (ipt.eq.7) then
        confun(ipt,1)=ffunc
        icntyp(ipt)=1
      end if


c     zero curvature for V at 0
c     EJW
c      ffunc2=0
c      do i=1,numfn
c         do j=1,numfn
c            ffunc2=ffunc2+x(i)*x(j)*xroot(j)**ll * xroot(i)**ll *
c     $           (xroot(i)*abs(xroot(i))-
c     $           xroot(j)*abs(xroot(j)))**2
c         enddo
c      enddo

c      if (ipt.eq.8) then
c         confun(ipt,1)=ffunc2
c         icntyp(ipt)=-2
c      end if
c      if (ipt.eq.9) then
c         confun(ipt,1)=-ffunc2
c         icntyp(ipt)=-2
c      end if


      return
      end
