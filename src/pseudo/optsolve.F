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
      subroutine optsolve(ixc)

c     *************************************************************************
c     This code finds an optimized pseudowavefunction following the paper
c     of Rappe et al. (PRB 41,1227 (1990)).
c     *************************************************************************
c     Now this routine can operate without conmax using the 'newmin' routine
c     which was coded by EJW and HK.     10/25/03 EJW

      implicit double precision (a-h,o-z)
      
#include "fortdim.h"

      parameter(maxitn=1000)

      logical lbes
      common /a/ a(numfn0)
      common /d/ d(numfn0,numfn0)
      common /e/ e(numfn0)
      common /f/ fint
      common /g/ gint
      common /roots/ xroot(numfn0)
      common /angm/ ll
      common /numfn/ numfn
      common /wavrc/ wavrc, slope, curvrc, indrc
      common /nnn/ nnn
      common /psdat/ nll,meth,npots
      common /grid/ h,r1,z,r(npdm),np

      common /rscore/ rscore(npdm),rdd(npdm),rddd(npdm),rscoretot(npdm),
     $     rsval(npdm)
      common /aorb/ ncore,nval,norb,nlm(n0),
     $     no(n0),lo(n0),nmax(n0),maxim
      common /adat/ wnl(n0),en(n0),so(n0),xion
      common /aval/ rcall(n0),rvap(n0),rnorm(n0),ibd(n0),etot
      common /wfn/ rnl(npdm,n0)

      common /bir/ bir(numfn0,npdm)
      common /frrv/ fr(npdm), rv(npdm)
      common /transum/ transumry

      common /totpot/ rvcore(npdm,n0),rvps(npdm,n0),rvcoul(npdm)
      
      common /scrconmax/ enc,tolc,lim,nstp,isw1,isw2,isw3
      
      common /psout/ wsumt(n0),sumc(n0),npsghost(n0)
      
c     internal common block 
      common /xlog/ ffunc,sumo,sumt,sumn
      
c     *************************************************************************
c     local variables
c     *************************************************************************

      dimension fr2(npdm),xguess(numfn0),v0(npdm)
      dimension iwork(7*7 + 7*numfn0 +3)
      dimension work(2*numfn0**2+4*7*numfn0+11*7+27*numfn0+13)
      dimension error(10),fun(1),pttbl(1,1),icntyp(7),confun(7,numfn0+1)

      do i = 1,np
        rv(i) = -z-z+rvcoul(i)
        if (i.gt.maxim) rv(i) = -xion-xion
        if (nnn.eq.1) rsval(i)=0.0
      enddo

      do j = 1,numfn
        do i = 1,np
          x = xroot(j) * r(i)
          bir(j,i) = besfn(x,ll)
        enddo
      enddo
 
      do i = 1,numfn0
        xguess(i) = 0.0
      enddo

      call guess(xguess)

cEJW added this for internal debugging (lets you load bessel functions)

      lbes = .false.
      inquire(file='BESSEL',exist=lbes)
      if (lbes) then
        open(unit=11,file='BESSEL',form='formatted')
        do i=1,numfn
         read(11,*,end=911) ijunk,xjunk,xguess(i)
        enddo
        call fnset(numfn,numgr,pttbl,iptb,indm,xguess,7,1,icntyp,confun)
        write(7,*) 'Using Bessel file'
        goto 912
 911    stop'Bessel file broken!'
      endif

      if (meth.eq.0) then
        call conmaxdriver(xguess)
      else if(meth.gt.0) then
        call newmin(numfn,xguess,sumt)
        call fnset(numfn,numgr,pttbl,iptb,indm,xguess,7,1,icntyp,confun)
      endif
 912  continue
      
      write(7,*) 'Bessel wavevectors and final coefficients'
      do i=1,numfn
        write(7,9012) i,xroot(i),xguess(i)
      enddo
 9012 format(1x,i3,2f20.10)
      sum2=0.0
      sum3=0.0
      do i = 1,numfn
         sum = 0.0
         do j = 1,numfn
            sum = sum + xguess(j) * d(i,j)
         enddo
         sum2 = sum2 + xguess(i)**2 * xroot(i)*abs(xroot(i)) * a(i)
         sum3 = sum3 + (sum + 2.0*e(i)) * xguess(i)
      enddo
      sum3=sum3+gint
      sum2=sum2-fint
      wsumt(nnn) = sumt*wnl(nnn)
      sumc(nnn) = sumt
      transumry=transumry+wsumt(nnn)

      write(7,*)
      write(7,9135) sum2
      write(7,9136) sum3
      write(7,9137) sum2-sum3
      write(7,9202) wnl(nnn)
      write(7,9102) sumt      
      write(7,9133) wsumt(nnn)*1000.0
      write(7,9134) wsumt(nnn)*13.6058*1000.0

c     *************************************************************************
c     Invert the wavefunction 
c     *************************************************************************

      call optinvert(fr,r,h,lo(nnn),en(nnn),indrc,np,maxim,xroot,xguess,
     $     numfn,xion,v0,ixc)

c     *************************************************************************
c     Step 5: Now put the potential into rvps and rnl and rsval.
c     *************************************************************************

      inode = 0
      maxim=min(maxim,np)
      if (maxim.lt.10) maxim=np

      do i = 1,np
        if (fr(i).lt.0.0.and.i.le.indrc) then
          inode = i
        endif
        fr2(i) = fr(i) * fr(i) * r(i) * r(i)
      enddo
      if (inode.gt.0) then
        write (7,*) 'Ang momemtum',ll,' has a node at r=',r(inode)
        stop
      endif

      tch = float(ll+ll+2)
      call radin(r,fr2,0,maxim,h,tch)
 
      tch = sqrt(tch)
      do i = 1,np
        fr(i) = fr(i)/tch
      enddo

      do i = 1,np
        fr2(i) = fr(i) * fr(i) * r(i) * r(i)
      enddo
      tch = float(ll+ll+2)
      call radin(r,fr2,0,maxim,h,tch)

      do i = 1,np
        rnl(i,nnn) = fr(i) * r(i)
        rsval(i) = rsval(i) + wnl(nnn) * rnl(i,nnn)**2

        if (ixc.ge.0) then
           if (i.gt.indrc+50) rvps(i,nnn) = -z-z+rvcoul(i)
           if (i.gt.maxim) rvps(i,nnn) = -xion-xion
        endif
      enddo
            
      tch = float(ll+ll+2)
      call radin(r,rsval,0,maxim,h,tch)
      
      write(7,*)
c      write(7,9200) tch
c      write(7,9202) wnl(nnn)

 9135 format(1x,'Total KE           (Ry) : ',3f16.10)
 9136 format(1x,'0  -> qc contribution   : ',3f16.10)
 9137 format(1x,'qc -> oo contribution   : ',f16.10,3x,
     $     '(=Convergence term)')
 9202 format(1x,'Occupancy of this state : ',f16.10)            
 9102 format(1x,'Convergence term (Ry/e) : ',f16.10)
 9133 format(1x,'Convergence error (mRy) : ',f16.10)
 9134 format(1x,'Convergence error (meV) : ',f16.10)

      return
      end
      
      
