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

      subroutine hfsolve(zeff,ixc,exccut,ipsp,ifc,iexit,irel,iprint)
     
c     *************************************************************************
c     determine the self-consistent solution
c     *************************************************************************
     
      implicit double precision (a-h,o-z)

#include "fortdim.h"
 
c     -------------------------------------------------------------------------
c     External (shared between C and Fortran) common blocks
c     -------------------------------------------------------------------------

      common /wfn/ rnl(npdm,n0)
      common /wfnrel/ rnla(npdm,n0),rnlb(npdm,n0)
      
      common /rscore/ rscore(npdm),rdd(npdm),rddd(npdm),rscoretot(npdm),
     $     rsval(npdm)
      common /totpot/ rvcore(npdm,n0),rvps(npdm,n0),rvcoul(npdm)

      common /consts/ etol,vtol,maxit,isoft
      common /grid/ h,r1,z,r(npdm),np

      common /filenames/ file_log

      common /box/ iboxstart(n0),iboxend(n0),boxheight(n0),numbox
      common /local/ nlghost(n0),iloc,idesign
      common /nlpot2/ inl,indrc(n0),IB,ID,IG

      common /ilogder/ ilogder
      common /logarith/ rphas,elogmax,elogmin,dlwf(npl0,n0)

      common /aorb/ ncore,nval,norb,nlm(n0),
     $     no(n0),lo(n0),nmax(n0),maxim
      common /adat/ wnl(n0),en(n0),so(n0),xion
      common /aval/ rcall(n0),rvap(n0),rnorm(n0),ibd(n0),etot
      common /edi/ e(n0,n0)
      common /ahf/ s(n0),az(n0),acc(n0),dpm(n0)

c     -------------------------------------------------------------------------
c     Internal (Fortran only) common blocks                
c     -------------------------------------------------------------------------
      common /nlpot1/ Flstar(npdm),phipsref(npdm,n0),phiguess(npdm,n0)
      common /partpot/ rvh(npdm),rvxc(npdm),rexc(npdm)
      common /nlcore/ rvloc(npdm)
      common /ipos/ ipos(n0),itermcount
      common /vdir/ vdir(npdm),pold(npdm)
      common /rold/ rnlold(npdm,n0),rpold(npdm,n0)
      common /eihom/ eihom(n0),rnlhom(npdm,n0)
      common /ind/ ind(n0)
      common /ameth/ meth(n0),imeth,nfails

      common /ighost/ ighost,ighost2
c     -------------------------------------------------------------------------
c     *************************************************************************
c     local variables
c     *************************************************************************
      dimension eold(n0)
      dimension rsatom(npdm),rsold(npdm)
      
      dimension wavea(npdm),icarray(30)
      dimension x(npdm),yr(npdm),fr(npdm),fr2(npdm)
      dimension lcount(10)
c     *************************************************************************
c     *************************************************************************
      character*80 file_log
      logical ilast,first

      vtl=vtol
      etl=etol*10.d0
      ilast=.false.
      iexit=0
      itermcount=0
      isoft=ipsp
      nig=0
      ipratt=6
      npratt = ipratt
      rvapmax=0.0
      ib=1
      ipr=0
      icorb = norb/4+2
      ighost=0
      ighost2=0
c      iprint=5

      open(unit=7,file=file_log,form='formatted',access='append')

      iflag=0
      nfails=0

      call hfcoef(iflag,ipsp)

c     EJW - adding this to keep allow to pick a hyrdogenic wavefunction 
c     starting guess for channels which don't have a psp-wfn guess 
      do i=1,10
         lcount(i)=0
      enddo

      do i=1,norb        
         en(i)=-en(i)
         meth(i)=1
         do j=1,norb
            e(i,j)=0.d0
            if (lo(i).eq.lo(j)) e(i,j)=1e-5
         enddo

         e(i,i)=en(i)
         
         nmax(i)=np
         dpm(i)=0.d0
         acc(i)=0.d0
         do j=1,np
            rvcore(j,i)=-rvcore(j,i)/2.d0
         enddo
      enddo

      if (ipsp.eq.0) then
         do j=1,norb
            do jj=1,norb
               e(j,jj)=0.d0
               if (lo(j).eq.lo(jj)) e(j,jj)=1e-5
            enddo
            e(j,j)=-en(j)
            do i = 1,np
               rnla(i,j) = 0
               rnlb(i,j) = 0
               rnl(i,j) = 0
            enddo
         enddo
      else
         do j=1,norb
            do i = 1,np
             rnl(i,j)=rnl(i,j)/sqrt(r(i))
            enddo
         enddo
      endif

      etotlast = 0
      ifock=0
      
      ss=0.d0
      do i=1,ncore
         imf=4*lo(i)+2
         s(i)=ss+imf/2.d0
         ss=ss+imf
      enddo
      do i=ncore+1,norb
         s(i)=ss+(wnl(i)-1.d0)/2
         ss=ss+wnl(i)
      enddo

      if (ipsp.eq.0) then
         do i=1,norb
            pnf=hnorm(no(i),lo(i),z-s(i))
            do k=1,np
               rnl(k,i)=pnf*hwf(no(i),lo(i),z-s(i),r(k))/sqrt(r(k))
            enddo
            az(i)=pnf*(2.d0*(z-s(i)/2)/no(i))**(lo(i)+1.d0)
         enddo
      else
         do i=1,norb
            lcount(lo(i)+1)=lcount(lo(i)+1)+1
            pnf=hnorm(no(i),lo(i),zeff-s(i))

            if (lcount(lo(i)+1).gt.1) then
               do k=1,np
                 rnl(k,i)=pnf*hwf(no(i),lo(i),zeff-s(i),r(k))/sqrt(r(k))
               enddo
            endif
            az(i)=1.d0
         enddo
      endif

      do i=1,norb
         if (no(i).eq.0.and.lo(i).eq.0) goto 33
         do j=1,i-1
            if (e(i,j).gt.0.d0) then
               rx=quadr(i,j,0)
               if (abs(rx).gt.1e-8) then
                  rx2=sqrt(1.d0-rx**2)
                  if (rnl(50,i)-rx*rnl(50,j).lt.0.d0) rx2=-rx2
                  do k=1,np
                     rnl(k,i)=(rnl(k,i)-rx*rnl(k,j))/rx2
                  enddo
               endif
            endif
         enddo
 33      continue
      enddo

      do niter = 1,maxit

C          if (ipsp.ne.0) then
C             do i=1,norb
C                do j=1,np
C                   if (r(j).lt.15) write(850+i,*) r(j),rvps(j,i)
C                   if (r(j).lt.15) write(860+i,*) r(j),
C      $                 x(j)/rnl(j,i),rnl(j,i)
C                   if (r(j).lt.15) write(870+i,*) r(j),yr(j)
C                   if (r(j).lt.15) write(370+i,*) r(j),rnl(j,i)
C                enddo
C             enddo
C             stop
C          endif

         call grange

c     ***********************************************************************
c     ***********************************************************************
         dp1= 0
         edmax = 0
         difmax = 0
         maxim = 0
         indey = 0
         do i = 1,np
            rsold(i)=rsatom(i)
            rsatom(i) = 0
         enddo

         dp1=0.d0
c         nfails=0
         do m=1,norb
            ed1=0
            eold(m)=en(m)
            ed2=e(m,m)

            call hfde(m,ipsp,zeff,dpm,iprint)
c            write(7,*) 'NFAILS', nfails
            if (nfails.gt.25) goto 911

            dp=dpm(m)
            if (dp1.ge.dp) goto 2
            dp1=dp
            ibad=m
 2          continue
            edif = abs((en(m) - eold(m))/eold(m))
            eold(m) = en(m)
            edmax = max(edif,edmax)
         enddo

c     ***********************************************************************
c      m loop
c     ***********************************************************************
         
         if (dp1.lt.vtl) goto 66

         dp1=0.d0
         do mm=1,6
            m=ibad
            ed1=0
            ed2=e(m,m)
            eold(m)=en(m)
            call hfde(m,ipsp,zeff,dpm,iprint)

            if (nfails.gt.25) goto 911
            dp1=0.d0
            do j=1,norb


               dp=dsqrt(wnl(j))*dpm(j)
               if (dp1.gt.dp) goto 5
               ibad=j
               dp1=dp
 5             continue

               edif = abs((en(m) - eold(m))/eold(m))
               eold(m) = en(m)
               edmax = max(edif,edmax)
            enddo

            if (dp1.lt.vtl) goto 66
         enddo
 66      continue

         call orthog(ilast)

         ebs=0.d0
         do m = 1,norb
            ebs = ebs + wnl(m) * en(m)
         enddo
         etot = 0.d0

         if (iprint.gt.1) write(7,8) nlm(ibad),dp1,vtl
         
 8       format(6x,'least self-consistent function is ',i3,
     $        ' weighted maximum change =',1pd10.2,1pd10.2)

         call hfenergy(ipsp,etot,ehart,exx,ekn)      

         evir=(etot-ekn)/ekn
         ekn=ekn*2.d0
         etot=etot*2.d0
         ehart=ehart*2.d0
         exx=exx*2.d0
         
         etot2=etot
         etot=-ebs-ehart-exx

c         etot2=-ebs-ehart-exx
         edf=etot-etot2

         ediff = etot - etotlast
         dvmax=dp1

         if (niter.eq.maxit) goto 911
         
c     911 = terminal error - exit scpot
c     ***********************************************************************
c     ***********************************************************************
         if (iprint.gt.0) then
            if (niter.eq.1) then
               write(7,*)
               write(7,191)
               write(7,491) niter,etot,ebs,ehart,edmax,dvmax
            else
               write(7,491) niter,etot,ebs,ehart,edmax,dvmax
            endif
            
            etotlast = etot
         
 191        format(' iter',7x,'Etot',13x,'Eeig',13x,'Ehart',4x
     $           ,'de_max',3x,'dv_max')
 490        format(1x,i3,2x,f15.7,1x,f15.7,1x,f15.7)
 491        format(1x,i3,2x,f15.7,1x,f15.7,1x,f15.7,1x,e8.2,1x,e8.2)
         endif
c         write(7,*) 'CONVV', dvmax,vtl
c         write(7,*) 'CONVE', edmax,etl
         if (dvmax.lt.vtl.and.edmax.lt.etl) goto 1000
         
c     ***********************************************************************
c     ***********************************************************************
         vtl=vtl*1.1d0
         etl=etl*1.1d0
      enddo

 1000 continue

      if (iprint.gt.0)  then
         write(7,*) 
         write(7,*) 'Hartree-Fock converged!'
      endif

      if (iprint.gt.1)  write(7,192)
      if (iprint.gt.1)  write(7,492) etot2,ekn,exx,evir,edf

 492  format(7x,f15.7,1x,f15.7,1x,f15.7,1x,f15.7,1x,f15.7)
 192  format(12x,'(Etot2)',13x,'Ekin',13x,'Exx',13x,'Evir',13x,'Edf')


c     EJW new ghost section
c     
      
c      do nloc=1,3
c            
c         do i=1,3
c            if (nloc.ne.i) then
c               ll=lo(i)
c               do j = 1,np
c                  fr(j) = -2.d0*rnl(j,i)**2*
c     $                 (rvcore(j,i)-rvcore(j,nloc))/r(j)
c
c               enddo
c               
c               tov = (float(ll+ll+2))
c               call radin(r,fr,0,np,h,tov)
c               xden = tov
c               
c               do j = 1,np
c                  fr(j) = -2.d0
c     $                 *fr(j) * (rvcore(j,i)-rvcore(j,nloc))/r(j)
c               enddo
c               
c               tov = (float(ll+ll+2))
c               call radin(r,fr,0,np,h,tov)
c               xnum = tov
c               
c               elkb = xnum/xden
c               
c               write(7,9111) nloc,i,elkb,sqrt(xnum),sqrt(xnum)/elkb
c 9111          format(2i5,1x,'KB energy : ',f10.6,2x
c     $              ,'KB strength: ',f10.6,2x,'KB cosine: ',f10.6,2x)
c            endif
c         enddo
c      enddo

      m=1
      n=2
      if (ipsp.eq.33) then
         write(7,*) 'E=',-en(m)
         iprint=5

         dpm(m)=0.d0
         ipsp=1
         first=.false.
         do k=1,np
            rsold(k)=rvcore(k,m)
            rvcore(k,m)=rvcore(k,n)
            rnl(k,m)=rnl(k,n)
         enddo
         ighost=1
         ighost2=1
         
         do k=1,1
            call hfmethd(m,ipsp,rsold,zeff,azz,fail,iprint)
            write(7,*) 'E=',-en(m)
            no(m)=no(m)+1

            call hfmethd(m,ipsp,rsold,zeff,azz,fail,iprint)
            write(7,*) 'E=',-en(m)
            no(m)=no(m)-1
         enddo
         

         stop
      endif


c      if (ipsp.ne.0) then
c         do i=1,3
c            call potl(yr,i)
c            call xch(x,i,1)
c            do j=1,np
c               write(600+i,8585) r(j),-2.d0*rvcore(j,i),2.d0*yr(j)+x(j)
c            enddo
c         enddo
c      endif
c 8585 format(4f20.10)
            
      ilast=.true.
      call orthog(ilast)

      rmmax=15.d0
      if (ipsp.eq.0) then
         do i=1,norb

            do j=1,np
               if ((abs(rnl(j,i)*sqrt(r(j))).lt.1e-8).
     $              and.(r(j).gt.15.0)) then
                  rmmax=r(j-1)
                  goto 788
               endif
            enddo
 788        continue
            if (iprint.gt.1) write(7,*) 'Max radius used for Veff:',
     $           i,rmmax
            call potl(yr,i)
            call xch(x,i,1)
            eps=1d-18
            do j=1,np
               if (r(j).gt.rmmax) goto 922
               xxx=x(j)/sqrt(r(j))
               rji=2.d0*yr(j)/r(j)
               rvps(j,i)=-2.d0*rvcore(j,i)+2.d0*yr(j)+x(j)
     $              /(rnl(j,i)+eps)
c     if (r(j).lt.15) write(850+i,*) r(j),rvps(j,i)
c     if (r(j).lt.15) write(860+i,*) r(j),x(j)/rnl(j,i),rnl(j,i)
c     if (r(j).lt.15) write(870+i,*) r(j),yr(j)
c     if (r(j).lt.15) write(370+i,*) r(j),rnl(j,i)
            enddo
 922        continue
         enddo
      endif

      do i=1,norb
         do j=1,np
            rnl(j,i)=rnl(j,i)*sqrt(r(j))
            rvcore(j,i)=-2.d0*rvcore(j,i)
         enddo
         en(i)=-en(i)
      enddo

      do i=1,np
         rsatom(i)=0.0
      enddo
      do nnn=1,ncore
         do i=1,np
            rsatom(i)=rsatom(i)+wnl(nnn)*rnl(i,nnn)**2
         enddo
      enddo
      call hrtree(np,h,r,rsatom,rvh)
      if (ipsp.eq.0) then
         open(unit=33,file="vcore",form="formatted")
         do i=1,np
            write(33,*) r(i),rvh(i)
         enddo
         close(33)
      endif


c     ***********************************************************************
c     ***********************************************************************

c     Flip wfn s.t. it is >0 as r->oo

      do i = ncore+1,norb
         do j=np,1,-1
            if (abs(rnl(j,i)).gt.1e-6) then
               nmax(i)=j
               goto 939
            endif
         enddo
 939     continue
         if (rnl(nmax(i)-1,i).lt.0.0) then
            do j = 1,np
               rnl(j,i) = -rnl(j,i)
            enddo
         endif
      enddo

c     compute valence density 
      do j = 1,np
         rsval(j) = 0
      enddo

      do i = ncore+1,norb
         do j = 1,nmax(i)
            rsval(j) = rsval(j) + wnl(i) * rnl(j,i)**2
         enddo
      enddo

c     compute core density

      if (ipsp.eq.0) then
         rsc=0.0d0
         do j = 1,np
            rscore(j) = 0
         enddo
         
         do i = 1,ncore
            do j = 1,nmax(i)
               rscore(j) = rscore(j) + wnl(i) * rnl(j,i)**2
            enddo
         enddo
         call radin(r,rscore,0,np,h,rsc)
      endif

c     loop over valence
      
      do i = ncore + 1,norb
         ic=i-ncore
c     compute icarray!! put OUT
         
         icarray(ic) = (log(rcall(ic)/r(1))/
     $        log(exp(1.0)))/h + 1
c     compute tail-norm (actually, here it is 1.d0-head norm)
         do k = 1,np
            wavea(k) = rnl(k,i) * rnl(k,i)
         enddo
         
         pow = 2 * lo(i) + 2
         call radin(r,wavea,0,icarray(ic),h,pow)
         rnorm(i) = 1.d0 - pow
      enddo
      
c     radius of outermost peak
      do i=ncore+1,norb
         rvapmax=0.0
         do j=maxim-1,1,-1
            if ((rnl(j,i).lt.rnl(j+1,i)).
     $           and.(abs(rnl(j,i)).gt.1e-6)) then
               rvap(i)=r(j+1)
               if (rvapmax.lt.rvap(i)) rvapmax=rvap(i)
               goto 944
            endif
         enddo
 944     continue
      enddo

c     ***********************************************************************
c     ***********************************************************************

      if (iprint.gt.0)  write(7,*) 

      if (niter.eq.1) then
         if (iprint.gt.0) write(7,700)
      else
         if (iprint.gt.0) write(7,500) niter
      endif

      if (iprint.gt.0) write(7,701) etot,-ebs
      if (iprint.gt.0) write(7,702) ehart,exx

c     write(7,701) etot,ebs,ehart+exx

      if (iprint.gt.0)  write (7,*)
      if (iprint.gt.0)  write (7,9222)

      do i = 1,ncore
         if (iprint.gt.0)  write (7,9223) nlm(i),wnl(i),en(i)
      enddo
      do i = ncore+1,norb
         if (iprint.gt.0)  write (7,9223) nlm(i),wnl(i),
     $        en(i),rnorm(i),rvap(i)
      enddo

 500  format(1x,i4,1x,'iterations: (units are Ry and bohr)')
 700  format(1x,'Converged in 1 iteration (probably reference state)')
 701  format(1x,'Energy: ',f16.8,4x,'Ebs: ',f16.8)
 702  format(1x,'Eh:     ',f16.8,4x,'Ex : ',f16.8)

c 701  format(1x,'Energy: ',f16.8,2x,'Eeig: ',f16.8,
c     $     2x,'Ehart+Exx: ',f16.8)

 9222 format(3x,'Orbital',4x,'Filling',7x,'Eigenvalues',4x,
     $     'Norm(rc->oo)',6x,'Peak')
 9223 format(4x,'|',i3,'>',5x,f6.3,5x,3f14.6)
 9224 format(4x,'|',i3,'>',5x,f6.3,5x,f14.6,6x,
     $     '---------',5x,'--------')

 9009 format(1x,'partial core radius : ', f10.4)
 9010 format(1x,'core charge         : ', f20.10)
 9011 format(1x,'partial core charge : ', f20.10)


      call flush(7)
      close(unit=7)

      return
      
 911  continue
c     911 = terminal error - exit scpot
      write(7,*) " !!ERROR!! Could not find self-consistent solution!"
      write(6,*) " !!ERROR!! Could not find self-consistent solution!"
      iexit=1
      
      call flush(7)
      close(unit=7)


 370  format(1x,i2,2i5,2f6.1,f10.4)
 390  format(//,' radial grid parameters',//,
     +     ' r(1) = .0 , r(2) =',e8.2,' , ... , r(',i5,') =',f6.2,
     +     /,' aa =',f5.2,'  bb =',f6.2,/)
      
 360  format(' input data for orbitals',//,
     +     '  i    n    l    s     j     occ',/)
      
 340  format(' nuclear charge             =',f10.6,/,
     +     ' number of core orbitals    =',i3,/,
     +     ' number of valence orbitals =',i3,/,
     +     ' electronic charge          =',f10.6,/,
     +     ' ionic charge               =',f10.6,//)
      
 332  format(' XC functional is GGA (Perdew-Burke-Ernzerhof)')
 334  format(' XC functional is LDA (Perdew-Zunger)')
 335  format(' XC functional is LDA (Perdew-Wang)')

      return
      
      end
      

