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
c change some of the 911 case to be less than 50 instead of less than 10
      subroutine dftsolve(zeff,ixc,exccut,ipsp,ifc,iexit,irel,
     $     irelxc,iprint)
     
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

      common /rpcc/ rpcc,rpccz

      common /aorb/ ncore,nval,norb,nlm(n0),
     $     no(n0),lo(n0),nmax(n0),maxim
      common /adat/ wnl(n0),en(n0),so(n0),xion
      common /aval/ rcall(n0),rvap(n0),rnorm(n0),ibd(n0),etot

      common /niter/ niter
      common /ahf/ s(n0),az(n0),acc(n0),dpm(n0)
c     -------------------------------------------------------------------------
c     Internal (Fortran only) common blocks                
c     -------------------------------------------------------------------------
      common /nlpot1/ Flstar(npdm),phipsref(npdm,n0),phiguess(npdm,n0)
      common /partpot/ rvh(npdm),rvxc(npdm),rexc(npdm),
     $     rvx(npdm),rex(npdm),rvc(npdm),rec(npdm),rvdft(npdm)
      common /nlcore/ rvloc(npdm)
      common /ipos/ ipos(n0),itermcount
      common /vdir/ vdir(npdm),pold(npdm)
      common /rold/ rnlold(npdm,n0),rpold(npdm,n0)
      common /ameth/ meth(n0),imeth,nfails
c     -------------------------------------------------------------------------

c     *************************************************************************
c     local variables
c     *************************************************************************
      dimension rvn(npdm),rvf(npdm),eold(n0),p(npdm)
      dimension rsatom(npdm),f(npdm),f2(npdm),f3(npdm),f4(npdm)
      dimension pder(npdm),ader(npdm)
      dimension dl(n0), rsatom2(npdm)
      dimension rsold(npdm),wnlfract(npdm)
      dimension wavea(npdm),icarray(n0)
      dimension emattest(n0,n0)
      dimension lcount(10)
c     rvpostloc(npdm),rvpreloc(npdm)
c     *************************************************************************
c     Section 1:  set up tolerance params and zero out arrays.
c     *************************************************************************
      character*80 file_log
      character*1 xc(0:3)
      
      xc(0)='s'
      xc(1)='p'
      xc(2)='d'
      xc(3)='f'

      iexit=0
      itermcount=0
      isoft=ipsp
      igh=0
      ipratt=6
      npratt = ipratt
      rvapmax=0.0

      open(unit=7,file=file_log,form='formatted',access='append')

c     User supplied eigenvalue guesses were never used in the SRL
c     routines ; disable for now.

      if (irel.eq.1) then
         do i=1,ncore
            en(i)=-z*z/(no(i)*no(i))
         enddo
         do i=1,nval
            ii=i+ncore
            if (ibd(i).ne.0) en(ii)=-z*z/(no(ii)*no(ii))
         enddo
      endif



c      if (ipsp.eq.1) then
c
c         do i=1,norb
c            write(7,*) i,en(i),no(i),lo(i),so(i),ibd(i)
c         enddo
c      endif

      write(7,*) '---------Eigenvalue guesses---------'      

      do i=1,norb
         ipos(i)=0
         write(7,9339) no(i),xc(lo(i)),en(i)
      enddo
 9339 format(1x,i1,a1,3x,f20.8)
     
      if (ipsp.eq.1) then
         nlc=nlm(iloc)/100
         llc=(nlm(iloc) - nlc * 100)/10
      endif

      if (inl.ne.0) then
         do i=1,np
            do j=1,norb
               phipsref(i,j)=rnl(i,j)
               phiguess(i,j)=rnl(i,j)
            enddo
         enddo
      endif
            
      if (ifc.le.12) then
         do j = 1,norb
            eold(j) = en(j)
            nmax(j) = 0
            if (irel.eq.0) then
               do i = 1,np
                  rnl(i,j) = 0
                  rsatom(i)=0
               enddo
            else
               do i = 1,np
                  rnla(i,j) = 0
                  rnlb(i,j) = 0
               enddo
            endif               
         enddo
      endif

      if (inl.ne.0) then
        do i=1,norb
          do j=1,np
            if (r(j).gt.rcall(i)) goto 993
          enddo
 993      continue
          indrc(i)=j
        enddo
      endif

      do i = 1,np         
        rvcoul(i) = rvps(i,1) - rvcore(i,1)
      enddo


      if (inl.ne.0) then
        write(7,*)"Performing non-local pseudopotential calculation"
      else
         if (ipsp.ne.0) then
       write(7,*) "Performing semi-local pseudopotential calculation"
         endif
      endif

c     *************************************************************************
c     If DNL calc, write out box info
c     *************************************************************************

      if (ipsp.eq.1.and.inl.ne.0) then
         write(7,9600) xc(llc)
      endif
 9600 format(1x,'Using the ',a1,' potential as the local potential')

      if (ipsp.eq.1.and.inl.ne.0.and.numbox.gt.0) then
         write(7,*)  
         write(7,9220)  
         write(7,9230) numbox
         rindmin = 1000000
         rindmax = 0
         do i=1,norb
            if (rcall(i).lt.rindmin) rindmin=rcall(i)
            if (rcall(i).gt.rindmax) rindmax=rcall(i)
         enddo
         ioutmin=0
         ioutmax=0
         do i=1,numbox
            write(7,9240) i,r(iboxstart(i)),r(iboxend(i)),boxheight(i) 
            if (r(iboxend(i)).gt.rindmin) ioutmin = 1
            if (r(iboxend(i)).gt.rindmax) ioutmax = 1
         enddo
         if (ioutmin.eq.1.and.ioutmax.eq.0) write(7,9250) 
         if (ioutmax.eq.1) write(7,9260) 
         write(7,9225)
      endif
      
 9220 format(1x,'--Augmentation operator info--')
 9225 format(1x,'------------------------------')
 9230 format(1x,'Number of functions: ',i3)
 9240 format(1x,'#',i2,1x,'range(a.u.): ',f10.6,' ----->',f10.6,3x,
     $     'size(Ry):',f10.6)
 9250 format(1x,'NOTE: One or more functions extend beyond the',
     $     ' minimum cutoff radius')
 9260 format(1x,'WARNING!: One or more functions extend beyond the',
     $     ' MAXIMUM cutoff radius')

c     *************************************************************************
c     Begin self-consistency loop.
c     *************************************************************************
      maxit=1000
      etotlast = 0
      ifock=0
      do niter = 1,maxit

c     ***********************************************************************
c     Section 2:  Initialization 
c     ***********************************************************************

c     unbound info
         if (niter.eq.1) then
            do i=1,norb
               ipos(i)=0
               nmax(i)=np
            enddo
         endif
         
         npratt = npratt + 1
         if (npratt.gt.ipratt) npratt = 0

c         do i = 1,np
c            write(80,*) rsatom(i)
c         enddo


         edmax = 0
         difmax = 0
         maxim = 0
         indey = 0
         do i = 1,np
            rsold(i)=rsatom(i)
            rsatom(i) = 0
         enddo
         
c     ***********************************************************************
c     Section 3:  call schsl for each orbital & compile charge density.
c     ***********************************************************************

         do m = 1,norb

            if (m.le.ncore) then
               ibound=1
            else
               ibound=ibd(m-ncore)               
            endif

            if (ifc.le.12.or.m.gt.ncore) then

c     if not a frozen core calc or if valence
               if (irel.eq.0) then
c     if non-rel
                  
                  if (inl.ne.0) then
c     if non-local
                     call applyaug(rvloc,rvps(1,iloc))
                     do i=1,np
                        Flstar(i) = (rvps(i,m) - rvloc(i))
     $                       *phipsref(i,m)/r(i)
                     enddo

                     call schsl (m,no(m),lo(m),en(m),nmax(m),
     $                    rvloc,p,igh,ibound)

                     if (itermcount.gt.50) then
                         goto 911
                     endif 
c     911 = terminal error - exit scpot
                     
                     call applyaug(rvloc,rvcore(1,iloc))
                     
                     do i = 1, np
                        phiguess(i,m) = p(i)
                     enddo
                     
                  else
c     if AE or local 
                     call schsl(m,no(m),lo(m),en(m),nmax(m),
     $                    rvps(1,m),p,igh,ibound)

                     if (itermcount.gt.50) then
                         write(7,*) "line 445"
                         goto 911
                     endif
c     911 = terminal error - exit scpot
                     
                  endif
                  
                  do i = 1,nmax(m)
                     rnl(i,m) = p(i)
                     rsatom(i) = rsatom(i) + wnl(m) * rnl(i,m)**2
                  enddo

                  do i = nmax(m) + 1,np
                     rnl(i,m) = 0
                  enddo                  
                  
                  indey = max(indey,nmax(m))
                  maxim = indey
                  
                  edif = abs((en(m) - eold(m))/eold(m))
                  eold(m) = en(m)
                  edmax = max(edif,edmax)
                  
               else
c     if relativistic
                  
                  call schslrel(m,rvps(1,m),rnla(1,m),
     $                 rnlb(1,m),z,rvcoul,en,nmax(m),ibound)
                  
                  if (itermcount.gt.50) then
                      goto 911
                  endif
c     911 = terminal error - exit scpot
                  
                  do i=1,np
                     rsatom(i)=rsatom(i)
     $                    +wnl(m)*(rnla(i,m)**2+rnlb(i,m)**2)
                  enddo
                  
                  do i = nmax(m) + 1,np
                     rnla(i,m) = 0
                     rnlb(i,m) = 0
                  enddo                  
                  
                  indey = max(indey,nmax(m))
                  maxim = indey
                  
c     hack for now
c     maxim=np
                  
                  edif = abs((en(m) - eold(m))/eold(m))
                  eold(m) = en(m)
                  edmax = max(edif,edmax)
                  
c     if relativistic
               endif
               
c     if not a frozen core calc or if valence 
            endif
         enddo

         if(maxim.eq.0) then
            goto 911
c     911 = terminal error - exit scpot
         endif

c     ***********************************************************************
c     Section 4:  mix density
c     ***********************************************************************

         alpha = 0.5
         if (niter.gt.3) then
            do i = 1,np
               rsatom(i) = alpha * rsatom(i) + (1.0-alpha) * rsold(i)
            enddo
         endif
            
c     ***********************************************************************
c     Section 5:  update self-consistent potentials and compute energy.
c     ***********************************************************************

         call hrtree(maxim,h,r,rsatom,rvh)
         
         if (rpcc.gt.1e-12.and.ipsp.ne.0) then
            do i=1,maxim
               rsatom2(i)=rsatom(i)+rscore(i)
            enddo
            call excorr(maxim,ixc,exccut,rsatom2,
     $           rvxc,rexc,rvx,rex,rvc,rec,irelxc)
            
            do i = 1,maxim
               f(i) = 0.5d0 * rvh(i)*rsatom(i)/r(i)
     $              + rexc(i)*rsatom2(i)/r(i)
     $              - (rvxc(i)*rsatom(i)+rvh(i)*rsatom(i))/r(i)
            enddo
         else
            call excorr(maxim,ixc,exccut,rsatom,
     $           rvxc,rexc,rvx,rex,rvc,rec,irelxc)
            
            do i = 1,maxim
               f(i) = (0.5d0 * rvh(i) + rexc(i) - (rvxc(i)+rvh(i)) ) 
     $              * rsatom(i)/r(i)
               f2(i) = 0.5d0 * rvh(i)* rsatom(i)/r(i)
               f3(i) = rexc(i)* rsatom(i)/r(i)
               f4(i) = rvxc(i)* rsatom(i)/r(i)
            enddo
         endif
         ehxc = 2
         call radin(r,f,0,maxim,h,ehxc)
         ehar = 2
         call radin(r,f2,0,maxim,h,ehar)
         exc = 2
         call radin(r,f3,0,maxim,h,exc)
         evxc = 2
         call radin(r,f4,0,maxim,h,evxc)
         ebs = 0
         do m = 1,norb
            ebs = ebs + wnl(m) * en(m)
         enddo
         etot = ebs + ehxc
         etot = ebs - ehar + exc - evxc

c     ***********************************************************************
c     Section 6:  Mix potentials
c     ***********************************************************************
               
         dvmax = 0
         xn2 = (zeff - xion) * 2
         do i = 1,maxim
            coul = rvh(i) + rvxc(i)  
            dv = abs(coul - rvcoul(i))
            if (dv.gt.dvmax) igdd=i
            dvmax = max(dv,dvmax)
            pr = 0.5
            if (npratt.gt.0) pr = pratt(rvn(i),rvf(i),rvcoul(i),coul)
            rvn(i) = rvcoul(i)
            rvf(i) = coul
            rvcoul(i)= pr * rvcoul(i) + (1.d0 - pr) * coul
         enddo
 9103    format(4f20.10)
         
         do i = maxim + 1,np
            rvcoul(i) = xn2
            rvn(i) = xn2
            rvf(i) = xn2
         enddo
         xmix = 0.9
         do i = 1,np
            do iorb = 1,norb
               rvps(i,iorb) = (rvcore(i,iorb) + rvcoul(i)) * xmix
     $              + (1.0 - xmix) * rvps(i,iorb)
            enddo
         enddo

         if (niter.eq.maxit) then
             goto 911
          endif
            
c     911 = terminal error - exit scpot
            
c     ***********************************************************************
c     Section 6: Check progress
c     ***********************************************************************
         ediff = etot - etotlast
         
         if (niter.eq.1) then
            write(7,*)
            write(7,191)
            write(7,491) niter,etot,ebs,ehxc,edmax,dvmax
         else
            write(7,491) niter,etot,ebs,ehxc,edmax,dvmax
         endif
         
         etotlast = etot

 191     format(' iter',7x,'Etot',13x,'Ebs',13x,'Ehxc',7x
     $        ,'de_max',3x,'dv_max')
 490     format(1x,i3,2x,f15.7,1x,f15.7,1x,f15.7)
 491     format(1x,i3,2x,f15.7,1x,f15.7,1x,f15.7,1x,e8.2,1x,e8.2)

         if (dvmax.lt.vtol.and.edmax.lt.etol) goto 1000
         
c     ***********************************************************************
c     Section 7:  predict trial eigenvals use 1st order perturb theory.
c     ***********************************************************************
         do i = 1,indey
            p(i) = (rvcoul(i) - rvn(i))/r(i)
         enddo
         do m = 1,norb
            if (m.le.ncore) then
               ibound=1
            else
               ibound=ibd(m-ncore)
            endif
            if (ifc.le.12.and.ibound.eq.1) then
               if (irel.eq.0) then
                  do i = 1,nmax(m)
                     f(i) = p(i) * rnl(i,m)**2
                  enddo
               else
                  do i = 1,nmax(m)
                     f(i) = p(i) * rnla(i,m)**2
                  enddo
               endif
               xl = 2 * lo(m) + 2
               
               if(nmax(m).le.10) then
                  goto 911
c     911 = terminal error - exit scpot
               endif

               call radin(r,f,0,nmax(m),h,xl)

               en(m) = en(m) + xl

               if (en(m).ge.0.0) then
                  a1 = lo(m) * (lo(m) + 1)
                  vmin = 1.0e6
                  do i = 1,np
                     vmin = min((rvps(i,m) + a1/r(i))/r(i),vmin)
                  enddo
                  en(m) = vmin * 0.5d0
                  write (7,505) nlm(m),en(m)               
                  
                  if (en(m).ne.en(m)) then
                     goto 911 
c     911 = terminal error - exit scpot
                  endif
               endif
            endif
         enddo

      enddo


 505  format(1x,'!WARNING! positive eigenvalue predicted',
     $     ' in scpot  nlm=',i4,5x,'new e =',e13.7)
     
 1000 continue

c     ***********************************************************************
c     Section 8:  logder
c     ***********************************************************************

      if (ilogder.eq.1) then
         if (ipsp.eq.0) then
            ic=0
            
            do m=ncore+1,norb
               ic=ic+1
               l = lo(m)
               dele = (elogmax-elogmin)/float(npl0-1)
               nr = 4
               isoft = 0
               do i = 1, np
                  pder(i) = 0.0
                  ader(i) = 0.0
               enddo
               do i = 1, npl0
                  e = elogmin + dele * float(i-1)
                  call logder(e,z,l,isoft,rvps(1,m),ader,r,np,h,
     $                 rphas,dl,npl0,nr,pder,wfld,0,1.0,1.0,1.0)
                  dlwf(i,ic) = wfld
               enddo
            enddo
         endif

         if (ipsp.eq.1) then
            if (inl.eq.0) then
               ic=0
               do m=1,nval
                  ic=ic+1
                  l = lo(m)
                  dele = (elogmax-elogmin)/float(npl0-1)
                  nr = 4
                  isoft = 1
                  do i = 1, np
                     pder(i) = 0.0
                     ader(i) = 0.0
                  enddo
                  do i = 1, npl0
                     e = elogmin + dele * float(i-1)
                     call logder(e,z,l,isoft,rvps(1,m),ader,r,np,h,
     $                    rphas,dl,npl0,nr,pder,wfld,0,1.0,1.0,1.0)
                     dlwf(i,m) = wfld
                  enddo
               enddo
            else
               do m=1,nval
                  call nllogd(m,lo(m))
               enddo
            endif
         endif
      endif

c     ***********************************************************************
c     Section 9:  Post-process
c     ***********************************************************************

c     Flip wfn s.t. it is >0 as r->oo

      do i = ncore+1,norb
         if (irel.eq.0) then
            do k=nmax(i),1,-1
               if (abs(rnl(k,i)).gt.1e-10) then
                  if (rnl(k,i).lt.0.0) then
                     do j = 1,np
                        rnl(j,i) = -rnl(j,i)
                     enddo
                  endif
                  goto 770
               endif
            enddo
 770        continue
         else
c            if (ibd(i-ncore).eq.0) then
c               icarray(i-ncore) = (log(rcall(i-ncore)/r(1))/
c     $              log(exp(1.0)))/h + 1
c               do k = icarray(i-ncore),np
c                  rnla(k,i)=0.0
c               enddo
c            endif
            do k=nmax(i),1,-1
               if (abs(rnla(k,i)).gt.1e-10) then
                  if (rnla(k,i).lt.0.0) then
                     do j = 1,np
                        rnla(j,i) = -rnla(j,i)
                     enddo
                  endif
                  goto 771
               endif
            enddo
 771        continue
            do k=nmax(i),1,-1
               if (abs(rnlb(k,i)).gt.1e-10) then
                  if (rnlb(k,i).lt.0.0) then
                     do j = 1,np
                        rnlb(j,i) = -rnlb(j,i)
                     enddo
                  endif
                  goto 772
               endif
            enddo
 772        continue
         endif
      enddo

c     compute valence density 
      do j = 1,np
         rsval(j) = 0
      enddo

      do i = ncore+1,norb
         if (irel.eq.0) then
            do j = 1,nmax(i)
               rsval(j) = rsval(j) + wnl(i) * rnl(j,i)**2
            enddo
         else
            do j = 1,nmax(i)
               rsval(j) = rsval(j) + wnl(i) * (rnla(j,i)**2
     $              +rnlb(j,i)**2)
               
            enddo
         endif
      enddo

c     compute core density

      if (ipsp.eq.0) then
         rsc=0.0d0
         do j = 1,np
            rscore(j) = 0
         enddo
         
         do i = 1,ncore
            if (irel.eq.0) then
               do j = 1,nmax(i)
                  rscore(j) = rscore(j) + wnl(i) * rnl(j,i)**2
               enddo
            else
               do j = 1,nmax(i)
                  rscore(j) = rscore(j) + wnl(i) *(rnla(j,i)**2
     $                 +rnlb(j,i)**2)
               enddo
            endif
         enddo
         call radin(r,rscore,0,np,h,rsc)
      endif

c     loop over valence
      if (irel.eq.0) then
         do i = ncore + 1,norb
            ic=i-ncore
c     compute icarray!! put OUT

            icarray(ic) = (log(rcall(ic)/r(1))/
     $           log(exp(1.0)))/h + 1
c     compute tail-norm (actually, here it is 1.d0-head norm)
            do k = 1,np
               wavea(k) = rnl(k,i) * rnl(k,i)
            enddo
            
            pow = 2 * lo(i) + 2

            call radin(r,wavea,0,icarray(ic),h,pow)
            
            if (ibd(ic).eq.1) then
               rnorm(i) = 1.d0 - pow
            else
c     special treatment for unbound states (norm from 0 to rc is 1)
c     may want to add bound radius <> rc?
c               do k = 1,icarray(ic)+20
               do k = 1,np
                  rnl(k,i)=rnl(k,i)/sqrt(pow)
               enddo
               do k = 1,icarray(ic)
                  wavea(k) = rnl(k,i) * rnl(k,i)
               enddo
               pow = 2 * lo(i) + 2
               call radin(r,wavea,0,icarray(ic),h,pow)            
               rnorm(i) = 1.0-pow
               
c     Should be zero right? Should we check?
            endif
         enddo                     
      else
         do i=ncore+1,norb
            ic=i-ncore
            icarray(ic) = (log(rcall(ic)/r(1))/
     $           log(exp(1.0)))/h + 1
            do k = 1,np
               wavea(k) = rnla(k,i) * rnla(k,i)
     $              +rnlb(k,i) * rnlb(k,i)
               rnl(k,i)=rnla(k,i)
            enddo
            pow = 2 * lo(i) + 2
            call radin(r,wavea,0,icarray(ic),h,pow)
            if (ibd(ic).eq.1) then
               rnorm(i)=1.d0-pow
            else
c               do k = 1,icarray(ic)+20
               do k = 1,np
                  rnla(k,i)=rnla(k,i)/sqrt(pow)
                  rnl(k,i)=rnla(k,i)
               enddo
               do k = 1,icarray(ic)
                  wavea(k) = rnla(k,i) * rnla(k,i)
               enddo
               pow = 2 * lo(i) + 2
               call radin(r,wavea,0,icarray(ic),h,pow)            
               rnorm(i) = 1.0-pow
            endif
         enddo
      endif

c      if (ipsp.eq.0) then
c         do i=ncore+1,norb
c            do j=1,np
c               if (r(j).gt.20) goto 943
c               write(80+i-ncore,*) r(j),rnl(j,i),rvps(j,i)
c            enddo
c 943        continue
c         enddo
c      endif


c     Flip wfn s.t. it is >0 as r->oo
c     HEY didn't we already do this?? for nmax though?
c      do i=1,norbrel
c        if (rnla(maxim,i).lt.0.0) then
c          do k=1,np
c            rnla(k,i)=-rnla(k,i)
c          enddo
c        endif
c      enddo

c     radius of outermost peak
      if (irel.eq.0) then
         do i=ncore+1,norb
            ic=i-ncore
            rvapmax=0.0
            if (ibd(ic).ne.0) then
               do j=maxim-1,1,-1
                  if (rnl(j,i).lt.rnl(j+1,i)) then
                     rvap(i)=r(j+1)
                     if (rvapmax.lt.rvap(i)) rvapmax=rvap(i)
                     goto 944
                  endif
               enddo
            else
               rvap(i)=0.0
            endif
 944        continue
         enddo
      else
         do i=ncore+1,norb
            rvapmax=0.0
            if (ibd(i-ncore).ne.0) then
               do j=maxim-1,1,-1
                  if (rnla(j,i).lt.rnla(j+1,i)) then
                     rvap(i)=r(j+1)
                     if (rvapmax.lt.rvap(i).
     $                    and.(abs(rnla(j,i)).gt.1e-6)) then
                        rvapmax=rvap(i)
                        goto 945
                     endif
                  endif
               enddo
            else
               rvap(i)=0.0
            endif
 945        continue
         enddo
      endif

c     ***********************************************************************
c     Section 10:  Write final INFO
c     ***********************************************************************

      write(7,*) 

      if (niter.eq.1) then
         write(7,700)
      else
         write(7,500) niter
      endif

      write(7,701) etot,ebs,ehxc
      write(7,702) ehar,exc,evxc

      write (7,*)
      write (7,9222)

      do i = 1,ncore
         if (irel.eq.0) then
            write (7,9223) nlm(i),wnl(i),en(i)
         else
            if (lo(i).eq.0) then
               write (7,9523) nlm(i),wnl(i),en(i)
            else
               if (so(i).lt.0) then
                  write (7,9323) nlm(i),wnl(i),en(i)
               else
                  write (7,9423) nlm(i),wnl(i),en(i)
               endif
            endif
         endif
      enddo
      do i = ncore+1,norb
         ic=i-ncore
         if (irel.eq.0) then
            if (ibd(ic).ne.0) then
               write (7,9223) nlm(i),wnl(i),en(i),rnorm(i),rvap(i)
            else
               write (7,9224) nlm(i),wnl(i),en(i)
            endif
         else
            if (lo(i).eq.0) then
               if (ibd(ic).ne.0) then
                  write (7,9523) nlm(i),wnl(i),en(i),rnorm(i),rvap(i)
               else
                  write (7,9524) nlm(i),wnl(i),en(i)
               endif
            else
               if (so(i).lt.0) then
                  if (ibd(ic).ne.0) then
                     write (7,9323) nlm(i),wnl(i),en(i),rnorm(i),rvap(i)
                  else
                     write (7,9324) nlm(i),wnl(i),en(i)
                  endif
               else
                  if (ibd(ic).ne.0) then
                     write (7,9423) nlm(i),wnl(i),en(i),rnorm(i),rvap(i)
                  else
                     write (7,9424) nlm(i),wnl(i),en(i)
                  endif
               endif
            endif
         endif
      enddo
      
      if (ipsp.ne.0.and.rpcc.gt.1e-12) then
         write(7,9009) rpcc
         call radin(r,rscore,0,np,h,rsc)
         write(7,9011) rsc
      endif

 500  format(1x,i4,1x,'iterations: (units are Ry and bohr)')
 700  format(1x,'Converged in 1 iteration (probably reference state)')
 701  format(1x,'Energy: ',f16.8,2x,'Ebs: ',f16.8,2x,'Ehxc: ',f16.8)
 702  format(1x,'Eh    : ',f16.8,2x,'Exc: ',f16.8,2x,'Evxc: ',f16.8)

 9222 format(3x,'Orbital',4x,'Filling',7x,'Eigenvalues',4x,
     $     'Norm(rc->oo)',6x,'Peak')
 9223 format(4x,'|',i3,'>',5x,f6.3,5x,3f14.6)
 9224 format(4x,'|',i3,'>',5x,f6.3,5x,f14.6,"*",5x,
     $     '---------',5x,'--------')
 9323 format(4x,'|',i3,'->',5x,f6.3,5x,3f14.6)
 9324 format(4x,'|',i3,'->',5x,f6.3,5x,f14.6,6x,
     $     '---------',5x,'--------')
 9423 format(4x,'|',i3,'+>',5x,f6.3,5x,3f14.6)
 9424 format(4x,'|',i3,'+>',5x,f6.3,5x,f14.6,6x,
     $     '---------',5x,'--------')
 9523 format(4x,'|',i3,' >',5x,f6.3,5x,3f14.6)
 9524 format(4x,'|',i3,' >',5x,f6.3,5x,f14.6,6x,
     $     '---------',5x,'--------')

 9009 format(1x,'partial core radius : ', f10.4)
 9010 format(1x,'core charge         : ', f20.10)
 9011 format(1x,'partial core charge : ', f20.10)

      if (inl.ne.0) then
         
         call applyaug(rvloc,rvps(1,iloc))
         call ghostnl (rvloc)
         call applyaug(rvloc,rvcore(1,iloc))
      endif
      call flush(7)
      close(unit=7)

      return

 911  continue
c     911 = terminal error - exit scpot
      write(7,*) " !!ERROR!! Could not find self-consistent solution!"
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
      
