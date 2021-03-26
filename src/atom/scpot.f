      subroutine scpot(zeff,ixc,ipsp)
     
c     *************************************************************************
c     determine the self-consistent solution
c     *************************************************************************
     
      implicit double precision (a-h,o-z)
      
#include "PARAMHFS"
 
c     -------------------------------------------------------------------------
c     External (shared between C and Fortran) common blocks
c     -------------------------------------------------------------------------
      common /atomic/ xion,rnl(npdm,n0),nlm(n0),wnl(n0),en(n0),norb
      common /consts/ etol,vtol,maxit,isoft
      common /grid/ h,r1,z,r(npdm),np
      common /box/ iboxstart(10),iboxend(10),boxheight(10),numbox
      common /ilogder/ ilogder
      common /ibound/ ibd(n0)
      common /local/ iloc,idesign
      common /logarith/ rphas,elogmax,elogmin,dlwf(npl0,n0)
      common /nlpot2/ inl,indrc(n0),IB,ID,IG
      common /npm/ ncores,nvales,iskip,maxip
      common /rpcc/ rpcc,rpccz
      common /rscore/ rscore(npdm),rdd(npdm),rddd(npdm),rscoretot(npdm)
      common /rcall/ rcall(n0)
      common /filenames/ file_log
      common /results/ rnorm(30), etot, lghost(30)
      common /nmax/ nmax(n0),maxim
      common /totpot/ rvcore(npdm,n0),rvps(npdm,n0),rvcoul(npdm)
      common /valden/ rsval(npdm)
c     -------------------------------------------------------------------------
 
c     -------------------------------------------------------------------------
c     Internal (Fortran only) common blocks                
c     -------------------------------------------------------------------------
      common /nlpot1/ Flstar(npdm),phipsref(npdm,n0),phiguess(npdm,n0)
      common /partpot/ rvh(npdm),rvxc(npdm),rexc(npdm)
      common /nnn/ nnn
      common /nlcore/ rvloc(npdm)
      common /ipos/ ipos(n0),itermcount
c     -------------------------------------------------------------------------

c     *************************************************************************
c     local variables
c     *************************************************************************

      dimension rvn(npdm),rvf(npdm),eold(n0),p(npdm)
      dimension rsatom(npdm),f(npdm)
      dimension pder(npdm),ader(npdm)
      dimension dl(npl0,4,3,2), rsatom2(npdm)
      dimension rsold(npdm)
      
      dimension wavea(npdm),icarray(30)
c     *************************************************************************
c     Section 1:  set up tolerance params and zero out arrays.
c     *************************************************************************
      character*80 file_log
      character*1 xc(0:3)

      xc(0)='s'
      xc(1)='p'
      xc(2)='d'
      xc(3)='f'

      open(unit=7,file=file_log,form='formatted',access='append')

      if (ipsp.eq.1) then
         nlc=nlm(iloc)/100
         llc=(nlm(iloc) - nlc * 100)/10
      endif

      isoft=ipsp
      ipratt=6
      npratt = ipratt
      nig=0

      if (inl.ne.0) then
        do i=1,norb
          do j=1,np
            if (r(j).gt.rcall(i)) goto 993
          enddo
 993      continue
          indrc(i)=j
        enddo
      endif

      pr = 0

      do i = 1, np
        do j = 1, norb
          phipsref(i,j) = rnl(i,j)
          phiguess(i,j) = phipsref(i,j)
        enddo
      enddo

      do j = 1,norb
        eold(j) = en(j)
        nmax(j) = 0
        do i = 1,np
          rnl(i,j) = 0
        enddo
      enddo
      do i = 1,np
        rvcoul(i) = rvps(i,1) - rvcore(i,1)
      enddo

c     *************************************************************************
c     If DNL calc, write out box info
c     *************************************************************************

      if (ipsp.eq.1) then
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

      etotlast = 0

      do niter = 1,maxit
         if (niter.eq.1) then
            do i=1,norb
               ipos(i)=0
            enddo
         endif

c     ***********************************************************************
c     Section 2:  Set up to call schsl.
c     ***********************************************************************
         
         npratt = npratt + 1
         if (npratt.gt.ipratt) npratt = 0
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

            n = nlm(m)/100
            l = nlm(m)/10 - 10 * n
            
            if (inl.ne.0) then
               
               call applyaug(rvloc,rvps(1,iloc))
               
               do i=1,np
                  Flstar(i) = (rvps(i,m) - rvloc(i))
     $                 *phipsref(i,m)/r(i)
               enddo
               
               call schsl (m,n,l,en(m),nmax(m),rvloc,p,nig)

               call applyaug(rvloc,rvcore(1,iloc))
               
            else
               call schsl(m,n,l,en(m),nmax(m),rvps(1,m),p,nig)
            endif

            do i = 1, np
               phiguess(i,m) = p(i)
            enddo

c     *********************************************************************
c     End of changes made by NJR
c     *********************************************************************

            do i = 1,nmax(m)
               rsatom(i) = rsatom(i) + wnl(m) * p(i)**2
               dif = abs(p(i) - rnl(i,m))
               if (dif.gt.difmax) then
                  difmax = dif
                  imax = i
                  nlmmax = nlm(m)
               endif
               rnl(i,m) = p(i)
            enddo

            do i = nmax(m) + 1,np
               rnl(i,m) = 0
            enddo
            indey = max(indey,nmax(m))
            maxim = indey
            if (ibd(m).ne.0) then 
               edif = abs((en(m) - eold(m))/eold(m))
            else
               edif = 0.0
            endif
            eold(m) = en(m)
            if (edif.gt.edmax) nlmed = nlm(m)
            edmax = max(edif,edmax)

 911        continue
         enddo

c     ***********************************************************************
c     This section aids in convergence of difficult atoms such as Cu.
c     Added 1/23/98 by AMR.
c     ***********************************************************************
         
         alpha = 0.5
         if (niter.gt.3) then
            do i = 1,np
               rsatom(i) = alpha * rsatom(i) + (1.0-alpha) * rsold(i)
            enddo
         endif
      
c     ***********************************************************************
c     Section 4:  update self-consistent potentials and compute energy.
c     ***********************************************************************

         call hrtree(maxim,h,r,rsatom,rvh)

         if (rpcc.gt.1e-12.and.ipsp.eq.1) then
            do i=1,maxim
               rsatom2(i)=rsatom(i)+rscore(i)
            enddo
            call excorr(maxim,ixc,r,rsatom2,rvxc,rexc)
            
            do i = 1,maxim
               f(i) = 0.5 * rvh(i)*rsatom(i)/r(i)
     $              + rexc(i)*rsatom2(i)/r(i)
     $              - (rvxc(i)*rsatom(i)+rvh(i)*rsatom(i))/r(i)
            enddo
         else
            call hrtree(maxim,h,r,rsatom,rvh)
            call excorr(maxim,ixc,r,rsatom,rvxc,rexc)
            
            do i = 1,maxim
               f(i) = (0.5 * rvh(i) + rexc(i) - (rvxc(i)+rvh(i))) 
     $              * rsatom(i)/r(i)
            enddo
         endif    
         
         ehxc = 2
         if(maxim.eq.0) then
            goto 1001
         endif

         call radin(r,f,0,maxim,h,ehxc)
         ebs = 0
         do m = 1,norb
            ebs = ebs + wnl(m) * en(m)
         enddo
         etot = ebs + ehxc
         dvmax = 0
         xn2 = (zeff - xion) * 2

         do i = 1,maxim
            coul = rvh(i) + rvxc(i)
            dv = abs(coul - rvcoul(i))
                       
            if (dv.gt.dvmax) then
               idvmax = i
            endif
 9221       format(i6,3f20.10)
            dvmax = max(dv,dvmax)
            pr = 0.5
            if (npratt.gt.0) pr = pratt(rvn(i),rvf(i),rvcoul(i),coul)
            rvn(i) = rvcoul(i)
            rvf(i) = coul
            rvcoul(i)= pr * rvcoul(i) + (1 - pr) * coul
            
         enddo
         
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

         if (niter.eq.maxit) goto 1001
         
c     ***********************************************************************
c     Check if time to exit
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
c     Section 5:  predict trial eigenvals use 1st order perturb theory.
c     ***********************************************************************

         do i = 1,indey
            p(i) = (rvcoul(i) - rvn(i))/r(i)
         enddo
         do m = 1,norb
            if (ibd(m).ne.0) then
               do i = 1,nmax(m)
                  f(i) = p(i) * rnl(i,m)**2
               enddo
               l = nlm(m)/10 - 10 * (nlm(m)/100)
               xl = 2 * l + 2
               
               if(nmax(m).le.10) then
                  goto 1001
               endif
               call radin(r,f,0,nmax(m),h,xl)
               
               en(m) = en(m) + xl
               if (en(m).ge.0.0) then
                  a1 = l * (l + 1)
                  vmin = 1.0e6
                  do i = 1,np
                     vmin = min((rvps(i,m) + a1/r(i))/r(i),vmin)
                  enddo
                  en(m) = vmin * 0.5
                  write (7,505) niter,nlm(m),en(m)               
 505              format('0****positive trial eigenvalue predicted in
     $                 scpot, niter=',i4,3x,'nlm=',i4,5x,
     $                 'new corrected e =',e13.7)
                  if (en(m).ne.en(m)) then
                     goto 1001
                  endif
               endif
            endif
         enddo
      enddo
 1000 continue

c     *********************************************************************
c     computation of logarithmic derivatives for an all-electron run.
c     NJR 5/29/97
c     *********************************************************************

      if (ilogder.eq.1) then
        if (ipsp.eq.0) then
          ic=0
          do m=ncores+1,norb
            ic=ic+1
            l = nlm(m)/10 - 10 * (nlm(m)/100)
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
     $             rphas,dl,npl0,nr,pder,wfld,0,1.0,1.0,1.0)
              dlwf(i,ic) = wfld
            enddo
          enddo
        endif
        if (ipsp.eq.1) then
          if (inl.eq.0) then
            ic=0
            do m=1,nvales
              ic=ic+1
              l = nlm(m)/10 - 10 * (nlm(m)/100)
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
     $               rphas,dl,npl0,nr,pder,wfld,0,1.0,1.0,1.0)
                dlwf(i,m) = wfld
              enddo
            enddo
          else

            do m=1,nvales
              n = nlm(m)/100
              l = nlm(m)/10 - 10 * n
              call nllogd(m,n,l,en(m),nmax(m))
            enddo
          endif
        endif
      endif

      do i = ncores+1,norb
         if (rnl(nmax(i)-1,i).lt.0.0) then
            do j = 1,np
               rnl(j,i) = -rnl(j,i)
            enddo
         endif
      enddo

      zv = z
      do i = 1,ncores
        zv = zv - wnl(i)
      enddo

      do j = 1,np
         rsval(j) = 0
      enddo
      do i = ncores+1,norb
         do j = 1,nmax(i)
            rsval(j) = rsval(j) + wnl(i) * rnl(j,i)**2
         enddo
      enddo

      if (ipsp.eq.0) then
         do j = 1,np
            rscore(j) = 0
         enddo
         do i = 1,ncores
            do j = 1,nmax(i)
               rscore(j) = rscore(j) + wnl(i) * rnl(j,i)**2
            enddo
         enddo
         call radin(r,rscore,0,np,h,rsc)
      endif

      write(7,*) 

      if (niter.eq.1) then
         write(7,700)
      else
         write(7,500) niter
      endif

      write(7,701) etot,ebs,ehxc

 500  format(1x,'After  ',i4,1x,'iterations...')
 700  format(1x,'Converged in 1 iteration (probably reference state)')
 701  format(1x,'Energy: ',f16.8,2x,'Ebs: ',f16.8,2x,'Ehxc: ',f16.8)

      do i = ncores + 1,norb
         icarray(i) = (log(rcall(i-ncores)/r(1))/
     $        log(exp(1.0)))/h + 1

         do k = 1,np
            wavea(k) = rnl(k,i) * rnl(k,i)
         enddo
         ll = nlm(i) - 100 * (nlm(i)/100)
         ll = ll/10
         pow = 2 * ll + 2
         call radin(r,wavea,0,icarray(i),h,pow)

         if (ibd(i).eq.1) then
            rnorm(i) = 1 - pow
         else
c     special treatment for unbound states (norm from 0 to rc is 1)
            do k = 1,np
               rnl(k,i)=rnl(k,i)/sqrt(pow)
            enddo
            do k = 1,icarray(i)
               wavea(k) = rnl(k,i) * rnl(k,i)
            enddo
            pow = 2 * ll + 2
            call radin(r,wavea,0,icarray(i),h,pow)            
            rnorm(i) = pow
         endif
      enddo
      
      do i=1,norb
        isign=1
        if (rnl(maxim,i).lt.0.0) then
          do k=1,np
            rnl(k,i)=-rnl(k,i)
          enddo
        endif
      enddo

      write (7,*)
      write (7,2002)

      do i = 1,ncores
         write (7,2001) nlm(i),wnl(i),en(i)
      enddo
      do i = ncores+1,norb
         write (7,2001) nlm(i),wnl(i),en(i),rnorm(i)
      enddo

 2002 format(3x,'Orbital',4x,'Filling',11x,'Eigenvalues',8x,
     $     'Norm(rc->oo)')
 2001 format(4x,'|',i3,'>',5x,f6.3,5x,2f19.10)


      if (ipsp.ne.0.and.rpcc.gt.1e-12) then
         write(7,9009) rpcc
         call radin(r,rscore,0,np,h,rsc)
         write(7,9011) rsc
      endif

 9009 format(1x,'partial core radius : ', f10.4)
 9010 format(1x,'core charge         : ', f20.10)
 9011 format(1x,'partial core charge : ', f20.10)



      if (inl.eq.1) then
         call applyaug(rvloc,rvps(1,iloc))
         call ghostnl (inl,iloc,rvloc)
         call applyaug(rvloc,rvcore(1,iloc))
      endif

      call flush(7)
      close(unit=7)

      return
      
 1001 continue
      write (7,601) niter
 601  format('****terminal error in scpot,niter=',i6,3x,
     $     'exceeds maxit or e=NAN')
      eigend=0.0
      call flush(7)
      close(unit=7)

      return
      
      end
      
