      subroutine interp(ifrl, itransf,ixc)
      
c     *************************************************************************
c     This code reads in the relativ. atomic information from fort.22. and 
c     processes it for the scalar-relativistic pseudopotential procedure.
c     The file *.ae is generated.
c     *************************************************************************

      implicit double precision (a-h,o-z)
      
      external val
      
#include "PARAMHFS"

      common /grid/ h,r1,z,r(npdm),np
      common /atomic/ xion,rnl(npdm,n0),nlm(n0),wnl(n0),en(n0),norb
      common /rcall/ rcall(nvale0)
      common /npm/ ncores,nvales,iskip,maxip
      common /rpcc/ rpcc,rpccz
      common /nll/nll    
      common /totpot/ rvcore(npdm,n0),rvps(npdm,n0),rvcoul(npdm)
      common /nmax/ nmax(nvale0),maxim
      common /lall/ lall(nvale0)
      common /rscore/ rscore(npdm),rdd(npdm),rddd(npdm),rscoretot(npdm)
      common /atom3/ rsvale(npdm)
      common /atom4/ v(npdm),rv(npdm),p(npdm),rsatom(npdm),
     $               c1(npdm),rvh(npdm),rvxc(npdm),rexc(npdm)
      common /relwf/ wfu(npdm,nvale0),wfd(npdm,nvale0)
      common /frlwf/ nlmp(nvale0),wnlp(nvale0),ev(nvale0*2),
     $     rnor(nvale0*2)
     
c     *************************************************************************
c     local variables
c     *************************************************************************

      dimension enp(nvale0),nmaxp(nvale0)
      dimension rin(npdm),wfin(npdm,nvale0*2),vin(npdm)
      dimension lo(nvale0*2),so(nvale0*2),no(nvale0*2)
      dimension wf2(npdm,nvale0*2),f(npdm),nu(nvale0)
      dimension nd(nvale0),cdc(npdm),rexct(npdm)
      dimension rvxct(npdm)
      dimension ibp(nvale0),indall(nvale0)
      dimension eigenu(100),eigend(100)
      dimension fe_old(nvale0),fe(nvale0), phom(npdm,nvale0)
      dimension rvscr(npdm)
      dimension xah(npdm),xdh(npdm),xch(npdm),xlh(npdm)
      dimension phom1(npdm),rsvale1(npdm),temp1(npdm),temp2(npdm)
      dimension ratio(nvale0),tempx(npdm),rclp(nvale0*2)
      
c     *************************************************************************
c     gjt: this parameter corresponds to the setting in f_src_rel/param.h
      parameter (nrmax=10000, norbp=40)
c     gjt: this is for communication between atm() and interp()      

      common /atmwave/ nr,numorb,rr(nrmax),aar(nrmax,30),bbr(nrmax,30),
     $     ccdc(nrmax),eev(norbp),llo(norbp),sso(norbp),
     $     nno(norbp)
      common /filenames/ file_log
      character*80 file_log
c     *************************************************************************
      
c     *************************************************************************
c     redirect stdout to the log file
c     *************************************************************************

      open(unit=7,file=file_log,form='formatted',access='append')

c     *************************************************************************
c     initialize some local arrays
c     *************************************************************************

      do i = 1,npdm
        rvxct(i) = 0.0
        rexct(i) = 0.0
        rvh(i) = 0.0
        rvxc(i) = 0.0
        rexc(i) = 0.0
        do j = 1,nvales
          rnl(i,j) = 0.0
        enddo
      enddo

      do iorb = 1,norb
        nlmp(iorb) = nlm(iorb)
        wnlp(iorb) = wnl(iorb)
        enp(iorb) = en(iorb)
      enddo
      ncores = norb - nvales
      zv = xion
      do iorb = ncores + 1,norb
        zv = zv + wnlp(iorb)
      enddo
      
c     Here we generate the new information.
c     Need:  maxim,en,nmax,rnl
      do i = 1,nvales
c         q1(i) = 0.0
         nlm(i) = nlmp(i+ncores)
         wnl(i) = wnlp(i+ncores)
      enddo
      
      irld = 0
      etot = 0.0
      ebs = 0.0
      ehxc = 0.0
      etr = 0.0
      ebsvr = 0.0
      ehxcvr = 0.0
      ecv = 0.0
      
c     *************************************************************************
c     Now we read in the information from the relativistic program output.
c     wfin holds r*wf and vin holds r*v, both on the rin grid.
c     *************************************************************************
      icount = 0
      numu = 0
      numd = 0
  9   continue      
c     gjt: either read AE info from file of expect through common block  
      if (itransf.gt.0) then
        read (22,*) nr, numorb
        do 10 j = 1,nr
          read (22,*) rin(j),wfin(j,icount+1),vin(j),cdc(j)
 10    continue
       read (22,*) ev(icount+1),lo(icount+1),so(icount+1),no(icount+1)
 121   format (4E32.24)
      else
        do j = 1,nr
          rin(j) = rr(j)
          wfin(j,icount+1) = aar(j,icount+1)
          cdc(j) = ccdc(j)
        enddo
        ev(icount+1) = eev(icount+1)
        lo(icount+1) = llo(icount+1)
        so(icount+1) = sso(icount+1)
        no(icount+1) = nno(icount+1)
      endif
            
      icount = icount + 1
c      write (7,*) 'icount / numorb = ', icount, numorb
      if (so(icount).lt.0.1.or.lo(icount).eq.0) then
        numd = numd + 1
        nd(numd) = icount
c        write (7,*) 'found j=l-1/2'
        nlmp(icount) = nlm(numd)
        rclp(icount) = rcall(numd)
        ll = nlm(numd) - 100 * (nlm(numd)/100)
        ll = ll/10
        wnlp(icount) = wnl(numd)*float(ll)/float(2*ll + 1)
      endif
      if (so(icount).gt.0.1.or.lo(icount).eq.0) then
        numu = numu + 1
        nu(numu) = icount
c        write (7,*) 'found j=l+1/2'
        nlmp(icount) = nlm(numu)
        rclp(icount) = rcall(numu)
        ll = nlm(numu) - 100 * (nlm(numu)/100)
        ll = ll/10
        wnlp(icount) = wnl(numu)*float(ll + 1)/float(2*ll + 1)
      endif
      if (icount.lt.numorb) goto 9
      
      
c     New section to include unfilled f orbital in core.  AMR 7/25/91
c      ifld = 0
c      do 901 i = 1,numd
c         if (lo(nd(i)).eq.3) ifld = i
c 901  continue
c      iflu = 0
c      do 902 i = 1,numu
c         if (lo(nu(i)).eq.3) iflu = i
c 902  continue
c      ii = nlmp(ncores)/100
c      ii = nlmp(ncores) - 100 * ii
c      ii = ii/10
c      if (ifld.ne.0) then 
c         if (ii.ne.3) then
c            write (7,*) 'make f orbital last core orbital in .hfsae'
c            stop
c         endif
c         numd = numd - 1
c         if (ifld.ne.numd+1) then
c            write (7,*) 'put valence wfs in order of increasing'
c            write (7,*) 'angular momentum.'
c            stop
c         endif
c         wg = wnlp(ncores) * 3.0/7.0
c         do 903 i = 1,nr
c            cdc(i) = cdc(i) + wg * wfin(i,nd(numd+1))**2
c 903     continue
c         numorb = numorb - 1
c      endif
c      if (iflu.ne.0) then 
c         if (ii.ne.3) then
c            write (7,*) 'make f orbital last core orbital in .hfsae'
c            stop
c         endif
c         numu = numu - 1
c         if (iflu.ne.numu+1) then
c            write (7,*) 'put valence wfs in order of increasing'
c            write (7,*) 'angular momentum.'
c            stop
c         endif
c         wg = wnlp(ncores) * 4.0/7.0
c         do 904 i = 1,nr
c            cdc(i) = cdc(i) + wg * wfin(i,nu(numu+1))**2
c 904     continue
c         numorb = numorb - 1
c      endif
c     More than 5 valence wavefunctions are available only if f
c     wavefunctions are included.  f pseudization does not yet exist,
c     so we do not process these f wavefunctions, thereby relegating
c     them to the core.  This will be improved in the future.  AMR 7/25/91
c     End of f processing section.



c      write (7,*) 'nu ',(nu(i),i=1,numu)
c      write (7,*) 'nd ',(nd(i),i=1,numd)
c      write (7,*)
c      write (7,*) 'nlmp ',(nlmp(i),i=1,numorb)
c      write (7,*) 'wnlp ',(wnlp(i),i=1,numorb)
      
      irel = 1
      if (numd.ne.numu) irel = 0
c      write (7,*) 'irel = ',irel
      
c     *************************************************************************
c     Now we create the new grid, and interpolate the wfs and pots onto it
c     *************************************************************************

      r(1) = r1/z**(1./3.)
      xmult = exp(h)
      do i = 2,np
        r(i) = r(i-1) * xmult
      enddo
      
      do i = 1,numorb
        jmax = nr
        jflag = 0
        do j = nr,1,-1
c          write (7,*) 'wfin ',j,rin(j),wfin(j,i),i
          f(j) = wfin(j,i)
          if (f(j).ne.0.0.and.jflag.eq.0) then
            jflag = 1
            jmax = j
          endif
        enddo
        do 12 j = 1,np
          oldval = -999.0
          iflag = 0
          if (r(j).gt.rin(jmax)) then
            wf2(j,i) = 0.0
            goto 12
          endif
          do m = 1,10
            fval = val(f,rin,nr,r(j),m)
            diff = oldval - fval
            if (abs(diff).lt.1.0e-8) then
              wf2(j,i) = fval
              iflag = 1
            endif
            oldval = fval
          enddo
c          write (7,*) 'wf2 ',j,r(j),wf2(j,i),i
          if (iflag.eq.0) then
            write (7,*) 'wf interpolation failed ',j,i
            write (7,*) 'nlm',nlm(i)
            stop
          endif
          if (abs(wf2(j,i)).gt.0.0) nmaxp(i) = j
 12     continue
      enddo

      do j = 1,nr
        f(j) = cdc(j)
c        write(900,*) rin(j),f(j)
      enddo

      do 26 j = 1,np
        oldval = -999.0
        iflag = 0
        if (r(j).gt.rin(nr-5)) then
          rscore(j) = 0.0
          goto 26
        endif

        do m = 1,20
          fval = val(f,rin,nr,r(j),m)
c          write (7,*) j,m,r(j),fval
          diff = oldval - fval
          if (abs(diff).lt.1.0e-7) then
            rscore(j) = fval
            if (rscore(j).lt.0.0) rscore(j) = 0.0
            iflag = 1
          endif
          oldval = fval
        enddo

        if (iflag.eq.0) then
          write (7,*) 'chd interpolation failed ',j
          stop
        endif

 26   continue
 
      pow = 2.0
      call radin(r,rscore,0,np,h,pow)
      xmul = (z-zv)/pow
      write (7,9000) z
      write (7,9010) zv
      do i = 1,np
        rscore(i) = rscore(i) * xmul
        rsvale(i) = 0.0
	rsvale1(i)=0.0
      enddo

      pow = 2.0
      call radin(r,rscore,0,np,h,pow)
      do i = 1,numorb
        xl = float(lo(i))
        if (so(i).gt.0.1.and.lo(i).ne.0) wsp = (xl+1.0)/(xl+xl+1.0)
        if (so(i).lt.0.1.and.lo(i).ne.0) wsp = xl/(xl+xl+1.0)
        if (lo(i).eq.0) wsp = 1.0
        ii = 0
        do k = 1,nvales
           nn=nlm(k)/100
          ll = nlm(k) - 100*(nlm(k)/100)
          ll = ll/10
          if (ll.eq.lo(i).and.nn.eq.no(i)) ii = k
        enddo
        if (ii.eq.0) then
          write (7,*) 'Didnt find weight '
          stop
        endif
        do j = 1,np
          f(j) = wf2(j,i)**2
        enddo
        pow = float(lo(i)*2+2)
        call radin(r,f,0,np,h,pow)
c        write (7,*) 'valence charge = ',i,pow

        psqrt = sqrt(pow)
        do j = 1,np
          wf2(j,i) = wf2(j,i)/psqrt
          rsvale(j) = rsvale(j) + f(j)/pow*wsp*wnl(ii)

          if (rsvale(j).lt.0.0) rsvale(j) = 0.0
        enddo

        pow = 2.0
        call radin(r,rsvale,0,np,h,pow)

c        write (7,*) 'total valence charge now = ',pow
      enddo
      write (7,9020) pow 

 9000 format(' Z atomic :',f10.6)
 9010 format(' Z valence :',f10.6)
 9020 format(' Total valence charge :',f20.10)

c     *************************************************************************
c     New section written by I. Grinberg for mapping relativistic
c     solutions to non-relativistic equations.
c     *************************************************************************
 
      ioutput = 0

      do jj = 1,1+irel
        iwrite = 9 + jj + irel
c       iwrite = 10 for irel = 0, iwrite = 11 and 12 for irel = 1      
c       iwrite = 11 for spin up, iwrite = 12 for spin down.
        do kk = 1,nvales
          en(kk) = ev(kk)
          if(iwrite.eq.11) then
            mm = nd(kk)
            eigend(kk)=ev(mm)
            write (7,*) 'DOWN nlm = ',nlm(kk),' Eigenvalue = ',
     $             eigend(kk)
            do j=1,np
              wfd(j,kk)=wf2(j,mm)
            enddo
          endif
          if(iwrite.eq.12) then 
            mm = nu(kk)
            eigenu(kk)=ev(mm)
            write (7,*) 'UP   nlm = ',nlm(kk),' Eigenvalue = ',
     $             eigenu(kk)
            do j=1,np
              wfu(j,kk)=wf2(j,mm)
            enddo
          endif
        enddo
      enddo
 
c     *************************************************************************
c     2 relativistic -> 1 non-relativistic starts here
c     *************************************************************************

      do i=1,nvales
        fe_old(i)=0.0
        ll = nlm(i) - 100 * (nlm(i)/100)
        ll = ll/10
        xt = float(4 * ll + 2)
        xu = float(2 * ll + 2)/xt
        xd = float(2 * ll)/xt
        do j=1,np
          if (rsvale(j).gt.0.0) maxim = j
        enddo
c       calculate charge density  and integrate it to get total valence charge
        do j=1,maxim
          phom1(j)=wfu(j,i)**2*xu+wfd(j,i)**2*xd
        enddo            
c       phom1 here has 2 powers of  r
        pow=2*ll+2
        call radin(r,phom1,0,maxim,h,pow)
        totalint=pow
        do ii = 1,np
          if (rcall(i).gt.r(ii)) indall(i) = ii
        enddo
c       calculate the within breakpoint charge, outside breakpoint charge     
        ibp(i) = indall(i) - 5
c       "breakpoint" for norm in tail region WAS indall-19 now indall-5.


cEJW    now compute the norm from rc->oo
        pow=2*ll+2
        call radin(r,phom1,0,indall(i),h,pow)
        rnor(i) = 1- pow
        write (7,*) 'norm rc to infinity = ',i,rnor(i)
cEJW

        pow=2*ll+2
        call radin(r,phom1,0,ibp(i),h,pow)
        powall = 1- pow
        write (7,*) 'norm rbp to infinity = ',i,powall

c       conmax insert
        if (ioutput.eq.1) then
          open(unit=56,file='AENORM',form='formatted')
          write(56,*) powall
        endif

        rcint=pow
        partint=abs(totalint-rcint)
        
c       assign the filling which must be preserve in the new nr psi and
c       get first guess for the nr wavefunc from the valence charge density
        fe_old(i)=partint
        do j=1,maxim
          phom(j,i)=sqrt(wfu(j,i)**2*xu+wfd(j,i)**2*xd)
        enddo
c       phom1 here has 1 power of  r              
c       overwrite guess wavefunc at breakpoint with arbitrary formula,
c       which still preserves the correct within breakpoint charge
        do j=1,ibp(i)
          rl=r(j)**(ll+1)
          chiro=phom(ibp(i),i)/r(ibp(i))**(ll+1)
          oneminus=(1-(r(j)/r(ibp(i)))**4)
          temp1(j)=(rl*chiro)**2
          temp2(j)=(rl*oneminus)**2
          tempx(j)=2*rl*rl*chiro*oneminus
        enddo
        pow=2*ll+2
        call radin(r,temp1,0,ibp(i),h,pow)
        temp_int1=pow
        pow=2*ll+2
        call radin(r,temp2,0,ibp(i),h,pow)
        temp_int2=pow
        pow=2*ll+2
        call radin(r,tempx,0,ibp(i),h,pow)
        temp_intx=pow
        const1=(-temp_intx+sqrt(temp_intx**2-4*temp_int2
     $        *(temp_int1-1+fe_old(i))))/2/temp_int2
         
        const2=(-temp_intx-sqrt(temp_intx**2-4*temp_int2
     $        *(temp_int1-1+fe_old(i))))/2/temp_int2
        c=const1
        if(const1.ge.0.0) c=const1
        if(const2.ge.0.0) c=const2
        do j=1,ibp(i)
          rl=r(j)**(ll+1)
          chiro=phom(ibp(i),i)/r(ibp(i))**(ll+1)
          oneminus=(1-(r(j)/r(ibp(i)))**4)
          phom(j,i)=rl*(chiro+c*oneminus)
        enddo
      enddo
      
      c=1.0
      do ii=1, np
        rsvale(ii)=0
      enddo
      do i=1, nvales
        ll = nlm(i) - 100 * (nlm(i)/100)
        ll = ll/10
        do ii=1, np
          phom1(ii)=0
        enddo
c       calculate the first guess for valence charge density
        do ii=1,np
          rsvale(ii)=rsvale(ii)+phom(ii,i)**2*wnl(i)
          phom1(ii)=phom1(ii)+phom(ii,i)**2
c          write(7,*) ii, rsvale(ii), phom1(ii)
        enddo
c        write (7,*) i, ii, phom(ii,i)
      enddo
      
c     *************************************************************************
c     iterative solution for the nr wavefunc with correct eigenvalues and 
c     within and outside breakpoint fillings 
c     *************************************************************************
      
      do ijk=1,100
c       calculate total atomic charge density and potentials
        do i=1,np
          rsatom(i)=rsvale(i)+rscore(i)
          rvxct(i) = 0.0
          rexct(i) = 0.0
          rvh(i) = 0.0
        enddo
        call excorr (np,ixc,r,rsatom,rvxct,rexct)
        call hrtree (np,h,r,rsatom,rvh)
        do i=1,np
          rvscr(i)=rvh(i)+rvxct(i)
          rv(i)=-2*z+rvscr(i)
        enddo
        
        do j =1,nvales
c         calculate eigenvalues, prepare for inward solve
          ll = nlm(j) - 100 * (nlm(j)/100)
          ll = ll/10
          xt = float(4 * ll + 2)
          xu = float(2 * ll + 2)/xt
          xd = float(2 * ll)/xt            
          x2=(ll+0.5)**2
          h2=h*h/12
          ei=eigenu(j)*xu+eigend(j)*xd
          do i = 1, np
            v(i) = ((rv(i) - ei * r(i)) * r(i) + x2) * h2
          enddo
c         do inward solve
          do i = 1,np
            if (rsvale(i).gt.0.0) maxim = i
            if (abs(phom(i,j)).gt.0.0) nmax(j) = i
          enddo
          if (nmax(j).gt.maxim) maxim = nmax(j)
          imax=maxim
          imat=ibp(j)+1
          imp1=imat+1
          xah(imat-1) = 1 - v(imat-1)
          xdh(imat-1) = -(2 + 10 * v(imat-1))
          xah(imat) = 1 - v(imat)
          xdh(imat) = -(2 + 10 * v(imat))
          xlh(imat) = xdh(imat) 
          xch(imat) = -(xah(imat-1)) * phom(imat-1,j)/sqrt(r(imat-1))
          do i = imp1,imax+1
            xah(i) = 1 - v(i)
            xch(i) = -xch(i-1) * xah(i-1)/xlh(i-1)
            xdh(i) = -(2 + 10 * v(i))
            xlh(i) = xdh(i) - xah(i) * xah(i-1)/xlh(i-1)
          enddo
          rootv = h * sqrt(v(imax)/h2)
          ap = exp(-rootv)
          phom(imax,j) = xch(imax)/(xlh(imax) + xah(imax+1) * ap)
          do i = imax-1,imat,-1
            phom(i,j) = (xch(i) - xah(i+1) * phom(i+1,j))/xlh(i)
          enddo
c         phom now has 1/2 powers of r, from ibp+1 to np, and 1 power of r
c         from 1 to ibp
        enddo
        
c       get new charge density
        sumff=0
        do i=1,nvales
          ll = nlm(i) - 100 * (nlm(i)/100)
          ll = ll/10
          fe(i)=0
          do ii=1,ibp(i)-1
            phom1(ii)=phom(ii,i)*phom(ii,i)
          enddo
          do ii=ibp(i)+1,np
            phom(ii,i)=phom(ii,i)*sqrt(r(ii))
          enddo
c         phom now has 1 power of r (1/2 +1/2 = 1), ibp point is 
c         untouched by the inward solve so it still has correct 1 power of r
          do ii=ibp(i),np
            phom1(ii)=phom(ii,i)*phom(ii,i)
          enddo
c         phom1 here has 2 powers of r (1+1=2)
c         get ratio between the filling past breakpoint for the solved 
c         wavefunc and the correct filling (integrated charge density)
          pow=2*ll+2
          call radin(r,phom1,0,np,h,pow)
          totalint=pow
          pow=2*ll+2
          call radin(r,phom1,0,ibp(i),h,pow)
          rcint=pow
          partint=abs(totalint-rcint)
          fe(i)=partint
          ratio(i)=fe_old(i)/fe(i)
          rtsqrt=sqrt(ratio(i))
c         multiply the wavefunction from breakpoint on by the sqrt of the ratio
c         to fix the filling 
c         (since ratio is obtained from the charge density which is psi**2)
          do ii=ibp(i),np
            phom(ii,i)=phom(ii,i)*rtsqrt
          enddo
          do ii=1,np
            phom1(ii)=phom(ii,i)*phom(ii,i)
          enddo
c         phom has 1 power of r, phom1 has 2 power of r
c         redo the within breakpoint part to make it continuous with the outside
c         breakpoint part
          do j=1,ibp(i)
            rl=r(j)**(ll+1)
            chiro=phom(ibp(i),i)/r(ibp(i))**(ll+1)
            oneminus=(1-(r(j)/r(ibp(i)))**4)
            temp1(j)=(rl*chiro)**2
            temp2(j)=(rl*oneminus)**2
            tempx(j)=2*rl*rl*chiro*oneminus
          enddo
          pow=2*ll+2
          call radin(r,temp1,0,ibp(i),h,pow)
          temp_int1=pow
          pow=2*ll+2
          call radin(r,temp2,0,ibp(i),h,pow)
          temp_int2=pow
          pow=2*ll+2
          call radin(r,tempx,0,ibp(i),h,pow)
          temp_intx=pow
          const1=(-temp_intx+sqrt(temp_intx**2-4*temp_int2
     $          *(temp_int1-1+fe_old(i))))/2/temp_int2
          const2=(-temp_intx-sqrt(temp_intx**2-4*temp_int2
     $          *(temp_int1-1+fe_old(i))))/2/temp_int2
          c=const1
          if(const1.ge.0.0) c=const1
          if(const2.ge.0.0) c=const2
          do ii=1,ibp(i)
            rl=r(ii)**(ll+1)
            chiro=phom(ibp(i),i)/r(ibp(i))**(ll+1)
            oneminus=(1-(r(ii)/r(ibp(i)))**4)
            phom(ii,i)=rl*(chiro+c*oneminus)
            phom1(ii)=phom(ii,i)**2
          enddo
        enddo
        
c       check if the ratio of old wavefunction past breakpoint filling 
c       to true filling is 1.0 if yes converged, if not continue with 
c       the iterative loop
        rttot=0.0
        do i=1,nvales
          rttot=rttot+abs((ratio(i)-1.00))
        enddo
        if(abs(rttot).lt.(1e-12)) then 
          iijjkk=1001
          goto 998
        endif
        do ii=1, np
          rsvale(ii)=0
        enddo
c       calculate new rsvale, phom has 1 power of r
        do i=1, nvales
          do ii=1,np
            rsvale(ii)=rsvale(ii)+phom(ii,i)**2*wnl(i)
          enddo
c         write (7,*) ijk, i, phom(1,i), wnl(i)
        enddo
      enddo
 
c     *************************************************************************
c     done with iterative solution and 
c     *************************************************************************
 
 998  write(7,*) "CONVERGED!",ijk,rttot
c     if converged use new wavefunctions to calculate all quantities that
c     will be outputed (rvcoul, rnl, etc.) 
      do j=1,nvales
        ll = nlm(j) - 100 * (nlm(j)/100)
        ll = ll/10
        xt = float(4 * ll + 2)
        xu = float(2 * ll + 2)/xt
        xd = float(2 * ll)/xt
        en(j)=eigenu(j)*xu+eigend(j)*xd
        do ii=1,np
          rnl(ii,j)=phom(ii,j)
c         rnl has 1 power of r, i.e. rnl proportional to psi*r
          if (abs(rnl(ii,kk)).gt.0.0) nmax(kk) = ii
        enddo
      enddo
      
      do i = 1,np
        rsatom(i) = rsvale(i) + rscore(i)
        if (rsatom(i).gt.0.0) maxim = i
      enddo
      pow = 2.0
      call radin(r,rsatom,0,maxim,h,pow)
      pow = 2.0
      call radin(r,rscore,0,maxim,h,pow)
      call excorr (maxim,ixc,r,rsatom,rvxct,rexct)
      call hrtree (maxim,h,r,rsatom,rvh)
      do  i = maxim+1,np
        rvh(i)  = z + z - xion - xion
      enddo
      xm = (z+z-xion-xion)/rvh(maxim)
      dif1 = abs(xm-1.0)
c     gjt: tol2 is defined in PARAMOPT
      if (dif1.gt.tol2) then
        write (7,*) 'Use more grid points. Potentials inaccurately '
        write (7,*) 'integrated from charge density. ',tol2
        write (7,*) dif1,(z+z-xion-xion),rvh(maxim),maxim 
      endif
      sum=0.0
      sum2=0.0
      do  i = 1,np
        rvh(i) = rvh(i) * xm
        rvcoul(i) = rvh(i) + rvxct(i)
        f1 = -z-z+rvcoul(i)
        do  j = 1,nvales
          rvcore(i,j) = -2 * z
          rvps(i,j) = f1
        enddo
      enddo
      tch=2.0
      
c     write out all the necessary information in 1 nr potential format
      im = maxim
      do i = 1,nvales
c       conmax insert
        if (ioutput.eq.1) then
          open(unit=55,file='AEEIG',form='formatted')
          write(55,*) en(i)
        endif

        write(7,*) 'nlm = ',nlm(i),' Averaged Eigenvalue = ',en(i)
      enddo
      
c     *************************************************************************
c     write the appropriate AE information to the appropriate binary file
c     *************************************************************************
 
      if (ifrl.eq.0) then
      
      else
         
c     frl output ; do norms and
c     ensure that the wave function converges to 0 from the positive side
         do i=1,numorb           
            do j=np,1,-1
               if (abs(wf2(j,i)).gt.1e-5) goto 1010
            enddo
 1010       continue
            if (wf2(j,i).lt.0.) then
               do j=np,1,-1
                  wf2(j,i) = -wf2(j,i)
               enddo
            endif
cEJW
            do j=1,np
               if (rclp(i).gt.r(j)) irc=j
               phom1(j)=wf2(j,i)**2
            enddo
            pow=float(lo(i)*2+2)
            call radin(r,phom1,0,irc,h,pow)
            rnor(i) = 1-pow
c            write(7,*) 'DEBUG: i,ll',i,irc,pow
cEJW

         enddo
         
      endif
      
c     *************************************************************************
c     finish up
c     *************************************************************************

      close(unit=7)
      end
