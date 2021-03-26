      subroutine ghostnl (inl,iloc,rvloc)
           
c     *************************************************************************
c     inl     unused
c     nloc    in
c     rvloc   in      array
c     *************************************************************************

c     *************************************************************************
c     Do the Gonze ghost analysis
c     *************************************************************************
c     This program finds the ground state e0l and first excited state e1l
c     of the local potential.
c     For nonlocal components, we calculate elkb.  The ghost theorem of
c     Xavier Gonze states the following:
c     For elkb>0, and eat>e1l the potential has a ghost below eat.
c     For elkb<0, and eat>el0 the potential has a ghost below eat.
c     Here eat is the reference eigenvalue of the nonlocal angular momentum.
c     *************************************************************************
      
      implicit double precision (a-h,o-z)
      
#include "PARAMHFS"
 
c     -------------------------------------------------------------------------
c     External (shared between C and Fortran) common blocks
c     -------------------------------------------------------------------------
      common /atomic/ xion,rnl(npdm,n0),nlm(n0),wnl(n0),en(n0),norb
      common /grid/ h,r1,z,r(npdm),np
      common /box/ iboxstart(10),iboxend(10),boxheight(10),numbox
      common /npm/ ncores,nvales,iskip,maxip
      common /results/ rnorm(30), etot, lghost(30)
      common /ibound/ ibd(n0)
c     -------------------------------------------------------------------------
 
c     -------------------------------------------------------------------------
c     Internal (Fortran only) common blocks                
c     -------------------------------------------------------------------------
      common /iterm/ iterm
      common /nlpot1/ Flstar(npdm),phipsref(npdm,n0),phiguess(npdm,n0)
      common /totpot/ rvcore(npdm,n0),rvps(npdm,n0),rvcoul(npdm)
      common /nmax/ nmax(nvale0),maxim
c     -------------------------------------------------------------------------
      
      dimension rvloc(npdm)
      
c     *************************************************************************
c     local variables
c     *************************************************************************
      
      dimension fr(npdm),p(npdm),gghost(npdm)
      dimension el0(n0),el1(n0),isemi(n0)
      character*1 xc(0:3)

      xc(0)='s'
      xc(1)='p'
      xc(2)='d'
      xc(3)='f'
      
      write(7,*) 
      write(7,*) '---Non-local ghost testing---'
      if (ibd(iloc).eq.0) then
         write(7,*) "!WARNING! Non-local ghost testing is not reliable"
         write(7,*) "!WARNING!   when the local pot is unbound!"
c         write(7,*) "Please change your choice of local potential, = "
c     $        , nlm(iloc)
c         stop
      endif

      do i = 1,np
        gghost(i) = 0
        flstar(i) = 0
      enddo
      
      ighost = 0
      ig=1

      nnn = nlm(iloc)/100
      lll = (nlm(iloc) - nnn * 100)/10
      
      write(7,9101) nnn,xc(lll)
      write(7,*)
      
      do k=1,nvales
        isemi(k)=0
      enddo
      do k=1,nvales
        do j=k+1,nvales
          n = nlm(k)/100
          l = nlm(k)/10 - 10 * n
          n2 = nlm(j)/100
          l2 = nlm(j)/10 - 10 * n2
          
          if (l2.eq.l) then
            isemi(j)=1
          endif
        enddo
      enddo

      do k = 1,nvales
        lghost(k) = 0
        if (ibd(k).eq.0) lghost(k)=0
        
        n = nlm(k)/100
        l = nlm(k)/10 - 10 * n

        if (k.eq.iloc.or.ibd(k).eq.0) goto 911

        write(7,9110) n,xc(l)

        isoft = 1                                                  
        n = l + 1
        if (isemi(k).eq.1) n=n+1
        ee = en(k)
        call schsl(iloc,n,l,ee,maxim,rvloc,p,ig)
        el0(k) = ee
        if (iterm.eq.1) el0(k) = 0.0

        n = l + 2
        if (isemi(k).eq.1) n=n+1
        call schsl(iloc,n,l,ee,maxim,rvloc,p,ig)
        el1(k) = ee
        if (iterm.eq.1) el1(k) = 0.0

        do i = 1,np
          fr(i) = phipsref(i,k)**2*(rvps(i,k)-rvloc(i))/r(i)
        enddo

        tov = (float(l+l+2))
        call radin(r,fr,0,maxim,h,tov)
        xden = tov

        do i = 1,np
          fr(i) = fr(i) * (rvps(i,k)-rvloc(i))/r(i)
        enddo

        tov = (float(l+l+2))
        call radin(r,fr,0,maxim,h,tov)
        xnum = tov

        if (xden.ne.0) then
           elkb = xnum/xden
           write(7,9111) elkb,sqrt(xnum),sqrt(xnum)/elkb
           write(7,9103) el0(k),el1(k),en(k)           

           if (elkb.gt.0.0.and.en(k).gt.el1(k)) then
              write(7,9020) n,xc(l),en(k),el1(k)
              lghost(k) = 1
              ighost = 1
           endif
           
           if (elkb.lt.0.0.and.en(k).gt.el0(k)) then
              write(7,9020) n,xc(l),en(k),el0(k)
              lghost(k) = 1
              ighost = 1
           endif
           
        else
           
           write (7,*)
     $          '...seems like this angular momentum is the local part!'
        endif

 911    continue
      enddo
 
c     *************************************************************************
c     report result
c     *************************************************************************

      write(7,*) '------------------------------'
      write(7,*)
      if (ighost.eq.0) then
        write (7,*) 'No ghosts present for designed local potential!'
      endif

 9300 format(1x,'No ghosts for local potential: ',i1,a1)
 9101 format(1x,'Local state: ',i1,a1)
 9110 format(1x,'Test  state: ',i1,a1)
 9111 format(1x,'KB energy : ',f10.6,2x
     $     ,'KB strength: ',f10.6,2x,'KB cosine: ',f10.6,2x)
 9103 format(1x,'el0       : ',f10.6,2x,
     $     'el1        : ',f10.6,2x,'eig      : ',f10.6)
 9020 format(1x,'***GHOST*** : ',i1,a1,
     $     1x,f10.6,2x,'Should be lower than',2x,f10.6)

      
      write(7,*)

      return
      end
