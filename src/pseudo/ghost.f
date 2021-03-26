      subroutine ghost
      implicit double precision (a-h,o-z)
      
#include "PARAMOPT"

c This program finds the ground state e0l and first excited state e1l
c of the local potential.
c For nonlocal components, we calculate elkb.  The ghost theorem of
c Xavier Gonze states the following:
c For elkb>0, and eat>e1l the potential has a ghost below eat.
c For elkb<0, and eat>el0 the potential has a ghost below eat.
c Here eat is the reference eigenvalue of the nonlocal angular momentum.

c -------------------------------------------------------------------------
c     External (shared between C and Fortran) common blocks
c -------------------------------------------------------------------------
      common /grid/ h,r1,z,r(npdm),np
      common /atomic/ xion,rnl(npdm,n0),nlm(n0),wnl(n0),en(n0),norb
      common /totpot/ rvcore(npdm,n0),rvps(npdm,n0),rvcoul(npdm)
      common /np/ ncores,nvales
      common /ibound/ ibd(n0)
      common /nmax/ nmax(nvale0),maxim
      common /convrpt/ wsumt(nvale0), sumc(nvale0), lghost(nvale0)
c -------------------------------------------------------------------------

c -------------------------------------------------------------------------
c     Internal (Fortran only) common blocks
c -------------------------------------------------------------------------
      common /el/ el0(nvale0),el1(nvale0)
      common /nloc/ nloc
      common /iterm/ iterm
c -------------------------------------------------------------------------
      
      dimension fr(npdm),p(npdm),g(npdm)
      character*1 xc(0:3)

      xc(0)='s'
      xc(1)='p'
      xc(2)='d'
      xc(3)='f'
      
      ig=0
      isoft=1
      inl=0

      do i = 1,np
         g(i) = 0
      enddo

      write(7,*) 
      write(7,*) '---Semilocal ghost testing---'
      
      do nloc = 1,nvales
         if (ibd(nloc).eq.0) goto 200
         lghost(nloc) = 0

         nnn = nlm(nloc)/100
         lll = (nlm(nloc) - nnn * 100)/10

         write(7,9101) nnn,xc(lll)

         do i = 1,nvales
            if (i.eq.nloc.or.ibd(nloc).eq.0.or.ibd(i).eq.0) goto 100

            write(7,*)
            
            nn = nlm(i)/100
            ll = (nlm(i) - nn * 100)/10
            im = maxim

            do j = 1,np
               fr(j) = rnl(j,i)**2*(rvps(j,i)-rvps(j,nloc))/r(j)
            enddo
            
            tov = (float(ll+ll+2))
            call radin(r,fr,0,im,h,tov)
            xden = tov
            
            do j = 1,np
               fr(j) = fr(j) * (rvps(j,i)-rvps(j,nloc))/r(j)
            enddo

            tov = (float(ll+ll+2))
            call radin(r,fr,0,im,h,tov)
            xnum = tov

            elkb = xnum/xden

            write(7,9110) nn,xc(ll)
            write(7,9111) elkb,sqrt(xnum),sqrt(xnum)/elkb

            ee = en(i)                                                 
            n = ll + 1            
            call schsl(nloc,n,ll,ee,im,rvps(1,nloc),p,ig)
            el0(i) = ee
            if (iterm.eq.1) el0(i) = 0.0
            n = ll + 2
            call schsl(nloc,n,ll,ee,im,rvps(1,nloc),p,ig)
            el1(i) = ee
            if (iterm.eq.1) el1(i) = 0.0

            write(7,9103) el0(i),el1(i),en(i)

            if (elkb.gt.0.0.and.en(i).gt.el1(i)) then
               write(7,9020) nn,xc(ll),en(i),el1(i)
               lghost(nloc) = 1
            endif

            if (elkb.lt.0.0.and.en(i).gt.el0(i)) then
               write(7,9020) nn,xc(ll),en(i),el0(i)
               lghost(nloc) = 1
            endif
 100        continue
         enddo
                  
         if (lghost(nloc).eq.0) then
            write (7,9300) nnn,xc(lll)
         endif
            
         write(7,*) '------------------------------'
         write(7,*)
 200     continue
      enddo

 9300 format(1x,'No ghosts for local potential: ',i1,a1)
 9101 format(1x,'Local state: ',i1,a1)
 9110 format(1x,'Test  state: ',i1,a1)
 9111 format(1x,'KB energy : ',f10.6,2x
     $     ,'KB strength: ',f10.6,2x,'KB cosine: ',f10.6,2x)
 9103 format(1x,'el0       : ',f10.6,2x,
     $     'el1        : ',f10.6,2x,'eig      : ',f10.6)
 9020 format(1x,'***GHOST*** : ',i1,a1,
     $     1x,f10.6,2x,'Should be lower than',2x,f10.6)

      return
      end


