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
      subroutine ghost
      implicit double precision (a-h,o-z)
      
#include "fortdim.h"

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
      common /totpot/ rvcore(npdm,n0),rvps(npdm,n0),rvcoul(npdm)

      common /rscore/ rscore(npdm),rdd(npdm),rddd(npdm),rscoretot(npdm),
     $     rsval(npdm)

      common /aorb/ ncore,nval,norb,nlm(n0),
     $     no(n0),lo(n0),nmax(n0),maxim
      common /adat/ wnl(n0),en(n0),so(n0),xion
      common /aval/ rcall(n0),rvap(n0),rnorm(n0),ibd(n0),etot

      common /wfn/ rnl(npdm,n0)
      common /psout/ wsumt(n0), sumc(n0), npsghost(n0)
c -------------------------------------------------------------------------

c -------------------------------------------------------------------------
c     Internal (Fortran only) common blocks
c -------------------------------------------------------------------------
      common /el/ el0(n0),el1(n0)
      common /nloc/ nloc
      common /iterm/ iterm
c -------------------------------------------------------------------------
      
      dimension fr(npdm),p(npdm),g(npdm)
      character*1 xc(0:3)

      xc(0)='s'
      xc(1)='p'
      xc(2)='d'
      xc(3)='f'
      
      ig=1
      isoft=1
      inl=0

      do i = 1,np
         g(i) = 0
      enddo

      write(7,*) 
      write(7,*) '---Semilocal ghost testing---'
      
      do nloc = 1,nval
         if (ibd(nloc).eq.0) goto 200
         npsghost(nloc) = 0

         nnn = nlm(nloc)/100
         lll = (nlm(nloc) - nnn * 100)/10

         write(7,9101) nnn,xc(lll)

         do i = 1,nval

            if (i.eq.nloc) goto 100

            write(7,*)
            
            nn = nlm(i)/100
            ll = (nlm(i) - nn * 100)/10
            im = maxim
            write(7,9110) nn,xc(ll)
            
            if (ibd(i).eq.0) then
               write(7,9220) 
               goto 100
            endif

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

            write(7,9111) elkb,sqrt(xnum),sqrt(xnum)/elkb

            ee = en(i)                                                 
            n = ll + 1            
            call schsl(i,n,ll,ee,im,rvps(1,nloc),p,ig,ibd(nloc))
            el0(i) = ee
            iterm0=iterm
            if (iterm0.eq.1) then
               el0(i) = 0.0
               write(7,9031)
            endif

            n = ll + 2
            call schsl(i,n,ll,ee,im,rvps(1,nloc),p,ig,ibd(nloc))
            el1(i) = ee
            iterm1=iterm
            if (iterm1.eq.1) then
               el1(i) = 0.0
               write(7,9032)
            endif
            
            write(7,9103) el0(i),el1(i),en(i)

c     Positive KB energy: e_loc0  < e_nl  < e_loc1
c     Will say no ghost if iterm0 = 0 and e_loc0 < e_nl 
c     and e_nl < e_loc1

            if (elkb.gt.0.0) then
               
c     if (iterm0.eq.0) then
               if (en(i).lt.el1(i).and.en(i).gt.el0(i)) then
                  write(7,9200) 
               else
                  write(7,9020) nn,xc(ll),en(i),el1(i)
                  npsghost(nloc) = 1
               endif
c     else
c     write(7,9300) 
c     if (npsghost(nloc).eq.0) npsghost(nloc) = -1
c     endif
               
            else
c     if (iterm0.eq.0) then
               if (en(i).lt.el0(i)) then
                  write(7,9210) 
               else
                  write(7,9020) nn,xc(ll),en(i),el0(i)
                  npsghost(nloc) = 1
               endif
c     else
c     write(7,9310) 
c     if (npsghost(nloc).eq.0) npsghost(nloc) = -1
            endif
c     endif
 100        continue
         enddo
         write(7,*)
         if (npsghost(nloc).eq.0) write(7,9400) nnn,xc(lll)
         if (npsghost(nloc).eq.1) write(7,9410) nnn,xc(lll)
         if (npsghost(nloc).eq.-1) write(7,9420) nnn,xc(lll)

         write(7,*) '------------------------------'
         write(7,*)
 200     continue
      enddo

 9400 format(1x,'No ghosts for local potential: ',i1,a1)
 9410 format(1x,' !WARNING! Ghosts for local potential: ',i1,a1)
 9420 format(1x,"  !NOTE! Ghost testing not conclusive for ",
     $     "local potential:  ",i1,a1)
 9101 format(1x,'Local state: ',i1,a1)
 9110 format(1x,'Test  state: ',i1,a1)
 9111 format(1x,'KB energy : ',f10.6,2x
     $     ,'KB strength: ',f10.6,2x,'KB cosine: ',f10.6,2x)
 9103 format(1x,'el0       : ',f10.6,2x,
     $     'el1        : ',f10.6,2x,'eig      : ',f10.6)
 9020 format(1x,'    !GHOST! : ',i1,a1,
     $     1x,f10.6,2x,'Should be lower than',2x,f10.6)
 9031 format(1x,'   !WARNING! No solution for ground state',
     $     ' of local potential.  Setting e=0.0')
 9032 format(1x,'   !WARNING! No solution for 1st excited state',
     $     ' of local potential.  Setting e=0.0')
 9220 format(1x,"   !WARNING! No ghost test for non-local",
     $     " unbound state")
 9200 format(1x,"No ghosts!  Ekb>0  and el0 < eig < el1")
 9210 format(1x,"No ghosts!  Ekb<0  and eig < el0")

 9300 format(1x,"   !WARNING! Ghost test unclear: ",
     $     "Must have two bound states for local pot when Ekb>0")
 9310 format(1x,"   !WARNING! Ghost test unclear: ",
     $     "Must have bound ground state for local pot when Ekb<0")
      

      return
      end


