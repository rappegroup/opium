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
      subroutine ghostnl(rvloc)
           
c     *************************************************************************
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
      
#include "fortdim.h"
 
c     -------------------------------------------------------------------------
c     External (shared between C and Fortran) common blocks
c     -------------------------------------------------------------------------

      common /aorb/ ncore,nval,norb,nlm(n0),
     $     no(n0),lo(n0),nmax(n0),maxim
      common /adat/ wnl(n0),en(n0),so(n0),xion
      common /aval/ rcall(n0),rvap(n0),rnorm(n0),ibd(n0),etot

      common /grid/ h,r1,z,r(npdm),np
      common /box/ iboxstart(n0),iboxend(n0),boxheight(n0),numbox
      common /local/ nlghost(n0),iloc,idesign
c     -------------------------------------------------------------------------
 
c     -------------------------------------------------------------------------
c     Internal (Fortran only) common blocks                
c     -------------------------------------------------------------------------
      common /iterm/ iterm
      common /nlpot1/ Flstar(npdm),phipsref(npdm,n0),phiguess(npdm,n0)
      common /totpot/ rvcore(npdm,n0),rvps(npdm,n0),rvcoul(npdm)
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
c      if (ibd(iloc).eq.0) then
c         write(7,*) " !WARNING! Non-local ghost testing is not reliable"
c         write(7,*) " !WARNING!   when the local pot is unbound!"
c         write(7,*) "Please change your choice of local potential, = "
c     $        , nlm(iloc)
c         stop
c      endif

      do i = 1,np
        gghost(i) = 0
        flstar(i) = 0
      enddo
      
      ighost = 0
      ig=1

      nnn = nlm(iloc)/100
      lll = (nlm(iloc) - nnn * 100)/10
      
      write(7,9101) nnn,xc(lll)
      
      do k=1,nval
        isemi(k)=0
      enddo
      do k=1,nval
        do j=k+1,nval
          n = nlm(k)/100
          l = nlm(k)/10 - 10 * n
          n2 = nlm(j)/100
          l2 = nlm(j)/10 - 10 * n2
          
          if (l2.eq.l.and.n2.ne.n) then
            isemi(j)=1
          endif
        enddo
      enddo

      do k = 1,nval
        nlghost(k) = 0

        n = nlm(k)/100
        l = nlm(k)/10 - 10 * n

c        if (k.eq.iloc.or.ibd(k).eq.0) goto 911
c        if (k.eq.iloc) goto 911

        write(7,*)        
        write(7,9110) n,xc(l)

        if (ibd(k).eq.0) then
           write(7,9220) 
           goto 911
        endif

        isoft = 1                                                  
        n = l + 1
        if (isemi(k).eq.1) n=n+1
        ee = en(k)
        call schsl(k,n,l,ee,maxim,rvloc,p,ig,ibd(k))

        el0(k) = ee
        iterm0=iterm
        if (iterm0.eq.1) then
           el0(k) = 0.0
           write(7,9031)
        endif

        n = l + 2
        if (isemi(k).eq.1) n=n+1
        call schsl(k,n,l,ee,maxim,rvloc,p,ig,ibd(k))

        el1(k) = ee
        iterm1=iterm
        if (iterm1.eq.1) then
           el1(k) = 0.0
           write(7,9032)
        endif

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
           
           if (elkb.gt.0.0) then

c     if (iterm1.eq.0.and.iterm0.eq.0) then
              
              if (en(k).lt.el1(k).and.en(k).gt.el0(k)) then
                 write(7,9200) 
              else
                 write(7,9020) n,xc(l),en(k),el1(k)
                 nlghost(k) = 1
                 ighost=1
              endif

c     else
c     write(7,9300) 
c     if (nlghost(k).eq.0) nlghost(k) = -1
c     if (ighost.eq.0) ighost=-1
c     endif
           else

c     if (iterm0.eq.0) then
              
              if (en(k).lt.el0(k)) then
                 write(7,9210) 
              else
                 write(7,9020) n,xc(l),en(k),el0(k)
                 nlghost(k) = 1
                 ighost=1
              endif

c     else
c     write(7,9310) 
c     if (nlghost(k).eq.0) nlghost(k) = -1
c     if (ighost.eq.0) ighost=-1

           endif

c     endif
           
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
         if (numbox.gt.0) then
            write (7,*) 'No ghosts present for designed local pot'
         else
            write (7,*) 'No ghosts present for local potential'
         endif
      else
         if(ighost.eq.-1) then
            write (7,*) "  !NOTE! Ghost testing not conclusive"
         else
            write (7,9555) 
         endif
      endif

 9555 format(1x,"!!ERROR!! There are ghosts!  Change one or more",
     $     " of the following...",/,
     $     " 1) the choice of local potential",/,
     $     " 2) the cut-off radii",/,
     $     " 3) the reference configuration")

 9400 format(1x,'No ghosts for local potential: ',i1,a1)
 9101 format(1x,'Local state: ',i1,a1)
 9110 format(1x,'Test  state: ',i1,a1)
 9111 format(1x,'KB energy : ',f10.6,2x
     $     ,'KB strength: ',f10.6,2x,'KB cosine: ',f10.6,2x)
 9103 format(1x,'el0       : ',f10.6,2x,
     $     'el1        : ',f10.6,2x,'eig      : ',f10.6)
 9020 format(1x,'    !GHOST! : ',i1,a1,
     $     1x,f10.6,2x,'Should be lower than',2x,f10.6)
 9031 format(1x,'!WARNING! No solution for ground state',
     $     ' of local potential.  Setting e=0.0')
 9032 format(1x,'!WARNING! No solution for 1st excited state',
     $     ' of local potential.  Setting e=0.0')
 9220 format(1x,"!NOTE! No ghost test for non-local",
     $     " unbound state")
 9200 format(1x,"No ghosts!  Ekb>0  and el0 < eig < el1")
 9210 format(1x,"No ghosts!  Ekb<0  and eig < el0")

 9300 format(1x,"!WARNING! Ghost test unclear: ",
     $     "Must have two bound states for local pot when Ekb>0")
 9310 format(1x,"!WARNING! Ghost test unclear: ",
     $     "Must have bound state for local pot for ghost test")
 
     
      write(7,*)

      return
      end
