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
      subroutine descreen(ixc,exccut)
      
      implicit double precision(a-h,o-z)
      
#include "fortdim.h"

c -------------------------------------------------------------------------
c     External (shared between C and Fortran) common blocks
c -------------------------------------------------------------------------
      common /grid/ h,r1,z,r(npdm),np
      common /rpcc/ rpcc,rpccz
      common /totpot/ rvcore(npdm,n0),rvps(npdm,n0),rvcoul(npdm)
      common /partpot/ rvh(npdm),rvxc(npdm),rexc(npdm),
     $     rvx(npdm),rex(npdm),rvc(npdm),rec(npdm),rvdft(npdm)
      common /rscore/ rscore(npdm),rdd(npdm),rddd(npdm),rscoretot(npdm),
     $     rsval(npdm)
      common /wfn/ rnl(npdm,n0)
      common /wavrc/ wavrc, slope, curvrc, indrc
      common /consts/ etol,vtol,maxit,isoft
      common /aorb/ ncore,nval,norb,nlm(n0),
     $     no(n0),lo(n0),nmax(n0),maxim
      common /adat/ wnl(n0),en(n0),so(n0),xion
      common /aval/ rcall(n0),rvap(n0),rnorm(n0),ibd(n0),etot
c -------------------------------------------------------------------------

c -------------------------------------------------------------------------
c     Internal (Fortran only) common blocks
c -------------------------------------------------------------------------
      common /atom4/ v(npdm),rv(npdm),p(npdm),rsatom(npdm)
c     $     c1(npdm),rvh(npdm),rvxc(npdm),rexc(npdm)
      common /iterm/ iterm
c -------------------------------------------------------------------------

      dimension rvxcv(npdm),rexcv(npdm),rvxct(npdm)
      dimension rexct(npdm),rvhc(npdm),rvhv(npdm),g(npdm)

      character*1 xc(0:3)

      xc(0)='s'
      xc(1)='p'
      xc(2)='d'
      xc(3)='f'

      pi=acos(-1.0)                                                   
      pi4=4.0*pi                                                      

      do i = 1,np
        g(i) = 0
      enddo

      isoft=1
      zero = 0.0
      iatom = 7
      charge = dfloat(2 * 1)

      write(7,*) '------------------------------'
      write(7,*) 'Descreening potential'

      call radin(r,rsval,0,maxim,h,charge)
      write (7,9010) charge

      rsc=0.0
         
      if (rpcc.gt.1e-12) then
         call radin(r,rscore,0,maxim,h,rsc)      
         write (7,9011) rsc
      else
         write (7,9011) rsc        
      endif 

 9010 format(1x,'valence charge          : ',f10.6)
 9011 format(1x,'core    charge          : ',f10.6)

c     *************************************************************************
c     section 1:  input/output from hfs files                             
c     *************************************************************************
      do i=1,maxim
        rsatom(i)=rsval(i)+rscore(i)
      enddo

      irelxc=0
      call excorr (maxim,ixc,exccut,rsval,rvxcv,rexcv,rvx,rex,rvc,
     $    rec,irelxc)
c      call excorr (maxim,ixc,exccut,rsval,rvxcv,rexcv,irelxc)
      call hrtree (maxim,h,r,rsval,rvhv)
      call excorr (maxim,ixc,exccut,rsatom,rvxct,rexct,rvx,rex,rvc,
     $    rec,irelxc)
c      call excorr (maxim,ixc,exccut,rsval,rvxcv,rexcv,irelxc)
c      call hrtree (maxim,h,r,rsval,rvhv)
c      call excorr (maxim,ixc,exccut,rsatom,rvxct,rexct,irelxc)
c      call excorr (maxim,ixc,exccut,rsatom,rvxct,rexct,irelxc)
      call hrtree (maxim,h,r,rscore,rvhc)
      
c     compute ionic pseudopotential by subtracting off screening terms.     

      wv = 0.         
      do i=1,nval
        wv = wv + wnl(i)
      enddo
      wc = z - xion - wv  
      zv = z - wc   
      zv2 = zv+zv     
      xn2 = zv2 - xion - xion   
      maxp1 = maxim + 1

      do i=1,maxim
        if (rpcc.gt.1e-12) then
          rvcoul(i) = rvhv(i)+rvxct(i)      
        else
          rvcoul(i) = rvhv(i)+rvxcv(i)      
        endif
                       
        do lp=1,nval      
          rvcore(i,lp) = rvps(i,lp)-rvcoul(i)    
        enddo
c       amr         rvcore(i,4) = -zv2      
      enddo

c     *************************************************************************
c     treat tail region separately (where wavefunctions are zero).
c     *************************************************************************

      do i = maxp1,np       
        do lp = 1,nval
          rvcore(i,lp) = -zv2   
          rvps(i,lp) = -xion-xion         
        enddo
        rvcoul(i)=rvps(i,1)-rvcore(i,1)
      enddo
 
c     *************************************************************************
c     solve schroedinger equation for each state in turn in following loop, 
c     *************************************************************************

      write(7,*) 
      write(7,*) '----Solving the Schrodinger equation',
     $           ' for all states----'
      ig=0
      do m=1,nval

        ee = en(m)         
        inl = 1  
        isoft = 1
        ebef = ee
        n = nlm(m)/100     
        l = nlm(m)/10 - 10*n         
        nn = l + 1
        im = maxim

        call schsl(m,nn,l,ee,im,rvps(1,m),p,ig,ibd(m))
        
        if (iterm.eq.0) then
          write(7,9111) n,xc(l),en(m),ee
 9111     format(1x,'State: ',i1,a1,2x,'AE eigenvalue = ', f10.6,2x,
     $         'PS eigenvalue = ', f10.6)
        endif            
        if (im.gt.maxim) maxim = im
        if (abs(ebef-ee).gt.0.1.or.iterm.eq.1) then
          write (7,*) 'This wavefunction has a node.  Try a '
          write (7,*) 'different starting guess for coefficients'
          write (7,*) ' or a different set of constraints and '
          write (7,*) 'weights.'
          stop
        endif

c     NOTE: the units are Ryd

      enddo
 
      return
      end
