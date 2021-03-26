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
      subroutine hfde(m,ipsp,zeff,dpm,iprint)
      implicit double precision (a-h,o-z)
#include "fortdim.h"

c -------------------------------------------------------------------------
c     External (shared between C and Fortran) common blocks
c -------------------------------------------------------------------------
      common /aorb/ ncore,nval,norb,nlm(n0),
     $     no(n0),lo(n0),nmax(n0),maxim
      common /adat/ wnl(n0),en(n0),so(n0),xion
      common /aval/ rcall(n0),rvap(n0),rnorm(n0),ibd(n0),etot
      common /wfn/ rnl(npdm,n0)
      common /consts/ etol,vtol,maxit,isoft
      common /grid/ h,r1,z,r(npdm),np
      common /edi/ e(n0,n0)
      common /ahf/ s(n0),az(n0),acc(n0)
      common /ipr/ ipr

c -------------------------------------------------------------------------
c     Internal (Fortran only) common blocks                
c -------------------------------------------------------------------------
      common /totpot/ rvcore(npdm,n0),rvps(npdm,n0),rvcoul(npdm)
      common /iterm/ iterm
      common /elim/ emore,eless,elim
      common /ipos/ ipos(n0),itermcount
      common /itoo/ ietoolo,ietoohi
      common /iskip/ iskip
      common /ameth/ meth(n0),imeth

      dimension p(npdm),dpm(n0),rsatom(npdm)
      
      logical fail,ilast

      ed1=0.d0
      ed2=e(m,m)
      imeth=max(1,meth(m))
      if (norb.eq.1) then
         meth(m)=2
         imeth=2
      endif

      call hfmethd(m,ipsp,p,zeff,azz,fail,iprint)

      if (fail) goto 25

 12   continue

      pnn=dsqrt(quad(m,nmax(m),p,p))
      do k=1,nmax(m)
         p(k)=p(k)/pnn
      enddo

c      write(7,*) 'azz 1',azz,pnn
      azz=azz/pnn

      if (imeth.eq.1) then
         if (abs(1.d0-ed2/e(m,m)).lt.0.005d0.and.
     $        dmax1(dabs(1.d0-pnn),dabs(1.d0/pnn - 1.d0)).
     $        gt.0.20d0) then
            meth(m)=2
            imeth=2
            goto 25
         endif
      endif

      if (ipr.ne.m) then
         acc(m)=0.75*acc(m)
      else
         ed2=ed2-e(m,m)
         if (ed1*ed2.gt.0) then
            acc(m)=0.75*acc(m)
         else
            acc(m)=(1.d0 + 3.d0*acc(m))/4.d0
         endif
      endif
      c=acc(m)
      cd=1.d0-c

c      write(7,*) 'ACC', m,acc(m),c,ipr

      if (acc(m).gt.0.99999) acc(m)=0.50

      dp=0.d0
      do j=1,nmax(m)
c      do j=1,np
         diff=rnl(j,m)-p(j)
         dp=max(dp,dabs(diff)*sqrt(r(j)))
         rnl(j,m)=p(j)+c*diff
      enddo
      do i=nmax(m)+1,np
         rnl(i,m)=0.d0
      enddo

c      write(7,*) 'before de adj', m,az(m),azz

      az(m)=cd*azz+c*az(m)
      azz=az(m)

      if (wnl(m).lt.1e-6) then
         dpw=dp/0.001
      else
         dpw=dp/sqrt(wnl(m))
      endif

      do j=1,norb
         if (abs(e(m,j)).ne.0.d0.and.m.ne.j) then
            if (dpm(j).lt.dsqrt(wnl(j))*dpw) then
               nnp=max(nmax(j),nmax(m))
               rx=quadr(m,j,0)
c               write(7,63) nlm(j),nlm(m),rx
               do jj=1,nnp
                  rnl(jj,m)=rnl(jj,m)-rx*rnl(jj,j)
               enddo
               azz=azz-rx*az(j)
c               write(7,*) 'orth', j,azz,az(j),rx
            endif
         endif
      enddo

 63   format(6x,'<',i3,'|',i3,'>=',1pd8.1)      

      if (abs(rx).gt.0) then
         rx=dsqrt(quadr(m,m,0))
         
         do jj=1,np
            rnl(jj,m)=rnl(jj,m)/rx
         enddo
         azz=azz/rx
      endif

      jjset=np
 687  if (dabs(rnl(jjset,m)).lt.1e-15) then
         rnl(jjset,m)=0.d0
         jjset=jjset-1
         goto 687
      endif
      nmax(m)=jjset

      if (azz.gt.0.d0) az(m)=max(azz,az(m)/2.d0)

c      write(7,*) 'after de adj', m,az(m),azz

      en(m)=e(m,m)

      if (iprint.gt.1) 
     $     write (7,17) nlm(m),e(m,m),az(m),pnn,dp,r(nmax(m))
 17   format(20x,i3,2f15.7,f12.7,1pd10.2,1pd10.2)            

      maxim=max(maxim,nmax(m))
      dpm(m)=dp
      if(ipr.eq.m) then
         ed1=ed2
      else
         ed1=ed2-e(m,m) 
      endif
      ipr=m
      do i=1,nmax(m)
         rsatom(i)=rsatom(i)+wnl(m)*rnl(i,m)**2
      enddo

      return

 25   continue
      write(7,*) "Solve for HF wfn failed, orthogonalizing..."
      if (norb.eq.1) goto 27
      call orthog(ilast)
      call grange
      
 27   continue
      call hfmethd(m,ipsp,p,zeff,azz,fail,iprint)
      if ( fail ) then
         if (imeth.eq.2) return
         imeth=2
         goto 27
      else
         goto 12
      endif
      
      return
      end
