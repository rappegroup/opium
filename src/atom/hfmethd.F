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
      subroutine hfmethd(m,ipsp,p,zeff,azz,fail,iprint)
      implicit double precision (a-h,o-z)
      
#include "fortdim.h"

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
c -------------------------------------------------------------------------
      common /totpot/ rvcore(npdm,n0),rvps(npdm,n0),rvcoul(npdm)
      common /iterm/ iterm
      common /elim/ emore,eless,elim
      common /ipos/ ipos(n0),itermcount
      common /itoo/ ietoolo,ietoohi
      common /iskip/ iskip
      common /ameth/ meth(n0),imeth,nfails

c -------------------------------------------------------------------------
      dimension p(npdm)
      logical v2,first,ilast,fail

      fail=.false.
      first=.true.

      if (ipsp.eq.0) then
         em=0.d0
         eu=((z-min(0.5*s(m),2.d0*s(m)))/no(m))**2
      else
         em=1.d-4
         eu=e(m,m)+1.d0
      endif

      fu=eu
      mk=0
      
 17   call hfint(m,ipsp,p,pp,zeff,ed,ntp,azd,deltae,first,eu,em,iprint)

      if (imeth.eq.3) goto 111

c     count the number of nodes
      node=no(m)-lo(m)-1

      diff1=p(40)-p(39)
      mm=nmax(m)
      dm=0.d0
      do j=40,mm
         dm=max(dm,dabs(p(j)))
      enddo
      ncc=0
      sign=0.d0
      do j=40,mm
         diff2=p(j+1)-p(j)
         if ( diff2*diff1.gt.0.d0.or.diff1.eq.0.d0) goto 22
         if (dabs(p(j))/dm .lt. 0.05d0) goto 22
         if (sign.ne.0.d0) goto 44
         mm=j
         goto 331
 44      if (p(j)*sign.gt.0.d0) goto 22
         ncc=ncc+1
 331     sign=p(j)
 22      diff1=diff2
      enddo

      if (iprint.gt.1) write (7,99) no(m),lo(m),ncc,r(mm),r(ntp),
     $     p(mm),ed,eu,em,deltae
 99   format(2x,2i3,' nc =',i3,' mn =',f10.6,' nj =',f10.6,' pde(mn) =',
     :     d10.2,' ed =',d10.2,' eu =',d10.2,' em =',d10.2,
     :     ' deltae =',d10.2)

      
      if (iabs(ncc-node).eq.1.and.dabs(deltae/ed).gt.0.02d0) goto 46

      if ((ncc-node).gt.0) then
         if (p(mm).lt.0.d0) then
            goto 11
         endif

         del=1.d0-em/ed
         em=ed
         if (del .lt. 0.05d0) fm = fm*((lo(m)+1.d0+ncc)/no(m))**2.5
         if (del .ge. 0.05d0) fm = ed*((lo(m)+1.d0+ncc)/no(m))**2.5
         if (fm .gt. eu) fm = 0.5*(eu + em)
         if (dabs(fm - ed) .lt. 0.00001d0) then
            go to 27
         endif
         ed = fm
         go to 33
      elseif((ncc-node).lt.0) then
         del=1.d0-ed/eu
         eu=ed
         if (del.lt.0.05d0) fu=fu*((lo(m)+1+ncc)/no(m))**2.5
         if (del.ge.0.05d0) fu=ed*((lo(m)+1+ncc)/no(m))**2.5
         if (fu.lt.em) fu=(eu+em)/2.d0
         if (dabs(fu-ed).lt.1e-4) goto 27
         ed=fu
         goto 33
      else
         v2=dabs(deltae).lt.1.d-3.or.dabs(deltae)/ed.lt.1d-5
         if (p(mm).lt.0.d0.and. .not.v2) goto 46
         if (p(mm).gt.0.d0) goto 51
         do j=1,np
            p(j)=-p(j)
         enddo
         pp=-2.d0-pp
 51      continue
      endif
      goto 111
         
         
 11   continue
      edp=ed
      ed=ed*((lo(m)+1+ncc)/float(no(m)))**2.5
      if (ed.ge.eu) ed=(eu+edp)/2.d0
      if (ed.le.em) ed=(em+edp)/2.d0

      
 33   mk=mk+1
      if (mk.gt.2*no(m).or.eu-em.lt.
     $     float(no(m))**(0.3333333)) goto 27
      first=.false.
      goto 17


 46   continue
      ed=ed-deltae

      if (ed.ge.em.and.ed.le.eu) goto 33
      edp=ed
      if (ncc-node.ne.0) ed=(ed+deltae)*((lo(m)+1+ncc)
     $     /float(no(m)))**2.5

      if (ed.ge.em.and.ed.le.eu) goto 33

      ed=edp+deltae*2.d0

      if (ed.ge.em.and.ed.le.eu) goto 33

      ed=ed-deltae
      deltae=deltae/2.d0
      goto 46

 27   continue
      if (iprint.gt.0) write (7,28) imeth,no(m),lo(m),ncc,mm,ed,em,eu
 28   format(1x,"method",i2,2x,"unable to solve equation for electron",
     $     2i2/10x,5hnc = ,i3,3x,5hnj = ,i3,3x,5hed = ,f10.6,3x,5hel = ,
     $     f10.6,3x,5heu = ,f10.6)
      call orthog(ilast)
      call grange

      fail=.true.
      nfails=nfails+1
      if (iprint.gt.0) write(7,*) "HF solver failed:  #fails=",nfails
      return
      
 111  continue
      
      azz=azd*(pp+1.d0)
      e(m,m)=ed
      en(m)=e(m,m)
      
      return
      end

