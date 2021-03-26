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
      subroutine nllogd(istate,lang)
      
      implicit double precision (a-h,o-z)
      
#include "fortdim.h"

c     -------------------------------------------------------------------------
c     External (shared between C and Fortran) common blocks
c     -------------------------------------------------------------------------
      common /grid/ h,r1,z,r(npdm),np
      common /box/ iboxstart(n0),iboxend(n0),boxheight(n0),numbox
      common /ilogder/ ilogder
      common /local/ nlghost(n0),iloc,idesign
      common /logarith/ rphas,elogmax,elogmin,dlwf(npl0,n0)
      common /nlpot2/ inl,indrc(n0),IB,ID,IG
      common /nlcore/ rvcorew(npdm)
      common /totpot/ rvcore(npdm,n0),rvps(npdm,n0),rvcoul(npdm)
      common /aval/ rcall(n0),rvap(n0),rnorm(n0),ibd(n0),etot
c     -------------------------------------------------------------------------
 
c     -------------------------------------------------------------------------
c     Internal (Fortran only) common blocks                
c     -------------------------------------------------------------------------
      common /nlpot1/ Flstar(npdm),phipsref(npdm,n0),phiguess(npdm,n0)
c     -------------------------------------------------------------------------

c     *************************************************************************
c     local variables
c     *************************************************************************

      dimension rvloc(npdm),g(npdm)
      dimension ader(npdm),phomder(npdm),pinhomdr(npdm),Flphider(npdm)
      dimension Flwder(npdm),Flxder(npdm),pder(npdm)
      dimension dl(n0)
      
      call applyaug(rvloc,rvps(1,iloc))

      do i=1,np
        Flstar(i) = (rvps(i,istate) - rvloc(i))
     $       *phipsref(i,istate)/r(i)
      enddo
      
      do i=1,np
        g(i)=-flstar(i)*sqrt(r(i))*r(i)
      enddo

      if (ilogder.eq.1) then
         if (rphas.gt.rcall(istate)) then
            nr = 4
            ider=0
            isoft=1
            dele = (elogmax-elogmin)/float(npl0-1)
            do i = 1,np
               ader(i) = 0.0
               phomder(i) = 0.0
               pinhomdr(i) = 0.0
            enddo
            do j = 1,npl0
               e=elogmin + dele *float(j-1)
               call logder(e,z,lang,isoft,rvloc,ader,r,np,h,rphas,dl,
     $              npl0,nr,phomder,wfld,ider,0,0,0)
               call logder(e,z,lang,isoft,rvloc,g,r,np,h,rphas,dl,
     $              npl0,nr,pinhomdr,wfld,ider,0,0,0)
               do i = 1, np
                  Flphider(i) = Flstar(i) * phipsref(i,istate)
                  Flwder(i) = Flstar(i) * phomder(i) * sqrt(r(i))
                  Flxder(i) = Flstar(i) * pinhomdr(i) * sqrt(r(i))
               enddo
               flinvder = 2 * lang + 2
               call radin (r, Flphider, 0, np, h, flinvder)
               if (abs(flinvder).le.1e-10) then
                  flder = 0.
               else
                  flder = 1/flinvder
               endif
               xtilder = 2 * lang + 2
               call radin (r, Flxder, 0, np, h, xtilder)
               xtilder = 1+flder*xtilder
               wtilder = 2 * lang + 2
               call radin (r, Flwder, 0, np, h, wtilder)
               wtilder=flder*wtilder
               if (inl.eq.0) then wtilder = 0
               do i = 1, np
                  pder(i) = (phomder(i) *xtilder - pinhomdr(i)*wtilder)
               enddo
               call cubint(istate,pder,rvloc,g,e,lang,dlwfinhm)
               dlwf(j,istate) = dlwfinhm
               write(600,*) istate,j,dlwf(j,istate)
            enddo
         else
            write(7,9333) rphas,istate,rcall(istate)
            write(7,*) 'Log deriv not computed!'
         endif
      endif
 9333 format(1x,'Log deriv radius < rc: r_log = ', f10.6,2x,
     $     'state,rc',i3,f10.6)
      
      return
      end

      subroutine applyaug(rvl,rvs)
      implicit double precision (a-h,o-z)

#include "fortdim.h"

c -------------------------------------------------------------------------
c     External (shared between C and Fortran) common blocks
c -------------------------------------------------------------------------
      common /grid/ h,r1,z,r(npdm),np
      common /box/ iboxstart(n0),iboxend(n0),boxheight(n0),numbox
      common /nlpot2/ inl,indrc(n0),IB,ID,IG

      dimension rvl(npdm),rvs(npdm)

c      IA     = iboxstart(1)
c      IB     = iboxend(1)
c      ADEPTH = boxheight(1)
c      IC     = iboxstart(2)
c      ID     = iboxend(2)
c      BDEPTH = boxheight(2)
c      IE     = iboxstart(3)
c      IG     = iboxend(3)
c      CDEPTH = boxheight(3)
c      IH     = iboxstart(4)
c      IJ     = iboxend(4)
c      DDEPTH = boxheight(4)
      
      do i = 1,np
        rvl(i) = rvs(i)
      enddo

c      do i = IA,IB
c         rvl(i) = rvl(i) + ADEPTH * r(i)
c      enddo
c      do i = IC,ID
c         rvl(i) = rvl(i) + BDEPTH * r(i)
c      enddo
c      do i = IE,IG
c         rvl(i) = rvl(i) + CDEPTH * r(i)
c      enddo
c      do i = IH,IJ
c         rvl(i) = rvl(i) + DDEPTH * r(i)
c      enddo

c     Replacement for the blocks commented out above; allows for more
c     than four boxes. src/main/parameter.c updated accordingly.
c     BST 09/08/10
      do i = 1,min(numbox,n0)
         do j = iboxstart(i),iboxend(i)
            rvl(j) = rvl(j) + boxheight(i) * r(j)
         enddo
      enddo
      
      return
      end

