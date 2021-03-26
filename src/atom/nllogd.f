      subroutine nllogd(istate,nqn,lang,ei,imax)
      
      implicit double precision (a-h,o-z)
      
#include "PARAMHFS"

c     -------------------------------------------------------------------------
c     External (shared between C and Fortran) common blocks
c     -------------------------------------------------------------------------
      common /grid/ h,r1,z,r(npdm),np
      common /box/ iboxstart(10),iboxend(10),boxheight(10),numbox
      common /ilogder/ ilogder
      common /local/ iloc,idesign
      common /logarith/ rphas,elogmax,elogmin,dlwf(npl0,n0)
      common /nlpot2/ inl,indrc(n0),IB,ID,IG
      common /rcall/ rcall(n0)
      common /nlcore/ rvcorew(npdm)
      common /npm/ ncores,nvales,iskip,maxip
      common /totpot/ rvcore(npdm,n0),rvps(npdm,n0),rvcoul(npdm)
c     -------------------------------------------------------------------------
 
c     -------------------------------------------------------------------------
c     Internal (Fortran only) common blocks                
c     -------------------------------------------------------------------------
      common /nlpot1/ Flstar(npdm),phiguess(npdm,n0)
c     -------------------------------------------------------------------------

c     *************************************************************************
c     local variables
c     *************************************************************************

      dimension rvloc(npdm),g(npdm),phipsref(npdm,n0)
      dimension ader(npdm),phomder(npdm),pinhomdr(npdm),Flphider(npdm)
      dimension Flwder(npdm),Flxder(npdm),pder(npdm)
      dimension dl(npl0,4,3,2)
      
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

#include "PARAMHFS"

c -------------------------------------------------------------------------
c     External (shared between C and Fortran) common blocks
c -------------------------------------------------------------------------
      common /grid/ h,r1,z,r(npdm),np
      common /box/ iboxstart(10),iboxend(10),boxheight(10),numbox
      common /nlpot2/ inl,indrc(n0),IB,ID,IG

      dimension rvl(npdm),rvs(npdm)

      IA     = iboxstart(1)
      IB     = iboxend(1)
      ADEPTH = boxheight(1)
      IC     = iboxstart(2)
      ID     = iboxend(2)
      BDEPTH = boxheight(2)
      IE     = iboxstart(3)
      IG     = iboxend(3)
      CDEPTH = boxheight(3)
      IH     = iboxstart(4)
      IJ     = iboxend(4)
      DDEPTH = boxheight(4)
      
      do i = 1,np
        rvl(i) = rvs(i)
      enddo

      do i = IA,IB
         rvl(i) = rvl(i) + ADEPTH * r(i)
      enddo
      do i = IC,ID
         rvl(i) = rvl(i) + BDEPTH * r(i)
      enddo
      do i = IE,IG
         rvl(i) = rvl(i) + CDEPTH * r(i)
      enddo
      do i = IH,IJ
         rvl(i) = rvl(i) + DDEPTH * r(i)
      enddo
      
      return
      end
