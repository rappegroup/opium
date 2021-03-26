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
      subroutine schsl(iorb,nqn,lang,ei,imax,rv,p,igh,ibound)
      implicit double precision (a-h,o-z)
      
#include "fortdim.h"

c -------------------------------------------------------------------------
c     External (shared between C and Fortran) common blocks
c -------------------------------------------------------------------------
      common /consts/ etol,vtol,maxit,isoft
      common /grid/ h,r1,z,r(npdm),np
      common /nlpot2/ inl,indrc(n0),IB,ID,IG

      common /etrial/ etrial(n0)
c -------------------------------------------------------------------------
c     Internal (Fortran only) common blocks                
c -------------------------------------------------------------------------
      common /iterm/ iterm
      common /nlpot1/ Flstar(npdm),phipsref(npdm,n0),phiguess(npdm,n0)
      common /elim/ emore,eless,elim
      common /ipos/ ipos(n0),itermcount
      common /itoo/ ietoolo,ietoohi
      common /iskip/ iskip
      logical bound
c -------------------------------------------------------------------------
      dimension rv(npdm), p(npdm), phom(npdm)
      dimension pinhom(npdm),wnlfract(npdm)

      iterm = 0
      iout = 7
      iclose = 0
      itmax = 100
      npm1 = np - 1
      x1 = lang * (lang + 1)
      x2 = (lang + 0.5)**2
      node = nqn - lang
      ncross = 0
      ietoohi = 0
      ietoolo = 0
      emore = -1.0e-6
      eless = 0
      de = 0
      iskipcount=0

      h2 = h * h/12.
      do i = 1,np
         wnlfract(i)=0.0
      enddo

      bound=.true.
      if (ibound.eq.0) then
         bound=.false.
      endif
      
      if (inl.ne.0) call nlcont(lang,phiguess(1,iorb),phipsref(1,iorb),
     $     flstar,wnlfract)

      if (bound) call coarseadj(nqn,lang,ei,rv,wnlfract)

      if (iterm.eq.1) then
         itermcount=itermcount+1
         write (7,9003) nqn,lang
c         if (itermcount.gt.5.and.igh.eq.0) then 
c            write(7,*) '!NOTE! : setting',nqn,lang,' to be unbound!'            
c            ei=etrial(iorb)
c            ibound=0
c         endif
         goto 912
      endif

      do iter = 1,itmax

         iskip=0

         if (bound.and.iclose.eq.0) then
            call newmatch(nqn,lang,ei,imax,imat,rv,wnlfract,iskip)
            if (iskip.eq.1) then 
               iskipcount=iskipcount+1
               if (iskipcount.eq.itmax) then 
                  write (7,9200) itmax
                  ei=0.0
                  goto 912
               endif
               goto 1
            endif
         else if (.not.bound) then
            imat=np-200
            imax=np-200
         endif            

         call outward(ei,lang,rv,p,phom,pinhom,
     $        flstar,phipsref(1,iorb),imat,xt,wt)


         if (bound) call cntnodes(nqn,lang,ei,p,imat,node,iskip,iclose)
         if (iterm.eq.1) goto 912
         if (iskip.eq.1) goto 1

         if (bound) call inward(ei,nqn,lang,rv,p,phom,pinhom,
     $        flstar,phipsref(1,iorb),imat,xt,wt,de,ap,iclose,imax)

         if (.not.bound) goto 911
         val = abs(de/ei)

         if (val.le.etol) goto 911
         
 1303    ep = ei + de
         if (ep.lt.emore) goto 1304
         de = de * 0.5
         val = val * 0.5

         if (val.gt.etol) goto 1303

         write (iout,9000) nqn,lang
         iterm=1
         ipos(iorb)=ipos(iorb)+1
c         if (igh.eq.0.and.ipos(iorb).gt.4) then 
c            write(7,*) '!NOTE! : setting',nqn,lang,' to be unbound!'
c            ei=etrial(iorb)
c            ibound=0
c         endif
         goto 912

 1304    ei = ep
         if (val.le.0.05) iclose = 1
c         write (iout,9100) iter,iorb,nqn,lang,ei,de,imat

 1       continue         
      enddo
      if (iter.eq.itmax) iterm=1
 911  continue

c      write (iout,9101) iter,iorb,nqn,lang,ei,de,imat

      do i = 1,imax
        p(i) = p(i)*sqrt(r(i))*1/sqrt(ap)
      enddo

      do i = imax+1,np
        p(i) = 0
      enddo

      return
 912  continue
      iterm=1
      itermcount=itermcount+1
c      write (iout,9102) iter,iorb,nqn,lang,ei,de,imat

 9000 format(1x,'!WARNING! : schsl found a positive eigenvalue.',
     $     ' n,l=:', 2i4)
 9003 format(1x,'!WARNING! : schsl found a positive eigenvalue.',
     $     ' n,l=:', 3i4)

 9100 format('Aiter:',i5,1x,'iorb:',i5,1x,'n:',i5,1x,'l:',
     $     i3,1x,'ee:',2f10.6,3x,'im  ',i5)
 9101 format('Biter:',i5,1x,'iorb:',i5,1x,'n:',i5,1x,'l:',
     $     i3,1x,'ee:',2f10.6,3x,'im  ',i5)
 9102 format('Citer:',i5,1x,'iorb:',i5,1x,'n:',i5,1x,'l:',
     $     i3,1x,'ee:',2f10.6,3x,'im  ',i5)

 9200 format(1x,'!WARNING! : classical turning point at r=0',
     $     i5,' times in a row. Solve FAILED!')

      return
      end
      
      subroutine nlcont(lang,pg,pref,fs,wnlf)
      implicit double precision (a-h,o-z)
      
#include "fortdim.h"
      common /grid/ h,r1,z,r(npdm),np
      dimension psvph(npdm),phvph(npdm),fs(npdm),
     $     pg(npdm),pref(npdm),wnlf(npdm)

      do i = 1,np
         psvph(i) = pg(i)*fs(i)
         phvph(i) = pref(i)*fs(i)
      enddo

      effnlnum = 2 * lang + 2
      effnlden = 2 * lang + 2

      call radin (r, psvph, 0, np, h, effnlnum)
      call radin (r, phvph, 0, np, h, effnlden)

      if (abs(effnlden).le.1e-10) then 
         eff=0.0
      else
        eff=effnlnum/effnlden
         do i = 1, np
            if (abs(pg(i)).gt.1e-10) then
               wnlf(i)=eff*fs(i)/pg(i)
            else
               wnlf(i)=0.0
            endif
         enddo
      endif

      return
      end

      subroutine coarseadj(nqn,lang,ei,rv,wnlf)
      implicit double precision (a-h,o-z)
      
#include "fortdim.h"
      common /grid/ h,r1,z,r(npdm),np
      common /elim/ emore,eless,elim
      common /iterm/ iterm
      dimension wnlf(npdm),rv(npdm),v(npdm)

      npm1=np-1
      x1 = lang * (lang + 1)
      iexit=0

      do i = 1, np
         v(i) = (rv(i) + x1/r(i))/r(i) + wnlf(i)
      enddo

      if (lang.ne.0) then
         elim = v(1)
         do i = 2,np
            elim = min(v(i),elim)
         enddo
         if (elim.ge.0.0) then
            write (7,9000)
            iterm=1
         endif
      else
         elim = -z * z/(float(nqn * nqn))
      endif
      eless = elim

      if (ei.lt.eless) ei = eless * 0.5
      if (ei.gt.v(npm1)) ei = 2.0 * v(npm1)
      if (ei.gt.emore) ei = emore * 1.1

 9000 format (1x,"!WARNING! : schsl stopped due since ",
     $     "v+l(l+1)/r^2 became positive")

      return
      end

      subroutine newmatch(nqn,lang,ei,imax,imat,rv,wnlf,iskip)
      implicit double precision (a-h,o-z)
      
#include "fortdim.h"
      common /consts/ etol,vtol,maxit,isoft
      common /grid/ h,r1,z,r(npdm),np
      common /nlpot2/ inl,indrc(n0),IB,ID,IG
      common /iterm/ iterm
      dimension rv(npdm),v(npdm),wnlf(npdm)

      iskip=0
      x2 = (lang + 0.5)**2
      h2 = h * h/12.     

      do i = 1,np
         v(i) = ((rv(i) - ei * r(i)) * r(i) + x2) * h2 
     $        +  h2 * wnlf(i) * r(i) * r(i)
      enddo
      
      do imat = np,1,-1
         if (v(imat).le.0.0) goto 213
      enddo
 213  continue

      if (imat.le.2) then
         ei = ei * 0.5
c     if (ei.ge.-etol.or.ncross.gt.node) then
         if (ei.ge.-etol) then
            write (7,9010) nqn,lang
            iterm=1
         endif
         iskip=1
         return
      endif
      imat = min(np - 5,imat)

      if (inl.ne.0) then
         indmax=0
         do i = 1,10
            if (indmax.lt.indrc(i)) indmax = indrc(i)
         enddo
         if (imat.lt.indmax) then
            imat = indmax
         endif
         if (imat.lt.IB) then
            imat = IB
         endif
         if (imat.lt.ID) then
            imat = ID
         endif
         if (imat.lt.IG) then
            imat = IG
         endif
      endif

 9010 format(1x,'!WARNING! : schsl found this state to be unbound:',
     $     ' n,l=',2i4,' --nm')

      return
      end
      
      subroutine outward(ei,lang,rv,p,phom,pinhom,fs,pref,imat,xt,wt)
      implicit double precision (a-h,o-z)
      
#include "fortdim.h"
      common /consts/ etol,vtol,maxit,isoft
      common /grid/ h,r1,z,r(npdm),np
      common /nlpot2/ inl,indrc(n0),IB,ID,IG
      
      dimension a(npdm),phom(npdm),pinhom(npdm),rv(npdm)
      dimension dl(n0),dld(5),g(npdm),fs(npdm),pref(npdm)
      dimension flph(npdm),flw(npdm),flx(npdm),p(npdm)

      do i =1, np
         a(i) = 0
         phom(i) = 0
         pinhom(i) = 0
      enddo
      nr=4
      ider=0
      call logder(ei,z,lang,isoft,rv,a,r,np,h,r(imat),dl,npl0,nr,
     $     phom,wfld,ider,de,qp,dld)
      
      if (inl.ne.0) then

         do i=1,np
            g(i)=-fs(i)*sqrt(r(i))*r(i)
         enddo
         call logder(ei,z,lang,isoft,rv,g,r,np,h,r(imat),dl,npl0,nr,
     $        pinhom,wfld,ider,de,qp,dld)
      endif

      do i = 1, imat
         flph(i) = fs(i) * pref(i)
         flw(i) = fs(i) * phom(i) * sqrt(r(i))
         flx(i) = fs(i) * pinhom(i) * sqrt(r(i))
      enddo
      do i = imat+1, np
         flph(i) = 0
         flw(i) = 0
         flx(i) = 0
      enddo
      
      flinv = 2 * lang + 2
      call radin (r, flph, 0, np, h, flinv)
      
      if (abs(flinv).le.1e-10) then
         fl = 0.
      else
         fl = 1/flinv
      endif
      
      xt = 2 * lang + 2
      call radin (r, flx, 0, np, h, xt)
      
      xt = 1+fl*xt
      
      wt = 2 * lang + 2
      call radin (r, flw, 0, np, h, wt)
      
      wt = fl*wt
      if (inl.eq.0) then wt = 0
      
      do i = 1, imat
         p(i) = phom(i) *xt - pinhom(i)*wt
      enddo
      
      return
      end
      
      subroutine cntnodes(nqn,lang,ei,p,imat,node,iskip,iclose)
      implicit double precision (a-h,o-z)
      
#include "fortdim.h"
      common /grid/ h,r1,z,r(npdm),np
      common /consts/ etol,vtol,maxit,isoft
      common /itoo/ ietoolo,ietoohi
      common /elim/ emore,eless,elim
      common /iterm/ iterm

      dimension p(npdm)

      iskip=0
      ncross = 1

      do i = 2, imat
         if (p(i)*p(i-1).le.0.0) ncross = ncross + 1
      enddo

      if (ncross.ne.node) then
         iclose = 0
         if (ncross.lt.node) then
            ietoolo = ietoolo + 1
            eless = max(ei,eless)
            if (ietoohi.eq.0) ei = 0.75 * ei
            if (ietoolo.gt.2.and.ei.ge.-etol) then
               write (7,9000) nqn,lang,ei
               iterm=1
            endif
         else
            ietoohi = ietoohi + 1
            emore = min(ei,emore)
            if (ietoolo.eq.0) ei = 1.25 * ei
            if (ietoohi.gt.2.and.ei.lt.elim) then
               write(7,9010) nqn,lang,ei
               iterm=1
            endif
         endif
         if (ietoolo * ietoohi.gt.0) ei = 0.5 * (emore + eless)
         iskip=1
      endif


 9000 format(1x,'!WARNING! : schsl found too few nodes',
     $     '   n,l=',2i4,2x,'e= ', e10.4)
 9010 format(1x,'!WARNING! : schsl found too many nodes',
     $     '   n,l=',2i4,2x,'e= ', e10.4)

      return
      end

      subroutine inward(ei,nqn,lang,rv,p,phom,pinhom,fs,pref,
     $     imat,xt,wt,de,ap,iclose,imax)
      implicit double precision (a-h,o-z)
      
#include "fortdim.h"
      common /consts/ etol,vtol,maxit,isoft
      common /grid/ h,r1,z,r(npdm),np
      common /nlpot2/ inl,indrc(n0),IB,ID,IG

      dimension rv(npdm),p(npdm)
      dimension rho(npdm)
      dimension v(npdm),pinhom(npdm),xai(npdm),xci(npdm),xdi(npdm)
      dimension xli(npdm),phom(npdm),xah(npdm),xch(npdm)
      dimension xdh(npdm),xlh(npdm),g(npdm),fs(npdm)

      x2 = (lang + 0.5)**2
      h2 = h * h/12.     
      npm1=np-1
      imp1 = imat + 1

      do i = 1, np
         v(i) = ((rv(i) - ei * r(i)) * r(i) + x2) * h2
         g(i)=-fs(i)*sqrt(r(i))*r(i)
      enddo

      if (iclose.eq.0) then
         do imax = npm1,1,-1
            val = v(imax)/h2 - x2
            if (val.le.750.) goto 752
         enddo
         write(7,9000) nqn,lang
         iterm=1
 752     continue
         if (imax.le.12) imax = 12
      endif
      
      xah(imat-1) = 1 - v(imat-1)
      xdh(imat-1) = -(2 + 10 * v(imat-1))
      xah(imat) = 1 - v(imat)
      xdh(imat) = -(2 + 10 * v(imat))
      xlh(imat) = xdh(imat) 
      xch(imat) = -(xah(imat-1)) * phom(imat-1)

      do i = imp1,imax+1
         xah(i) = 1 - v(i)
         xch(i) = -xch(i-1) * xah(i-1)/xlh(i-1)
         xdh(i) = -(2 + 10 * v(i))
         xlh(i) = xdh(i) - xah(i) * xah(i-1)/xlh(i-1)
      enddo
      rootv = h * sqrt(abs(v(imax)/h2))
      ap = exp(-rootv)
      phom(imax) = xch(imax)/(xlh(imax) + xah(imax+1) * ap)
      
      do i = imax-1,imat,-1
         phom(i) = (xch(i) - xah(i+1) * phom(i+1))/xlh(i)
      enddo

      xai(imat-1) = 1 - v(imat-1)
      xdi(imat-1) = -(2 + 10 * v(imat-1))
      xai(imat) = 1 - v(imat)
      xdi(imat) = -(2 + 10 * v(imat))
      xli(imat) = xdi(imat) 
      xci(imat) = -(xai(imat-1)) * pinhom(imat-1) +
     $     h2*(g(imat+1)+10*g(imat)+g(imat-1))

      do i = imp1,imax+1
         xai(i) = 1 - v(i)
         xci(i) = h2*(g(i+1)+10*g(i)+g(i-1))-xci(i-1) * 
     $        xai(i-1)/xli(i-1)
         xdi(i) = -(2 + 10 * v(i))
         xli(i) = xdi(i) - xai(i) * xai(i-1)/xli(i-1)
      enddo
      pinhom(imax) = xci(imax)/(xli(imax) + xai(imax+1) * ap)

      do i = imax-1,imat,-1
         pinhom(i) = (xci(i) - xai(i+1) * pinhom(i+1))/xli(i)
      enddo
      
      do i = imat, np
         p(i) = phom(i)*xt - pinhom(i)*wt
      enddo
      
      do i = 1, np
         rho(i) = p(i) * p(i)
      enddo
      ap = 2 * lang + 1
      
      call radin (r,rho, 1, imax, h, ap)
      
      ap = ap + 0.5 * p(imax) * p(imax) * r(imax)/rootv
      
      xah(imat-2) = 1 - v(imat-2)
      res = xah(imat-2)*p(imat-2)+xdh(imat-1)*p(imat-1)+xah(imat)*
     $     p(imat)
      res1= h2*(wt)*(g(imat-2)+10*g(imat-1)+g(imat))
      de = (-p(imat-1) * (res+res1))/(h*ap)

 9000 format(1x,'!WARNING! : schsl found (v-e)*r^2 to be too low',
     $     ' for n,l=',2i4)
      
      return
      end





 
