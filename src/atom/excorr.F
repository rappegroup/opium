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
      subroutine excorr(maxim,ixc,exccut,rsatom,
     $                  rvxc,rexc,rvx,rex,rvc,rec,
     $                  irelxc)
      implicit double precision (a-h,o-z)
      
#include "fortdim.h"
      parameter(thrd=1.d0/3.d0,thrd2=2.d0*thrd)
      parameter(pi32=29.608813203268075856503472999628d0)
      parameter(pi=3.1415926535897932384626433832795d0)
      parameter(alpha=1.91915829267751300662482032624669d0)

      common /grid/ h,r1,z,r(npdm),np

      dimension rsatom(npdm),rvxc(npdm),rexc(npdm)
      dimension agde(npdm),agdc(npdm),ggagde(npdm),ggagdc(npdm),
     $     rlape(npdm),rlapc(npdm),dens(npdm),
     $     rex(npdm),rec(npdm),rvx(npdm),rvc(npdm),tau(npdm)
      
      dimension bb(3,3),vcc(3)
      dimension ipvt(3)

c     This routine is now the driver for all other XC routines - EJW

c     irel - 0: no rel corrections, 1: use corrections

c     ixc - 0 : PZ LDA 81
c     ixc - 1 : PW LDA 92
c     ixc - 2 : PBE GGA 96
c     ixc - 3 : PW GGA 91
c     ixc - 4 : WC GGA 06
c     ixc - 5 : PBEsol GGA 08
c     ixc - 6 : VWN LDA V (the CA formulation)
c     ixc - 7 : PBE0 hybrid functional 

      ictype=0
      if (exccut.lt.0) ictype=1

c     ictype=0  : smoothing at origin is done by slow crossover to LDA
c     ictype=1  : smoothing at origin is done by polynomial fit when descreening ONLY
      
c      if ((ixc.gt.7).or.(ixc.lt.0)) then
      if ((ixc.gt.7)) then
         write(7,*) 'ixc out of range - stop'
         write(6,*) 'ixc out of range - stop'
         stop
      endif

      iset=0
      do 10 i=1,maxim
         rvxc(i)=0.0
         rexc(i)=0.0
         rvx(i)=0.0
         rvc(i)=0.0
         rex(i)=0.0
         rec(i)=0.0

         dens(i) = rsatom(i)/(4.0*pi*r(i)*r(i))
         if (iset.eq.0) then
            if (dens(i).lt.1e-16.and.r(i).gt.5.d0) then
               inxc = i-1
               iset = 1
            endif
         endif
 10   continue
      if (iset.eq.0) inxc=maxim-1

      if (ixc.ge.2.and.ixc.le.7) then
         call grad(dens,agde,
     $        ggagde,rlape,agdc,ggagdc,rlapc,tau,inxc)
         if (irelxc.eq.1) then
            a1l=2.21259
            a2l=.669152
            b1l=1.32998	
            b2l=0.794803
            
            a1t=3.48754
            a2t=0.218599
            b1t=1.15417
            b2t=0.015802
            
         endif
      endif

      zero=0.d0
      do i=1,inxc
         rho=rsatom(i)/(4.0*pi*r(i)*r(i))
         fk=(pi32*rho)**thrd
         rs=alpha/fk

c     add the regular smoothing for all GGA atoms by default

c     This should already be done ; why repeat?
c         if (ixc.ge.2.and.ixc.lt.5.and.exccut.lt.-10) exccut=1e-3

c     saying exccut.gt.1e-10 is ok since we want to skip
c     this if-then if it is negative anyway
         rwgt=1.d0
         if (r(i).lt.exccut.and.exccut.gt.1e-10) then
            rexp=-20.d0/exccut**2
            rwgt=exp(rexp*(exccut-r(i))**2)
         endif

         if (ixc.eq.1.or.ixc.eq.0.or.ixc.eq.6) then
c     This is PW92 LDA X which is the same LDA X for PZ X LDA and VWN X LDA
            call exchpbe(rho,zero,zero,zero,0,1,ex,vx,fx,fs,fss)
            if (irelxc.eq.1) then
               beta = 0.0140D0/rs
               sb = sqrt(1+beta*beta)
               alb = log(beta+sb)
               ex=ex*(1.D0 - 1.5D0 
     $              * ((beta*sb-alb) / beta**2)**2)
               vx=vx*(-0.5D0 + 1.5D0 * alb / (beta*sb))
            endif
         endif

         if (ixc.ge.2.and.ixc.le.7) then

            if (ixc.eq.2.or.ixc.eq.7) then
c     This is PBEGGA X
               call exchpbe(rho,agde(i),ggagde(i),rlape(i),1,1,
     $              ex,vx,fx,fs,fss)

c     This is PW92LDA X - just for terms needed in rel correction
               call exchpbe(rho,zero,zero,zero,0,1,
     $              exlda,vxlda,f1,f2,f3)
            endif

            if(ixc.eq.3) then
c     This is PW91GGA X
               call exchpw91(rho,agde(i),ggagde(i),rlape(i),
     $              ex,vx,fx,fs,fss)

c     This is PW92LDA X - just for terms needed in rel correction
               call exchpbe(rho,zero,zero,zero,0,1,
     $              exlda,vxlda,f1,f2,f3)
            endif

            if(ixc.eq.4) then
c     This is WCGGA X
               call exchpbewc(rho,agde(i),ggagde(i),rlape(i),1,1,
     $              ex,vx,fx,fs,fss)

c     This is WCLDA X - just for terms needed in rel correction
               call exchpbewc(rho,zero,zero,zero,0,1,
     $              exlda,vxlda,f1,f2,f3)
            endif

            if (ixc.eq.5) then
c     This is PBEsolGGA X
               call exchpbesol(rho,agde(i),ggagde(i),rlape(i),1,1,
     $              ex,vx,fx,fs,fss)
               vx=vx*rwgt
               
c     This is PW92LDA X - just for terms needed in rel correction
               call exchpbesol(rho,zero,zero,zero,0,1,
     $              exlda,vxlda,f1,f2,f3)
            endif
            
            if (irelxc.eq.1) then
               beta = 0.0140D0/rs
               sb = sqrt(1+beta*beta)
               alb = log(beta+sb)
               phi0=(1.D0 - 1.5D0 
     $              * ((beta*sb-alb) / beta**2)**2)

               chighl=(1+a1l*beta*beta+a2l*beta**4)
               clowl=(1+b1l*beta*beta+b2l*beta**4)
               phil=chighl/clowl
               
               dhighl=2*a1l*beta+4*a2l*beta**3
               dlowl=2*b1l*beta+4*b2l*beta**3
               sqlowl=(1+b1l*beta*beta+b2l*beta**4)**2
               dphil=(clowl*dhighl-chighl*dlowl)/sqlowl
               
               chight=(a1t*beta*beta+a2t*beta**4)
               clowt=(1+b1t*beta*beta+b2t*beta**4)
               phit=chight/clowt
               
               dhight=2*a1t*beta+4*a2t*beta**3
               dlowt=2*b1t*beta+4*b2t*beta**3
               sqlowt=(1+b1t*beta*beta+b2t*beta**4)**2
               dphit=(clowt*dhight-chight*dlowt)/sqlowt
               
               phi2=phil+phit
               dphilt=dphil+dphit

               s=agde(i)
               u=tau(i)
               v=rlape(i)

               if (r(i).lt.1e-4) then
                  u=0.0
                  v=0.0
               endif
               g=fx-1.d0
               dg=fs/2.d0
               d2g=fss/(4.d0*s)
               
               exlsdr = exlda*phi0
               ex = exlsdr+(ex-exlda)*phi2
               
               vxpr = vxlda*(-0.5D0 + 1.5D0 * alb / (beta*sb))
               vx=vxpr+vxlda*(0.25*beta*dphilt*(g-2*s**2*dg) 
     $              +phi2*(g-1.5*v*dg-1.5*d2g*u))
            endif
         endif
         rexc(i)=2.d0*r(i)*(ex*rwgt+(1.d0-rwgt)*exlda)
         rvxc(i)=2.d0*r(i)*(vx*rwgt+(1.d0-rwgt)*vxlda)
         rex(i)=rexc(i) 
         rvx(i)=rvxc(i) 

         if (ixc.eq.0) then
c     This is PZ LDA C
            call corpz(rs,ec,vc)
         endif
         if (ixc.eq.6) then
c     This is VWN3 LDA C
            call corvwn(rs,ec,vc,dens(i))
         endif
         if (ixc.eq.1) then
c     This is PW92 LDA C
            call corlsd(rs,zero,ec,vc,vc,ecrs,eczet,alfc)
         endif
         if (ixc.eq.2.or.ixc.eq.4.or.ixc.eq.7) then
c     This is PBE GGA C - used for WCGGA also
            call corpbe(rs,zero,agdc(i),ggagdc(i),rlapc(i),zero,1,1,
     $           ec,vc,vc,hc,dvc,dvc)
            ec=ec+hc*rwgt
            vc=vc+dvc*rwgt
         endif            
         if (ixc.eq.3) then
c     This is PW91 GGA C 
            gg=1.d0
            call corlsd(rs,zero,ec,vc,vc,ecrs,eczet,alfc)
            call corpw91(rs,zero,gg,ec,ecrs,eczet,agdc(i),
     $           ggagdc(i),rlapc(i),zero,hc,dvc,dvc)
            ec=ec+hc*rwgt
            vc=vc+dvc*rwgt
         endif            
         if (ixc.eq.5) then
c     This is PBEsol GGA C
            call corpbesol(rs,zero,agdc(i),ggagdc(i),rlapc(i),zero,1,1,
     $                     ec,vc,vc,hc,dvc,dvc)
            ec=ec+hc*rwgt
            vc=vc+dvc*rwgt
         endif

         rec(i)=r(i)*ec*2.d0 
         rvc(i)=r(i)*vc*2.d0 

         rexc(i)=rexc(i)+r(i)*ec*2.d0 
         rvxc(i)=rvxc(i)+r(i)*vc*2.d0 
      enddo
      if (ictype.eq.1) then
         do i=1,maxim
            if (r(i).gt.-exccut) then 
               icut=i
               goto 922
            endif
         enddo
 922     continue

         !smoothen 
         call smoothen(rvxc,icut,r)
         call smoothen(rexc,icut,r)
         call smoothen(rvx,icut,r)
         call smoothen(rex,icut,r)
         call smoothen(rvc,icut,r)
         call smoothen(rec,icut,r)

      endif

      return
      end



      subroutine grad(dens,agde,ggagde,rlape,agdc,ggagdc,
     $     rlapc,tau,inxc)
      implicit double precision(a-h,o-z)
      
#include "fortdim.h"
      
      parameter(o3=1.0/3.0,pi=3.141592653589793238462643,pi2=pi*pi)

c     -------------------------------------------------------------------------
c     External (shared between C and Fortran) common blocks
c     -------------------------------------------------------------------------
      common /grid/ h,r1,z,r(npdm),np
c     -------------------------------------------------------------------------

      dimension dens(npdm),agrd(npdm),grd(npdm)
      dimension agde(npdm),agdc(npdm),ggagde(npdm),ggagdc(npdm),
     $     rlape(npdm),rlapc(npdm),rgrd(npdm),tau(npdm)

c     prepare the various gradient terms

      do i=1,npdm
         agrd(i)=0.0
         grd(i)=0.0
         agde(i)=0.0
         agdc(i)=0.0
         ggagde(i)=0.0
         ggagdc(i)=0.0
         rlap=0.0
         rlape(i)=0.0
         rlapc(i)=0.0
         tau(i)=0.0
      enddo

      do i=inxc+1,npdm
         dens(i)=0.0
      enddo
      
      do i=2,inxc-1
         fk=(3*pi2*dens(i))**o3
         fs=sqrt(4*fk/pi)
         grd(i)=(dens(i+1)-dens(i-1))/2.0
         grd(i)=grd(i)*(1.0/(h*r(i)))
         rgrd(i)=r(i)*grd(i)+dens(i)
         agrd(i)=abs(grd(i)) 
         agde(i)=(1.0/(2*fk*dens(i))) * agrd(i)
         agdc(i)=(1.0/(2*fs*dens(i))) * agrd(i)
      enddo

      rgrd(1)=2*rgrd(2)-rgrd(3)
      grd(1)=2*grd(2)-grd(3)
      agrd(1)=2*agrd(2)-agrd(3)
      agde(1)=2*agde(2)-agde(3)
      agdc(1)=2*agdc(2)-agdc(3)

      rgrd(inxc)=2*rgrd(inxc-1)-rgrd(inxc-2)
      grd(inxc)=2*grd(inxc-1)-grd(inxc-2)
      agrd(inxc)=2*agrd(inxc-1)-agrd(inxc-2)
      agde(inxc)=2*agde(inxc-1)-agde(inxc-2)
      agdc(inxc)=2*agdc(inxc-1)-agdc(inxc-2)
      
      do i=2,inxc-1
         fk=(3*pi2*dens(i))**o3
         fs=sqrt(4*fk/pi)
         ggagd=(agrd(i+1)-agrd(i-1))/2.0
         ggagd=(1.0/(h*(r(i)))) * ggagd

         ggagde(i)=ggagd * (grd(i)/(dens(i) 
     $        * dens(i)*(2*fk)**3))
         ggagdc(i)=ggagd * (grd(i)/(dens(i) 
     $        * dens(i)*(2*fs)**3))

         rgrd2=(rgrd(i+1)-rgrd(i-1))/2.0
         rlap=rgrd2/(h*r(i)*r(i))
         rlape(i)=rlap/(dens(i)*(2*fk)**2)
         rlapc(i)=rlap/(dens(i)*(2*fs)**2)
      enddo
      rlape(1)=2*rlape(2)-rlape(3)
      rlapc(1)=2*rlapc(2)-rlapc(3)
      ggagde(1)=2*ggagde(2)-ggagde(3)
      ggagdc(1)=2*ggagdc(2)-ggagdc(3)


      rlape(inxc)=-2*rlape(inxc-1)-rlape(inxc-2)
      rlapc(inxc)=-2*rlapc(inxc-1)-rlapc(inxc-2)
      ggagde(inxc)=-2*ggagde(inxc-1)-ggagde(inxc-2)
      ggagdc(inxc)=-2*ggagdc(inxc-1)-ggagdc(inxc-2)


      do i=4,inxc-3
         fkplus=(3*pi2*dens(i+1))**o3
         fkminus=(3*pi2*dens(i-1))**o3
         fk=(3*pi2*dens(i))**o3
         splus=abs(dens(i+2)-dens(i))/(4*h*r(i+1)*dens(i+1)*fkplus)
         sminus=abs(dens(i)-dens(i-2))/(4*h*r(i-1)*dens(i-1)*fkminus)
         dels2=(splus**2 - sminus**2)/(2*h*r(i))
         deln=(dens(i+1)-dens(i-1))/(2*h*r(i))
         tau(i)=dels2*deln/(4*dens(i)*fk**2)
      enddo

 922  continue

      return
      end

      subroutine smoothen(rv,icut,r)
        implicit double precision(a-h,o-z)
#include "fortdim.h"
        !common /grid/ h,r1,z,r(npdm),np

        dimension bb(3,3),vcc(3)
        dimension ipvt(3)
        dimension rv(npdm),r(npdm)

        rvx1=rv(icut-1)/r(icut-1)
        rvx2=rv(icut)/r(icut)
        rvx3=rv(icut+1)/r(icut+1)
        rr1=r(icut-1)
        rr2=r(icut)
        rr3=r(icut+1)

        bb(1,1)=rr1**2
        bb(1,2)=rr1
        bb(1,3)=1.d0

        bb(2,1)=rr2**2
        bb(2,2)=rr2
        bb(2,3)=1.d0

        bb(3,1)=rr3**2
        bb(3,2)=rr3
        bb(3,3)=1.d0

        vcc(1)=rvx1
        vcc(2)=rvx2
        vcc(3)=rvx3

        call dgesv(3,1,bb,3,ipvt,vcc,3,ierr)
        if (ierr.ne.0) stop

        a=vcc(1)
        b=vcc(2)
        c=vcc(3)

        do i=1,icut
          rv(i)=r(i)*(a*r(i)**2+b*r(i)+c)
        enddo


        return
        end

