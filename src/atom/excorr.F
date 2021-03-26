c
c Copyright (c) 1998-2005 The OPIUM Group
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
      subroutine excorr(maxim,ixc,r,rsatom,rvxc,rexc)
      implicit double precision (a-h,o-z)
      
#include "fortdim.h"

      parameter(um=0.2195149727645171,uk=0.8040,ul=um/uk)
      parameter(o3=1.0/3.0,third4=4.0/3.0,thrd=1.0/3.0)
      parameter(thrdm=-thrd,thrd2=2.0*thrd,sixthm=thrdm/2.0)
      parameter(gam=0.5198420997897463295344212145565)
      parameter(gamma=0.03109069086965489503494086371273)
      parameter(bet=0.06672455060314922,delt=bet/gamma)


      dimension r(npdm),rsatom(npdm),rvxc(npdm),rexc(npdm)
      dimension agde(npdm),agdc(npdm),ggagde(npdm),ggagdc(npdm),
     $     rlape(npdm),rlapc(npdm),dens(npdm)

      data  pi/3.141592653589793238462643/

    5 continue
      xl= (18./pi**2)**o3
      xconst = xl / (3.**o3)
c     Get density from rsatom (=4*Pi*r*r*rho)
c     
c     Since we are taking many derivatives, we set inxc to be
c     the grid point at which the dens < 1e-18.  If not, the
c     gradient terms may overflow.

      iset=0
      do 10 i=1,maxim
         rvxc(i)=0.0
         rexc(i)=0.0

         dens(i) = rsatom(i)/(4.0*pi*r(i)*r(i))
         if (iset.eq.0) then
            if (dens(i).lt.1e-18.and.i.gt.200) then
               inxc = i
               iset = 1
            endif
         endif
         rvxci   =  xconst * (r(i) * rsatom(i))**o3
         rvxc(i) =  - max(rvxci,1.e-25)
            
 10   continue

c     rvxc(i) is the pure exchange potential * r

      b11 = 1.0529
      b22 = 0.3334
      aa = 0.0622
      bb =-0.0960
      cc = 0.0040
      dd =-0.0232
      ga =-0.2846
      t11 = 7./6.
      t22 = (4./3.)*b22      
      if (ixc.eq.2) call grad(dens,agde,
     $     ggagde,rlape,agdc,ggagdc,rlapc,inxc)

      do 510 i=1,inxc
         
         if (ixc.eq.2) then
            ex=0.75*rvxc(i)
            fxpbe=1.0+uk-uk/(1.0+ul*agde(i)**2)
            fs=2.0*uk*ul/((1.0+ul*agde(i)**2)**2)
            fss=-4.0*ul*agde(i)*fs/(1.0+ul*agde(i)**2)
            vxcon=(third4*fxpbe
     $           -(ggagde(i)-third4*agde(i)**3)
     $           *fss-rlape(i)*fs)
            vxpbe=ex*vxcon
         else
            fxpbe=1.0
            vxpbe=rvxc(i)
         endif
         rr = r(i)
         rs = -xl*rr/rvxc(i)

         if (ixc.eq.0) then

c     The Perdew-Zunger correlation functional
c     Low density correlation
            
            if (rs.lt.1) goto 506
            root = sqrt(rs) * b11
            xd = 1. + root + b22 * rs
            ec = ga/xd
            uc = (1.+t11*root+t22*rs)*ec/xd
            goto 509
            
c     High density correlation
            
 506        dl = log(rs)
            t33 = cc*rs
            t44 = dd*rs
            t55 = t33*dl
            ec = aa*dl +bb + t55 + t44
            uc = ec - (aa+ t55 + t44 + t33)*o3
            
         else

c     This is the Perdew-Wang LSD correlation functional
c     which is used with the PBE GGA enhancement factores

            a=0.031091
            a1=0.21370
            b1=7.5957
            b2=3.5876
            b3=1.6382
            b4=0.49294
            rts=sqrt(rs)
            
            q0=-2.0*a*(1.0+a1*rts*rts)
            q1=2.0*a*rts*(b1+rts*(b2+rts*(b3+b4*rts)))
            q2=log(1.0+1.0/q1)
            ec=q0*q2
            q3=a*(b1/rts+2.0*b2+rts*(3.0*b3+4.0*b4*rts))
            rdec=-2.0*a*a1*q2-q0*q3/(q1*(1.0+q1))
            uc=ec-rs*rdec/3.0
            ec=ec*2.0
            uc=uc*2.0
         endif
         
 509     continue
         
         if (ixc.eq.2) then
            suc=uc
            sec=ec
            ec=ec/2
            uc=uc/2
            pon=max(1e-8,-(ec)/gamma)
            b=delt/(exp(pon)-1.0)
            b2=b*b
            t=agdc(i)
            t2=t*t
            t4=t2*t2
            t6=t4*t2
            rs2=rs*rs
            rs3=rs2*rs
            q4=1.0+b*t2
            q5=1.0+b*t2+b2*t4
            hc=(bet/delt)*log(1.0+delt*q4*t2/q5)
            rsthrd=rs/3.0
            fac=delt/b+1.0
            ecrs=(3.0/rs)*(ec-uc)
            bec=b2*fac/bet
            q8=q5*q5+delt*q4*q5*t2
            q9=1.0+2.0*b*t2
            hb=-bet*b*t6*(2.0+b*t2)/q8
            hrs=-rsthrd*hb*bec*ecrs
            fact0=2.0*delt-6.0*b
            fact1=q5*q9+q4*q9*q9
            hbt=2.0*bet*t4*((q4*q5*fact0-delt*fact1)/q8)/q8
            hrst=rsthrd*t2*hbt*bec*ecrs
            ht=2.0*bet*q9/q8
            fact2=q4*q5+b*t2*(q4*q9+q5)
            fact3=2.0*b*q5*q9+delt*fact2
            htt=4.0*bet*t*(2.0*b/q8-(q9*fact3/q8)/q8)
            uu=ggagdc(i)
            vv=rlapc(i)

            ec=(ec+hc)*2.0
            ucu=uc
            uc=hc+hrs+hrst+t2*ht/6.0+7.0*t2*t*htt/6.0
            uc=uc-uu*htt-vv*ht
            uc=(uc+ucu)*2.0
         endif

         rexc(i) = fxpbe*rvxc(i) * 0.75 + ec * rr
         rvxc(i) = (vxpbe+uc*rr)

 510  continue

      return
      end
      
      subroutine grad(dens,agde,ggagde,rlape,agdc,ggagdc,
     $     rlapc,inxc)
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
     $     rlape(npdm),rlapc(npdm),rgrd(npdm)

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
      enddo

      do i=inxc,npdm
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

c      do i=1,inxc
c        if (r(i).gt.rxccut) goto 922
c        agde(i)=agde(i)*rg
c        agdc(i)=agdc(i)*rg
c        rlape(i)=rlape(i)*rg
c        rlapc(i)=rlapc(i)*rg
c        ggagde(i)=ggagde(i)*rg
c        ggagdc(i)=ggagdc(i)*rg
c      enddo
 922  continue

      return
      end
