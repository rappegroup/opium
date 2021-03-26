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
      subroutine tmpot(ixc)
      implicit double precision (a-h,o-z)
      
#include "fortdim.h"
      common /grid/ h,r1,z,r(npdm),np

      common /angm/ ll
      common /nnn/ nnn
      common /rscore/ rscore(npdm),rdd(npdm),rddd(npdm),rscoretot(npdm),
     $     rsval(npdm)
      common /aorb/ ncore,nval,norb,nlm(n0),
     $     no(n0),lo(n0),nmax(n0),maxim
      common /adat/ wnl(n0),en(n0),so(n0),xion
      common /aval/ rcall(n0),rvap(n0),rnorm(n0),ibd(n0),etot
      common /wfn/ rnl(npdm,n0)

      common /numfn/ numfn
      common /wavrc/ wavrc, slope, curvrc, indrc
      common /totpot/ rvcore(npdm,n0),rvps(npdm,n0),rvcoul(npdm)
      
      dimension rpsi(npdm),rtemp(npdm),vtemp(npdm)

      pi=acos(-1.0)

      astep = 0.2
      idir = 1
      iter = 1
      itmax = 10000

      llp=ll+1

c     Save R_nl in rpsi
      do i = 1,np
         rpsi(i) = rnl(i,nnn)*r(i)
      enddo

c     Save core_charge (headnorm)
      do i=1,indrc
         rtemp(i)=(rnl(i,nnn)*r(i))**2
      enddo

      rrnorm=2*ll+2
      call radin(r,rtemp,0,indrc,h,rrnorm)
      headnorm=rrnorm

c     rtemp is tail part of rho
      do i=indrc+1,np
         rtemp(i)=(rnl(i,nnn)*r(i))**2
      enddo

c     ------------------------------------------------------------
c     Fit rpsi to TM form
      
      delta=0.0
      rc=r(indrc)
      rc1=rc
      rc2=rc*rc
      rc3=rc2*rc
      rc4=rc3*rc
      rc5=rc4*rc
      rc6=rc5*rc
      rc8=rc6*rc2
      rc10=rc8*rc2
      rc12=rc10*rc2
      rr=rc**(llp)


c     curvrc = d2  R_nl / dr2
c     curvrc2= d2 (r*R_nl) / dr2

      curvrc2 = curvrc - 2.0/rc*slope + float(ll*(ll+1))/rc/rc*wavrc
      wavrc=wavrc*rc

c     factor of 0.5 => Ha units
      ei=0.5*en(nnn)
      do i = 1,np
         vtemp(i) = 0.5*rvps(i,nnn)/r(i)
      enddo

      rvps0=vtemp(indrc)
      rvps2=(16*(vtemp(indrc-1)+vtemp(indrc+1))-(30*rvps0+
     $     vtemp(indrc-2)+vtemp(indrc+2)))/12.0d0
      rvps1=(vtemp(indrc-2)+8*(vtemp(indrc+1)-vtemp(indrc-1))-
     $     vtemp(indrc+2))/12.0d0
      rvps2=(rvps2-h*rvps1)/(h*rc)**2
      rvps1=rvps1/(h*rc)

      psi0=rpsi(indrc)
      psi2=(16*(rpsi(indrc-1)+rpsi(indrc+1))-(30*psi0+
     $     rpsi(indrc-2)+rpsi(indrc+2)))/12.0d0
      psi1=(rpsi(indrc-2)+8*(rpsi(indrc+1)-rpsi(indrc-1))-
     $     rpsi(indrc+2))/12.0d0

      psi2=(psi2-h*psi1)/(h*rc)**2
      psi1=psi1/(h*rc)

      slope=psi1
      curvrc2=psi2

      vir2=rvps2
      vir1=rvps1
      vir0=rvps0

      eq1=log(wavrc/rr)      
      eq2=slope/wavrc-llp/rc
      eq3=2*rvps0-2*ei-(2*llp/rc)*eq2-eq2*eq2
      eq4=2*rvps1+2*LLP*eq2/rc2-2*llp*eq3/rc-2*eq2*eq3
      eq5=2*(rvps2+llp*(2*(eq3-eq2/rc1)/rc1-eq4)/rc1-eq3*eq3-eq2*eq4)

      c0=0.0d0
      c2=0.0d0
      
      c4n =-(5760*c0 - 5760*eq1 + 2295*eq2*rc+ 1920*c2*rc2-375*eq3*rc2
     $     + 30*eq4*rc3 - eq5*rc4)/(384.*rc4)
      
      c6n = -(-3840*c0 + 3840*eq1 - 1761*eq2*rc - 960*c2*rc2 + 
     $     321*eq3*rc2 - 28*eq4*rc3 + eq5*rc4)/(96.*rc6)
      
      c8n = -(2880*c0 - 2880*eq1 + 1395*eq2*rc + 640*c2*rc2 - 
     $     275  *eq3*rc2 + 26*eq4*rc3 - eq5*rc4)/(64.*rc8)
      
      c10n =-(-2304*c0 + 2304*eq1 - 1149*eq2*rc - 480*c2*rc2 + 
     $     237*eq3*rc2 - 24*eq4*rc3 + eq5*rc4)/(96.*rc10)
 
      c12n =-(1920*c0 - 1920*eq1 + 975*eq2*rc + 384*c2*rc2 - 
     $     207  *eq3*rc2 + 22*eq4*rc3 - eq5*rc4)/(384.*rc12)


      max_tm_iter=20

      c2=0.0
      c0=0.0
      c2_shift=0.25d0
      eqc2_old=0.0d0

c     iterations over c2 begin here 

      do iter_c2=1,max_tm_iter

         c4=c4n-(5760*c0 + 1920*c2*rc2)/(384.*rc4)
         c6=c6n-(-3840*c0- 960*c2*rc2)/(96.*rc6)
         c8=c8n-(2880*c0+ 640*c2*rc2)/(64.*rc8)
         c10=c10n-(-2304*c0-480*c2*rc2)/(96.*rc10)
         c12=c12n-(1920*c0+384*c2*rc2)/(384.*rc12)
         
         eqc0_old=0.0d0
         c0_shift=0.25d0

c     iterations over c0 begin here 

         do iter_c0=1,max_tm_iter
            
c     Calc head norm
            do i = 1,indrc
               ri2=r(i)*r(i)
               ri4=ri2*ri2
               ri6=ri4*ri2
               ri8=ri6*ri2
               ri10=ri8*ri2
               ri12=ri10*ri2
               temp=c2*ri2+c4*ri4+c6*ri6+c8*ri8+c10*ri10+c12*ri12
               vtemp(i) = r(i)**llp * exp(temp)
            enddo

            rrnorm = float(ll+ll+2)
            do i=1,indrc
               rtemp(i)=(vtemp(i))**2
            enddo
            
            call radin(r,rtemp,0,indrc,h,rrnorm)

            eqc0=log(headnorm/rrnorm)-2*c0
            
            if (abs(eqc0) .gt. 1e-10) then
               c0_shift=-eqc0/((eqc0-eqc0_old)/c0_shift)
               eqc0_old=eqc0
               c0=c0+c0_shift
            else
               goto 18
            endif

            c4=c4n-(5760*c0 + 1920*c2*rc2)/(384.*rc4)
            c6=c6n-(-3840*c0- 960*c2*rc2)/(96.*rc6)
            c8=c8n-(2880*c0+ 640*c2*rc2)/(64.*rc8)
            c10=c10n-(-2304*c0-480*c2*rc2)/(96.*rc10)
            c12=c12n-(1920*c0+384*c2*rc2)/(384.*rc12)
            
            if (iter_c0 .eq. max_tm_iter) then
               write (*,*) "Couldn't find c0 in TM tmpot.f"
               write (7,*) "Couldn't find c0 in TM tmpot.f"
               stop
            endif
         enddo
         
 18      continue
         
         eqc2=c2*c2+c4*(2*ll+5)
         if (abs(eqc2) .gt. 1e-10) then
            c2_shift=-eqc2/((eqc2-eqc2_old)/c2_shift)
            eqc2_old=eqc2
            c2=c2+c2_shift
         else
            goto 19
         endif

         c4=c4n-(5760*c0 + 1920*c2*rc2)/(384.*rc4)
         c6=c6n-(-3840*c0- 960*c2*rc2)/(96.*rc6)
         c8=c8n-(2880*c0+ 640*c2*rc2)/(64.*rc8)
         c10=c10n-(-2304*c0-480*c2*rc2)/(96.*rc10)
         c12=c12n-(1920*c0+384*c2*rc2)/(384.*rc12)
         
         if (iter_c2 .eq. max_tm_iter) then
            write (*,*) "Couldn't find c2 in TM tmpot.f"
            write (7,*) "Couldn't find c2 in TM tmpot.f"
            stop
         endif
      enddo
      
 19   continue

      c4=c4n-(5760*c0 + 1920*c2*rc2)/(384.*rc4)
      c6=c6n-(-3840*c0- 960*c2*rc2)/(96.*rc6)
      c8=c8n-(2880*c0+ 640*c2*rc2)/(64.*rc8)
      c10=c10n-(-2304*c0-480*c2*rc2)/(96.*rc10)
      c12=c12n-(1920*c0+384*c2*rc2)/(384.*rc12)
      

      write(7,*) "Troullier-Martins has converged"

      do i = 1,indrc
         ri2=r(i)*r(i)
         ri4=ri2*ri2
         ri6=ri4*ri2
         ri8=ri6*ri2
         ri10=ri8*ri2
         ri12=ri10*ri2
         temp=c0+c2*ri2+c4*ri4+c6*ri6+c8*ri8+c10*ri10+c12*ri12
         rpsi(i) = r(i)**llp * exp(temp)
      enddo

c     ============================================================
c     rvps is V_scr^PS(r)

      ell=ll

      cc0=c0
      cc1=c2
      cc2=c4
      cc3=c6
      cc4=c8
      cc5=c10
      cc6=c12

      do i=1,indrc
         ri=r(i)
         ri2=r(i)*r(i)
         ri4=ri2*ri2
         ri5=ri4*ri
         ri6=ri4*ri2
         ri7=ri6*ri
         ri8=ri6*ri2
         ri10=ri8*ri2
         ri12=ri10*ri2
         
         polynom1=r(i)*(2*cc1+ri2*(4*cc2+ri2
     $        *(6*cc3+ri2*(8*cc4+ri2*(10*cc5+ri2*12*cc6)))))
         polynom2=2*cc1+ri2*(12*cc2+ri2*(30*cc3+ri2
     $        *(56*cc4+ri2*(90*cc5+ri2*132*cc6))))
         
         pp1= 2*(c2 + 2*c4*ri2 + 3*c6*ri4 + 4*c8*ri6 
     $        + 5*c10*ri8 + 6*c12*ri10)*ri

         pp2=2*c2+12*c4*ri2 + 30*c6*ri4 + 56*c8*ri6 + 
     $        90*c10*ri8 + 132*c12*ri10
         temp=ei+llp*pp1/ri+0.5*(pp2+pp1*pp1)
         rvps(i,nnn)=(ri*temp)*2

      enddo

c     put r * R_nl into rnl
      do i = 1,np
         rnl(i,nnn) = rpsi(i) 
         rsval(i)=rsval(i) + wnl(nnn)*rpsi(i)**2
         if (ixc.ge.0) then
            if (i.gt.indrc+50) rvps(i,nnn) = -z-z+rvcoul(i)
            if (i.gt.maxim) rvps(i,nnn) = -xion-xion
         endif
      enddo

      return
      end
