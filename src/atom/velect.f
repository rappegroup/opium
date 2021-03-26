      subroutine velect(iter,iconv,ixc,
     +     cdd,cdu,cdc,
     +     viod,viou,vid,viu,vod,vou,
     +     etot,ev,ek,ep)
      implicit double precision(a-h,o-z)
c     
c     generate the electronic output potential from
c     the electron charge density.
c     the ionic part is added in dsolv[12].

#include "param.h"

      parameter (nrmax3  = 3 * nrmax)
      parameter(eps=1e-18)

      common /rxccut/ rxccut
      common /reli/ norb,ncore,nval,no(norbp),lo(norbp)
      common /reld/ zcore,zval,zo(norbp),so(norbp)
      common /rgrid/ a,b,r(nrmax),rab(nrmax)
      common /nrgrid/ nr

      dimension 
     +     cdd(nrmax),cdu(nrmax),cdc(nrmax),
     +     viod(lmax,nrmax),viou(lmax,nrmax),vid(nrmax),viu(nrmax),
     $     vod(nrmax),vou(nrmax),
     +     etot(10),ev(norb),ek(norb),ep(norb)
      dimension y(nrmax),yp(nrmax),ypp(nrmax),w(nrmax3),
     +     s1(nrmax),s2(nrmax)
      dimension cdr(10000)

      pi = 4*atan(1.D0)

      if (iter.eq.0) then
c     set up initial charge density.
c     cdd and cdu  =  2 pi r**2 rho(r)
         aa = 0.5
         do i=1,nr
            cdd(i) = zel*aa**3*exp(-aa*r(i))*r(i)**2/4
            cdu(i) = cdd(i)
         enddo
      endif


c     fit cd/r by splines
      ex_ilya=0.0
      y(1) = 0.D0
      do 100 i=2,nr
         y(i) = (cdd(i)+cdu(i))/r(i)
 100  continue
      isx = 0
      a1 = 0.D0
      an = 0.D0
      b1 = 0.D0
      bn = 0.D0
      call splift(r,y,yp,ypp,nr,w,ierr,isx,a1,b1,an,bn)

c     compute the integrals of cd/r and cd from
c     r(1)=0 to r(i)

      xlo = 0.D0
      call spliq(r,y,yp,ypp,nr,xlo,r,nr,s2,ierr)
      do 110 i=1,nr
         ypp(i) = r(i)*ypp(i) + 2*yp(i)
         yp(i)  = r(i)*yp(i)  + y(i)
         y(i)   = r(i)*y(i)
  110 continue
      call spliq(r,y,yp,ypp,nr,xlo,r,nr,s1,ierr)

c     check normalization

      delta = zval+zcore - s1(nr)
      if (iconv .eq. 1 .and. abs(delta) .gt. 1.D-5) then
         write(nout,120) delta
  120    format(/,' warning *** unnormalized charge density in velect',
     +        'delta e =',e10.2,/)
      end if

c     compute new hartree potential

      do 130 i=2,nr
         vod(i) = 2 * (s1(i)/r(i) + s2(nr) - s2(i))
         vou(i) = vod(i)
  130 continue

c     compute hartree contribution to total energy
      
      if (iconv .eq. 1) then
         ehart = 0.D0
         ll = 4
         do 140 i=2,nr
            ehart = ehart + ll * (cdd(i)+cdu(i)) * vod(i) * rab(i)
            if (ifcore .ge. 2) then
               ehart = ehart + ll * cdc(i) * vod(i) * rab(i)
            end if
            ll = 6 - ll
  140    continue
         ehart = ehart / 6
      end if


c     add exchange and correlation

      trd = 1.D0/3.D0
      ftrd = 4*trd
      tftm = 2**ftrd-2
      a0 = (4/(9*pi))**trd

c     set x-alpha

      alp = 2 * trd
      vxc = 0.D0
      vc  = 0.D0
      exc = 0.D0
      ec  = 0.D0
      pi = 4*atan(1.D0)
      do i=1,nr
	cdr(i)=cdu(i)+cdd(i)
	cdr(i)=cdr(i)/(4*pi*r(i)**2+eps)
      enddo

c     start loop
c           
      ll = 4
      do 170 i=2,nr
         cdsum = cdd(i) + cdu(i)
c     if (ifcore .ge. 1) cdsum = cdsum + cdc(i)
         if (cdsum .le. 0.D0) goto 170
         rs = (3*r(i)**2/cdsum)**trd
         z = 0.D0
         fz = 0.D0
         fzp = 0.D0

c     exchange (only use (xa))

         vxp = -3*alp/(pi*a0*rs)
         exxp = 3*vxp/4
         beta = 0.0140D0/rs
         sb = sqrt(1+beta*beta)
         alb = log(beta+sb)
         vxp = vxp * (-0.5D0 + 1.5D0 * alb / (beta*sb))
         exxp = exxp * (1.D0 - 1.5D0 * ((beta*sb-alb) / beta**2)**2)
         vxf = 2**trd*vxp
         exf = 2**trd*exxp

         if (ixc.eq.0) then
            
c     correlation
c     ceperly-alder (ca)
c     
c     The Perdew-Zunger parameterization is used.
c     See Phys. Rev. B 23 5075 (1981).
c     
            if (rs .gt. 1.D0) then
               te = 1.D0+(7.D0/6.D0)*1.0529D0*sqrt(rs)
     +              + (4.D0/3.D0)*0.3334D0*rs
               be = 1.D0+1.0529D0*sqrt(rs)+0.3334D0*rs
               ecp = -0.2846D0/be
               vcp = -0.2846D0*te/be**2
               te = 1.D0+(7.D0/6.D0)*1.3981D0*sqrt(rs)
     +              + (4.D0/3.D0)*0.2611D0*rs
               be = 1.D0+1.3981D0*sqrt(rs)+0.2611D0*rs
               ecf = -0.1686D0/be
               vcf = -0.1686D0*te/be**2
            else
               ecp = +2*((0.0311D0+0.0020D0*rs)*log(rs)
     +              - 0.048D0-0.0116D0*rs)
               vcp = +2*((0.0311D0+2.D0/3.D0*0.0020D0*rs)*log(rs)
     +              - (0.048D0+0.0311D0/3.D0)
     +              - (2.D0/3.D0*0.0116D0+0.0020D0/3.D0)*rs)
               ecf = +2*((0.01555D0+0.0007D0*rs)*log(rs)
     +              - 0.0269D0-0.0048D0*rs)
               vcf = +2*((0.01555D0+2.D0/3.D0*0.0007D0*rs)*log(rs)
     +              - (0.0269D0+0.01555D0/3.D0)
     +              - (2.D0/3.D0*0.0048D0+0.0007D0/3.D0)*rs)
            end if
         else
           
            if(i.lt.2) goto 170

            exlsd=0.0
            fxpbe=0.0
            
            call pbe(cdr,cdu,r,nr,b,exxp,vxp,rab,
     $           ecp,vcp,i,ixc,fxpbe,exlsd,vxlsd,rholap)

            if (ixc.eq.2) then
               expbe=exxp
               vxpbe=vxp
               beta = 0.0140D0/rs
               sb = sqrt(1+beta*beta)
               alb = log(beta+sb)
               
               a1l=2.21259
               a2l=.669152
               b1l=1.32998	
               b2l=0.794803
               
               a1t=3.48754
               a2t=0.218599
               b1t=1.15417
               b2t=0.015802
               
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
               
               philt=phil+phit
               dphilt=dhil+dphit
               pi2=pi*pi
               
               if (i.ge.4) then
                  fk=(3*pi2*cdr(i))**(1.0/3.0)
                  fkplus=(3*pi2*cdr(i+1))**(1.0/3.0)
                  fkminus=(3*pi2*cdr(i-1))**(1.0/3.0)
                  s=abs(cdr(i+1)-cdr(i-1))
                  s=s/(2*rab(i))/(2*cdr(i)*fk)
                  if(cdr(i+2).gt.0.0) then
                     splus=abs(cdr(i+2)-cdr(i))
                     splus=splus/(2*rab(i+1))/(2*cdr(i+1)*fkplus)
                  endif
                  sminus=abs(cdr(i)-cdr(i-2))
                  sminus=sminus/(2*rab(i-1))/(2*cdr(i-1)*fkminus)
                  dels2=(splus**2-sminus**2)/(2*rab(i))
                  deln=(cdr(i+1)-cdr(i-1))/(2*rab(i))
                  
                  fnu=rholap/(4*cdr(i)*fk**2)
                  
                  tau=dels2*deln/(4*cdr(i)*fk**2)		
                  fkappa=0.804
                  fmu=0.21951
                  
                  dg=fmu/(1+fmu*s**2/fkappa)**2
                  d2g=-2.0*fmu**2/fkappa/(1+fmu*s**2/fkappa)**3
                  g=fkappa-fkappa/(1+fmu*s**2/fkappa)
               endif
               
               exlsdr=exlsd*(1.D0 - 1.5D0 *((beta*sb-alb)/beta**2)**2)
               vxpr = vxlsd * (-0.5D0 + 1.5D0 * alb / (beta*sb))
               vxprgga=vxpr+vxlsd*(0.25*beta*dphilt*(g-2*s**2*dg)
     $              +philt*(g-1.5*fnu*dg-1.5*tau*d2g))
               exxp = exlsdr+(expbe-exlsd)*philt		
               vxp=vxprgga     
            endif
         endif
         vxcp = vxp  + vcp
         vxcf = vxf + vcf
         vxcd = vxcp
         vxcu = vxcp
         excp = exxp + ecp
         excf = exf + ecf
         vcd = vcp
         vcu = vcp
         exct = excp
         ect = ecp
         
         exf=0.0
         ecf=0.0
         vxf=0.0
         vcf=0.0
         excf=0.0
         vxcf=0.0
         
         vxcd = vxcd + fz*(vxcf-vxcp) + (1-z)*fzp*(excf-excp)
         vxcu = vxcu + fz*(vxcf-vxcp) - (1+z)*fzp*(excf-excp)
         vcd = vcd + fz*(vcf-vcp) + (1-z)*fzp*(ecf-ecp)
         vcu = vcu + fz*(vcf-vcp) - (1+z)*fzp*(ecf-ecp)
         exct = exct + fz*(excf-excp)
         ect = ect + fz*(ecf-ecp)
	 
         vod(i) = vod(i) + vxcd
         vou(i) = vou(i) + vxcu
	 
         vxc = vxc + ll * (cdd(i)*vxcd + cdu(i)*vxcu) * rab(i)
         vc  = vc  + ll * (cdd(i)*vcd  + cdu(i)*vcu ) * rab(i)
         exc = exc + ll * cdsum * exct * rab(i)
         ec  = ec  + ll * cdsum * ect  * rab(i)
	
         ll = 6 - ll
 170  continue
	
      etot(4) = ehart
      etot(5) = vxc / 3
      etot(6) = (3*vc - 4*ec) / 3
      etot(7) = exc / 3

      vod(1) = vod(2) - (vod(3)-vod(2))*r(2)/(r(3)-r(2))
      vou(1) = vou(2) - (vou(3)-vou(2))*r(2)/(r(3)-r(2))

      return
      end



