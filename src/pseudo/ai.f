      subroutine ai
      
c     *************************************************************************
c     This subroutine computes the definite integral from 0 to rc of
c     r**2 * (jl(qr))**2 dr analytically.
c     Now if q<0 integrand is r**2 * (il(abs(q)r))**2 , where
c     il(x) = i**-l * jl(ix)
c     *************************************************************************

      implicit double precision(a-h,o-z)
      
#include "PARAMOPT"

      common /a/ a(numfn0)
      common /cuts/ qc,rc
      common /roots/ xroot(numfn0)
      common /angm/ ll
      common /numfn/ numfn
      common /bs/ bs(numfn0)
      common /bd/ bd(numfn0)
      
      do i = 1,numfn
         x = xroot(i) * rc
         xs = 1.0
         if (xroot(i).lt.0.0) xs = -1.0
         bs(i) = besfn(x,ll)
         bd(i) = besder(x,ll)
         t1 = 1.0 - xs * dfloat(ll * (ll + 1))/x/x
         a(i) = (t1 * bs(i)**2 + xs * bd(i) * bs(i)/abs(x) + 
     $        xs * bd(i)**2) * rc**3/2.0
      enddo
      
      return
      end
      
