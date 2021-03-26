      subroutine biq
      
c     *************************************************************************
c     This subroutine calculates the definite integral from 0 to rc of
c     r**2 * jl(xroot(i) * r) * jl(qq(j) * r) dr analytically.
c     *************************************************************************

      implicit double precision(a-h,o-z)
      
#include "PARAMOPT"

      common /quads/ qq(nquad0)
      common /b/ b(numfn0,nquad0)
      common /cuts/ qc,rc
      common /roots/ xroot(numfn0)
      common /angm/ ll
      common /numfn/ numfn
      common /bs/ bs(numfn0)
      common /bd/ bd(numfn0)
      
      pi = 3.14159265358979323
      con1 = sqrt(2.0/pi) * rc**3
      
      do j = 1,nquad0
        x2 = qq(j) * rc
        bd2 = besder(x2,ll)
        bs2 = besfn(x2,ll)
        do i = 1,numfn
          x1 = xroot(i) * rc
          xs = 1.0
          if (xroot(i).lt.0.0) xs = -1.0
          t1 = bs2 * bd(i) * abs(x1) - bs(i) * bd2 * x2
          b(i,j) = t1 * con1/(x2 * x2 - xs * x1 * x1)
        enddo
      enddo
 
      return
      end
      
