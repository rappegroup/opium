      subroutine ei
      implicit double precision (a-h,o-z)

c This code calculates the integral from 0 to qc of
c q**4 * bi(q) * c(q) using gaussian quadrature.

#include "PARAMOPT"
c -------------------------------------------------------------------------
c     Internal (Fortran only) common blocks
c -------------------------------------------------------------------------
      common/gauss/ xquad(nquad0),wquad(nquad0)
      common/quads/ qq(nquad0)
      common/cuts/ qc,rc
      common /b/ b(numfn0,nquad0)
      common /c/ c(nquad0)
      common /e/ e(numfn0)
      common /numfn/ numfn
c -------------------------------------------------------------------------

      do i = 1,numfn
         sum = 0.0
         do j = 1,nquad0
            sum = sum + b(i,j) * c(j) * qq(j)**4 * wquad(j)
         enddo
         e(i) = sum * qc/2.0
      enddo

      return
      end
