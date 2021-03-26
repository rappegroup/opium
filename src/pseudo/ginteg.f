      subroutine ginteg
      implicit double precision (a-h,o-z)
c This code calculates the integral from 0 to qc of
c q**4 * c(q)**2 using gaussian quadrature.

#include "PARAMOPT"

c -------------------------------------------------------------------------
c     Internal (Fortran only) common blocks
c -------------------------------------------------------------------------
      common/gauss/ xquad(nquad0),wquad(nquad0)
      common/quads/ qq(nquad0)
      common/cuts/ qc,rc
      common /c/ c(nquad0)
      common /g/ gint
c -------------------------------------------------------------------------

      sum = 0.0

      do j = 1,nquad0
         sum = sum + c(j)**2 * qq(j)**4 * wquad(j)
      enddo

      gint = sum*qc/2.0

      return
      end
