      subroutine dij
c This subroutine calculates the definite integral from 0 to qc of
c q**4 * b(i,...) * b(j,...) dq by gaussian quadrature.
c That is why b(i,...) is actually a 2-d array stored for each 
c b(i,...) at all the gaussian quadrature points.

      implicit double precision(a-h,o-z)
      
#include "PARAMOPT"

c -------------------------------------------------------------------------
c     Internal (Fortran only) common blocks
c -------------------------------------------------------------------------
      common/gauss/ xquad(nquad0),wquad(nquad0)
      common/quads/ qq(nquad0)
      common/b/ b(numfn0,nquad0)
      common/cuts/ qc,rc
      common/d/ d(numfn0,numfn0)
      common /numfn/ numfn
c -------------------------------------------------------------------------

      do i = 1,numfn
         do j = i,numfn
            sum = 0.0
            do k = 1,nquad0
               sum = sum + wquad(k) * b(i,k) * b(j,k) * qq(k)**4
            enddo
            d(i,j) = sum * qc/2.0
            d(j,i) = d(i,j)
         enddo
      enddo

      return
      end

