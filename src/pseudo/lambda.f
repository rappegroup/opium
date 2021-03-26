      subroutine lambda(rc,numfn)
      implicit double precision(a-h,o-z)

#include "PARAMOPT"

      common /rlm/ rlmat(3,numfn0),rlvec(numfn0)

      dimension rc(numfn0)

      rc(numfn-1)=0.0
      rc(numfn)=0.0

      alpha = rlmat(2,numfn)*rlmat(1,numfn-1)
     $     - rlmat(2,numfn-1)*rlmat(1,numfn)
      beta  = (-rlmat(2,numfn)*rlvec(1)+rlmat(1,numfn)*rlvec(2))

      do i=1,numfn-2
         rc(numfn-1) = rc(numfn-1) + rc(i)*(rlmat(1,numfn)*rlmat(2,i)
     $        - rlmat(2,numfn)*rlmat(1,i))
      enddo

      rc(numfn-1)=(rc(numfn-1)+beta)/alpha

      do i=1,numfn-1
         rc(numfn) = rc(numfn) + rlmat(1,i)*rc(i)
      enddo
      rc(numfn)=(-rlvec(1)-rc(numfn))/rlmat(1,numfn)

      end

               
