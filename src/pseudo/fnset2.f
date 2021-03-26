      subroutine fnset2(ipr,rc,numfn)
      implicit double precision(a-h,o-z)

#include "PARAMOPT"

      common /rke/ rkmat(numfn0,numfn0),rkvec(numfn0)
      common /rn/ rnmat(numfn0,numfn0),rnvec(numfn0)
      common /rlm/ rlmat(3,numfn0),rlvec(numfn0)
      common /rconst/ rkcon,rncon,rlcon(3)
      common /re/ rkin,rnorm
      common /emat/ emat(3,3),rfl(numfn0)
      common /rf/ afor(3),rlam(3)
      common /rb/ rb(numfn0,numfn0)
      common /ncon/ ncon

      dimension rc(numfn0)

      rkin=0
      rnorm=0
      do i=1,numfn
         do j=1,numfn
            rkin=rkin+rc(i)*rkmat(i,j)*rc(j) 
            rnorm=rnorm+rc(i)*rnmat(i,j)*rc(j)
         enddo
         rkin=rkin+rkvec(i)*rc(i)
      enddo

      do i=1,ncon
         rlcon(i)=0.0
         do j=1,numfn
            rlcon(i)=rlcon(i)+rlmat(i,j)*rc(j)
         enddo
      enddo

      if (ipr.ne.0) then
c         write(7,9004)
c         write(7,9005) (rc(i),i=1,numfn)
c         write(7,9004)
         write(7,9000) rkin+rkcon
         write(7,9001) rnorm+rncon
         write(7,9002) rlcon(1)+rlvec(1)
         write(7,9003) rlcon(2)+rlvec(2)
         if (ncon.eq.3) write(7,9007) rlcon(3)+rlvec(3)
c         write(7,9004)
      endif
 
 9000 format(1x,'Resid KE (Ry)     :    ',f16.10)
 9001 format(1x,'Norm error        :    ',e10.3)
 9002 format(1x,'Continuity error  :    ',e10.3)
 9003 format(1x,'Curvature error   :    ',e10.3)
 9007 format(1x,'Sum       error   :    ',2f10.6)
 9005 format(1x,'Coeff:    ',20f10.6)
 9006 format(1x,'Coeff SUM         :    ',f10.6)
 9004 format(1x,'-------------------------------------',
     $     '--------------------------')
      
      return
      end
