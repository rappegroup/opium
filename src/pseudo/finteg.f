      subroutine finteg
      implicit double precision (a-h,o-z)
c Here we integrate  phi(r)*delsq(phi(r)) from rc to infinity.
c We do this in one interval (relatively low accuracy) since it is
c a costly calculation and is relatively unimportant to the calc.
      
#include "PARAMOPT"
c -------------------------------------------------------------------------
c     External (shared between C and Fortran) common blocks
c -------------------------------------------------------------------------
      common /grid/ h,r1,z,r(npdm),np
      common /atomic/ xion,rnl(npdm,n0),nlm(n0),wnl(n0),en(n0),norb
      common /np/ ncores,nvales
      common /nmax/ nmax(nvale0),maxim
c -------------------------------------------------------------------------

c -------------------------------------------------------------------------
c     Internal (Fortran only) common blocks
c -------------------------------------------------------------------------
      common/gauss/ xquad(nquad0),wquad(nquad0)
      common /nnn/ nnn
      common/cuts/ qc,rc
      common/angm/ ll
      common /f/ fint
c -------------------------------------------------------------------------
      dimension xf(nquad0),f(npdm)

      nord = 10
      do i = 1,np
         f(i) = rnl(i,nnn)
      enddo

      nint = 1
      width = (r(nmax(nnn)) - rc)/float(nint)
      zsum = 0.0

      do i = 1,nint
         xlend = float(i-1) * width + rc
         do j = 1,nquad0
            xf(j) = (xquad(j) + 1.0)*width/2.0 + xlend

            vv = val(f,r,np,xf(j),nord)
            v2 = val2(f,r,np,xf(j),nord)
            v3 = val3(f,r,np,xf(j),nord)

            xl = float(ll * (ll + 1))
            dl = v3+2.0*v2/xf(j)-xl*vv/xf(j)**2
            zsum=zsum+vv*dl*xf(j)**2*wquad(j)
         enddo
      enddo

      zsum = zsum * width/2.0
      fint = zsum

      return
      end
