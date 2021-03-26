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
      subroutine btrans(zeff,name1,name2)
      implicit double precision (a-h,o-z)
      
#include "fortdim.h"

      common /aorb/ ncore,nval,norb,nlm(n0),
     $     no(n0),lo(n0),nmax(n0),maxim
      common /adat/ wnl(n0),en(n0),so(n0),xion
      common /aval/ rcall(n0),rvap(n0),rnorm(n0),ibd(n0),etot
      common /wfn/ rnl(npdm,n0)
      common /grid/ h,r1,z,r(npdm),np
      common /rscore/ rscore(npdm),rdd(npdm),rddd(npdm),rscoretot(npdm),
     $     rsval(npdm)
      common /filenames/ file_log
      common /totpot/ rvcore(npdm,n0),rvps(npdm,n0),rvcoul(npdm)
      common /nlcore/ rvloc(npdm)
      common /psdat/ nll

      character*80 name1
      character*80 name2

      dimension pl(1000),rtemp(npdm)

      open(unit=70,file=name1,form='formatted')
      open(unit=60,file=name2,form='formatted')

      rewind(70)
      rewind(60)

      pi=acos(-1.0)                                                   
      pi4=4.0*pi                                                      
      rnq=sqrt(2.0/pi)

      do m=1,nval
         l = nlm(m)/10 - 10 * (nlm(m)/100)
         qmax=20.d0
         qp=500.d0
         dql=qmax/qp

         do j=1,nint(qp)
            ql=dql*dble(j-1)                                             

            do k=1,maxim                                                
               qr=r(k)*ql                                             
               rtemp(k) = rnl(k,m)*r(k) * besfn(qr,l)
            enddo
            pow=dble(l+1)**2
            call radin(r,rtemp,0,maxim,h,pow)
            pl(j)=pow * rnq
            write(70,*) ql,pl(j)
         enddo
         write(70,*) '@'
      enddo

      do m=1,nll
         l = nlm(m)/10 - 10 * (nlm(m)/100)
         qmax=100.d0
         qp=1000.d0
         dql=qmax/qp
         do j=1,nint(qp)
            ql=dql*dble(j-1)                                             
            do k=1,maxim                                                
               qr=r(k)*ql                                             
               rtemp(k) = (rvcore(k,m)*r(k)+zeff*2.0*r(k)) 
     $              * besfn(qr,l)
            enddo
            pow=dble(l+1)**2
            call radin(r,rtemp,0,maxim,h,pow)
            pl(j)=pow * rnq    
            write(60,*) ql,pl(j)
         enddo
         write(60,*) '@'
      enddo

      do j=1,nint(qp)                                                    
         ql=dql*dble(j-1)                                             
         do k=1,maxim                                                
            qr=r(k)*ql                                             
            rtemp(k) = (rvloc(k)*r(k)+zeff*2.0*r(k)) 
     $           * besfn(qr,0)
         enddo
         pow=1.0
         call radin(r,rtemp,0,maxim,h,pow)
         pl(j)=pow * rnq
         write(60,*) ql,pl(j)
      enddo
      write(60,*) '@'

      close(70)

      if (rpcc.gt.1e-6) then

         do j=1,nint(qp)
            ql=dql*dble(j-1)                                             
            do k=1,maxim                                                
               qr=r(k)*ql                                             
               rtemp(k) = rscore(k) * besfn(qr,0)
            enddo
            pow=1.0
            call radin(r,rtemp,0,maxim,h,pow)
            pl(j)=pow * pi4    
            write(60,*) ql,pl(j)
         enddo
      endif

      close(60)
      return
      end
