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
      subroutine kcomp(irel,name,rkstor,ikstor)
      implicit double precision (a-h,o-z)
      
#include "fortdim.h"

      common /grid/ h,r1,z,r(npdm),np
      common /filenames/ file_log
      common /totpot/ rvcore(npdm,n0),rvps(npdm,n0),rvcoul(npdm)
      common /aorb/ ncore,nval,norb,nlm(n0),
     $     no(n0),lo(n0),nmax(n0),maxim
      common /adat/ wnl(n0),en(n0),so(n0),xion
      common /aval/ rcall(n0),rvap(n0),rnorm(n0),ibd(n0),etot
      common /wfn/ rnl(npdm,n0)

      character*80 name
      character*80 file_log

      dimension pl(5000),rtemp(npdm),rke(n0,5000),rkstor(n0,10)
      dimension ikstor(n0,10),ill(n0)

      open(unit=7,file=file_log,form='formatted',access='append')

      open(unit=71,file=name,form='formatted')

      rewind(71)

      pi=acos(-1.0)                                                   
      pi4=4.0*pi                                                      
      rnq=sqrt(2.0/pi)

      write(7,9393)
      write(7,9394) 


      do m=1,nval
         ill(m)=0
      enddo

      do m=1,nval
         l = nlm(m)/10 - 10 * (nlm(m)/100)
         if (ill(m).gt.0) goto 911
         ill(m)=ill(m)+1
         
         if (irel.eq.1) then
            if (so(m).gt.0) then
               write(7,9001) nlm(m)
            else
               write(7,9002) nlm(m)
            endif
         else
            write(7,9000) nlm(m)
         endif

         write(7,*) "     Ecut(Ry)  KE/electron(Ry)",
     $        "         KE error/electron (mRy,meV)"

         rks=0.0
         mxk=np
         qmax=20
         qp=500
         iqp=int(qp)
         dqq=qmax/qp

         do j=1,iqp+1
            rke(m,j)=0.0
         enddo

         do k=np,1,-1
            if (rnl(k,m).gt.1e-9) then 
               mxk=k
               goto 912
            endif
         enddo
 912     continue

         do j=1,iqp+1                                                    
            qq=dqq*dble(j-1)                                             

            do k=1,mxk
               qr=r(k)*qq                                             
               rtemp(k) = rnl(k,m)*r(k) * besfn(qr,l)
            enddo
            pow=dble(l+1)**2
            call radin(r,rtemp,0,mxk,h,pow)
            pl(j)=pow *rnq 

            rks=rks+pl(j)*pl(j)*qq*qq*qq*qq*dqq
            rke(m,j)=rks
         enddo
         write(71,*) '@'

         do kk=1,4
            ikstor(m,kk)=-1
            rkstor(m,kk)=0.0
         enddo
         rkkmlast=0.0
         do j=1,iqp                                                    
            qq=dqq*dble(j-1)                            
            rkk=rks-rke(m,j)
            rkkmlast=rkkm
            rkkm=rkk*1000*13.6057
            if ((rkkm.lt.1000).and.(rkkmlast.gt.1000)) then
               write(7,9020) nint(qq*qq),rke(m,j),rkk*1000,rkkm
               ikstor(m,1)=nint(qq*qq)
               rkstor(m,1)=rkkm
            endif
            if ((rkkm.lt.100).and.(rkkmlast.gt.100)) then
               write(7,9020) nint(qq*qq),rke(m,j),rkk*1000,rkkm
               ikstor(m,2)=nint(qq*qq)
               rkstor(m,2)=rkkm
            endif
            if ((rkkm.lt.10).and.(rkkmlast.gt.10)) then
               write(7,9020) nint(qq*qq),rke(m,j),rkk*1000,rkkm
               ikstor(m,3)=nint(qq*qq)
               rkstor(m,3)=rkkm
            endif
            if ((rkkm.lt.1).and.(rkkmlast.gt.1)) then
               write(7,9020) nint(qq*qq),rke(m,j),rkk*1000,rkkm
               ikstor(m,4)=nint(qq*qq)
               rkstor(m,4)=rkkm
            endif
            write(71,*) qq**2,(rks-rke(m,j))*13.6057
         enddo

 911     continue
      enddo

      close(7)
      close(71)

 9393 format(1x,70('-'))
 9394 format(1x,20x,"Kinetic energy convergence")
 9000 format(1x,5x,30('='),'|',i3,'>',30('='),10x)
 9001 format(1x,5x,30('='),'|',i3,'+>',30('='),10x)
 9002 format(1x,5x,30('='),'|',i3,'->',30('='),10x)
 9020 format(i10,f18.9,10x,f15.6,f15.6)
      return
      end
