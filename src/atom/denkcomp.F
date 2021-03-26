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
      subroutine denkcomp(name,rkstor,ikstor)
      implicit double precision (a-h,o-z)
      
#include "fortdim.h"


      common /grid/ h,r1,z,r(npdm),np
      common /rscore/ rscore(npdm),rdd(npdm),rddd(npdm),
     $     rscoretot(npdm),rsval(npdm)


      character*80 name
      character*80 file_log

      dimension pl(5000),rtemp(npdm),rke(5000),rkstor(3)
      dimension ikstor(3)

      open(unit=7,file=file_log,form='formatted',access='append')

      open(unit=71,file=name,form='formatted')

      rewind(71)

      pi=acos(-1.0)                                                   
      pi4=4.0*pi                                                      
      rnq=sqrt(2.0/pi)


      write(7,*) '-----------------------------------------------------'
      write(7,*) '              PCC KE convergence  '

      write(7,*) "     Ecut(Ry)      dKE(mRy)         dKE(meV) ",
     $     "              KE(Ry)" 

      rks=0.0
      
      qmax=20
      iqp=500
      dqq=qmax/iqp
      
      rsc=0.0
      call radin(r,rscore,0,np,h,rsc)

      do j=1,iqp+1
         rke(j)=0.0
      enddo
      
      do j=1,iqp+1                                                    
         qq=dqq*dble(j-1)                                             
         
         do k=1,maxim                                                
            qr=r(k)*qq                                             
            rtemp(k) = sqrt(rscore(k))*r(k) * besfn(qr,0)
         enddo
         pow=0
         call radin(r,rtemp,0,maxim,h,pow)
         pl(j)=pow *rnq
         
         rks=rks+pl(j)*pl(j)*dqq*qq*qq*qq*qq
         rke(j)=rks
      enddo
      write(71,*) '@'
      
      rkkmlast=0.0

         rkkmlast=0.0
         do j=1,iqp                                                    
            qq=dqq*dble(j-1)                            
            rkk=rks-rke(j)
            rkkmlast=rkkm
            rkkm=rkk*1000*13.6057
            if ((rkkm.lt.100).and.(rkkmlast.gt.100)) then
               write(7,9020) nint(qq*qq),rkk,rkkm,rke(j)
               ikstor(1)=nint(qq*qq)
               rkstor(1)=rkkm
            endif
            if ((rkkm.lt.10).and.(rkkmlast.gt.10)) then
               write(7,9020) nint(qq*qq),rkk,rkkm,rke(j)
               ikstor(2)=nint(qq*qq)
               rkstor(2)=rkkm
            endif
            if ((rkkm.lt.1).and.(rkkmlast.gt.1)) then
               write(7,9020) nint(qq*qq),rkk,rkkm,rke(j)
               ikstor(3)=nint(qq*qq)
               rkstor(3)=rkkm
            endif

            write(71,*) qq,(rks-rke(j))*13.6057
         enddo

      
 911  continue

      write(7,*) '-----------------------------------------------------'
      
      close(7)
      close(71)
 9000 format(1x,10x,'==========|',i3,'>==========',10x)
 9020 format(i10,f18.9,f18.9,f22.9)
      return
      end
