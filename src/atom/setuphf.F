
c------------------------------------------------------------
c
c The following code is adapted from:
c
c ADEA_v1_0: Single configuration Dirac-Fock atom code.
c
c by A.L. Ankudinov, S.I. Zabinsky, J.J. Rehr
c
c Published in Computer Physics Communications 98 (1996) 359.
c
c------------------------------------------------------------

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
      subroutine setuphf(epss,ka,scc,ipsp,zeff)
     
c     *************************************************************************
c     determine the self-consistent solution
c     *************************************************************************
     
      implicit double precision (a-h,o-z)

#include "fortdim.h"

      common /aorb/ ncore,nval,norb,nlm(n0),
     $     no(n0),lo(n0),nmax(n0),maxim
      common /adat/ wnl(n0),en(n0),so(n0),xion
      common /aval/ rcall(n0),rvap(n0),rnorm(n0),ibd(n0),etot
      common /grid/ h,r1,z,r(npdm),np
      common /totpot/ rvcore(npdm,n0),rvps(npdm,n0),rvcoul(npdm)
      common /wfnrel/ rnla(npdm,n0),rnlb(npdm,n0)

      common /hfrel/ eg(npdm),ep(npdm)
      common /hfrel2/ bg(npdm,n0),bp(npdm,n0)
      common /hfrel3/ cg(npdm,n0),cp(npdm,n0)
      common /hfrel4/ fl(n0)
      common /hfrel5/ ceg(10),cep(10),ag(10),ap(10),av(10)
      common /hfrel6/ gg(npdm),gp(npdm)
      common /mulabk/ afgk(30,30,0:30)
      common /scrhf1/ eps(2*npdm),nre(N0),ipl    
      common /cl/ cl

      dimension epss(npdm),ka(n0),scc(n0),anoy(10)

      character*80 file_log


      ndor=10
c     inmuat
      do i=1,np
         epss(i)=0.d0
      enddo

      ch=-xion-1
c      write(7,*) xion,ch,zeff
      ipl=0
      do i=1,norb
         ka(i)=-nint(2*so(i)*(lo(i)+so(i)+0.5))
c         if ((so(i).lt.0.1).and.lo(i).ne.0) ka(i)=-lo(i)
         nre(i)=-1
         llq=abs(ka(i))
         l=llq+llq
         if (ka(i).lt.0) llq=llq-1
         if (llq.lt.0.or.llq.gt.no(i).or.llq.gt.3) stop 'llq'
         nmax(i)=np
         scc(i)=0.1
         if (wnl(i).lt.l) nre(i)=1
         
         do j=1,i-1
            if (ka(j).eq.ka(i)) then
               if (nre(j).gt.0.or.nre(i).gt.0) ipl=ipl+1
            endif
         enddo
         if (ipsp.eq.0) en(i)=-z*z/(2.d0*no(i)*no(i))
         if (ipsp.ne.0) en(i)=en(i)/2.d0
      enddo

c     wfirdf
      r1=r(1)
c     nucdev
      do i=1,ndor
         anoy(i)=0.d0
      enddo
      anoy(1)=-z
      a=(zeff/cl)**2
      do i=1,norb
         fl(i)=sqrt(ka(i)**2 - a)
      enddo
      
      do k=1,norb
         do i=1,np
            rvps(i,k)=(dentfa(r(i),zeff,ch)+rvcore(i,k)/(2.d0*r(i)))/cl
            eg(i)=0.d0
            ep(i)=0.d0
         enddo
      enddo

      

      do i=1,10
         ceg(i)=0.d0
         cep(i)=0.d0
         av(i)=anoy(i)/cl
      enddo
      av(2)=av(2)+dentfa(r(1),zeff,ch)/cl
      testy=1.d-5
      rap1=1.d2
      rap2=1.d1
      test1=testy/rap1
      test2=testy/rap2
      b=test1
      ainf=test1
      nuc=1

c     soldir      
      do j=1,norb

         bg(1,j)=1.d0
         i=no(j)-abs(ka(j))
         if (ka(j).lt.0) i=i-1
         if (mod(i,2).eq.0) bg(1,j)=-bg(1,j)
         if (ka(j).ge.0) then
            bp(1,j)=bg(1,j)*cl*(ka(j)+fl(j))/zeff 
                                  
         else
            bp(1,j)=bg(1,j)*zeff/(ka(j)-fl(j))/cl 
         endif
         method=0 
         call soldir
     $        (j,no(j),lo(j),en(j), fl(j),bg(1,j),bp(1,j),rvps(1,j), 
     $        ainf,test1,test2, ka(j),nmax(j),ifail,method)

c         write(7,'(4i5,4f10.6)') j,ka(j),no(j),lo(j),so(j),wnl(j),en(j)
         do i=1,10
            bg(i,j)=ag(i)
            bp(i,j)=ap(i) 
         enddo
         do i=1,np
            cg(i,j)=gg(i)
            cp(i,j)=gp(i)
         enddo
         nem=0   
      enddo                        

      call hfcoef_rel(ka)

      return
      end




