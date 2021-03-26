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
c
      subroutine hfcoef_rel(ka)
     
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
      common/scrhf1/eps(2*npdm),nre(N0),ipl    

      dimension ka(n0)

      character*80 file_log

c     muatco
      
      do i=1,30                                                             
         do j=1,30                                                             
            do k=0,3                                                              
               afgk(i,j,k)=0.0d 00                                                       
            enddo
         enddo
      enddo

 601  do i=1,norb                                                           
         if (wnl(i).lt.1e-7) wnl(i)=1e-7
         isumi=wnl(i)
         pisumi=isumi
         dsumi=wnl(i)-isumi
         li= abs(ka(i))*2-1                                                    

         if (lo(i).ne.0) then
            if (ka(i).gt.0) rwgt=float(lo(i))/(2.d0*lo(i)+1.d0)
            if (ka(i).lt.0) rwgt=(lo(i)+1.d0)/(2.d0*lo(i)+1.d0)
            pisumi=rwgt* nint((wnl(i)/rwgt))
            dsumi=rwgt* (wnl(i)/rwgt - nint(wnl(i)/rwgt))
            if (abs(dsumi).lt.1d-10) dsumi=0.d0
         endif
         
         do j=1,i                                                           
         if (wnl(j).lt.1e-7) wnl(j)=1e-7
            isumj=wnl(j)
            pisumj=isumj
            dsumj=wnl(j)-isumj


            lj= abs(ka(j))*2-1                                                 
            if (lo(j).ne.0) then
               if (ka(j).gt.0) rwgt=float(lo(j))/(2.d0*lo(j)+1.d0)
               if (ka(j).lt.0) rwgt=(lo(j)+1.d0)/(2.d0*lo(j)+1.d0)
               pisumj=rwgt* nint(wnl(j)/rwgt)
               dsumj=rwgt* (wnl(j)/rwgt - nint(wnl(j)/rwgt))
               if (abs(dsumj).lt.1d-10) dsumj=0.d0
            endif

            kmax=(li+lj)/2                                                      
            kmin= abs(li-lj)/2                                                  

            if ((ka(i)*ka(j)).lt.0) kmin=kmin+1                               
c calculate a_k(i,j)                                                            
            m=0                                                                 
            if (j.eq.i) then
               m=1                                                     
            endif
            afgk(j,i,0)=afgk(j,i,0)+wnl(i)*(wnl(j)-m)                        

            if (j.ne.i) then
c               c=wnl(j)
               if (dsumi.ne.0.0.and.dsumj.ne.0.0)
     $              afgk(j,i,0)=(dsumi*(pisumi+1.d0)*pisumj+
     $              dsumj*(pisumj+1.d0)*pisumi)
            else
c               c=wnl(i)-1.d0
               if (dsumi.ne.0.0)
     $              afgk(j,i,0)=(pisumi*(wnl(i)+dsumi-1.d0))

            endif

c calculate b_k(i,j)                                                            
            b=afgk(j,i,0)                                                       
            if (j.eq.i) then                                                    
               a=li                                                             
               b=-b*(a+1.0)/a                                               
               kmin = kmin+2                                                    
            endif                                                               
            do k = kmin, kmax,2 
               afgk(i,j,k/2)=b*(cwig3j(li,k*2,lj,1,0,2)**2)                     


c               write(6,'(3i5,3f10.6)') i,j,k,b,
c     $              cwig3j(li,k*2,lj,1,0,2)**2,
c     $              b*(cwig3j(li,k*2,lj,1,0,2)**2)                    
            enddo
         enddo
      enddo

      return
      end




