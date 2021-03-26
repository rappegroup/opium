
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
      subroutine potrdf (rvn,rvs,ia,kap,method)
                                                    
c        this programm uses akeato(bkeato),aprdev,multrk,yzkrdf                 
      implicit double precision (a-h,o-z)                                       
#include "fortdim.h"
      common/scrhf1/eps(npdm*2),nre(N0),ipl                                        

      common /grid/ h,r1,z,r(npdm),np
      common /aorb/ ncore,nval,norb,nlm(n0),
     $     no(n0),lo(n0),nmax(n0),maxim
      common /adat/ wnl(n0),en(n0),so(n0),xion
      common /aval/ rcall(n0),rvap(n0),rnorm(n0),ibd(n0),etot

      common /hfrel/ eg(npdm),ep(npdm)
      common /hfrel2/ bg(npdm,n0),bp(npdm,n0)
      common /hfrel3/ cg(npdm,n0),cp(npdm,n0)
      common /hfrel4/ fl(n0)
      common /hfrel5/ ceg(10),cep(10),ag(10),ap(10),av(10)
      common /hfrel6/ gg(npdm),gp(npdm)
      common /mulabk/ afgk(30,30,0:30)                                                        
      common /cl/ cl
      common /xg/ xg(npdm,n0)

      dimension bgj(10),bpj(10),rjunk(npdm)                                                 
      dimension at(npdm),bt(npdm),kap(n0)
      dimension rvs(npdm),rvn(npdm),anoy(10)

      do i=1,10
         anoy(i)=0.d0
      enddo
      anoy(1)=-z
                                                                                
      do 9 i=1,10                                                             
         cep(i)=0.0d 00                                                         
         ceg(i)=0.0d 00                                                         
 9       av(i)=anoy(i)                                                          
      do 11 i=1,np                                                            
         at(i)=0.0d 00                                                          
         bt(i)=0.0d 00                                                          
         ep(i)=0.0d 00                                                          
         eg(i)=0.0d 00                                                          
 11      rvs(i)=0.0d 00                                                          
                                                                                
c     coulomb terms                                                             
      jia=2* abs(kap(ia))-1                                                     
      k=0                                                                       
 21   do 25 i=1,np                                                            
 25   gg(i)=0.0d 00                                                             
      do 31 i=1,10                                                            
 31   ag(i)=0.0d 00                                                             
      max0=0                                                                    
      do 51 j=1,norb                                                            
         do 33 i = 1,10                                                         
            bgj(i) = bg(i,j)                                                    
 33         bpj(i) = bp(i,j)                                                    
         m=2* abs(kap(j))-1                                                     
         if (k.gt.m) go to 51     

          a=akeato(ia,j,k)/wnl(ia)                                              
         
c         write(6,8333) ia,j,k,a*wnl(ia),wnl(ia),akeato(ia,j,k)
 8333    format('hart ', 3i5,3f10.6)

         if (a.eq.0.0d 00) go to 51                                             
         m=nmax(j)                                                              
         do 35 i=1,m                                                            
            gg(i)=gg(i)+a*(cg(i,j)*cg(i,j)+cp(i,j)*cp(i,j))                     
 35      continue
         n=2* abs(kap(j))-k                                                     
         l=10+2-n                                                             
         if (l.le.0) go to 51                                                   
         do 41 i=1,l                                                            
            m=n-2+i                                                             
 41         ag(m)=ag(m)+a*(aprdev(bgj,bgj,i)+                                   
     1            aprdev(bpj,bpj,i))                                            
 51   max0= max(max0,nmax(j))                                                
      call yzkrdf (0,max0,k,cg,cp)                                                    
      do 61 i=1,10                                                            
         l=k+i+3                                                                
         if (l.gt.10) go to 61                                                
         av(l)=av(l)-ag(i)                                                      
 61   continue                                                                  
      do 81 i=1,np                                                            
         rvs(i)=rvs(i)+gg(i)
         rjunk(i)=cg(i,ia)
 81   continue
c     power of r?

      k=k+2                                                                     
      if (k.le.10) av(k)=av(k)+ap(1)                                          
      if (k.lt.jia) go to 21                                                    

c     exchange terms                                                            
      if (method.eq.0) go to 411                                                
      do 201 j=1,norb                                                           

         if (j-ia) 105,201,105                                                  
 105     max0=nmax(j)                                                           
         jj=2* abs(kap(j))-1                                                    
         kma=(jj+jia)/2                                                         
         k= abs(jj-kma)                                                         
         if ((kap(j)*kap(ia)).lt.0) k=k+1                                       
                                                                                
 111     a=bkeato(j,ia,k)/wnl(ia)  
c         write(8,8334) ia,j,k,a*wnl(ia),wnl(ia),a
 8334    format(' xc  ', 3i5,3f10.6)
         if (a.eq.0.0d 00) go to 151                                            
         call yzkrdf (j,ia,k,cg,cp)                                                   
         do 121 i=1,max0                                                        
            eg(i)=eg(i)+a*gg(i)*cg(i,j)                                         
            ep(i)=ep(i)+a*gg(i)*cp(i,j)                                         
 121     continue
         n=k+1+ abs(kap(j))- abs(kap(ia))                                       
         if (n.gt.10) go to 141                                               
         do 135 i=n,10                                                        
            ceg(i)=ceg(i)+bg(i+1-n,j)*a*ap(1)                                   
 135        cep(i)=cep(i)+bp(i+1-n,j)*a*ap(1)                                   
 141     i=2* abs(kap(j))+1                                                     
         if (i.gt.10) go to 151                                               
         do 143 i = 1,10                                                        
            bgj(i) = bg(i,j)                                                    
 143        bpj(i) = bp(i,j)                                                    
         do 145 n=i,10                                                        
            ceg(n)=ceg(n)-a*aprdev(ag,bgj,n+1-i)                                
 145        cep(n)=cep(n)-a*aprdev(ag,bpj,n+1-i)                                
 151     k=k+2                                                                  
         if (k.le.kma) go to 111                                                
 201  continue                                                                  

 411  if (ipl.eq.0) go to 511                                                   
      do 481 j=1,norb                                                        
         if (kap(j).ne.kap(ia).or.j.eq.ia) go to 481                            
         if (nre(j).lt.0.and.nre(ia).lt.0) go to 481                            
         m= max(j,ia)                                                           
         i= min(j,ia)+((m-1)*(m-2))/2                                           
         a=eps(i)*wnl(j)                                                       
         max0=nmax(j)                                                           
         do 461 i=1,max0                                                        
            at(i)=at(i)+a*cg(i,j)                                               
 461        bt(i)=bt(i)+a*cp(i,j)                                               
         do 471 i=1,10                                                        
            ceg(i)=ceg(i)+bg(i,j)*a                                             
 471        cep(i)=cep(i)+bp(i,j)*a                                             
 481  continue                                                                  
                                                                                
c addition of nuclear potential and division of potentials and                  
c       their development limits by speed of light                              

 511  do 527 i=1,10                                                           
         av(i)=av(i)/cl                                                         
         cep(i)=cep(i)/cl                                                       
 527     ceg(i)=ceg(i)/cl                                                       


      do 531 i=1,np                                                           

         rvs(i)=(rvs(i)/r(i)+rvn(i)/(2.d0*r(i)))/cl
         ep(i)=(ep(i)+bt(i)*r(i))/cl                                           
         eg(i)=(eg(i)+at(i)*r(i))/cl                                           

         xg(i,ia)=0.d0
         if (abs(cg(i,ia)).gt.1e-15) then
            xg(i,ia)=-2.d0*eg(i)*cl
         endif

c     EJW - important stuff ->

c         x3old = 0.0
c         diff3o = 1.0e-8
c         do m = 3,10
c            x3 = val3(rjunk,r,np,r(i),m)
c            diff3 = abs(x3-x3old)
c            if (diff3.lt.diff3o) ttpsi = x3
c            x3old = x3
c            diff3o = diff3
c         enddo

c        write(870+ia,*) r(i),ttpsi

c         if (r(i).gt.0.5) write(880+ia,9222) r(i),rvs(i)*cg(i,ia)*cl,
c     $        -en(ia)*cg(i,ia),-0.5*ttpsi,cg(i,ia)*kap(ia)
c     $        *(kap(ia)+1)/(2.d0*r(i)**2), - eg(i)*cl/r(i)
c 9222    format(8f20.10)

c         if (r(i).lt.4) write(850+ia,*) r(i),rvs(i)*cl
c     $        +eg(i)*cl/(r(i)*cg(i,ia))

 531  continue
      return                                                                    
      end                                                                       

