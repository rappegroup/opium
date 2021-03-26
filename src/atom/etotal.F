
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

      subroutine etotal(kap,eatom,eh,ex)
c combined from original subroutines tabfgk,tabbre,tabrat.                      
c kap quantique  number "kappa"                                                 
c xnel occupation of  orbitales (can be fractional)                             
c en one-electron energies                                                      
c fdrirk function calculating radial integrals rk                               
c akeato angular coefficient for integrals  fk, for the                         
c integrals fk(i;i) gives angular coefficients multiplied by 2                  
c bkeato angular coefficient for integrals  gk                                  
c coul ener(1) direct coulomb interaction                                       
c ech  ener(2) exchange coulomb interaction                                     
c        * average value of the breit hamiltonian *                             
c fdrocc function of the orbitals' occupations.                                 
c bkmrdf is a programm to calculate angular coefficients                        
c ema ener(3) magnetic energy                                                   
c ere ener(4) retardation term                                                  
c        sous programmes utilises akeato,bkeato                                 
c        fdrocc fdrirk bkmrdf                                                   
      
      implicit double precision (a-h,o-z)                                       
#include "fortdim.h"
      common /hfrel/ eg(npdm),ep(npdm)
      common /hfrel2/ bg(npdm,n0),bp(npdm,n0)
      common /hfrel3/ cg(npdm,n0),cp(npdm,n0)
      common /hfrel4/ fl(n0)
      common /hfrel5/ ceg(10),cep(10),ag(10),ap(10),av(10)
      common /hfrel6/ gg(npdm),gp(npdm)
      common/scrhf1/eps(2*npdm),nre(N0),ipl    
      common /aorb/ ncore,nval,norb,nlm(n0),
     $     no(n0),lo(n0),nmax(n0),maxim
      common /adat/ wnl(n0),en(n0),so(n0),xion
      common /aval/ rcall(n0),rvap(n0),rnorm(n0),ibd(n0),etot

      dimension kap(N0)
c      common/itescf/testy,rap(2),teste,nz,norb,norbsc                           
      dimension mk(12),ener(4)                                                  
      dimension cer(17),mbi(9),mii(9),mjj(9)                                    
      common/tabre/cmag(3),cret(3)                                              
      common/inelma/nem                                                         
      iprint=5
      do 10 i = 1,4                                                             
 10   ener(i)=0.0d 00                                                           
      iv=0                                                                      
c       fk  integrales                                                          
      do 40 i=1,norb                                                            
         l= abs(kap(i))-1                                                       
         do 40 j=1,i                                                            
            a=1.0d 00                                                           
            if (j.eq.i) a=a+a                                                   
            m= abs(kap(j))-1                                                    
            kmi=2* min(l,m)                                                     
            k=0                                                                 
 20         iv=iv+1                                                             
            cer(iv)=fdrirk(i,i,j,j,k,kap)                                           
            ener(1)=ener(1)+cer(iv)*akeato(i,j,k)/a                             
c            write(7,*) 'hart', i,j,k,kap(j),iv,ener(1),
c     $           cer(iv),akeato(i,j,k),a
            mk(iv)=k                                                            
            if (iv.lt.3) go to 30                                               
            iv=0                                                                
 30         k=k+2                                                               
            if (k.le.kmi) go to 20                                              
 40   continue                                                                  
      iv=0                                                                      
      if (norb.gt.1) then                                                       
c       gk  integrales                                                          
      do 70 i=2,norb                                                            
         i1=i-1                                                                 
         do 70 j=1,i1                                                           
            l= abs(kap(i))                                                      
            m= abs(kap(j))                                                      
            k= abs(l-m)                                                         
            if ((kap(i)*kap(j)).lt.0) k=k+1                                     
            kmi=l+m-1                                                           
 50         iv=iv+1                                                             
            cer(iv)=fdrirk(i,j,i,j,k,kap)                                           
             ener(2) = ener(2) -cer(iv)*bkeato(i,j,k)                           
            mk(iv)=k                                                            
            if (iv.lt.3) go to 60                                               
            iv=0                                                                
 60         k=k+2                                                               
            if (k.le.kmi) go to 50                                              
 70   continue                                                                  
      endif                                                                     
c                                                                               
      nem=1                                                                     
c       direct  integrales                                                      
      ik=0                                                                      
      do 140 j=1,norb                                                           
         jj=2* abs(kap(j))-1                                                    
         do 140 i=1,j                                                           
            ji=2* abs(kap(i))-1                                                 
            k=1                                                                 
            kma= min(ji,jj)                                                     
 110        ik=ik+1                                                             
            mbi(ik)=k                                                           
            mii(ik)=i                                                           
            mjj(ik)=j                                                           
            cer(ik)=fdrirk(j,j,i,i,k,kap)                                           
            if (i.ne.j) go to 120                                               
            call bkmrdf (j,j,k,kap)                                                 
            ener(3)=ener(3)+(cmag(1)+cmag(2)+cmag(3))*cer(ik)*                  
     1              fdmocc(j,j,kap)/2.0d 00                                         
 120        if (ik.lt.3) go to 130                                              
            ik=0                                                                
 130        k=k+2                                                               
            if (k.le.kma) go to 110                                             
 140  continue                                                                  
      if (norb.gt.1) then                                                       
c       exchange  integrales                                                    
      do 201 j=2,norb                                                           
         lj= abs(kap(j))                                                        
         na=-1                                                                  
         if (kap(j).gt.0) go to 121                                             
         na=-na                                                                 
         lj=lj-1                                                                
 121     jp=j-1                                                                 
         do 201 l=1,jp                                                          
            ll= abs(kap(l))                                                     
            nb=-1                                                               
            if (kap(l).gt.0) go to 131                                          
            nb=-nb                                                              
            ll=ll-1                                                             
 131        b=fdmocc(j,l,kap)                                                       
            nm1= abs(lj+na-ll)                                                  
            nmp1=ll+lj+nb                                                       
            nmm1=ll+lj+na                                                       
            np1= abs(ll+nb-lj)                                                  
            k= min(nm1,np1)                                                     
            kma=max(nmp1,nmm1)                                                  
            if (mod(k+ll+lj,2).eq.0) k=k+1                                      
            nb= abs(kap(j))+ abs(kap(l))                                        
 141        call bkmrdf (j,l,k,kap)                                                 
            do 151 i=1,3                                                        
 151           cer(i)=0.0d 00                                                   
            if (nb.le.k.and.kap(l).lt.0.and.kap(j).gt.0) go to 161              
            cer(1)=fdrirk(l,j,l,j,k,kap)                                            
            cer(2)=fdrirk(0,0,j,l,k,kap)                                            
 161        if (nb.le.k.and.kap(l).gt.0.and.kap(j).lt.0) go to 171              
            cer(3)=fdrirk(j,l,j,l,k,kap)                                            
            if (cer(2).ne.0.0d 00) go to 171                                    
            cer(2)=fdrirk(0,0,l,j,k,kap)                                            
 171        do 185 i=1,3                                                        
               ener(3) =ener(3) +cmag(i)*cer(i)*b                               
 185           ener(4) =ener(4) +cret(i)*cer(i)*b                               
            k=k+2                                                               
            if (k.le.kma) go to 141                                             
 201  continue                                                                  
      endif                                                                     
                                                                                
c       total   energy                                                          
c      eatom = -(ener(1)+ener(2))
c      eatom = -(ener(1)+ener(2))+ener(3)+ener(4)                                

      eh=ener(1)
      ex=ener(2)
      eatom=0.d0
      do 212 j=1,norb                                                           
 212     eatom = eatom + en(j)*wnl(j)                                          
c      if (iprint .ge. 1)  write(16,'(a,f15.6)') 'etot',eatom                  
c      write(7,'(a,f15.6)') 'etot',eatom                                    
c      call wlog(slog)                                                           
c      do 215 i=1,4                                                              
c         if (iprint .ge. 1) write(16,'(a4,1pd18.7)') iner(i),ener(i)            
c         write(7,'(a4,1pd18.7)') iner(i),ener(i)                             
c 215  continue   
c     call wlog(slog)                                                        
      return                                                                    
      end                                                                       

      subroutine bkmrdf(i,j,k,kap)                                                 
c     angular coefficients for the breit term                                   
c i and j are the numbers of orbitals                                           
c k is the value of k in uk(1,2)                                                
c        this programm uses cwig3j                                              
c coefficients for magnetic interaction  are in cmag                            
c and those for retarded term are in cret                                       
c the order correspond to -1 0 and +1                                           
                                                                                
      implicit double precision (a-h,o-z)                                       
#include "fortdim.h"
c      common/ratom1/xnel(N0),en(N0),scc(N0),scw(N0),sce(N0),                    
c     1nq(N0),kap(N0),nmax(N0)                                                   
      common/tabre/cmag(3),cret(3)                                              
      dimension kap(n0)

      do 12 l=1,3                                                               
        cmag(l)=0.0d 00                                                         
 12     cret(l)=0.0d 00                                                         
      ji=2* abs(kap(i))-1                                                       
      jj=2* abs(kap(j))-1                                                       
      kam=kap(j)-kap(i)                                                         
      l=k-1                                                                     
      do 51 m=1,3                                                               
         if (l.lt.0) go to 51                                                   
         a=cwig3j(ji,jj,l+l,-1,1,2)**2                                          
         if (a.eq.0.0d 00) go to 51                                             
         c=l+l+1                                                                
         if (m-2) 14,16,17                                                      
 14      cm=(kam+k)**2                                                          
         cz=kam*kam-k*k                                                         
         cp=(k-kam)**2                                                          
         n=k                                                                    
 15      l1=l+1                                                                 
         am=(kam-l)*(kam+l1)/c                                                  
         az=(kam*kam+l*l1)/c                                                    
         ap=(l+kam)*(kam-l1)/c                                                  
         d=n*(k+k+1)                                                            
         go to 31                                                               
                                                                                
 16      d=k*(k+1)                                                              
         cm=(kap(i)+kap(j))**2                                                  
         cz=cm                                                                  
         cp=cm                                                                  
         go to 41                                                               
                                                                                
 17      cm=(kam-l)**2                                                          
         cz=kam*kam-l*l                                                         
         cp=(kam+l)**2                                                          
         n=l                                                                    
         c=-c                                                                   
         go to 15                                                               
                                                                                
 31      c= abs(c)*d                                                            
         if (c.ne.0.0d 00) c=n/c                                                
         cret(1)=cret(1)+a*(am-c*cm)                                            
         cret(2)=cret(2)+(a+a)*(az-c*cz)                                        
         cret(3)=cret(3)+a*(ap-c*cp)                                            
 41      if (d.eq.0.0d 00) go to 51                                             
         a=a/d                                                                  
         cmag(1)=cmag(1)+cm*a                                                   
         cmag(2)=cmag(2)+cz*(a+a)                                               
         cmag(3)=cmag(3)+cp*a                                                   
 51      l=l+1                                                                  
      return                                                                    
      end                                                                       


      double precision function fdmocc (i,j,kap)                                    
c     product of the occupation numbers of the orbitals i and j                 
                                                                                
      implicit double precision (a-h,o-z)                                       
#include "fortdim.h"
c      common/ratom1/xnel(N0),en(N0),scc(N0),scw(N0),sce(N0),                    
c     1nq(N0),kap(N0),nmax(N0)
      common /adat/ wnl(n0),en(n0),so(n0),xion
      dimension kap(n0)
                                                                                
      if (j.eq.i) then                                                          
         fdmocc=wnl(i)*(wnl(j)-1)                                             
         a=2* abs(kap(i))                                                       
         fdmocc=fdmocc*a/(a-1.0)                                                
      else                                                                      
         fdmocc=wnl(i)*wnl(j)                                                 
      endif                                                                     
      return                                                                    
      end                                                                       
