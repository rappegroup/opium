
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
      subroutine ortdat (ia,kap)                                                    
c        * orthogonalization by the schmidt procedure*                          
c the ia orbital is orthogonalized toa all orbitals of the same                 
c symmetry if ia is positive, otherwise all orbitals of the same                
c symmetry are orthogonalized                                                   
c        this program uses dsordf                                               
                                                                                
      implicit double precision (a-h,o-z)                                       
#include "fortdim.h"
      common /hfrel/ eg(npdm),ep(npdm)
      common /hfrel2/ bg(npdm,n0),bp(npdm,n0)
      common /hfrel3/ cg(npdm,n0),cp(npdm,n0)
      common /hfrel4/ fl(n0)
      common /hfrel5/ ceg(10),cep(10),ag(10),ap(10),av(10)
      common /hfrel6/ gg(npdm),gp(npdm)
      common /grid/ h,r1,z,r(npdm),np
      common /aorb/ ncore,nval,norb,nlm(n0),
     $     no(n0),lo(n0),nmax(n0),maxim
      common /adat/ wnl(n0),en(n0),so(n0),xion
      common /aval/ rcall(n0),rvap(n0),rnorm(n0),ibd(n0),etot
      dimension kap(n0)
      
      ndor=10
      m=norb                                                                    
      l= max(ia,1)                                                              
      if (ia.gt.0) go to 11                                                     
 5    m=l                                                                       
      l=l+1                                                                     
      if (l.gt.norb) go to 999                                                  
 11   do 15 i=1,np                                                            
         gg(i)=0.0d 00                                                          
 15      gp(i)=0.0d 00                                                          
         maxl=nmax(l)                                                              
      do 21 i=1,maxl                                                            
         gg(i)=cg(i,l)                                                          
 21      gp(i)=cp(i,l)                                                          
      do 25 i=1,ndor                                                            
         ag(i)=bg(i,l)                                                          
 25      ap(i)=bp(i,l)                                                          

      do 51 j=1,m                                                               
         if (j.eq.l.or.kap(j).ne.kap(l)) go to 51                               
         max0=nmax(j)                                                           
         a=dsordf (j,j,0,3,fl(j))                                               
         do 41 i=1,max0                                                         
            gg(i)=gg(i)-a*cg(i,j)                                               
 41         gp(i)=gp(i)-a*cp(i,j)                                               
         do 45 i=1,ndor                                                         
            ag(i)=ag(i)-a*bg(i,j)                                               
 45         ap(i)=ap(i)-a*bp(i,j)                                               
         maxl= max(maxl,max0)                                                   
 51   continue                                                                  
      max0= maxl                                                                
      nmax(l)=max0                                                              
      a=dsordf (l,max0,0,4,fl(l))                                               
      a= sqrt(a)                                                                
      do 71 i=1,max0                                                            
         cg(i,l)=gg(i)/a                                                        
 71      cp(i,l)=gp(i)/a                                                        
         write(7,*) 'in ortdat', ag(1),ap(1),a
      do 75 i=1,ndor                                                            
         bg(i,l)=ag(i)/a                                                        
 75      bp(i,l)=ap(i)/a                                                        
      if (ia.le.0) go to 5                                                      
 999  return                                                                    
      end                                                                       

