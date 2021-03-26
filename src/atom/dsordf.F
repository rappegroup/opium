
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
      double precision function dsordf (i,j,n,jnd,a)                            
c              * calculation of diff. integrals*                                
c        integration by simpson method of the   hg*(r**n)                       
c        hg(l)=cg(l,i)*cg(l,j)+cp(l,i)*cp(l,j)  if jnd=1                        
c        hg=expression above multiplied by  dg  if jnd=-1                       
c        hg(l)=cg(l,i)*cp(l,j)                  if jnd=2                        
c        hg=expression above multiplied by  dg  if jnd=-2                       
c        hg(l)=dg(l)*cg(l,i)+dp(l)*cp(l,j)      if jnd=3                        
c        hg(l)=dg(l)*dg(l)+dp(l)*dp(l)          if jnd=4                        
c        hg is constructed by calling program   if jnd>=5                       
c                  cg(l,i)  large component of the orbital i                    
c                  cp(l,j)  small component of the orbital j                    
c        a is such that dg,dp or hg following the case                          
c        behave at the origin as cte*r**a                                       
c        the integration is made as far as r(j) for jnd>3                      
c                                                                               
c        the development limits at the origin (used for calculation             
c        of integral form 0 to r(1) ) of functions dg,dp and hg are            
c        supposed to be in blocks ag,ap and chg respectively                    
c        this program utilises   aprdev                                         
c                                                                               
      implicit double precision (a-h,o-z)                                       
#include "fortdim.h"
      common /hfrel/ eg(npdm),ep(npdm)
      common /hfrel2/ bg(npdm,n0),bp(npdm,n0)
      common /hfrel3/ cg(npdm,n0),cp(npdm,n0)
      common /hfrel4/ fl(n0)
      common /grid/ h,r1,z,r(npdm),np
      common /hfrel5/ ceg(10),cep(10),ag(10),ap(10),av(10)
      common /hfrel6/ gg(npdm),gp(npdm)

      common /aorb/ ncore,nval,norb,nlm(n0),
     $     no(n0),lo(n0),nmax(n0),maxim
      common /adat/ wnl(n0),en(n0),so(n0),xion
      common /aval/ rcall(n0),rvap(n0),rnorm(n0),ibd(n0),etot

      
      dimension hg(npdm),chg(10)                                                 
      dimension bgi(10),bgj(10),bpi(10),bpj(10)                                 


c      write(7,*) i,dg(100),dp(100)
c        construction of the array hg                                           
      if (jnd.le.3) go to 11                                                    
      max0=j
      b=a                                                                       
      go to 101                                                                 
 11   max0= min(nmax(i),nmax(j)) 
      do  15 l= 1,10                                                          
        bgi(l) = bg(l,i)                                                        
        bgj(l) = bg(l,j)                                                        
        bpi(l) = bp(l,i)                                                        
 15     bpj(l) = bp(l,j)  
      if ( abs(jnd)-2) 21,55,101                                                
 21   do 31 l=1,max0                                                            
 31      hg(l)=cg(l,i)*cg(l,j)+cp(l,i)*cp(l,j)                                  
      do 45 l=1,10                                                            
 45   chg(l)=aprdev(bgi,bgj,l)+aprdev(bpi,bpj,l)                                
      go to 81                                                                  
                                                                                
 55   do 61 l=1,max0                                                            
         hg(l)=cg(l,i)*cp(l,j)                                                  
 61   continue
      do 71 l=1,10                                                            
 71      chg(l)=aprdev(bgi,bpj,l)                                               
 81   b=fl(i)+fl(j)                                                             
      if (jnd.gt.0) go to 301                                                   
                                                                                
      do 85 l=1,max0                                                            
 85      hg(l)=hg(l)*gg(l)                                                      
      do 87 l=1,10                                                            
 87      ap(l)=chg(l)                                                           
      b=b+a                                                                     
      do 95 l=1,10                                                            
 95      chg(l)=aprdev(ap,ag,l)                                                 
      go to 301                                                                 
                                                                                
 101  if (jnd-4) 201,111,301                                                    
 111  do 121 l=1,max0                                                           
         hg(l)=gg(l)*gg(l)+gp(l)*gp(l)                                          
 121     continue

      b=b+b                                                                     
      do 131 l=1,10                                                           
 131     chg(l)=aprdev(ag,ag,l)+aprdev(ap,ap,l)                                 
      go to 301                                                                 
                                                                                
 201  do 221 l=1,max0                                                           
 221     hg(l)=gg(l)*cg(l,i)+gp(l)*cp(l,j)                                      
      b=a+fl(i)                                                                 
      do 241 l=1,10                                                           
 241     chg(l)=aprdev(bgi,ag,l)+aprdev(bpj,ap,l)                               
                                                                                
c        integration of the hg                                                  
 301  dsordf=0.0d 00                                                            
      io=n+1                                                                    
      do 305 l=1,max0                                                           
 305     hg(l)=hg(l)*(r(l)**io)                                                
      do 311 l=2,max0,2                                                         
         dsordf=dsordf+hg(l)+hg(l)+hg(l+1)                                      
 311     continue
      dsordf=h*(dsordf+dsordf+hg(1)-hg(max0))/3.0d 00                          
c        integral from 0 to r(1)                                               
      b=b+n                                                                     
      do 331 l=1,10                                                           
         b=b+1.0d 00                                                            
 331     dsordf=dsordf+chg(l)*(r(1)**b)/b                                      
      return                                                                    
      end                                                                       
      
