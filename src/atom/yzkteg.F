
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
      subroutine yzkteg (f,af,g,ag,dr,ap,k,nd,npp,idim)                        
c calculation of yk(r)=zk(r)+ r**(k+1) * integral from r to                     
c   infinity of  f(u) * u**(-k-1)                                               
c zk(r) = r**(-k) * integral from 0 to r of f(u) * u**k                         
                                                                                
c at the origin f(r)=sum from i=1 to nd of af(i)*r**(ap+i-1)                    
c dr tabulation points   h exponential step                                     
c np number of tabulation points for f                                          
c idim dimension of the blocks f,g and dr                                       
                                                                                
c at the origin yk=cte*r**(k+1)-developement limit                              
c the constant for yk lies in ap                                                
c output functions yk and zk lie in f and g, and their                          
c development coefficients at the origin in af and ag.                          
                                                                                
c integration from point to point by a 4 points method.                         
c integral from r to r+h = h*(-f(r-h)+13*f(r)+13*f(r+h)-f(r+h+h))/24            
                                                                                
      implicit double precision (a-h,o-z)                                       
#include "fortdim.h"
      common /grid/ h,r1,z,r(npdm),np
      dimension f(npdm),af(10),g(npdm),ag(10),dr(npdm)                             
                                                                                
c    initialisation and development coefficients of yk                          
      npp= min(npp,idim-2)                                                        
      b=ap                                                                      
      ap=0.0d 00                                                                
      g(1)=0.0d 00                                                              
      g(2)=0.0d 00                                                              
      do 15 i=1,nd                                                              
         b=b+1.0d 00                                                            
         ag(i)=af(i)/(b+k)                                                      
         if (af(i).ne.0.0d 00) then                                             
            c=r(1)**b                                                          
            g(1)=g(1)+ag(i)*c                                                   
            g(2)=g(2)+ag(i)*(r(2)**b)                                          
            af(i)=(k+k+1)*ag(i)/(b-k-1)                                         
            ap=ap+af(i)*c                                                       
         endif                                                                  
 15   continue                                                                  
      do 21 i=1,npp                                                              
 21   f(i)=f(i)*r(i)                                                           
      np1=npp+1                                                                  
      f(np1)=0.0d 00                                                            
      f(np1+1)=0.0d 00                                                          
                                                                                
c     calcualation of zk                                                        
      eh= exp(h)                                                                
      e=eh**(-k)                                                                
      b=h/2.4d+01                                                               
      c=1.3d+01*b                                                               
      ee=e*e*b                                                                  
      b=b/e                                                                     
      do 51 i=3,np1                                                             
 51   g(i)=g(i-1)*e+(c*(f(i)+f(i-1)*e)-(f(i-2)*ee+f(i+1)*b))                    
                                                                                
c     calcualation of yk                                                        
      f(npp)=g(npp)                                                               
      do 61 i=np1,idim                                                          
 61   f(i)=f(i-1)*e                                                             
      i=k+k+1                                                                   
      b=i*b*eh                                                                  
      ee=i*ee/(eh*eh)                                                           
      e=e/eh                                                                    
      c=i*c                                                                     
      do 71  i=npp-1,2,-1                                                        
 71   f(i)=f(i+1)*e+(c*(g(i)+g(i+1)*e)-(g(i+2)*ee+g(i-1)*b))                    
      ee=e*e                                                                    
      c=8.0d 00*c/1.3d+01                                                       
      f(1)=f(3)*ee+c*(g(3)*ee+4.0d 00*e*g(2)+g(1))                              
      ap=(ap+f(1))/(r(1)**(k+1))                                               
      return                                                                    
      end                                                                       
cREF. IN COMP. PHYS. COMMUN. 98 (1996) 359                               ADEA0000
c                                                                        ADEA****
