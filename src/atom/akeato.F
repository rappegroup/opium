
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
      double precision function akeato (i,j,k)                                  
c     angular coefficient by the direct coulomb integral fk                     
c     for orbitals i and j                                                      
                                                                                
      implicit double precision (a-h,o-z)                                       
      common/mulabk/afgk(30,30,0:3)                                                 
                                                                                
c     afgk angular coefficients by integrales fk and gk                         
c        coefficient of integral fk(i;j) is in  afgk(min,max)                   
c        and that of integral gk(i;j) is in  afgk(max,min)                      
c        max=max(i,j) min=min(i,j)                                              
                                                                                
      if (i .le. j) then                                                        
         akeato=afgk(i,j,k/2)                                                   
      else                                                                      
         akeato=afgk(j,i,k/2)                                                   
      endif                                                                     
      return                                                                    
                                                                                
      entry bkeato (i,j,k)                                                      
c angular coefficient at the exchange coulomb integral gk                       
                                                                                
      bkeato=0.0d 00                                                            
      if (i .lt. j) then                                                        
         bkeato=afgk(j,i,k/2)                                                   
      elseif (i.gt.j) then                                                      
         bkeato=afgk(i,j,k/2)                                                   
      endif                                                                     
      return                                                                    
      end                                                                       

