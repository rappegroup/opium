
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
      subroutine norm(b,hp,dr,gg,gp,ag,ap,method,hx,ndor,                       
     1 gpmat,fl,max0,mat)                                                       
c     calculate norm b. this part of original code was used twice,              
c     causing  difficult block structure. so it was rearranged into             
c     separate subroutine. ala                                                  
      implicit double precision (a-h, o-z)                                      
#include "fortdim.h"
      dimension hp(npdm),dr(npdm),gg(npdm),
     $     gp(npdm),ag(10),ap(10)                   
                                                                                
      b=0.0d 00                                                                 
      do 311 i=1,max0                                                           
 311  hp(i)=dr(i)*(gg(i)*gg(i)+gp(i)*gp(i))                                     
      if (method.ne.1) go to 315                                                
      hp(mat)=hp(mat)+dr(mat)*(gpmat**2-gp(mat)**2)/2.0d 00                     
 315  do 321 i=2,max0,2                                                         
         b=b+hp(i)+hp(i)+hp(i+1)                                                            
 321  continue
c      write(7,*) 'b1',b,max0,hp(400)
      b=hx*(b+b+hp(1)-hp(max0))/3.0d 00                                         
      do 325 i=1,ndor                                                           
         g=fl+fl+i                                                              
         g=(dr(1)**g)/g                                                         
         do 325 j=1,i                                                                       
            b=b+ag(j)*g*ag(i+1-j)+ap(j)*g*ap(i+1-j)                                            
 325     continue
c         write(7,*) 'b2',b,ndor
      return                                                                    
      end                                                                       

