
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
      subroutine cofcon (a,b,p,q)                                               
c         acceleration of the convergence in the iterative process              
c b is the part of final iteration n is a function of the error (p)             
c (p) at iteration n and the error (q) at the iteration n-1.                    
c if the product p*q is positive  b is increased by 0.1                         
c                        zero b is unchanged                                    
c                        negative b is decreased by 0.1                         
c b is between 0.1 and 0.9                                                      
c                a = 1. - b                                                     
c             ** at the end makes q=p                                           
c                                                                               
      implicit double precision (a-h,o-z)                                       
                                                                                
      if (p*q)  11,31,21                                                        
 11   if (b .ge. 0.2) b = b - 0.1
      go to 31                                                                  
                                                                                
 21   if (b .le. 0.9) b = b + 0.1                                               
                                                                                
 31   a = 1.0 - b                                                               
      q=p           
      
      return                                                                    
      end                                                                       

