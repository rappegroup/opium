
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
      double precision function aprdev (a,b,l)                                  
c     the result of this function is the coefficient of the term of             
c     power for the product of two polynomes, whose coefficients are            
c     in rows a and b                                                           
                                                                                
      implicit double precision (a-h,o-z)                                       
      dimension a(10),b(10)                                                     
                                                                                
      aprdev=0.0d 00                                                            
      do 11 m=1,l                                                               
 11      aprdev=aprdev+a(m)*b(l+1-m)                                            
      return                                                                    
      end                                                                       

