
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
      double precision function dentfa (dr,dz,ch)                               
c     analitical approximation of potential is created for electrons in         
c     thomas-fermi model for atom or free ion. dr distance from nucleus         
c     with charge dz                                                            
c        ch=ionicity = number of electrons-dz-1                                 
      implicit double precision (a-h,o-z)                                       
                                                                                
      dentfa=0.0d 00                                                            
      if ((dz+ch).lt.1.0d-04) return                                            
      w=dr*(dz+ch)**(1./3.)                                                     
      w=sqrt(w/0.8853)                                                          
      t=w*(0.60112*w+1.81061)+1.                                                
      w=w*(w*(w*(w*(0.04793*w+0.21465)+0.77112)+1.39515)+1.81061)+1             
      dentfa=(dz+ch)*(1.0d 00-(t/w)**2)/dr                                      
      return                                                                    
      end                                                                       

