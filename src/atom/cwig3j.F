
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
      double precision function cwig3j (j1,j2,j3,m1,m2,ient)                    
c        wigner 3j coefficient for integers  (ient=1)                           
c                            or semiintegers (ient=2)                           
c        other arguments should be multiplied by ient                           
                                                                                
      implicit double precision (a-h,o-z)                                       
      save                                                                      
      character*512 slog                                                        
      dimension al(32),m(12)                                                    
      data ini/1/,idim/31/                                                      
c     idim-1 is the largest argument of factorial in calculations               

      m3=-m1-m2                                                                 
      if (ini) 1,21,1                                                           
c        initialisation of the log's of the factorials                          
 1    ini=0                                                                     
      al(1)=0.0d 00                                                             
      do 11 i=1,idim                                                            
         b=i                                                                    
 11      al(i+1)=al(i)+ log(b)                                                  
 21   cwig3j=0.0d 00                                                            
      if (((ient-1)*(ient-2)).ne.0) go to 101                                   
      ii=ient+ient                                                              
c        test triangular inequalities, parity and maximum values of m           
      if (( abs(m1)+ abs(m2)).eq.0.and.mod(j1+j2+j3,ii).ne.0) go to 99          
      m(1)=j1+j2-j3                                                             
      m(2)=j2+j3-j1                                                             
      m(3)=j3+j1-j2                                                             
      m(4)=j1+m1                                                                
      m(5)=j1-m1                                                                
      m(6)=j2+m2                                                                
      m(7)=j2-m2                                                                
      m(8)=j3+m3                                                                
      m(9)=j3-m3                                                                
      m(10)=j1+j2+j3+ient                                                       
      m(11)=j2-j3-m1                                                            
      m(12)=j1-j3+m2                                                            
      do 41 i=1,12                                                              
         if (i.gt.10) go to 31                                                  
         if (m(i).lt.0) go to 99                                                
 31      if (mod(m(i),ient).ne.0) go to 101                                     
         m(i)=m(i)/ient                                                         
         if (m(i).gt.idim) go to 101                                            
 41   continue                                                                  
                                                                                
c  calculate 3j coefficient                                                     
      max0= max(m(11),m(12),0)+1                                                
      min0= min(m(1),m(5),m(6))+1                                               
      isig=1                                                                    
      if (mod(max0-1,2).ne.0) isig=-isig                                        
      c=-al(m(10)+1)                                                            
      do 61 i=1,9                                                               
 61   c=c+al(m(i)+1)                                                            
      c=c/2.0d 00                                                               
      do 71 i=max0,min0                                                         
      j=2-i                                                                     
      b=al(i)+al(j+m(1))+al(j+m(5))+al(j+m(6))+al(i-m(11))+al(i-m(12))          
      cwig3j=cwig3j+isig* exp(c-b)                                              
 71   isig=-isig                                                                
      if (mod(j1-j2-m3,ii).ne.0) cwig3j=-cwig3j                                 


 99   continue
c      write(7,9222) cwig3j,j1,j2,j3,m1,m2,ient
 9222 format('CWIG',f10.6,6i5)

      return
 101     write(slog,'(a,6i5)') 'error in cwig3j ',j1,j2,j3,m1,m2,ient           
      stop                                                                      
      end                                                                       
