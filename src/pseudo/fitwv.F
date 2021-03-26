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
      subroutine fitwv
      
c     *************************************************************************
c     this routine determines the grid index closest to rc and the
c     value, slope and del^2 of the wavefunction at that grid point
c     by fitting to polynomials of increasing order
c     *************************************************************************

      implicit double precision (a-h,o-z)    
                                     
#include "fortdim.h"

      common /grid/ h,r1,z,r(npdm),np

      common /aorb/ ncore,nval,norb,nlm(n0),
     $     no(n0),lo(n0),nmax(n0),maxim
      common /adat/ wnl(n0),en(n0),so(n0),xion
      common /aval/ rcall(n0),rvap(n0),rnorm(n0),ibd(n0),etot

      common /wfn/ rnl(npdm,n0)

      common /cuts/ qc,rc
      common /angm/ ll

      common /wavrc/ wavrc, slope, curvrc, indrc
      common /nnn/ nnn

      dimension psitem(npdm)                                            

c     *************************************************************************
c     find where rc falls in the r array                                  
c     *************************************************************************

      two = 2.0
      indrc = 1                                                         
      do i = 1,np                                                    
         rnl(i,nnn)=rnl(i,nnn)/r(i)
         psitem(i) = rnl(i,nnn)                                           
         if (rc.gt.r(i)) indrc = i
      enddo

      if (two*rc.gt.r(indrc)+r(indrc+1)) indrc = indrc + 1              

      write(7,9004) indrc
      write(7,9005) rc
      write(7,9006) r(indrc) 
      write(7,*)
      rc = r(indrc)                                                     
      
c     *************************************************************************
c     determine slope, curvature, etc. in order to fit rnl()
c     *************************************************************************
c     the derivatives in this section use the logarithmic grid properties.  

      slope = (rnl(indrc+1,nnn) - rnl(indrc-1,nnn))/two/h/r(indrc)

      temp =  rnl(indrc+1,nnn) + rnl(indrc-1,nnn) - two * rnl(indrc,nnn)      
      temp = temp/h/h/r(indrc)/r(indrc)                                 
      curvrc = temp - slope/r(indrc)
      
      wavrc = rnl(indrc,nnn)
      
      x2old = slope
      x3old = curvrc
      diff2o = 1.0e-10
      diff3o = 1.0e-10
      
      do m = 3,10
        x2 = val2(psitem,r,np,rc,m)                                    
        x3 = val3(psitem,r,np,rc,m)                                    
        diff2 = abs(x2-x2old)
        diff3 = abs(x3-x3old)
        if (diff2.lt.diff2o) slope = x2
        if (diff3.lt.diff3o) curvrc = x3
        x2old = x2                                                     
        x3old = x3                                                     
        diff2o = diff2
        diff3o = diff3
      enddo
                                                           
      write(7,9010) rnl(indrc,nnn)
      write(7,9011) slope
      write(7,9012) curvrc
      curvrc = curvrc + 2.0/rc*slope - float(ll*(ll+1))/rc/rc*wavrc

      write(7,9013) curvrc

      if (rnl(indrc,nnn).ne.0.0d0) then                                   
         write(7,9014) slope/rnl(indrc,nnn)            
      endif                                                             

 9004 format(1x,'point nearest rc :',i5)
 9005 format(1x,'rc               :',f10.6)
 9006 format(1x,'actual rc        :',f10.6)

 9010 format(1x,'Psi(rc)          : ',f16.10) 
 9011 format(1x,'Slope at rc      : ',f16.10)
 9012 format(1x,'Curvature at rc  : ',f16.10)
 9013 format(1x,'Del^2 Psi at rc  : ',f16.10)
 9014 format(1x,'Log Deriv at rc  : ',f16.10)

      return                                                            
      end                                                               


c     #########################################################################


      function val2(f,x,n,z,m)                                          
      
c     *************************************************************************
c     interpolate by the lagrange formula up to order m.  the function      
c     f is defined at the points f(x(i)), and f and x are dimensioned       
c     for n.  the function's derivative at abscissa z is returned.      
c     *************************************************************************

      implicit double precision(a-h,o-z)                                          

c#include "fortdim.h"      
      
c step 1:  find where z is located in the x array.                      
      dimension f(1), x(1)                                              
      ind = 0                                                           
      do 1 i = 1,n                                                      
         if (x(i).lt.z) ind = i                                         
 1    continue                                                          
      if (ind.gt.n-m/2-1) ind = n - m/2 - 1                             
      if (ind.lt.m/2) ind = m/2                                         
      ind = ind - m/2 + 1                                               
      val2 = 0.0d0                                                      
      do 2 i = ind,ind+m-1                                              
         prod = f(i)                                                    
         xi = x(i)                                                      
         do 3 j = ind,ind+m-1                                           
            if (j.ne.i) then                                            
               prod = prod/(xi - x(j))                                  
            endif                                                       
 3       continue                                                       
         do 4 k = ind,ind+m-1                                           
            prodj = 1.0d0                                               
            if (k.eq.i) goto 4                                          
            do 5 j = ind,ind+m-1                                        
               if (j.ne.i.and.j.ne.k) then                              
                  prodj = prodj * (z - x(j))                            
               endif                                                    
 5          continue                                                    
            val2 = val2 + prod * prodj                                  
 4       continue                                                       
 2    continue                                                          
      return                                                            
      end                                                               
      
      
c     #########################################################################

      
      function val3(f,x,n,z,m)                                          
      
c     *************************************************************************
c     interpolate by the lagrange formula up to order m.  the function      
c     f is defined at the points f(x(i)), and f and x are dimensioned       
c     for n.  the function's second derivative at abscissa z is returned.
c     *************************************************************************

      implicit double precision(a-h,o-z)                                          

c#include "fortdim.h"
      
c step 1:  find where z is located in the x array.                      
      dimension f(1), x(1)                                              
      ind = 0                                                           
      do 1 i = 1,n                                                      
         if (x(i).lt.z) ind = i                                         
 1    continue                                                          
      if (ind.gt.n-m/2-1) ind = n - m/2 - 1                             
      if (ind.lt.m/2) ind = m/2                                         
      ind = ind - m/2 + 1                                               
      val3 = 0.0d0                                                      
      do 2 i = ind,ind+m-1                                              
         prod = f(i)                                                    
         xi = x(i)                                                      
         do 3 j = ind,ind+m-1                                           
            if (j.ne.i) then                                            
               prod = prod/(xi - x(j))                                  
            endif                                                       
 3       continue                                                       
         do 4 k = ind,ind+m-1                                           
            if (k.eq.i) goto 4                                          
            do 6 l = ind,ind+m-1                                        
               if (l.eq.k.or.l.eq.i) goto 6                             
               prodj = 1.0d0                                            
               do 5 j = ind,ind+m-1                                     
                  if (j.ne.i.and.j.ne.k.and.j.ne.l) then                
                     prodj = prodj * (z - x(j))                         
                  endif                                                 
 5             continue                                                 
               val3 = val3 + prod * prodj                               
 6          continue                                                    
 4       continue                                                       
 2    continue                                                          
      return                                                            
      end                                                               
      
      
c     #########################################################################

      
      function val(f,x,n,z,m)                                           

c     *************************************************************************
c     interpolate by the lagrange formula up to order m.  the function      
c     f is defined at the points f(x(i)), and f and x are have n values
c     tabulated.  the function's value at abscissa z is returned.  
c     *************************************************************************

      implicit double precision(a-h,o-z)  

c#include "fortdim.h"
      
c step 1:  find where z is located in the x array.                      
      dimension f(1), x(1)                                              
      ind = 0                                                           
      do 1 i = 1,n                                                      
         if (x(i).lt.z) ind = i                                         
 1    continue                                                          
      if (ind.gt.n-m/2-1) ind = n - m/2 - 1                             
      if (ind.lt.m/2) ind = m/2                                         
      ind = ind - m/2 + 1                                               
      val = 0.0d0                                                       
      der1 = 0.0d0                                                      
      der2 = 0.0d0                                                      
      do 2 i = ind,ind+m-1                                              
         prod = f(i)                                                    
         xi = x(i)                                                      
         prodj = 1.0d0                                                  
         do 3 j = ind,ind+m-1                                           
            if (j.ne.i) then                                            
               prod = prod/(xi - x(j))                                  
               prodj = prodj * (z - x(j))                               
            endif                                                       
 3       continue                                                       
         val = val + prod * prodj                                       
 2    continue                                                          

      return                                                            
      end                                                               
