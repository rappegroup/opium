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
      subroutine intdir(gg,gp,ag,ap,ggmat,gpmat,en,dv,av,
     $     fl,agi,api,ainf,test1,test2,imax,mat,emin,esup,ka,imm)        
c            solution of the inhomogenios dirac equation                        
c gg gp initially exchage terms, at the time of return - wave functions         
c ag and ap development coefficients of  gg and gp                              
c ggmat gpmat  values at the matching point for the inward integration          
c en one-electron energy                                                        
c fl power of the first development term at the origin                          
c agi (api) initial values of the first development coefficients                
c at the origin of a large (small) component                                    
c ainf initial value for large component at point dr(max0)                      
c   - at the end of tabulation of gg gp                                         
      implicit double precision (a-h,o-z)                                       
#include "fortdim.h"  
      save                                                                      
      common /grid/ h,r1,z,r(npdm),np
      common /cl/ cl

      dimension gg(npdm),gp(npdm),ag(10),ap(10),
     $     coc(5),cop(5),dg(5),dp(5),av(10)
      dimension dv(npdm)

      data cop/2.51d+02,-1.274d+03,2.616d+03,-2.774d+03,1.901d+03/,             
     $     coc/-1.9d+01,1.06d+02,-2.64d+02,6.46d+02,2.51d+02/,                       
     $     cmixn/4.73d+02/,cmixd/5.02d+02/,hxd/7.2d+02/,npi/5/,icall/0/              
                                                                                
c numerical method is a 5-point predictor-corrector method                      
c predicted value    p(n) = y(n-1) + c * somme de i=1,5 cop(i)*y'(n-i)          
c corrected value    c(n) = y(n-1) + c * somme de i=1,4 coc(i)*y'(n-i)          
c                                  + coc(5)*p'(n)                               
c final value        y(n) = cmix*c(n) + (1.-cmix)*p(n)                          
c                           cmix=cmixn/cmixd                                    

c      write(7,*) 'in intdir', test1
      ndor=10
      fk=ka
      ccl=cl+cl
      ell=fk*(fk+1.0d 00)/ccl                                                   
      if (icall.eq.0) then                                                      
         icall=1                                                                
         c=cmixn/cmixd                                                          
         a=1.0d 00-c                                                            
         cmc=c*coc(5)                                                           
         f=coc(1)                                                               
         do 1 j=2,npi                                                           
            g=coc(j)                                                            
            coc(j)=c*f+a*cop(j)                                                 
 1          f=g                                                                 
         coc(1)=c*cop(1)                                                        
      endif                                                                     
      c=h/hxd                                                                  
      ec=en/cl                                                                  
      ag(1)=agi                                                                 
      ap(1)=api                                                                 
c      do i=1,10
c         write(6,*) 'av',av(i)
c      enddo
c      write(7,*) 'imm=',imm
      if (imm) 81,15,26                                                         
c      search for the second sign change point                                  
 15   continue
      mat=npi                                                                
      j=1                                                                       
 16   continue

      mat=mat+2                                                                 
      if (mat.ge.imax) then    
         write(7,*) 'mat,imax',mat,imax
c     i had trouble with screened k-hole for la, for f-electrons.                 
c     below i still define matching point if one electron energy is               
c     not less than -1ev. ala, january 1995                                       
         if (ec .gt. -0.0003) then                                           
            mat = imax - 12                                                     
            go to 25                                                          
         endif                                                               
         stop 'problem'
      endif   
      f=dv(mat)+ell/(r(mat)*r(mat))                                        
      f=(f-ec)*j                                                             
c      write(7,*) 'ff',f,dv(mat),ec
      if (f) 25,25,16                                                        
 25   continue
c      write(7,*) 'At 25'
      j=-j                                                                   
      if (j.lt.0) go to 16                                                      
      if (mat .ge. imax-npi) mat=imax-12                                            
c      write(7,*) 'out of loop', mat,r(mat)
c     initial values for the outward integration                                
 26   do 35 j=2,ndor                                                            
         k=j-1                                                                  
         a=fl+fk+k                                                              
         b=fl-fk+k                                                              
         ep=a*b+av(1)*av(1)                                                     
         f=(ec+ccl)*ap(k)+ap(j)                                                 
         g=ec*ag(k)+ag(j)                                                       
         do 31 i=1,k          
c            write(7,*) 'hey',i,j,av(i+1),ap(j-i)                                  
            f=f-av(i+1)*ap(j-i)                                                 
 31         g=g-av(i+1)*ag(j-i)                                                 
                                                                                
         ag(j)=(b*f+av(1)*g)/ep                                                 
         ap(j)=(av(1)*f-a*g)/ep                                                 
 35      continue
      do 41 i=1,npi                                                             
         gg(i)=0.0d 00                                                          
         gp(i)=0.0d 00                                                          
         dg(i)=0.0d 00                                                          
         dp(i)=0.0d 00                                                          
         do 41 j=1,ndor                                                         
c            write(7,*) 'ndor',j,a,ag(j),ap(j)
            a=fl+j-1                                                            
            b=r(i)**a                                                          
            a=a*b*c                                                             
            gg(i)=gg(i)+b*ag(j)                                                 
            gp(i)=gp(i)+b*ap(j)                                                 
            dg(i)=dg(i)+a*ag(j)                                                 
 41         dp(i)=dp(i)+a*ap(j)                                                 
      i=npi                                                                     
      k=1                                                                       
      ggmat=gg(mat)                                                             
      gpmat=gp(mat)                                                             

c      write(7,*) 'after the start'
c     integration of the inhomogenious system                                   
 51   cmcc=cmc*c                                                                
                                                                                
c      write(7,*) 'At 55',dg(1),dp(1)
 55   continue                                                                  
         a=gg(i)+dg(1)*cop(1)                                                   
         b=gp(i)+dp(1)*cop(1)                                                   
         i=i+k                                                                  
         ep=gp(i)                                                               
         eg=gg(i)                                                               
         gg(i)=a-dg(1)*coc(1)                                                   
         gp(i)=b-dp(1)*coc(1)                                                   
         do 61 j=2,npi                                                          
            a=a+dg(j)*cop(j)                                                    
            b=b+dp(j)*cop(j)                                                    
            gg(i)=gg(i)+dg(j)*coc(j)                                            
            gp(i)=gp(i)+dp(j)*coc(j)                                            
            dg(j-1)=dg(j)                                                       
            dp(j-1)=dp(j)                                                       
c            if ((i+10.gt.mat).or.(i-10.lt.npi))
c     $           write(7,*) i,j,gg(i),dg(j),cop(j),coc(j)
 61         continue
         f=(ec-dv(i))*r(i)                                                     
         g=f+ccl*r(i)                                                          
         gg(i)=gg(i)+cmcc*(g*b-fk*a+ep)                                         
         gp(i)=gp(i)+cmcc*(fk*b-f*a-eg)                                         
         dg(npi)=c*(g*gp(i)-fk*gg(i)+ep)                                        
         dp(npi)=c*(fk*gp(i)-f*gg(i)-eg)                                        
      if (i.ne.mat) go to 55    

c      write(7,*) 'gg',ec,dv(mat),r(mat),ccl
                                                                                
      if (k.lt.0) go to 999                                                     
      a=ggmat                                                                   
      ggmat=gg(mat)                                                             
      gg(mat)=a                                                                 
      a=gpmat                                                                   
      gpmat=gp(mat)                                                             
      gp(mat)=a                                                                 
      if (imm.ne.0) go to 81                                                    
                                                                                
c     initial values for inward integration                                     
c      write(7,*) 'init a test1 ggmat,ainf',a,test1,abs(ggmat),ainf
      a=test1* abs(ggmat)                                                       
      if (ainf.gt.a) ainf=a                                                     
c      imax=np+2                                                            
      imax=np+1
      if (mod(np+1,2).eq.0) imax=np
 73   a=7.0d+02/cl                                                              
 75   imax=imax-2                                                               
         if ((imax+1).le.(mat+npi)) then                                        
            stop 'problem2'
c            numerr=138021                                                       
c          *the last tabulation point is too close to the matching point        
c            return                                                              
         endif                                                                  
      if (((dv(imax)-ec)*r(imax)*r(imax)).gt.a) go to 75                      
                                                                                
 81   continue
      c=-c                                                                      
      a=- sqrt(-ec*(ccl+ec))                                                    
      if ((a*r(imax)).lt.-1.7d+02) go to 73                                    
      b=a/(ccl+ec)                                                              
      f=ainf/ exp(a*r(imax))                                                   
      if (f.eq.0.0d 00) f=1.0d 00                                               
      do 91 i=1,npi                                                             
c         write(7,*) 'at 91',f,ainf,imax,a,b
         j=imax+1-i                                                             
         gg(j)=f* exp(a*r(j))                                                  
         gp(j)=b*gg(j)                                                          
         dg(i)=a*r(j)*gg(j)*c                                                  
 91      dp(i)=b*dg(i)                                                          
      i=imax-npi+1                                                              
      k=-1                                                                      
      go to 51                                                                  
      
c      write(7,*) 'end test1', test1
                                                                    
 999  return                                                                    
      end                                                                       

