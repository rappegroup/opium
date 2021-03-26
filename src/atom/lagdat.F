
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
      subroutine lagdat (ia,iex,kap)                                                
c        * non diagonal lagrange parameteres *                                  
c     lagrange parameters involving orbital ia if ia is positive                
c     all lagrange parameters are calculated if ia is negative or zero          
c     contribution of the exchange terms is omitted if iex=0                    
c        this program uses akeato(bkeato) fdrirk multrk                         
                                                                                
      implicit double precision (a-h,o-z)                                       
#include "fortdim.h"
      common /hfrel/ eg(npdm),ep(npdm)
      common /hfrel2/ bg(npdm,n0),bp(npdm,n0)
      common /hfrel3/ cg(npdm,n0),cp(npdm,n0)
      common /hfrel4/ fl(n0)
      common /hfrel5/ ceg(10),cep(10),ag(10),ap(10),av(10)
      common /grid/ h,r1,z,r(npdm),np
      common /aorb/ ncore,nval,norb,nlm(n0),
     $     no(n0),lo(n0),nmax(n0),maxim
      common /adat/ wnl(n0),en(n0),so(n0),xion
      common /aval/ rcall(n0),rvap(n0),rnorm(n0),ibd(n0),etot
      common /hfrel6/ gg(npdm),gp(npdm)

      common/scrhf1/eps(2*npdm),nre(N0),ipl                                        
      dimension kap(n0)
      dimension ni(2),nj(2)                                                     

      norbsc=norb
      i1= max(ia,1)                                                             
      idep=1                                                                    
      if (ia.gt.0) go to 15                                                     
 11   idep=i1+1                                                                 
 15   ni(1)=i1                                                                  
      nj(2)=i1                                                                  
      ji1=2* abs(kap(i1))-1                                                     
      do 201 i2=idep,norbsc                                                     
         if (i2.eq.i1.or.kap(i2).ne.kap(i1)) go to 201                          
         if (nre(i1).lt.0.and.nre(i2).lt.0) go to 201                           
c    the following line was included to handle the case of single               
c    electron in 2 s-shells                                                     
c    probably need to use schmidt orthogonalization in this case                
         if (wnl(i1).eq.wnl(i2)) go to 201                                    
         ni(2)=i2                                                               
         nj(1)=i2                                                               
         d=0.0d 00                                                              
         do 101 l=1,norbsc                                                      
            k=0                                                                 
            jjl=2* abs(kap(l))-1                                                
            kma= min(ji1,jjl)                                                   
 41         a=akeato(l,i1,k)/wnl(i1)                                           
            b=a-akeato(l,i2,k)/wnl(i2)                                         
            c=b                                                                 
            if (a.ne.0.0d 00) c=c/a                                             
            if ( abs(c).lt.1.0d-07) go to 51                                    
            d=d+b*fdrirk(l,l,i1,i2,k,kap)                                           
 51         k=k+2                                                               
            if (k.le.kma) go to 41                                              
            if (iex.eq.0) go to 101                                             
            kma=(ji1+jjl)/2                                                     
            k= abs(jjl-kma)                                                     
            if ((kap(i1)*kap(l)).lt.0) k=k+1                                    
 61         a=bkeato(l,i2,k)/wnl(i2)                                           
            b=a-bkeato(l,i1,k)/wnl(i1)                                         
            c=b                                                                 
            if (a.ne.0.0d 00) c=c/a                                             
            if ( abs(c).lt.1.0d-07) go to 71                                    
            d=d+b*fdrirk(i1,l,i2,l,k,kap)                                           
 71         k=k+2                                                               
            if (k.le.kma) go to 61                                              
 101     continue                                                               
         i= min(i1,i2)                                                          
         j= max(i1,i2)                                                          
         eps(i+((j-1)*(j-2))/2)=d/(wnl(i2)-wnl(i1))                           
 201  continue                                                                  
      if (ia.gt.0) go to 999                                                    
      i1=i1+1                                                                   
      if (i1.lt.norbsc) go to 11                                                
 999  return                                                                    
      end                                                                       
