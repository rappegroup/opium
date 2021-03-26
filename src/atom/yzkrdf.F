
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
      subroutine yzkrdf (i,j,k,cg,cp)                                                 
c       * calculate  function yk *                                              
c yk = r * integral of f(s)*uk(r,s)                                             
c uk(r,s) = rinf**k/rsup**(k+1)   rinf=min(r,s)   rsup=max(r,s)                 
c f(s)=cg(s,i)*cg(s,j)+cp(s,i)*cp(s,j)      if nem=0                            
c f(s)=cg(s,i)*cp(s,j)                      if nem is non zero                  
c f(s) is constructed by the calling programm  if i < or =0                     
c in the last case a function f (lies in the block dg) is supposedly            
c tabulated untill point dr(j), and its' devlopment coefficients                
c at the origin are in ag and the power in r of the first term is k+2           
                                                                                
c the output functions yk and zk are in the blocks dp and dg.                   
c at the origin  yk = cte * r**(k+1) - developement limit,                      
c cte lies in ap(1) and development coefficients in ag.                         
c        this programm uses aprdev and yzkteg                                   
                                                                                
      implicit double precision (a-h,o-z)                                       
#include "fortdim.h"
c      common cg(NP,N0),cp(NP,N0),bg(10,N0),bp(10,N0),fl(N0),ibgp              
c      common/comdir/cl,dz,dg(NP),ag(10),dp(NP),ap(10),bidcom(783)             
c      common/ratom1/xnel(N0),en(N0),scc(N0),scw(N0),sce(N0),                    
c     $     nq(N0),kap(N0),nmax(N0)                                                   
c      common/tabtes/hx,dr(NP),test1,test2,ndor,npp,nes,method,idim              
      common /grid/ h,r1,z,r(npdm),np
      common /aorb/ ncore,nval,norb,nlm(n0),
     $     no(n0),lo(n0),nmax(n0),maxim
      common /adat/ wnl(n0),en(n0),so(n0),xion
      common /aval/ rcall(n0),rvap(n0),rnorm(n0),ibd(n0),etot

      common /hfrel/ eg(npdm),ep(npdm)
      common /hfrel2/ bg(npdm,n0),bp(npdm,n0)
c      common /hfrel3/ cg(npdm,n0),cp(npdm,n0)
      common /hfrel4/ fl(n0)
      common /hfrel5/ ceg(10),cep(10),ag(10),ap(10),av(10)
      common /hfrel6/ gg(npdm),gp(npdm)

      common/inelma/nem                                                         

      dimension chg(10)                                                         
      dimension bgi(10),bgj(10),bpi(10),bpj(10)                                 
      dimension cg(npdm,n0),cp(npdm,n0)

      ndor=10
      if (i.le.0) go to 51                                                      
c     construction of the function f                                            
      do  5 l= 1,10                                                           
        bgi(l) = bg(l,i)                                                        
        bgj(l) = bg(l,j)                                                        
        bpi(l) = bp(l,i)                                                        
  5     bpj(l) = bp(l,j)                                                        
      id= min(nmax(i),nmax(j))                                                  
      ap(1)=fl(i)+fl(j)                                                         
      if (nem.ne.0) go to 31                                                    
      do 11 l=1,id                                                              
 11   gg(l)=cg(l,i)*cg(l,j)+cp(l,i)*cp(l,j)                                     
      do 21 l=1,ndor                                                            
 21   ag(l)=aprdev(bgi,bgj,l)+aprdev(bpi,bpj,l)                                 
      go to 55                                                                  
                                                                                
 31   do 35 l=1,id                                                              
 35   gg(l)=cg(l,i)*cp(l,j)                                                     
      do 41 l=1,ndor                                                            
 41   ag(l)=aprdev(bgi,bpj,l)                                                   
      go to 55                                                                  
                                                                                
 51   ap(1)=k+2                                                                 
      id=j                                                                      
 55   continue

      call yzkteg (gg,ag,gp,chg,dr,ap(1),k,10,id,np)                     

      return                                                                    
      end                                                                       

