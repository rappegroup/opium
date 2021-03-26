
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
      double precision function fdrirk (i,j,l,m,k,kap)                              
c                       * calculate radial integrales rk *                      
c        rk = integral of f(r) * uk(r,s) * g(s)                                 
c uk(r,s) = rinf**k / rsup**(k+1)    rinf=min(r,s)   rsup=max(r,s)              
c        if nem=0  f(.)=cg(.,i)*cg(.,j)+cp(.,i)*cp(.,j)                         
c                  g(.)=cg(.,l)*cg(.,m)+cp(.,l)*cp(.,m)                         
c        if nem non zero f(.)=cg(.,i)*cp(.,j)                                   
c                        g(.)=cg(.,l)*cp(.,m)                                   
c                  cg (cp) large (small) componenents of the orbitales          
c moreover if nem > or =0 the integration is made from 0 to infinity,           
c and otherwise from 0 to r.                                                    
c        this programm uses yzkrdf and dsordf                                   
                                                                                
      implicit double precision (a-h,o-z)                                       
#include "fortdim.h"
      common /hfrel/ eg(npdm),ep(npdm)
      common /hfrel6/ gg(npdm),gp(npdm)
      common /hfrel3/ cg(npdm,n0),cp(npdm,n0)
      common /hfrel5/ ceg(10),cep(10),ag(10),ap(10),av(10)
      dimension hg(npdm),kap(n0)                                                         
      common/inelma/nem                                                         
      save                                                                      
                                                                                
      fdrirk=0.0d 00                                                            
      if (i.le.0.or.j.le.0) go to 201                                           
      call yzkrdf (i,j,k,cg,cp)                                                       
      nn= abs(kap(i))+ abs(kap(j))                                              
      nn=max(nn-k,1)                                                            
      a=k+1                                                                     
      do 21 n=1,10                                                            
 21   hg(n)=0.0d 00                                                             
      do 31 n=1,10                                                            
         if (nn.gt.10) go to 31                                               
         hg(nn)=-ag(n)                                                          
 31      nn=nn+1                                                                
      do 41 n=1,10                                                            
 41      ag(n)=hg(n)                                                            
      ag(1)=ag(1)+ap(1)                                                         
                                                                                
 201  if (l.le.0.or.m.le.0) return                                              
      n=-1                                                                      
      if (nem.ne.0) n=-2                                                        
      fdrirk=dsordf(l,m,-1,n,a)                                                 
      return                                                                    
      end                                                                       
