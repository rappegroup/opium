
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

      subroutine soldir(iorb,nqn,lang,en,fl,agi,api,dv,
     $     ainf,test1,test2,kap,imax,ifail,method)  
      implicit double precision (a-h,o-z)
      
#include "fortdim.h"

c -------------------------------------------------------------------------
c     External (shared between C and Fortran) common blocks
c -------------------------------------------------------------------------
      common /consts/ etol,vtol,maxit,isoft
      common /grid/ h,r1,z,r(npdm),np

c -------------------------------------------------------------------------
c     Internal (Fortran only) common blocks                
c -------------------------------------------------------------------------
      common /iterm/ iterm
      common /elim/ emore,eless,elim
      common /ipos/ ipos(n0),itermcount
      common /itoo/ ietoolo,ietoohi
      common /iskip/ iskip
      common /vdir/ vdir(npdm),pold(npdm)
      common /rold/ rnlold(npdm,n0),rpold(npdm,n0)
      common /wfn/ rnl(npdm,n0)
      common /cl/ cl

      logical bound
c -------------------------------------------------------------------------

      common /hfrel/ eg(npdm),ep(npdm)
      common /hfrel3/ cg(npdm,n0),cp(npdm,n0)
      common /hfrel5/ ceg(10),cep(10),ag(10),ap(10),av(10)
      common /hfrel6/ gg(npdm),gp(npdm)

      dimension hg(npdm),hp(npdm)
      dimension dv(npdm),bg(npdm),bp(npdm)
      dimension agh(10),aph(10),bgh(10),bph(10)

      ndor=10
      enav=1.0d 00                                                              
      ainf= abs(ainf)                                                           
      ccl=cl+cl                                                                 
      iex=method                                                                
      nes=150
      ah=0.d0

c      write(7,*) 'entered soldir',iex,method,agi,api
      if (method.le.0) method=1                                                 
c notice that below iex=0,1 and method=1,2 only.                                
c this was used to simplify block structure of program. ala 11/22/94            
      fk=kap                                                                    
      if (av(1).lt.0.0d 00.and.kap.gt.0) api=-agi*(fk+fl)/av(1)                 
      if (av(1).lt.0.0d 00.and.kap.lt.0) api=-agi*av(1)/(fk-fl)                 
c      write(7,*) 'av',av(1),api,agi,fl,fk
      ell=fk*(fk+1.0d 00)/ccl                                                   
      node=nqn- abs(kap)                                                         
      if (kap.lt.0) node=node+1                                                 
      emin=0.0                                                                  
      do 91 i=1,imax                                                              
         a=(ell/(r(i)*r(i))+dv(i))*cl                                         
         if (a.lt.emin) emin=a                                                  
 91   continue                                                                  

c      write(7,*) 'hi',emin,a,imax
      if (emin .ge. 0.0) then                                                   
c         write(6,*) kap,iorb,nqn
c         stop 'emin ge 0'
         goto 899
      endif                                                                     
      if (en.lt.emin) en=emin*0.9d 00                                           
      edep=en                                                                   
                                                                                
 101  continue
c      write(7,*) 'at 101'
      numerr=0                                                                  
      test=test1                                                                
      if (method.gt.1) test=test2                                               
      einf=1.0d 00                                                              
      esup=emin                                                                 
      en=edep                                                                   
      ies=0                                                                     
      nd=0                                                                      
 105  continue
c      write(7,*) 'at 105'
      jes=0                                                                     
 106  continue
c      write(7,*) 'at 106'
      modmat=0                                                                  
      imm=0                                                                     
      if ( abs((enav-en)/en).lt.1.0d-01) imm=1                                  
      enav=en                                                                   

c      write(7,*) 'entered soldir2',iex,method
c     integration of the inhomogenious system                                   
 107  continue
c      write(7,*) 'at 107'
      do 111 i=1,imax                                                           
         gg(i)=eg(i)                                                            
 111     gp(i)=ep(i)                                                            
      do 115 i=2,ndor                                                           
         ag(i)=ceg(i-1)                                                         
 115     ap(i)=cep(i-1)                                                         
c         write(7,*) 'before first int test1=',en
      call intdir (gg,gp,ag,ap,ggmat,gpmat,en,dv,av,
     $        fl,agi,api,ainf,test1,test2,imax,mat,emin,esup,kap,imm)             
c      write(7,*) 'after first int iex, test1=',en

      if (iex.ne.0) go to 141                                                   
                                                                                
c     match large component for the homogenios system(method=0)                 
      a=ggmat/gg(mat)                                                           
      do 135 i=mat,imax                                                         
         gg(i)=a*gg(i)                                                          
 135     gp(i)=a*gp(i)                                                          
      j=mat                                                                     
      go to 215                                                                 
                                                                                
c     integration of the homogenios system                                      
 141  continue
c      write(7,*) 'at 141'
      do 151 i=1,imax                                                           
            hg(i)=0.0d 00                                                       
 151     hp(i)=0.0d 00                                                          
      do 155 i=1,ndor                                                           
         agh(i)=0.0d 00                                                         
 155     aph(i)=0.0d 00                                                         
c         write(7,*) 'imm,method',imm,method
      imm=1                                                                     
      if (method.eq.1) imm=-1                                                   
c         write(7,*) 'imm,method',imm,method
c         write(7,*) 'before second int',mat
      call intdir (hg,hp,agh,aph,hgmat,hpmat,en,dv,av,
     $        fl,agi,api,ainf,test1,test2,imax,mat,emin,esup,kap,imm)           
c       write(7,*) 'after second int iex',mat                           

c     match the large component for inhomogenious system(method=1)              
      a=gg(mat)-ggmat                                                           
      if (method.lt.2) then                                                     
         b=-a/hg(mat)                                                           
      else                                                                      
         b=gp(mat)-gpmat                                                        
         ah=hpmat*hg(mat)-hgmat*hp(mat)                                         
         if (abs(ah).lt.1d-8) go to 263                                           
         c=(b*hg(mat)-a*hp(mat))/ah                                             
         b=(b*hgmat-a*hpmat)/ah                                                 
         do 165 i=1,ndor                                                        
            ag(i)=ag(i)+c*agh(i)                                                
 165        ap(i)=ap(i)+c*aph(i)                                                
         j=mat-1                                                                
         do 168 i=1,j                                                           
            gg(i)=gg(i)+c*hg(i)                                                 
 168        gp(i)=gp(i)+c*hp(i)                                                 
      endif                                                                     
c      write(7,*) 'before dash', b,a,hg(mat),gg(mat),ggmat
      do 173 i=mat,imax                                                         
         gg(i)=gg(i)+b*hg(i)                                                    
 173     gp(i)=gp(i)+b*hp(i)                                                    

c     ------------------------------------------------------------
c         write(7,*) 'method',method
      if (method.ge.2) then                                                     
c        integration of the system derived from disagreement in energy          
         do 175 i=2,ndor                                                        
            bgh(i)=ag(i-1)/cl                                                   
 175        bph(i)=ap(i-1)/cl                                                   
         do 177 i=1,imax                                                        
            bg(i)=gg(i)*r(i)/cl                                                
 177        bp(i)=gp(i)*r(i)/cl  
c         write(7,*) 'before third int',iex
         call intdir (bg,bp,bgh,bph,bgmat,bpmat,en,dv,av,
     $        fl,agi,api,ainf,test1,test2,imax,mat,emin,esup,kap,imm)        
c         write(7,*) 'after third int iex',iex
                                                                                
c        match both components for inhomogenious system (method=2)              
         f=bg(mat)-bgmat                                                        
         g=bp(mat)-bpmat                                                        
         a=(g*hg(mat)-f*hp(mat))/ah                                             
         g=(g*hgmat-f*hpmat)/ah                                                 
         do 181 i=1,j                                                           
            bg(i)=bg(i)+a*hg(i)                                                 
 181        bp(i)=bp(i)+a*hp(i)                                                 
         do 182 i=1,ndor                                                        
            bgh(i)=bgh(i)+a*agh(i)                                              
 182        bph(i)=bph(i)+a*aph(i)                                              
         do 183 i=mat,imax                                                      
            bg(i)=bg(i)+g*hg(i)                                                 
 183        bp(i)=bp(i)+g*hp(i)                                                 
c        calculate the norm                                                     
         call norm(b,hp,r,gg,gp,ag,ap,method,hx,ndor,                          
     1     gpmat,fl,imax,mat)                                                   
         
c        correction to the energy (method=2)                                    
         do 186 i=1,imax                                                        
 186     hg(i)=(gg(i)*bg(i)+gp(i)*bp(i))*r(i)                                  
         ah=0.0d 00                                                             
         c=0.0d 00                                                              
         do 187 i=2,imax,2                                                      
 187     ah=ah+hg(i)+hg(i)+hg(i+1)                                              
         ah=hx*(ah+ah+hg(1)-hg(imax))/3.0d 00+hg(1)/(fl+fl+1.0d 00)             
         f=(1.0d 00-b)/(ah+ah)                                                  
         c=1.0d 00-b                                                            
         do 191 i=1,imax                                                        
            gg(i)=gg(i)+f*bg(i)                                                 
 191        gp(i)=gp(i)+f*bp(i)                                                 
         do 195 i=1,ndor                                                        
            ag(i)=ag(i)+f*bgh(i)                                                
 195        ap(i)=ap(i)+f*bph(i)                                                
      endif                                                                     
c     ------------------------------------------------------------
                                                                                
c     search for the maximum of the modulus of large component                  
      a=0.0d 00                                                                 
      bgh(1)=b                                                                  
      bph(1)=ah                                                                 
      do 211 i=1,imax                                                           
         g=gg(i)*gg(i)                                                          
         if (g.le.a) go to 211                                                  
         a=g                                                                    
         j=i                                                                    
 211  continue                                                                  
      if (j.gt.mat .and. modmat.eq.0) then                                      
         modmat=1                                                               
         mat=j                                                                  
         if (mod(mat,2).eq.0) mat=mat+1                                         
         imm=1                                                                  
         if (mat.lt.(imax-10)) go to 107                                        
                                                                                
         mat=imax-12                                                            
         j=mat                                                                  
         if (mod(mat,2).eq.0) mat=mat+1                                         
c         write(7,*) "warning mat,imax",mat,imax
      endif                                                                     
c this case can happen due to bad starting point in scf procedure.              
c ignore this warning unless you are getting it at final norb calls of          
c soldir.  redirected by ala 11/21/94.                                          
c     numerr=220021                                                             
c * impossible matching point                                                   
c     go to 899                                                                 
                                                                                
c compute number of nodes                                                       
 215  nd=1                                                                      
      j= max(j,mat)
      do 231 i=2,j                                                              
         if (gg(i-1).eq.0.0d 00) go to 231                                      
         if ((gg(i)/gg(i-1)).le.0.0d 00) nd=nd+1                                
 231  continue                                                                  

c      write(7,*) 'nodes',nd,node,en

      if (nd-node) 251,305,261                                                  
 251  esup=en                                                                   
      if (einf.lt.0.0d 00) go to 271                                            
      en=en*8.0d-01                                                             
      if ( abs(en).gt.test1) go to 285                                          
      numerr=238031                                                             
c    *zero energy                                                               
      go to 899                                                                 
                                                                                
 261  einf=en     
c      write(7,*) 'eee',einf,en,esup,emin

      if (esup.gt.emin) go to 271                                               
 263  en=en*1.2d 00                                                             

      if (en.gt.emin) go to 285                                                 
      numerr=245041                                                             
c    *energy is lower than the minimum of apparent potential                    
      go to 899                                                                 
                                                                                
 271  if ( abs(einf-esup).gt.test1) go to 281                                   
      numerr=249051                                                             
c    *the upper and lower limits of energy are identical                        
      go to 899                                                                 
                                                                                
 281  en=(einf+esup)/2.0d 00                                                    
                       
 285  jes=jes+1                                                                 


      if (jes.le.nes) go to 106                                                 
      
c *number of attempts to find good number of nodes is over the limit            
c this case can happen due to bad starting point in scf procedure.              
c ignore this warning unless you are getting it at final norb calls of          
c soldir	                                                                       
c      write(7,*) "warning jes>nes"
      ifail=1                                                                   
c    *redirected by ala 11/21/94.                                               
c     numerr=255061                                                             
c     go to 899                                                                 
                                                                                
c     calculation of the norm                                                   
 305  continue
c      write(7,*) 'before norm', hp(100),gp(100)
      call norm(b,hp,r,gg,gp,ag,ap,method,h,ndor,                             
     1     gpmat,fl,imax,mat)                                                   
c         write(7,*) 'NORM=',b,gpmat,gp(mat),gg(mat)

      if (method.eq.1) then                                                     
c        correction to the energy (method=1)                                    
         c=gpmat-gp(mat)                                                        
         f=gg(mat)*c*cl/b                                                       
         if (gpmat.ne.0.0d 00) c=c/gpmat                                        
      endif                                                                     
                                                                                
      en=en+f                                                                   
      g= abs(f/(en-f))  

c      write(7,*) 'adj1',en
c      write(7,*) 'adj2',f,g
c      write(7,*) 'adj3',b,c
c      write(7,*) 'adj4',mat
c      write(7,*) 'adj5',gpmat,gp(mat),gg(mat)
c      write(7,*) 'c test esup einf',c,test,esup,einf

 371  if ((en.ge.0 .or. g.gt.2.0d-01) .or.                                      
     1 (abs(c).gt.test .and. (en.lt.esup.or.en.gt.einf))) then                  
c        try smaller step in enrgy under above conditions                       

         f=f/2.0d 00                                                            
         g=g/2.0d 00                                                            
         en=en-f                                                                

c         write(7,*) 'adjusted f',f,g,en

         if (g.gt.test1) go to 371                                              
c         write(7,*) 'g > test1'
         go to 899                                                              
      endif                                                                     

c      write(7,*) 'after norm', b,en,ies

c      write(7,*) '2nd comp',abs(c),test
      if ( abs(c).gt.test)  then                                                
         ies=ies+1                                                              
c         write(7,*) ies,nes
         if (ies.le.nes) go to 105                                              
         ifail=1                                                                
c         write(7,*) "warning ies=nes"
c     everything is fine unless you are getting this message                    
c     on the latest stage selfconsistent process.                               
c     just stopped trying to match lower component                              
c     because number of trials exceeded limit.                                  
c     lines below were commented out.  ala 11/18/94                             
      endif                                                                     

c     divide by a square root of the norm, and test the sign of w.f.            
      b= sqrt(b)                                                                
      c=b                                                                       
      if ((ag(1)*agi).lt.0.0d 00.or.(ap(1)*api).lt.0.0d 00) c=-c                
      do 711 i=1,ndor                                                           
         ag(i)=ag(i)/c                                                          
 711     ap(i)=ap(i)/c                                                          
      if ((gg(1)*agi).lt.0.0d 00.or.(gp(1)*api).lt.0.0d 00) b=-b                
      do 721 i=1,imax                                                           
         gg(i)=gg(i)/b                                                          
 721     gp(i)=gp(i)/b                                                          
      if (imax.ge.np) return                                                    
      j=imax+1                                                                  
      do 741 i=j,imax                                                            
         gg(i)=0.0d 00                                                          
 741     gp(i)=0.0d 00                                                          
c     if everything o'k , exit is here.                                                  
c         write(7,*) 'returning',ifail
      return                                                                    

 899  if (iex.eq.0 .or. method.eq.2) go to 999                                  
      method=method+1                                                           
      go to 101                            

 999  continue
      return
      end                                                                       

