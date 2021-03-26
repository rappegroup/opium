c
c Copyright (c) 1998-2008 The OPIUM Group
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
      subroutine writerecpot(ifp, ifp_param, psmeth)
      
c -- Perform KB transformation and output the info in the 
c -- CASTEP .recpot pseudopotential file format
c -- KR Aug 2004.  Adapted from writepwf.f
      
      implicit double precision (a-h,o-z)
      
#include "fortdim.h"

c      common /atomic/ xion,rnl(npdm,n0),nlm(n0),wnl(n0),en(n0),norb
      common /grid/ h,r1,z,r(npdm),np
      common /en/ ien(10),inv                                           
      common /dql/ dql
      common /maxx/ maxx(5)
      common /flqeta/ eta(5),flq(maxflq,5),nflq
      common /nlcore/ rvloc(npdm)
      common /pwf/ zeff, rpcc, spacing, first,
     &             nval2, igrid, inl, ist1, ncut, ilocal,ilocalind
c      common /nmax/ nmax(n0),maxim

      common /rscore/ rscore(npdm),rdd(npdm),rddd(npdm),rscoretot(npdm)
      common /totpot/ rvcore(npdm,n0),rvps(npdm,n0),rvcoul(npdm)

      common /aorb/ ncore,nval,norb,nlm(n0),
     $     no(n0),lo(n0),nmax(n0),maxim
      common /adat/ wnl(n0),en(n0),so(n0),xion
      common /aval/ rcall(n0),rvap(n0),rnorm(n0),ibd(n0),etot
      common /wfn/ rnl(npdm,n0)
      common /psdat/ nll

c      common /np/ ncores,nvales
c      common /npp/ npots

      common /optparam/ qcl(10),nbl(10)

      dimension vl(npspt0),bb(npdm)
      dimension xr(npdm),bb3(npdm),qcore(npdm)
      dimension ill(4)
      dimension psp(npspt0)
      dimension rrvloc(npdm)

      character*1 psmeth

      common /filenames/ file_log
      character*80 file_log

      data bohrad,ryd / 0.529177249d0,13.6056981d0 /
      
c     *************************************************************************
c     redirect stdout to the log file
c     *************************************************************************

      open(unit=7,file=file_log,form='formatted',access='append')

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c -- some hardcoded things      
      
      pi=3.14159265358979323846d0
      pi4=4.d0*pi                                                      

c -- Info originally out of *.crysft file
c     ist1: is the method (2=regular, 3=local pot only, 
c                         NOTHING ELSE supported)
c     inv: number of intervals for the piecewise intergration
c                              of the potential
c     EJW: ilocal is the 'l' value of the local pot
c          ilocalind is the valence index of the local pot
c          use ilocalind + 1 for fortran

c -- hardcoded settings gjt moved out of runNL to where they are used (here):

      inv  = 1
      mq = 2000
      dql = 0.05d0*bohrad
      dqnl = dql
      nflq = mq+1

c     Get non-local projectors 
      if(ist1.eq.2) call klby                                           

c     Take care of long range tail of potential

      do i = maxim+1,np
         do j = 1,nval
            rvcore(i,j) = -zeff-zeff
         enddo
      enddo

c     Write real space potentials to output file
      do k=1,np
         rrvloc(k) = rvloc(k)/r(k)
c     NOTE: redefinition of rvloc no factors of r  
         
         do i=1,nval 
            rvcore(k,i)=rvcore(k,i)/r(k)                                    
c     NOTE: redefinition of rvcore no factors of r 
         enddo

      enddo

c     Decide what the local potential is
c      if (inl.ne.0) then
c         iloc = 0
c         goto 170
c      endif
c      if (ist1.eq.2) then                                                
c         iloc=ilocalind+1                                                    
c         goto 170                                                       
c      endif                                                          
c      if (ist1.eq.3) then                                               
c         if (ilocal.ne.0) then                                            
c            iloc = ilocalind+1                                               
c            goto 170                                                    
c         end if                                                      
c         iloc = 1                                                     
c         goto 211                                                       
c      end if                                                         
c     Wrong method (ist1 <> 1 or 2)
c      write(7,180)                                                      
c      stop                                                              
c  170 continue 

      if (inl.eq.1) then
         write(7,185)
         ilocal=-67
c     this allows the local+box to be local!
      else
         if (ilocal.eq.0) write(7,190) 
         if (ilocal.eq.1) write(7,200)                                      
         if (ilocal.eq.2) write(7,210)                                      
      endif

c 211  continue 

c Compute a likely cutoff
      if( psmeth .eq. 'o' ) then
        qcmax = 0.0
        do j=1,nval
          qcmax = max(qcl(j),qcmax)
        end do   
      else
c Can we make a better guess for Kerker? KR        
        qcmax = 7.07
      end if
      ecutev = qcmax**2*ryd
ccccc
c     Now put local potential into rvlr

      do k=1,np                                                
         rx=r(k)

         if (inl.ne.0) then
            rvloc(k) = rrvloc(k)
         else
            rvloc(k)=rvcore(k,ilocalind+1)  
         endif

      enddo

c     Calculate fourier components of vlocal(r)                         
c     Get q grid
      do k=1,np                                                
         xr(k)=r(k)                                                     
         bb(k)=rvloc(k)*r(k)**2+zeff*2.d0*r(k)                            
      enddo

      pow=1.0e0                                                         
      if (zeff.eq.0.0e0) pow = 2.0e0                                      
c     compute q=0 term
      call radina(xr,bb,0,np,h,pow,inv,ien) 

      sumg0 = pow * pi4                                                   
      
      do j=1,mq                                                    
         ql=dql*dble(j)                                             
         pow=1.0e0                                                      
         if (zeff.eq.0.0e0) pow=2.0e0                                     
         do k=1,np                                                
            dot=xr(k)*ql                                             
            bb3(k) = bb(k) * sin(dot)/dot                               
         enddo
         
         call radina(xr,bb3,0,np,h,pow,inv,ien) 
         vl(j)=pow * pi4                                                

c     NOTE: the units are Ryd
      enddo
      
c     Now, change untis and write to .recpot
      npspts=mq+1
      gcut = float(npspts-1)*dql


c      close(unit=7)

      if (rpcc.gt.1e-12) then
        call writeparam(ifp,ifp_param,ecutev,6)
      else
        call writeparam(ifp,ifp_param,ecutev,5)
      end if
c      open(unit=7,file=file_log,form='formatted',access='append')      

      call nwrite2(ifp, 1, "%18.8f", gcut/bohrad)
      psp(1)=sumg0*ryd*(bohrad**3)                                     
      do n=2,npspts
         psp(n)=(vl(n-1)-2.0d0*pi4*zeff/(dql*(n-1))**2)*ryd*(bohrad**3)
      enddo

c      call nclear(ifp)
c -- use the C I/O call to output the local potential coefficients
      do i=1,npspts
        call nwrite2(ifp, 3, "%20.10f", psp(i))
      enddo
      call nclear(ifp)
        

c     for kleinman-bylander non-local psp                               
c     change units
      do m=1,nll                                                      
         do n=1,nflq                                                    
            flq(n,m)=flq(n,m)*ryd*sqrt((bohrad)**3)                           
         enddo
         eta(m)=eta(m)*ryd                                                 
      enddo

c -- New lines that omit the 0's but still produce the BH format

      do i=1,4
         ill(i)=0
      enddo

      do k=1,nll
         if (k.ne.ilocal+1) then
           call iwrite2(ifp, 1, "%d",k-1)
           call nwrite2(ifp, 1, "%20.10f",eta(k))
           do i=1,nflq
               call nwrite2(ifp, 3, "%20.10f", flq(i,k))
           enddo
         endif
      enddo

      call iwrite2(ifp, 1, "%d", 1000)

c
c NLCC
c
      if (rpcc.gt.1e-12) then

c     Convert to recip space

        do j=1,mq                                                    
          ql=dql*dble(j-1)                                             
          pow=2.0d0                                                      
          if( j .le. 1 ) then
            do k=1,np                                                
              bb3(k) = rscore(k)
            enddo
          else
            do k=1,np                                                
              dot=xr(k)*ql                                             
              bb3(k) = rscore(k) * sin(dot)/dot                               
            enddo
          end if
          call radina(xr,bb3,0,np,h,pow,inv,ien) 
          qcore(j)=pow
          
        enddo
        do i=1,npspts
          call nwrite2(ifp, 3, "%20.10f", qcore(i))
        enddo
      end if
      call nclear(ifp)

      
 100  format(//,' npoint=',i5,2x,'atomic number=',i3,2x,                
     $     'valence charge=',i3,/)                                           

 180  format(' sorry, partner, but ist1 must equal 2, or 3.')         

 185  format(1x,'local potential has been designed according to njr')
 190  format(1x,'local potential is s')                                 
 200  format(1x,'local potential is p')                                 
 210  format(1x,'local potential is d')                                 

      close(unit=7)

      end
      
