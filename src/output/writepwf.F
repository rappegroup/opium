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
      subroutine writepwf(ifp, symbol)
      
c -- Perform KB transformation and output the info in the 
c -- PWF pseudopotential file format
      
      implicit double precision (a-h,o-z)
      
#include "fortdim.h"

      common /grid/ h,r1,z,r(npdm),np
      common /en/ ien(10),inv                                           
      common /dql/ dql
      common /maxx/ maxx(5)
      common /flqeta/ eta(5),flq(maxflq,5),nflq
      common /nlcore/ rvloc(npdm)
      common /pwf/ zeff, rpcc, spacing, first,
     &             nval2, igrid, inl, ist1, ncut, ilocal,ilocalind

      common /rscore/ rscore(npdm),rdd(npdm),rddd(npdm),rscoretot(npdm)
      common /totpot/ rvcore(npdm,n0),rvps(npdm,n0),rvcoul(npdm)
      common /aorb/ ncore,nval,norb,nlm(n0),
     $     no(n0),lo(n0),nmax(n0),maxim
      common /adat/ wnl(n0),en(n0),so(n0),xion
      common /aval/ rcall(n0),rvap(n0),rnorm(n0),ibd(n0),etot
      common /wfn/ rnl(npdm,n0)
      common /psdat/ nll

      dimension vl(npspt0),bb(npdm)
      dimension xr(npdm),bb3(npdm)
      dimension dflq(maxflq,15),ill(4)
      dimension psp(npspt0),dpsp(npspt0)
      dimension rrvloc(npdm)

      character*1 xrpcc
      character*5 xproj

      character*8 xform
      character*(*) symbol

      common /filenames/ file_log
      character*80 file_log
      character*80 namepwf

      data bohrad,ryd / 0.529177,13.6057 /
      
c     *************************************************************************
c     redirect stdout to the log file
c     *************************************************************************

      open(unit=7,file=file_log,form='formatted',access='append')

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c -- some hardcoded things      
      
      pi=acos(-1.0)                                                   
      pi4=4.0*pi                                                      

c -- Info originally out of *.crysft file
c     ist1: is the method (2=regular, 3=local pot only, 
c                         NOTHING ELSE supported)
c     inv: number of intervals for the piecewise intergration
c                              of the potential

c -- hardcoded settings gjt moved out of runNL to where they are used (here):

      inv  = 1
      mq = 4000
      dql = 0.020000
      dqnl = 0.075000
      nflq=501

c -- original writepwf comes here ---

c     EJW: ilocal is the 'l' value of the local pot
c          ilocalind is the valence index of the local pot
c          use ilocalind + 1 for fortran

      if (rpcc.gt.1e-12) then
         xrpcc = 'y'
      else
         xrpcc = 'n'
      endif

      xproj="     "

      if (nll.eq.1) then
         iproj = 0
         xproj = "    "
      elseif (nll.eq.2) then
         if (inl.eq.1) then
            iproj = 2
            xproj = "sp  "
         else
            iproj = 1
            if (ilocal.eq.0) then
               xproj="p   "
            else
               xproj="s   "
            endif
         endif
      elseif (nll.eq.3) then
         if (inl.eq.1) then
            iproj = 3
            xproj = "spd "
         else
            iproj = 2
            if (ilocal.eq.0) then
               xproj="pd  "
            elseif (ilocal.eq.1) then
               xproj="sd  "
            else
               xproj="sp  "
            endif
         endif
      elseif (nll.eq.4) then
         if (inl.eq.1) then
            iproj = 4
            xproj = "spdf"
         else
            iproj = 3
            if (ilocal.eq.0) then
               xproj="pdf "
            elseif (ilocal.eq.1) then
               xproj="sdf "
            elseif (ilocal.eq.2) then
               xproj="spf "
            else
               xproj="spd "
            endif
         endif
      endif
      
c -- number of |nlm> AO (including m-multiplicity!)
      numAOs = iproj**2
      
c -- hardcoded settings for now...
      aorad = 5.0

c -- now call the new C routine that writes the header info      
      call cwritepwf(ifp, symbol, zeff, NPSPT0, NFLQ, iproj, nval,
     $       xproj, xrpcc, np, spacing, first, numAOs, aorad)

c -- end of original writepwf comes here ---

c -- original pspfft stuff

c     Take care of long range tail of potential

      do i = maxim+1,np
         do j = 1,nll
            rvcore(i,j) = -zeff-zeff
         enddo
      enddo

c     Get non-local projectors 
      if(ist1.eq.2) call klby                                           

c     Write real space potentials to output file
      do k=1,np
         rrvloc(k) = rvloc(k)/r(k)
c     NOTE: redefinition of rvloc no factors of r  
         
         do i=1,nll 
            rvcore(k,i)=rvcore(k,i)/r(k)                                    
c     NOTE: redefinition of rvcore no factors of r 
         enddo

      enddo

      if (inl.eq.1) then
         write(7,185)
         ilocal=-67
c     this allows the local+box to be local!
      else
         if (ilocal.eq.0) write(7,190) 
         if (ilocal.eq.1) write(7,200)                                      
         if (ilocal.eq.2) write(7,210)                                      
      endif

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
         bb(k)=rvloc(k)*r(k)**2+zeff*2.0*r(k)                            
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
      
c     Now, change untis and write to .dta
      npspts=mq+2
      gcut = float(npspts-2)*dql
      psp(1)=sumg0*ryd*(bohrad**3)                                     
      do n=2,npspts-1
         psp(n)=vl(n-1)*ryd*(bohrad**3)                                     
      enddo
      do n=2,npspts-2                                                 
         dpsp(n)=(psp(n+1)-psp(n-1))*bohrad/(2.0*dql)                    
      enddo
c     linearly interpolate dpsp at endpoints  1,npspts-1 see note at top   
      dpsp(1)=2.0*dpsp(2)-dpsp(3)                                     
      dpsp(npspts-1)=2.0*dpsp(npspts-2)-dpsp(npspts-3)                
      psp(npspts)=gcut/bohrad                                           
      dpsp(npspts)=dql/bohrad                                           


c -- use the C I/O call to output the local potential coefficients
c
c      call nclear(ifp)
      do i=1,npspts
         call nwrite2(ifp, 4, "%18.8f", psp(i))
      enddo
      call nclear(ifp)
      do i=1,npspts
         call nwrite2(ifp, 4, "%18.8f", dpsp(i))
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

      do m=1,nll                                                      
         do n=2,nflq-1                                                  
            dflq(n,m)=(flq(n+1,m)-flq(n-1,m))*bohrad/(dql*2.0)                
         enddo
         if (nll.eq.1) then                                                  
            dflq(1,m)=2.0*dflq(2,m)-dflq(3,m)                          
         else                                                         
            dflq(1,m)=0.0                                              
         endif                                                        
         dflq(nflq,m)=2.0*dflq(nflq-1,m)-dflq(nflq-2,m)                  
      enddo

c -- New lines that omit the 0's but still produce the BH format

      do i=1,4
         ill(i)=0
      enddo

c      xform="%20.10f"
      call nclear(ifp)
      do k=1,nll
         if (k.ne.ilocal+1) then
            do i=1,nflq
               call nwrite2(ifp, 3, "%20.10f", flq(i,k))
            enddo
         endif
      enddo

      call nclear(ifp)
      do k=1,nll
         if (k.ne.ilocal+1) then
            do i=1,nflq
               call nwrite2(ifp,  3, "%20.10f", dflq(i,k))
            enddo
         endif
      enddo
      call nclear(ifp)
      
      do k=1,nll
         if (k.ne.ilocal+1) then
            call nwrite2(ifp, 3, "%20.10f", eta(k))
         endif
      enddo
      call nclear(ifp)
      
 100  format(//,' npoint=',i5,2x,'atomic number=',i3,2x,                
     $     'valence charge=',i3,/)                                           

 180  format(' sorry, partner, but ist1 must equal 2, or 3.')         

 185  format(1x,'local potential has been designed according to njr')
 190  format(1x,'local potential is s')                                 
 200  format(1x,'local potential is p')                                 
 210  format(1x,'local potential is d')                                 


      if (rpcc.gt.1e-12) then
         do i=1,np
          call nwrite2(ifp, 2, "%26.16e",r(i))
          call nwrite2(ifp, 2, "%26.16e",rscore(i))
        enddo
      endif

      call nclear(ifp)


      close(unit=7)

      end
      
