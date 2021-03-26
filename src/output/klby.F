c
c Copyright (c) 1998-2005 The OPIUM Group
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
      subroutine klby                                                   
      implicit double precision (a-h,o-z)
                                             
#include "fortdim.h"

c      common /atomic2/ rnl(npdm,5)                                      
c      common /pseudo2/ r(npdm),rvps(5,npdm),rvlr(npdm),
c     &                 npoint,iz,izv,idi,mcut
      common /atomic/ xion,rnl(npdm,n0),nlm(n0),wnl(n0),en(n0),norb
      common /grid/ h,r1,z,r(npdm),np
      common /dql/ dql
      common /maxx/ maxx(5)
      common /en/ ien(10),inv                                           
      common /flqeta/ eta(5),flq(maxflq,5),nflq
      common /np/ ncores,nvales      
      common /nlcore/ rvloc(npdm)
      common /pwf/ zeff, rpcc, spacing, first,
     &             nval, igrid, inl, ist1, ncut, ilocal,ilocalind
      common /totpot/ rvcore(npdm,n0),rvps(npdm,n0),rvcoul(npdm)
      common /nmax/ nmax(n0),maxim
      dimension delcor(npdm),f(npdm),f3(npdm),ill(4)                      

c -gjt-------------------------------------------------------------
c klby is now using besfn(), which is the sqherical Bessel function
c -----------------------------------------------------------------

c ------ input parameters ---------------------------------
c 
c rnl (i,j):  radial wave function j on grid i
c r(i):       real space grid
c rvcore(j,i):  descreened pseudopot j on grid i
c rvlr(i):    - not used here -
c npoint:     - not used here - 
c iz:         - not used here - 
c izv:        - not used here - 
c nvales:        number of angular momenta
c mcut:       - not used here - 
c h:          r(i+1)/r(i)
c dql:        ??? hardcoded in pspfft.f
c maxx(j):    number of grid points for projector j
c ien(i):     dummy array for radina()???
c inv:        number of intervals for radina()
c eta(j):     KB denominator for projector j
c flq(j,i):   KB factors for projector j, PW state i
c rvloc(i):   local potential in case of designed local part
c 
c zeff:     - not used here - effective Z
c rpcc:     - not used here - partical core radius
c spacing:  - not used here - spacing of the grid
c first:    - not used here - first grid point
c nval:     - not used here - number of valence orbitals
c igrid:    - not used here - number of grid points
c inl:      =0 no designed, =1 designed local part
c ist1:     - not used here - local pot flag
c ncut:     - not used here - grid point cutoff
c ilocal:   local potential angular momentum index
c 
c ---------------------------------------------------------

c -- print and check the number of projectors

      write(7,10) nval                                                   
      if (nval.le.1) then                                                 
         write(7,2223)                                                 
         stop                                                          
      endif                                                         

c -- set up some variables used later on:
      
      gcut=dql*float(nflq-1)

c -- main loop over projectors
      do i=1,4
         ill(i)=0
      enddo

      do k=1,nval
         do lp=1,nvales
            nnj=nlm(lp)/100
            ljp=(nlm(lp) - nnj * 100)/10 + 1
            if ((ill(ljp).eq.0).and.(ljp.eq.k)) then
               ill(ljp)=ill(ljp)+1
               
               im=nmax(lp)

c -- -- get the local potential into the right format
        
               do j=1,im
                  if (inl.ne.0) then
                     delcor(j) = rvcore(j,lp)-rvloc(j)
                  else
                     delcor(j)=rvcore(j,lp)-rvcore(j,ilocalind+1)
                  endif

                  f3(j)=delcor(j)*rnl(j,lp)**2                    
               enddo
               
c     -- -- calculate the KB denominator: eta
               
               fx = dble(k+k+1)
               call radina(r,f3,-1,im,h,fx,inv,ien)

               eta(k) = fx
               write(7,998) k-1,eta(k)
               
c     -- -- KB factor loop
               
               do iq=1,nflq
                  
                  q = dql*(iq-1)
                  
c     -- ---- use besfn() to do Fourier transform of the potential
                  
                  do i=1,im        
                                             
                     qr = r(i)*q
                     f(i) = delcor(i)*rnl(i,lp)*besfn(qr,k-1)
                  enddo
                  
                  fx = dble(k+k)
                  
                  call radina(r,f,0,im,h,fx,inv,ien)
                  
                  flq(iq,k) = fx
                  
               enddo
            
            endif   
         enddo
      enddo

 10   format(1x,'number of angular momentum states',i5)                 
 2223 format(1x,'not two or more angular momentum components')          
 998  format(1x,'for l=',i3,' eta equals',d16.8)
 999  format(1x,'as written,kleinb can only compute for nonlocal s,p,d')

      return
      end
