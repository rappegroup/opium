      subroutine descreen(ixc)
      
      implicit double precision(a-h,o-z)
      
#include "PARAMOPT"

c -------------------------------------------------------------------------
c     External (shared between C and Fortran) common blocks
c -------------------------------------------------------------------------
      common /grid/ h,r1,z,r(npdm),np
      common /atomic/ xion,rnl(npdm,n0),nlm(n0),wnl(n0),en(n0),norb
      common /rpcc/ rpcc,rpccz
      common /totpot/ rvcore(npdm,n0),rvps(npdm,n0),rvcoul(npdm)
      common /np/ ncores,nvales
      common /nmax/ nmax(nvale0),maxim
      common /rscore/ rscore(npdm),rdd(npdm),rddd(npdm),rscoretot(npdm)
      common /wavrc/ wavrc, slope, curvrc, indrc
      common /consts/ etol,vtol,maxit,isoft
c -------------------------------------------------------------------------

c -------------------------------------------------------------------------
c     Internal (Fortran only) common blocks
c -------------------------------------------------------------------------
      common /atom4/ v(npdm),rv(npdm),p(npdm),rsatom(npdm),
     $     c1(npdm),rvh(npdm),rvxc(npdm),rexc(npdm)
      common /atom3/ rsvale(npdm)
      common /iterm/ iterm
c -------------------------------------------------------------------------

      dimension rvxcv(npdm),rexcv(npdm),rvxct(npdm)
      dimension rexct(npdm),rvhc(npdm),rvhv(npdm),g(npdm)

      character*1 xc(0:3)

      xc(0)='s'
      xc(1)='p'
      xc(2)='d'
      xc(3)='f'

      pi=acos(-1.0)                                                   
      pi4=4.0*pi                                                      

      do i = 1,np
        g(i) = 0
      enddo

      isoft=1
      zero = 0.0
      iatom = 7
      charge = dfloat(2 * 1)

      write(7,*) '------------------------------'
      write(7,*) 'Descreening potential'

      call radin(r,rsvale,0,maxim,h,charge)
      write (7,9010) charge

      rsc=0.0
         
      if (rpcc.gt.1e-12) then
         call radin(r,rscore,0,maxim,h,rsc)      
         write (7,9011) rsc
      else
         write (7,9011) rsc        
      endif 

 9010 format(1x,'valence charge          : ',f10.6)
 9011 format(1x,'core    charge          : ',f10.6)

c     *************************************************************************
c     section 1:  input/output from hfs files                             
c     *************************************************************************

      do i=1,maxim
        rsatom(i)=rsvale(i)+rscore(i)
      enddo
      
      call excorr (maxim,ixc,r,rsvale,rvxcv,rexcv)        
      call hrtree (maxim,h,r,rsvale,rvhv)
      call excorr (maxim,ixc,r,rsatom,rvxct,rexct)        
      call hrtree (maxim,h,r,rscore,rvhc)
      
c     compute ionic pseudopotential by subtracting off screening terms.     

      wv = 0.         
      do i=1,nvales
        wv = wv + wnl(i)
      enddo
      wc = z - xion - wv  
      zv = z - wc   
      zv2 = zv+zv     
      xn2 = zv2 - xion - xion   
      maxp1 = maxim + 1

      do i=1,maxim
        if (rpcc.gt.1e-12) then
          rvcoul(i) = rvhv(i)+rvxct(i)      
        else
          rvcoul(i) = rvhv(i)+rvxcv(i)      
        endif
                       
        do lp=1,nvales      
          rvcore(i,lp) = rvps(i,lp)-rvcoul(i)    
        enddo
c       amr         rvcore(i,4) = -zv2      
      enddo

c     *************************************************************************
c     treat tail region separately (where wavefunctions are zero).
c     *************************************************************************

      do i = maxp1,np       
        do lp = 1,nvales
          rvcore(i,lp) = -zv2   
          rvps(i,lp) = -xion-xion         
        enddo
        rvcoul(i)=rvps(i,1)-rvcore(i,1)
      enddo
 
c     *************************************************************************
c     solve schroedinger equation for each state in turn in following loop, 
c     *************************************************************************

      write(7,*) 
      write(7,*) '----Solving the Schodinger equation',
     $           ' for all states----'
      ig=0
      do m=1,nvales

        ee = en(m)         
        inl = 1  
        isoft = 1
        ebef = ee
        n = nlm(m)/100     
        l = nlm(m)/10 - 10*n         
        nn = l + 1
        im = maxim

        call schsl(m,nn,l,ee,im,rvps(1,m),p,ig)
        
        if (iterm.eq.0) then
          write(7,9111) n,xc(l),ee,en(m)
 9111     format(1x,'State: ',i1,a1,2x,'AE eigenvalue = ', f10.6,2x,
     $         'PS eigenvalue = ', f10.6)
        endif            
        if (im.gt.maxim) maxim = im
        if (abs(ebef-ee).gt.0.1.or.iterm.eq.1) then
          write (7,*) 'This wavefunction has a node.  Try a '
          write (7,*) 'different starting guess for coefficients'
          write (7,*) ' or a different set of constraints and '
          write (7,*) 'weights.'
          stop
        endif

c     NOTE: the units are Ryd

      enddo
 
      return
      end
