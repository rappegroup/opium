
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
c
      subroutine dfsolve(zeff,ixc,exccut,ipsp,ifc,iexit,irel,iprint)
     
c     *************************************************************************
c     determine the self-consistent solution
c     *************************************************************************
     
      implicit double precision (a-h,o-z)

#include "fortdim.h"
 
c     -------------------------------------------------------------------------
c     External (shared between C and Fortran) common blocks
c     -------------------------------------------------------------------------

      common /wfnrel/ rnla(npdm,n0),rnlb(npdm,n0)
      
      common /rscore/ rscore(npdm),rdd(npdm),rddd(npdm),rscoretot(npdm),
     $     rsval(npdm)
      common /totpot/ rvcore(npdm,n0),rvps(npdm,n0),rvcoul(npdm)

      common /consts/ etol,vtol,maxit,isoft
      common /grid/ h,r1,z,r(npdm),np

      common /filenames/ file_log

      common /aorb/ ncore,nval,norb,nlm(n0),
     $     no(n0),lo(n0),nmax(n0),maxim
      common /adat/ wnl(n0),en(n0),so(n0),xion
      common /aval/ rcall(n0),rvap(n0),rnorm(n0),ibd(n0),etot
      common /wfn/ rnl(npdm,n0)
c     -------------------------------------------------------------------------
c     Internal (Fortran only) common blocks                
c     -------------------------------------------------------------------------
c     relhf
      common /hfrel/ eg(npdm),ep(npdm)
      common /hfrel2/ bg(npdm,n0),bp(npdm,n0)
      common /hfrel3/ cg(npdm,n0),cp(npdm,n0)
      common /hfrel4/ fl(n0)
      common /hfrel5/ ceg(10),cep(10),ag(10),ap(10),av(10)
      common /hfrel6/ gg(npdm),gp(npdm)
      common/inelma/nem                                                         
      common/scrhf1/eps(2*npdm),nre(N0),ipl    
      common /cl/ cl

c     *************************************************************************
c     local variables
c     *************************************************************************
      dimension sce(n0),scw(n0)
      dimension wavea(npdm),icarray(30)
      dimension epss(npdm),ka(n0),scc(n0)
      dimension rhart(npdm),xg(npdm)
c     *************************************************************************
c     *************************************************************************
      character*80 file_log

      cl=1.370373d+02

c     this is important! -EJW
      nem=0

      iexit=0
      itermcount=0
      isoft=ipsp
      nig=0
      ipratt=6
      npratt = ipratt
      rvapmax=0.0
      ib=1
      ndor=10
      jfail=0
      numerr=0

      open(unit=7,file=file_log,form='formatted',access='append')

      do j=1,norb
c         if (wnl(j).lt.1e-8) then
c            write(7,*) 'Occupation of all states must be > 1e-8'
c            write(7,*) ' for HF exchange -- STOP'
c            write(6,*) 'Occupation of all states must be > 1e-8'
c            write(6,*) ' for HF exchange -- STOP'
c            stop
c         endif
         do i = 1,np
            rnla(i,j) = 0
            rnlb(i,j) = 0
            rnl(i,j) = 0
         enddo
      enddo

      call setuphf(epss,ka,scc,ipsp,zeff)
      
      niter=500
      teste=etol
      testy=vtol
c      teste=5.d-6
c      testy=1.d-5
      j=1
      ind=1
      nter=0
      do i=1,norb
         scw(i)=0
      enddo
      test1=teste*100.d0
      test2=teste*10.d0
      netir=abs(niter)*norb

      if (iprint.gt.0) write(7,*) 'starting DF iterations'
      if (iprint.gt.0) write(7,191) 
 101  iort=0
      nter=nter+1

      if (niter.ge.0) goto 105
 104  call ortdat(j,ka)
 105  method=1

      if ((nre(j).gt.0).and.ipl.ne.0) call lagdat(j,1,ka)
      call potrdf(rvcore(1,j),rvps(1,j),j,ka,method)

      e=en(j)
      npp=nmax(j)
      idim=npp
      
      ifail=0

      ainf=cg(nmax(j),j)
      
      call soldir(j,no(j),lo(j),en(j),
     $     fl(j),bg(1,j),bp(1,j),rvps(1,j),
     $     ainf,test1,test2,
     $     ka(j),nmax(j),ifail,method)

      if ((ifail.ne.0).and.(jfail.eq.0)) jfail=j
      if ((jfail.eq.j).and.(ifail.eq.0)) jfail=0

      if (numerr.eq.0) goto 111
      iort=1
      goto 104

 111  sce(j)=abs((e-en(j))/en(j))

      pr=0.d0
      do i=1,nmax(j)
         w=cg(i,j)-gg(i)
         if (abs(w).lt.abs(pr)) goto 115
         pr=w
         a=cg(i,j)
         b=gg(i)
 115     w=cp(i,j)-gp(i)
         if (abs(w).lt.abs(pr)) goto 121
         pr=w
         a=cp(i,j)
         b=gp(i)
 121     continue
      enddo

      if (iprint.gt.0) write(7,491) 
     $     nter,j,en(j)*2.d0,sce(j)*2.d0,pr*2.d0

 491  format(1x,i6,2x,i6,f18.8,14x,e10.3,8x,e10.3)
 191  format(' iter',7x,'orb',13x,'eig',18x,
     $     'de_max',12x,'dv_max')

c     accel conv
      b=scc(j)
      call cofcon(a,b,pr,scw(j))
      scc(j)=b
      
      do i=1,nmax(j)
         gg(i)=b*gg(i)+a*cg(i,j)
         gp(i)=b*gp(i)+a*cp(i,j)
      enddo
      do i=1,ndor
         ag(i)=b*ag(i)+a*bg(i,j)
         ap(i)=b*ap(i)+a*bp(i,j)
      enddo
      mm=0

      a=dsordf(j,nmax(j),mm,4,fl(j))
      a=sqrt(a)
      
      do i=1,npp
         cp(i,j)=gp(i)/a
         cg(i,j)=gg(i)/a
      enddo
      do i=1,ndor
         bg(i,j)=ag(i)/a
         bp(i,j)=ap(i)/a
      enddo

      if ((nter.lt.norb).or.(ind.lt.0.and.j.lt.norb)) then
         j=j+1
         goto 451
      endif
      j=j+1
      pr=0.d0
      do i=1,norb
         w=abs(scw(i))
         if (w.gt.pr) then
            pr=w
            j=i
         endif
      enddo
      if (j.gt.norb) j=1
      if (pr.gt.testy) goto 421
      pr=0.d0
      do i=1,norb
         w=abs(sce(i))
         if (w.gt.pr) then
            pr=w
            j=i
         endif
      enddo
      if (pr.gt.teste) goto 421
      if (ind.lt.0) goto 999
      ind=-1
      j=1
      goto 451
 421  ind=1
 451  continue

      if (nter.le.netir) goto 101
      stop "out of its"
      
 999  continue
      if (numerr.eq.0) then         
         if (jfail.ne.0) then
            stop 'lower component not matched'
         endif
         
      endif

      call etotal(ka,ebs,eh,ex)

      
      write(7,*) "Overlap Integrals: {nlk|nlk}"
      zero=0.d0
      do 351 i=1,norb-1                                                         
         do 331 j=i+1,norb                                                      
            if (ka(j).ne.ka(i)) go to 331                                     
            at=dsordf(i,j,0,1, zero)                                         
            if (ka(j).gt.0) then
               write(7,2091)  no(i),lo(i),iabs(ka(i)),
     $              no(j),lo(j),iabs(ka(j)),at                 
            else
               write(7,2092)  no(i),lo(i),iabs(ka(i)),
     $              no(j),lo(j),iabs(ka(j)),at                 
            endif
 331     continue                                                               
 351  continue                                                                  
 2091 format(6x,'{',2i1,'+',i1,'|',2i1,'+',i1,'}=',1pd8.1)
 2092 format(6x,'{',2i1,'-',i1,'|',2i1,'-',i1,'}=',1pd8.1)


      do i=1,norb
         en(i)=en(i)*2.d0
         do k=1,np
            rnla(k,i)=cg(k,i)
            rnlb(k,i)=cp(k,i)
         enddo
      enddo
      if (iprint.gt.0) write(7,*) 'Dirac-Fock converged!'

      if (ipsp.eq.0) then
         method=1     
         maxim=0
         do i=1,norb
            maxim=max(maxim,nmax(i))
            call potrdf2(rvcore(1,i),rvps(1,i),i,ka,method,
     $           rnla,rnlb,rhart,xg,ipsp)
            do j=1,np
               rnl(j,i)=sqrt(rnla(j,i)**2+rnlb(j,i)**2)
               rvps(j,i)=rvcore(j,i)+2.d0*rhart(j)+xg(j)
            enddo
         enddo
      else
         do i=1,norb
            do j=1,np
               rnl(j,i)=rnla(j,i)
            enddo
         enddo
      endif
c     ***********************************************************************
c     ***********************************************************************

c     Flip wfn s.t. it is >0 as r->oo

      do i = ncore+1,norb
         do j=np,1,-1
            if (abs(rnlb(j,i)).gt.1e-10) then
               nmax(i)=j
               goto 939
            endif
         enddo
 939     continue

         if (rnlb(nmax(i)-1,i).lt.0.0) then
            do j = 1,np
               rnlb(j,i) = -rnlb(j,i)
            enddo
         endif

         do j=np,1,-1
            if (abs(rnla(j,i)).gt.1e-10) then
               nmax(i)=j
               goto 938
            endif
         enddo
 938     continue
         if (rnla(nmax(i)-1,i).lt.0.0) then
            do j = 1,np
               rnla(j,i) = -rnla(j,i)
            enddo
         endif
      enddo
      
c     compute valence density 
      do j = 1,np
         rsval(j) = 0
      enddo

      do i = ncore+1,norb
         if (irel.eq.0) then
            do j = 1,nmax(i)
               rsval(j) = rsval(j) + wnl(i) * rnl(j,i)**2
            enddo
         else
            do j = 1,nmax(i)
               rsval(j) = rsval(j) + wnl(i) * (rnla(j,i)**2
     $              +rnlb(j,i)**2)
               
            enddo
         endif
      enddo

c     compute core density

      if (ipsp.eq.0) then
         rsc=0.0d0
         do j = 1,np
            rscore(j) = 0
         enddo
         
         do i = 1,ncore
            if (irel.eq.0) then
               do j = 1,nmax(i)
                  rscore(j) = rscore(j) + wnl(i) * rnl(j,i)**2
               enddo
            else
               do j = 1,nmax(i)
                  rscore(j) = rscore(j) + wnl(i) *(rnla(j,i)**2
     $                 +rnlb(j,i)**2)
               enddo
            endif
         enddo
         call radin(r,rscore,0,np,h,rsc)
      endif

c     loop over valence
      do i=ncore+1,norb
         ic=i-ncore
         icarray(ic) = (log(rcall(ic)/r(1))/
     $        log(exp(1.0)))/h + 1
         do k = 1,np
            wavea(k) = rnla(k,i) * rnla(k,i)
     $           +rnlb(k,i) * rnlb(k,i)
         enddo
         pow = 2 * lo(i) + 2
         call radin(r,wavea,0,icarray(ic),h,pow)
         rnorm(i)=1.d0-pow
      enddo

c     radius of outermost peak
      do i=ncore+1,norb
         rvapmax=0.0
         do j=maxim-1,1,-1
            if (rnla(j,i).lt.rnla(j+1,i)) then
               rvap(i)=r(j+1)
               if (rvapmax.lt.rvap(i).
     $              and.(abs(rnla(j,i)).gt.1e-6)) then
                  rvapmax=rvap(i)
                  goto 945
               endif
            endif
         enddo
 945     continue
      enddo

c     ***********************************************************************
c     ***********************************************************************
      
      ebs=ebs*2.d0
      eh=-eh*2.d0
      ex=ex*2.d0
      
      etot=ebs+eh-ex
      if (iprint.gt.0) write(7,*) 

      if (niter.eq.1) then
         if (iprint.gt.0) write(7,700)
      else
         if (iprint.gt.0) write(7,500) nter
      endif

      if (iprint.gt.0) write(7,701) etot,ebs
      if (iprint.gt.0) write(7,702) -eh,ex

      if (iprint.gt.0) write (7,*)
      if (iprint.gt.0) write (7,9222)

      do i = 1,ncore
         if (iprint.gt.0) then 
            if (lo(i).eq.0) then
               write (7,9523) nlm(i),wnl(i),en(i)
            else
               if (so(i).lt.0) then
                  write (7,9323) nlm(i),wnl(i),en(i)
               else
                  write (7,9423) nlm(i),wnl(i),en(i)
               endif
            endif
         endif
      enddo
      
      do i = ncore+1,norb
         if (iprint.gt.0) then
            if (lo(i).eq.0) then
               write (7,9523) nlm(i),wnl(i),en(i),rnorm(i),rvap(i)
            else
               if (so(i).lt.0) then
                  write (7,9323) nlm(i),wnl(i),en(i),rnorm(i),rvap(i)
               else
                  write (7,9423) nlm(i),wnl(i),en(i),rnorm(i),rvap(i)
               endif
            endif
         endif
      enddo
      
      do i=1,norb
         do j=1,np
            cg(j,i)=rnla(j,i)
         enddo
      enddo

 500  format(1x,i4,1x,'iterations: (units are Ry and bohr)')
 700  format(1x,'Converged in 1 iteration (probably reference state)')
 701  format(1x,'Energy: ',f16.8,4x,'Ebs: ',f16.8)
 702  format(1x,'Eh:     ',f16.8,4x,'Ex : ',f16.8)

 9222 format(3x,'Orbital',4x,'Filling',7x,'Eigenvalues',4x,
     $     'Norm(rc->oo)',6x,'Peak')

 9223 format(4x,'|',i3,'>',5x,f6.3,5x,3f14.6)
 9224 format(4x,'|',i3,'>',5x,f6.3,5x,f14.6,"*",5x,
     $     '---------',5x,'--------')
 9323 format(4x,'|',i3,'->',5x,f6.3,5x,3f14.6)
 9324 format(4x,'|',i3,'->',5x,f6.3,5x,f14.6,6x,
     $     '---------',5x,'--------')
 9423 format(4x,'|',i3,'+>',5x,f6.3,5x,3f14.6)
 9424 format(4x,'|',i3,'+>',5x,f6.3,5x,f14.6,6x,
     $     '---------',5x,'--------')
 9523 format(4x,'|',i3,' >',5x,f6.3,5x,3f14.6)
 9524 format(4x,'|',i3,' >',5x,f6.3,5x,f14.6,6x,
     $     '---------',5x,'--------')

 9009 format(1x,'partial core radius : ', f10.4)
 9010 format(1x,'core charge         : ', f20.10)
 9011 format(1x,'partial core charge : ', f20.10)

      call flush(7)
      close(unit=7)

      return
      
 911  continue
c     911 = terminal error - exit scpot
      write(7,*) " !!ERROR!! Could not find self-consistent solution!"
      iexit=1
      
      call flush(7)
      close(unit=7)


 370  format(1x,i2,2i5,2f6.1,f10.4)
 390  format(//,' radial grid parameters',//,
     +     ' r(1) = .0 , r(2) =',e8.2,' , ... , r(',i5,') =',f6.2,
     +     /,' aa =',f5.2,'  bb =',f6.2,/)
      
 360  format(' input data for orbitals',//,
     +     '  i    n    l    s     j     occ',/)
      
 340  format(' nuclear charge             =',f10.6,/,
     +     ' number of core orbitals    =',i3,/,
     +     ' number of valence orbitals =',i3,/,
     +     ' electronic charge          =',f10.6,/,
     +     ' ionic charge               =',f10.6,//)
      
 332  format(' XC functional is GGA (Perdew-Burke-Ernzerhof)')
 334  format(' XC functional is LDA (Perdew-Zunger)')
 335  format(' XC functional is LDA (Perdew-Wang)')

      return
      
      end



