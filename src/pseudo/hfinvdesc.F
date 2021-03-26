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
      subroutine hfinvdesc(ixc,irel,zeff)
      
      implicit double precision(a-h,o-z)
      
#include "fortdim.h"

      common /grid/ h,r1,z,r(npdm),np      
      common /consts/ etol,vtol,maxit,isoft
      common /lall/ lall(n0)
      common /angm/ ll
      common /cuts/ qc,rc
      common /nnn/ nnn
      common /numfn/ numfn

      common /aorb/ ncore,nval,norb,nlm(n0),
     $     no(n0),lo(n0),nmax(n0),maxim
      common /adat/ wnl(n0),en(n0),so(n0),xion
      common /aval/ rcall(n0),rvap(n0),rnorm(n0),ibd(n0),etot
      common /wfn/ rnl(npdm,n0)
      common /totpot/ rvcore(npdm,n0),rvps(npdm,n0),rvcoul(npdm)
      common /partpot/ rvh(npdm),rvxc(npdm),rexc(npdm),
     $     rvx(npdm),rex(npdm),rvc(npdm),rec(npdm),rvdft(npdm)
      common /nlall/ nlall(n0)
      common /filenames/ file_log      
      common /inelma/ nem                                                         
      common /psout/ wsumt(n0),sumc(n0),npsghost(n0)
      common /hybrid/ exxw, ecw

      character*80 file_log
      
      dimension xg(npdm)
      dimension yr(npdm),x(npdm),ka(n0),rjunk(npdm),rhart(npdm)
      dimension rnljunk(npdm,n0),esave(n0)
      dimension rvnew(npdm,n0),rvnew2(npdm,n0),rnnew(npdm,n0)
      dimension ennew(n0),wnlnew(n0)
      dimension rsatom(npdm)

      character*1 xc(0:3)

      xc(0)='s'
      xc(1)='p'
      xc(2)='d'
      xc(3)='f'


      rmmax=15.d0
      ipsp=1
      norb=nval
      iprint=0



c     EJW
      nem=0
c     EJW this is for dirac-fock potrdf2.F

      do i=1,np
         rjunk(i)=0.d0
         do nnn=1,nval
            rnljunk(i,nnn)=0.d0
         enddo
      enddo

      if (irel.eq.0) then
         call hfcoef(iflag,ipsp)
         do nnn=1,nval
            do i=1,np
               rnl(i,nnn)=rnl(i,nnn)/sqrt(r(i))
            enddo
         enddo              
      else
         do nnn=1,nval
            ka(nnn) = -(lo(nnn)+1)
            if (so(nnn) .lt. 0.1D0 
     $           .and. lo(nnn) .ne. 0) ka(nnn)=lo(nnn)
         enddo
         call hfcoef_rel(ka)
      endif

      do nnn=1,nval
         
         do j=1,np
            if ((abs(rnl(j,nnn)*sqrt(r(j))).lt.1e-8).
     $           and.(r(j).gt.5.0)) then
               rmmax=r(j-5)
               nmv=j-10
               goto 788
            endif
         enddo
 788     continue
         if (iprint.gt. -1) write(7,*) 'Max radius used for Veff:',
     $        j,rmmax

         if (irel.eq.0) then

            !get charge density
             do m=1,norb
               do i = 1,np
                  rsatom(i) = 0
               enddo
            enddo

            maxim = 0 
            do m=1,norb
               do i = 1,np
                  rsatom(i) = rsatom(i) + wnl(m) * rnl(i,m)**2 * r(i)
               enddo
               maxim = max(maxim,nmax(m))
            enddo

            !set up dft potentials

            do i=1,np
               rvxc(i) = 0
               rvx(i) = 0
               rvc(i) = 0
               rvdft(i) = 0
            enddo

            call hrtree(np,h,r,rsatom,rvh)

            call excorr(np,ixc,exccut,rsatom,
     $          rvxc,rexc,rvx,rex,rvc,rec,irel)


            do i=1,np
              rvdft(i) = (1.0-exxw)*(rvh(i)+rvx(i)) + rvc(i)*ecw
            enddo

            !done setting up dft potentials

            call potl(yr,nnn)
            call xch(x,nnn,1)

            !debug: print out potentials
c            write(*,*) "in psp, maxim=", maxim," nval=",nval," nnn=",nnn
c            open(unit=91,file="rvps_psp.dat")
c            open(unit=92,file="yr_psp.dat")
c            open(unit=93,file="xrnl_psp.dat")
c            open(unit=94,file="rvdft_psp.dat")
c            open(unit=95,file="rvh_psp.dat")
c            open(unit=96,file="rvx_psp.dat")
c            open(unit=97,file="rvc_psp.dat")
c            do j=1,np
c               write(91,*) rvps(j,nnn)
c               write(92,*) 2.d0*yr(j)
c               write(93,*) x(j)/rnl(j,nnn)
c               write(94,*) rvdft(j)
c               write(95,*) rvh(j)
c               write(96,*) rvx(j)
c               write(97,*) rvc(j)
c            enddo
c
c            close(unit=91)
c            close(unit=92)
c            close(unit=93)
c            close(unit=94)
c            close(unit=95)
c            close(unit=96)
c            close(unit=97)

            
            do j=1,np
               rps=rnl(j,nnn)*sqrt(r(j))
               if (r(j).lt.rmmax) then
                  xxx=x(j)/sqrt(r(j))
                  rji=2.d0*yr(j)/r(j)
                  rvcore(j,nnn)=(rvps(j,nnn)-r(j)*(rji+xxx/rps))
               else
                  rvcore(j,nnn)=-2.d0*zeff
                  rvps(j,nnn)=-xion-xion-2
               endif
            enddo
         else
            m=1
            do j=1,np
               rjunk(j)=0.d0
               do mm=1,nval
                  rnljunk(j,mm)=0.d0
               enddo
            enddo
            
            call potrdf2(rjunk,rjunk,nnn,ka,m,
     $           rnl,rnljunk,rhart,xg,ipsp)

            do j=1,np
               if (r(j).lt.rmmax) then
                  rvcore(j,nnn)=(rvps(j,nnn)-2.d0*rhart(j)-xg(j))
               else
                  rvcore(j,nnn)=-2.d0*zeff
                  rvps(j,nnn)=xion-xion-2
               endif
            enddo
         endif
      enddo
 
c      if (irel.eq.1) then
c         iexit=0
c         ifc=0
c         exccut=0.d0
c         do i=1,norb
c            no(i)=lo(i)+1
c            esave(i)=en(i)
c         enddo
c         close(7)
c         call dfsolve(zeff,ixc,exccut,ipsp,ifc,iexit,irel,iprint)
c         
c         open(unit=7,file=file_log,form='formatted',access='append')
c         write(7,*) "Results of DF psp solve..."
c         do k=1,nval
c            xst=" "
c            if (ka(k).gt.0) xst="*"
c         
c            write(7,9112) no(k),xc(lo(k)),xst,esave(k),en(k)
c         enddo
c      endif


      if (irel.eq.1) then
         nc=0
         lolast=-444
         do nnn=1,nval
            write(7,*) "LO",nnn,no(nnn),lo(nnn)
            if (lo(nnn).eq.0) then
               nc=nc+1
               do j=1,np
                  rvnew(j,nc)=rvcore(j,nnn)
                  rvnew2(j,nc)=rvps(j,nnn)
                  rnnew(j,nc)=rnl(j,nnn)
               enddo
               ennew(nc)=en(nnn)
               wnlnew(nc)=wnl(nnn)
               lo(nc)=lo(nnn)
               no(nc)=lo(nc)+1
               wsumt(nc)=wsumt(nnn)
               sumc(nc)=sumc(nnn)
            endif

            xt=2.d0*lo(nnn)+1.d0
            xd=float(lo(nnn))/xt
            xu=float(lo(nnn)+1)/xt
            
            if ((lo(nnn).ne.0).and.(lo(nnn).ne.lolast)) then
               nc=nc+1
               do j=1,np
                  rvnew(j,nc)=rvcore(j,nnn)*xd
                  rvnew2(j,nc)=rvps(j,nnn)*xd
                  rnnew(j,nc)=rnl(j,nnn)*xd
               enddo
               ennew(nc)=en(nnn)*xd
               wnlnew(nc)=wnl(nnn)
               wsumt(nc)=wsumt(nnn)*xd
               sumc(nc)=sumc(nnn)*xd

            endif
            
            if ((lo(nnn).ne.0).and.(lo(nnn).eq.lolast)) then
               do j=1,np
                  rvnew(j,nc)=rvnew(j,nc)+rvcore(j,nnn)*xu
                  rvnew2(j,nc)=rvnew2(j,nc)+rvps(j,nnn)*xu
                  rnnew(j,nc)=rnnew(j,nc)+rnl(j,nnn)*xu
               enddo
               ennew(nc)=ennew(nc)+en(nnn)*xu
               wnlnew(nc)=wnlnew(nc)+wnl(nnn)
               lo(nc)=lo(nnn)
               no(nc)=lo(nc)+1
               wsumt(nc)=wsumt(nc)+wsumt(nnn)*xu
               sumc(nc)=sumc(nc)+sumc(nnn)*xu
            endif
            lolast=lo(nnn)
         enddo
         nval=nc
         norb=nc
         do k=1,nval
            en(k)=ennew(k)
            esave(k)=ennew(k)
            wnl(k)=wnlnew(k)
            nlm(k)=no(k)*100+lo(k)*10
            do j=1,np
               rvcore(j,k)=rvnew(j,k)
               rvps(j,k)=rvnew2(j,k)
               rnl(j,k)=rnnew(j,k)
            enddo
         enddo
      else
         do k=1,nval
            esave(k)=en(k)
            no(k)=lo(k)+1
            nlm(k)=no(k)*100+lo(k)*10
         enddo
      endif

      if (irel.eq.0) then
         do nnn=1,nval
            do j=1,np
               rnl(j,nnn)=rnl(j,nnn)*sqrt(r(j))
            enddo
         enddo
      endif

      write(7,*) 
      write(7,*) "Running HF solver for raw psp..."
      irel=0
      iexit=0
      ifc=0
      exccut=0.d0
      iprint=1

      close(7)

      call hfsolve(zeff,ixc,exccut,ipsp,ifc,iexit,irel,iprint)

      open(unit=7,file=file_log,form='formatted',access='append')
      write(7,*) "Results of HF psp solve..."
      do k=1,nval
         write(7,9111) no(k),xc(lo(k)),esave(k),en(k)
      enddo


 9111 format(1x,'State: ',i1,a1,2x,'AE eigenvalue = ', f10.6,2x,
     $     'PS eigenvalue = ', f10.6)
      
 9112 format(1x,'State: ',i1,a1,a1,2x,'AE eigenvalue = ', f10.6,2x,
     $     'PS eigenvalue = ', f10.6)

 9393 format(i5,2f20.10)
 6060 format(4f25.15)
 9000 format(1x,'--------------------Pseudizing state: ',
     $       '|',i3,'>',3x,'--------------------')
 9001 format(1x,'# basis functions        : ',i4)
 9002 format(1x,'rc                       : ',f8.4)
 9003 format(1x,'qc                       : ',f8.4)

 9011 format(1x,'Total Convergence Error: ',f10.6,1x,
     $       'mRy',3x,f10.6,1x,'meV') 
      
      return
      end
