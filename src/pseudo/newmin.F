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
      subroutine newmin(numfn,rx,rkinfin)
      implicit double precision(a-h,o-z)

#include "fortdim.h"

      parameter(ncon0=3)
      parameter(lwork=numfn0*numfn0)
      common /angm/ ll
      common /rke/ rKmat(numfn0,numfn0),rKvec(numfn0)
      common /rn/ rNmat(numfn0,numfn0)
      common /rlm/ rlmat(3,numfn0),rlvec(numfn0)
      common /rconst/ rkcon,rncon,rlcon(3)
      common /re/ rkin,rnorm
      common /ncon/ ncon

      dimension rx(numfn0),ry(numfn0),rz(numfn0)
      dimension rltemp(ncon0,numfn0)

      dimension rfor(numfn0),grad(numfn0)
      dimension rwork(lwork)
      dimension VNV(numfn0,numfn0),
     $     rMi(numfn0,numfn0),rKMisN(numfn0)
      dimension VKV(numfn0,numfn0),rMisN(numfn0)
      dimension sN(numfn0),sK(numfn0),rJ(numfn0,numfn0),
     $     rJi(numfn0,numfn0),raV(numfn0)
      dimension Umat(ncon0,ncon0),Vmat(numfn0,numfn0),Sv(ncon0),
     $     Vb(numfn0),Vsum(numfn0),VSz(numfn0)
      dimension rKtil(numfn0,numfn0),rkvtil(numfn0)

      dimension ryt(numfn0),rforold(numfn0)

      ncon=2
      iprint=1

      write(7,*) 
      write(7,*) 'Starting KE minimization...'

      call lambda(rx,numfn)

c      rlvec(3)=0.0
c      do i=1,numfn
c         rlmat(3,i)=1.0
c         rlvec(3)=rlvec(3)-rx(i)
c      enddo
c      if (ll.eq.0.and.-rlvec(3).lt.0.1) ncon=3

c      call fnset2(iprint,rx,numfn)

      do i=1,ncon
         do j=1,numfn
            rltemp(i,j)=rlmat(i,j)
         enddo
      enddo

      call dgesvd('A','A',ncon,numfn,rltemp,
     $     ncon0,Sv,Umat,ncon0,Vmat,numfn0,
     $     rwork,lwork,info)
      if (info.ne.0) then
         write(7,*) 'newmin SVD error: INFO=',info
         stop
      endif

c     Vmat is the actually the transpose of V
      do i=1,numfn
         ry(i)=0.0
         do j=1,numfn
            ry(i)=ry(i)+Vmat(i,j)*rx(j)
         enddo
      enddo
      
      do i=1,numfn
         do j=1,numfn
            VNV(i,j)=0.0
            VKV(i,j)=0.0
            do n=1,numfn
               do m=1,numfn
                  VNV(i,j)=VNV(i,j)+Vmat(i,n)*Vmat(j,m)*rNmat(n,m)
                  VKV(i,j)=VKV(i,j)+Vmat(i,n)*Vmat(j,m)*rKmat(n,m)
               enddo
            enddo
         enddo
      enddo

      sNs=0.0
      sKs=0.0
      do i=1,ncon
         do j=1,ncon
            sNs=sNs+VNV(i,j)*ry(i)*ry(j)
            sKs=sKs+VKV(i,j)*ry(i)*ry(j)
         enddo
      enddo

      do i=ncon+1,numfn
         sN(i)=0.0
         sK(i)=0.0
         do j=1,ncon
            sN(i)=sN(i)+2.0*ry(j)*VNV(i,j)
            sK(i)=sK(i)+2.0*ry(j)*VKV(i,j)
         enddo
      enddo

      yNy=0.0
      RNy=0.0
      yKy=0.0
      RKy=0.0
      do i=ncon+1,numfn
         do j=ncon+1,numfn
            yNy=yNy+VNV(i,j)*ry(i)*ry(j)
            yKy=yKy+VKV(i,j)*ry(i)*ry(j)
         enddo
         RNy=RNy+sN(i)*ry(i)
         RKy=RKy+sK(i)*ry(i)
      enddo

      do i=1,numfn
         Vb(i)=0.0
         Vsum(i)=0.0
         do j=1,numfn
            Vb(i)=Vb(i)+Vmat(i,j)*rKvec(j)
            Vsum(i)=Vsum(i)+Vmat(i,j)
         enddo
      enddo

      Vbs=0.0
      VSs=0.0
      do i=1,ncon
         Vbs=Vbs+Vb(i)*ry(i)
         VSs=VSs+Vsum(i)*ry(i)
      enddo

      Vby=0.0
      Vsy=0.0
      do i=ncon+1,numfn
         Vby=Vby+Vb(i)*ry(i)
         Vsy=Vsy+Vsum(i)*ry(i)
      enddo

c      write(7,*) 'YSUM1', vsy+vss

c     ------------------------------------------------

      do i=1,numfn-ncon
         do j=1,numfn-ncon
            rJ(i,j)=VNV(i+ncon,j+ncon)
         enddo
      enddo

      call dpotrf('U',numfn-ncon,rJ,numfn0,info)

      if (info.ne.0) then
         write(7,*) 'newmin M Cholesky INFO=',INFO
         stop
      endif

      do i=1,numfn-ncon
         do j=i+1,numfn-ncon
            rJ(j,i)=0.0
         enddo
      enddo

      do i=1,numfn-ncon
         do j=1,numfn-ncon
            rMi(i,j)=rJ(i,j)
            rJi(i,j)=rJ(i,j)
         enddo
      enddo
      
      call dpotri('U',numfn-ncon,rMi,numfn0,info)
      if (info.ne.0) then
         write(7,*) 'newmin M inversion INFO=',INFO
         stop
      endif

      call dtrtri('U','N',numfn-ncon,rJi,numfn0,info)
      if (info.ne.0) then
         write(7,*) 'newmin J inversion INFO=',INFO
         stop
      endif

      do i=1,numfn-ncon
         do j=i+1,numfn-ncon
            rMi(j,i)=rMi(i,j)
         enddo
      enddo

      do i=ncon+1,numfn
         rMisN(i)=0.0
         do j=ncon+1,numfn
            rMisN(i)=rMisN(i)+rMi(i-ncon,j-ncon)*sN(j)
         enddo
         rMisN(i)=rMisN(i)/2.0
      enddo
      
c     y=y+M'b/2
      do i=ncon+1,numfn
         ry(i)=ry(i)+rMisN(i)
      enddo
      
      rNrem=0.0
      do i=ncon+1,numfn
         rNrem=rNrem+sN(i)*rMisN(i)
      enddo
      rNrem=rNrem/2.0

      yMy=0.0
      do i=ncon+1,numfn
         do j=ncon+1,numfn
            yMy=yMy+VNV(i,j)*ry(i)*ry(j)
         enddo
      enddo
      
      do i=ncon+1,numfn       
         rKMisN(i)=0.0
         do j=ncon+1,numfn       
            rKMisN(i)=rKMisN(i)+2.0*VKV(i,j)*rMisN(j)
         enddo
      enddo

      sVKV=0.0
      do i=ncon+1,numfn
         do j=ncon+1,numfn
            sVKV=sVKV-rMisN(i)*rMisN(j)*VKV(i,j)
         enddo
         sVKV=sVKV-(sK(i)+Vb(i)-rKMisN(i))*rMisN(i)
      enddo

      raVy=0.0
      vsy=0.0
      vsy2=0.0
      do i=ncon+1,numfn
         raV(i)=sK(i)+Vb(i)-rKMisN(i)
         raVy=raVy+raV(i)*ry(i)
         vsy=vsy+vsum(i)*(ry(i)-rMiSN(i))
         vsy2=vsy2-vsum(i)*rMiSN(i)
      enddo

      yKy=0.0
      do i=ncon+1,numfn
         do j=ncon+1,numfn
            yKy=yKy+ry(i)*ry(j)*VKV(i,j)
         enddo
      enddo

c      write(7,*) 'YSUM', vsy+vss
c      write(7,*) 'YKY+aY+c:',yky+raVy+sKs+Vbs+rkcon+sVKV
c      write(7,*) 'YNY+d   :',yMy-rNrem+sNs+rncon

c     -------------------------------------------------------
      rkconstant=sKs+Vbs+sVKV + rkcon
      rnormconstant=sNs-rNrem + rncon

      do i=ncon+1,numfn
         rz(i)=0.0
         do j=ncon+1,numfn
            rz(i)=rz(i)+ry(j)*rJ(i-ncon,j-ncon)
         enddo
      enddo

c     rz is the new vector of coeffs

c     Lets fix the norm

      rnn=0.0
      do i=ncon+1,numfn
         rnn=rnn+rz(i)**2
      enddo
      rmm=sqrt(abs(rnormconstant)/rnn)
      do i=ncon+1,numfn
         rz(i)=rz(i)*rmm
      enddo

c     this is how to get the norm
      rnn=0.0
      do i=ncon+1,numfn
         rnn=rnn+rz(i)**2
      enddo

c     transform Kmat -> Kmat_tilde
c     rz' * Kmat_tilde * rz = yKy
      rkk=0.0
      do i=ncon+1,numfn
         do j=ncon+1,numfn
            rKtil(i,j)=0.0
            do n=ncon+1,numfn
               do m=ncon+1,numfn
                  rKtil(i,j)=rKtil(i,j)+VKV(n,m)
     $                 *rJi(n-ncon,i-ncon)*rJi(m-ncon,j-ncon)
               enddo
            enddo
            rkk=rkk+rKtil(i,j)*rz(i)*rz(j)
         enddo
      enddo

c     transform Kvec -> kvec_tilde
c     kvec_tilde * rz = ay
      rkv=0.0
      do i=ncon+1,numfn
         rkvtil(i)=0.0
         VSz(i)=0.0
         do j=ncon+1,numfn
            rkvtil(i)=rkvtil(i)+raV(j)*rJi(j-ncon,i-ncon)
            VSz(i)=VSz(i)+Vsum(j)*rJi(j-ncon,i-ncon)
         enddo
         rkv=rkv+rkvtil(i)*rz(i)
      enddo

c      write(7,*) 'KE=',rkk+rkv+rkconstant
c      write(7,*) 'Norm=',rnn+rnormconstant

c     ----------------------------------------
c     Ok, now lets find the force

      write(7,9045)
      call flush(7)
      do ntim=1,1500

         do i=ncon+1,numfn
            rforold(i)=rfor(i)
c            rlam(i)=-(rkk+rkv)/rnormconstant
         enddo

         do i=ncon+1,numfn
            rfor(i)=0.0
            do j=ncon+1,numfn
               rfor(i)=rfor(i)+2.0*rKtil(i,j)*rz(j)
            enddo
            rfor(i)=rfor(i)+rkvtil(i)
         enddo
         
         modv=10
         if (ntim.eq.1.or.(ntim/modv)*modv.eq.ntim) then
            do i=ncon+1,numfn
               grad(i)=-rfor(i)
            enddo
         else
            rtop=0.0
            rbot=0.0
            do i=ncon+1,numfn
               rtop=rtop+rfor(i)**2
               rbot=rbot+rforold(i)**2
            enddo
            do i=ncon+1,numfn
               grad(i)=rfor(i)+(rtop/rbot)*grad(i)
            enddo
         endif
         rtop=0.0
         rbot=0.0
         do i=ncon+1,numfn
            rtop=rtop+grad(i)*rz(i)
            rbot=rbot+rz(i)*rz(i)
         enddo

         do i=ncon+1,numfn
            grad(i)=grad(i)-(rtop/rbot)*rz(i)
         enddo
                  
         rov=0.0
         do i=ncon+1,numfn
            rov=rov+grad(i)*rz(i)
         enddo
         
         zkz=0.0
         gkg=0.0
         ag=0.0
         az=0.0
         zkg=0.0
         do i=ncon+1,numfn
            do j=ncon+1,numfn
               zkz=zkz+rz(i)*rKtil(i,j)*rz(j)
               gkg=gkg+grad(i)*rKtil(i,j)*grad(j)
               zkg=zkg+rz(i)*rKtil(i,j)*grad(j)
            enddo
            az=az+rkvtil(i)*rz(i)
            ag=ag+rkvtil(i)*grad(i)
         enddo
         
         rnz=0.0
         rng=0.0
         do i=ncon+1,numfn
            rnz=rnz+rz(i)**2
            rng=rng+grad(i)**2
         enddo
         rcs=sqrt(rnz/rng)
         
         rtop=(ag+2.0*zkg)*rcs
         rbot=2*(gkg*rcs*rcs-zkz)-az
         rth=-rtop/rbot
 
         do i=ncon+1,numfn
            rz(i)=rz(i)*cos(rth)+grad(i)*sin(rth)*rcs
         enddo
         
         rnn=0.0
         do i=ncon+1,numfn
            rnn=rnn+rz(i)**2
         enddo
         
         rkk=0.0
         do i=ncon+1,numfn
            do j=ncon+1,numfn
               rkk=rkk+rKtil(i,j)*rz(i)*rz(j)
            enddo
         enddo
         
         rkv=0.0
         do i=ncon+1,numfn
            rkv=rkv+rkvtil(i)*rz(i)
         enddo

         rsum=0.0
         do i=ncon+1,numfn
            rsum=rsum+VSz(i)*rz(i)
         enddo
         rsum=rsum+Vss+vsy2

c      write(7,9030) rkk+rkv+rkconstant,rnn+rnormconstant
         if ((mod(ntim-1,10).eq.0).or.(ntim.lt.5)) then
            write(7,9040) ntim,rth,rtop,rbot/2,rkk+rkv+rkconstant,rsum
         endif
         if (abs(rtop).lt.1e-3.or.rkk+rkv+rkconstant.lt.1e-5.
     $        or.rsum.lt.0.01) goto 911
      enddo
      
      write(7,*) 'After 500 CG steps, the state is still not converged!'
      write(7,*) 'Use CONMAX'
      stop

 911  continue
      write(7,*)
      write(7,9025) ntim
      write(7,9020) rkk+rkv+rkconstant
      write(7,9030) rsum

 9025 format(1x,'# steps: ',i5)
 9020 format(1x,'Final KE error:',f16.10)
 9030 format(1x,'Sum of coeffs :',f16.10)
 9040 format(1x,i5,5f12.6)
 9045 format(1x,' step #    theta       slope       ',
     $     'curv       KEresid    coeffsum')
c     ----------------------------------------
c     Lets see if we know how to z -> y -> x

      do i=ncon+1,numfn
         ryt(i)=0.0
         do j=ncon+1,numfn
            k=i-ncon
            kk=j-ncon
            ryt(i)=ryt(i)+rz(j)*rJi(k,kk)
         enddo
      enddo

      do i=1,ncon
         ryt(i)=ry(i)
      enddo

      do i=1,numfn
         rx(i)=0.0
         do j=1,numfn
            rx(i)=rx(i)+Vmat(j,i)*(ryt(j)-rMisN(j))
         enddo
      enddo

      call fnset2(iprint,rx,numfn)
      
      return
      end

