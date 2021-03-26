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
      subroutine hfsmooth(iqpopt,rlocalr,ixc,tole)
      
c     *************************************************************************
c     smooth the hartree-fock potential
c     *************************************************************************

      implicit double precision(a-h,o-z)
      
#include "fortdim.h"
      common /consts/ etol,vtol,maxit,isoft
      common /grid/ h,r1,z,r(npdm),np      

      common /filenames/ file_log      
      common /cuts/ qc,rc
      common /nnn/ nnn
      common /angm/ ll
      common /aorb/ ncore,nval,norb,nlm(n0),
     $     no(n0),lo(n0),nmax(n0),maxim
      common /adat/ wnl(n0),en(n0),so(n0),xion
      common /aval/ rcall(n0),rvap(n0),rnorm(n0),ibd(n0),etot
      common /wfn/ rnl(npdm,n0)
      common /totpot/ rvcore(npdm,n0),rvps(npdm,n0),rvcoul(npdm)

      character*80 file_log

      dimension rgam(npdm),rnlo(npdm,n0),rvcoresave(npdm,n0)
      dimension fgam(npdm),wavea(npdm),waveb(npdm)
      dimension qar(n0),par(n0),rvh(npdm)
      dimension ilocalr(n0),zion(npdm),esave(n0)
      dimension rhead(n0),rhead0(n0),rnlsave(npdm,n0)
      dimension rlocalr(n0),qq(n0),pp(n0)

      open(unit=7,file=file_log,form='formatted',access='append')

      iprint=0

      tolpp=1e-3
      tolqq=1e-3

      ivcore=0
      open(unit=33,file="vcore",form="formatted")
      do i=1,np
         read(33,*,end=944) rt,rvh(i)
      enddo
      close(33)
      ivcore=1
 944  continue

      do i=1,np
         wavea(i)=0.0
         waveb(i)=0.0
         rgam(i)=0.0
      enddo
      
      ipsp=1      
      zeff=xion
      do nnn=1,nval
         zeff=zeff+wnl(nnn)
      enddo
      
      if (ivcore.eq.1) then
         do j=1,np
            zion(j)=-rvh(j)/2.d0+z
         enddo
      else
         do j=1,np
            zion(j)=zeff
         enddo
      endif

      do korb=1,nval
         if (rlocalr(korb).lt.1e-6) rlocalr(korb)=rcall(korb)
         do i=1,np
            rnlo(i,korb)=rnl(i,korb)
            rvcoresave(i,korb)=rvcore(i,korb)
         enddo
         esave(korb)=en(korb)
      enddo

      write (7,*) 
      write (7,*) 'Entering HF smoother'
c      write (7,*) 'HF options:  iqpopt,rlocalR',
c     $     iqpopt,(rlocalr(j),j=1,nval)

c     ------------------------------------------------------------------------
c     Determine localization radius
c
c     uncomment the next line to hardwire that rloc=rc
c      rlocalr=0

c     first find indrc
      rmaxrc=0.0
      do korb=1,nval
         indrc = 1                                                         
         do i = 1,np                                                    
            if (r(i).gt.rcall(korb)) goto 919
         enddo
 919     continue
         indrc=i

c     if rlocalr>0 then we have some other rloc
         if (rlocalr(korb).gt.0) then
            ilocalr(korb) = 1 
            do i = 1,np                                                    
               if (r(i).gt.rlocalr(korb)) goto 918
            enddo
 918        continue
            ilocalr(korb)=i
         else
            ilocalr(korb)=indrc
         endif
         write(7,*) 'localization radius: ',korb,ilocalr(korb),
     $        r(ilocalr(korb))

         rmaxrc=max(rmaxrc,r(ilocalr(korb)))
         alp=16.0**2/rmaxrc**2

      enddo
c     ------------------------------------------------------------------------
c     compute head norms for the original vion

      call headnorm(ilocalr,rnl,rhead0)
      if (iqpopt.gt.0) then
c     ------------------------------------------------------------
c     apply the p=0 and q=0 correction

         do korb=1,nval
            qq(korb)=0.0
            pp(korb)=0.0
            
            call setgam(iqpopt,ilocalr(korb),fgam)
            call applygamma(korb,ilocalr(korb),qq(korb),pp(korb),alp,
     $           zion,fgam,rvcoresave,rvcore)            
            
         enddo

         write(7,*) "QP results using q=0,p=0 adjustment"

         close(7)
         iprint=0
         call hfsolve(zeff,ixc,exccut,ipsp,ifc,iexit,irel,iprint)
         open(unit=7,file=file_log,form='formatted',access='append')
         
         do korb=1,nval
            write(7,9393) nlm(korb),en(korb),en(korb)-esave(korb)
         enddo

c         stop
c     ------------------------------------------------------------

         write(7,*)
         write(7,9292)

         mloop=0
         qq(korb)=1e-8
c         qq=1e-6
 990     continue
         tolpp=max(tolpp*0.9,1e-7)
         tolqq=max(tolqq*0.9,1e-7)
         
c         write(7,*) 'QP : new tol',tolpp,tolqq
         do korb=1,nval
            call setgam(iqpopt,ilocalr(korb),fgam)            

            pp(korb)=0.0
            del=0.0
            loop_qq=0
c     =q=q=q=q=q=q=q=q=q=q=q=q=q=q=q=q=q=q=q=q=q=q=q=q=q=q=q=q=q=q=q=q=q=q=q=q
c     ------------------------------------------------------------------------
 976        continue
            loop_qq=loop_qq+1
            if (mod(loop_qq,50).eq.0) then
               if (loop_qq.eq.500) stop 'QQ'
c               tolqq=tolpp/2.0
            endif

            loop_pp=0
c     =p=p=p=p=p=p=p=p=p=p=p=p=p=p=p=p=p=p=p=p=p=p=p=p=p=p=p=p=p=p=p=p=p=p=p=p
c     ------------------------------------------------------------------------
 977        continue
            loop_pp=loop_pp+1
            if (loop_pp.gt.500) stop 'PP'

            call applygamma(korb,ilocalr(korb),qq(korb),pp(korb),alp,
     $           zion,fgam,rvcoresave,rvcore)
            
            call flush(7)
            close(7)
            call hfsolve(zeff,ixc,exccut,ipsp,ifc,iexit,irel,iprint)
            open(unit=7,file=file_log,form='formatted',access='append')

            if (iexit.ne.0) then
               stop
            endif

            ppold=pp(korb)
            call computep(korb,ilocalr,qq(korb),fgam,rnlo,rnl,pp(korb))

c            write(7,*) 'QP  pp',korb,pp,loop_pp,loop_qq

            if (abs(ppold-pp(korb)).lt.tolpp) goto 958
            goto 977
c     ------------------------------------------------------------------------
c     =p=p=p=p=p=p=p=p=p=p=p=p=p=p=p=p=p=p=p=p=p=p=p=p=p=p=p=p=p=p=p=p=p=p=p=p
 958        continue
c     p is found, now update q
            
            write(7,9335) nlm(korb),mloop,loop_qq,loop_pp,
     $           pp(korb),qq(korb),(((en(k)-esave(k))),k=1,nval)

            del=abs(en(korb)-esave(korb))
            
            if (loop_qq.eq.1) then
               qqdel=1e-8
c               qqdel=1e-1
            else
               qqdel=-del*qqdel/(del-delold)
            endif
            
            delold=del
c            write(7,*) "QP  qq ", qq,qqdel

            qq(korb)=qq(korb)+qqdel

            if (del.lt.tolqq) goto 959
            goto 976
c     ------------------------------------------------------------------------
c     =q=q=q=q=q=q=q=q=q=q=q=q=q=q=q=q=q=q=q=q=q=q=q=q=q=q=q=q=q=q=q=q=q=q=q=q
 959        continue
c     q is found, save q&p and goto next orb
            qar(korb)=qq(korb)
            par(korb)=pp(korb)
         enddo   
c     korb loop
c     all orbs have been individually satisfied, check for global convergence
         
         dd=0.0
         do i=1,nval
            dd=dd+abs(en(i)-esave(i))            
         enddo
c         write(7,*) 'QP at end of loop:', dd,tole
         if ((dd.gt.tole)) then 
            mloop=mloop+1
            if (mloop.gt.500) stop 'EE'
c     loop over orbs again
            goto 990
         endif
c     Done!
c     ------------------------------------------------------------------------
         write(7,*) "QP results using optimal q,error (mRy), conv q & p"
         call flush(7)
         close(7)
         call hfsolve(zeff,ixc,exccut,ipsp,ifc,iexit,irel,iprint)
         open(unit=7,file=file_log,form='formatted',access='append')

         call headnorm(ilocalr,rnl,rhead)
                 
         etsum=0.0
         do korb=1,nval
            write(7,9394) nlm(korb),rhead(korb)-rhead0(korb)
            etsum=etsum+(rhead(korb)-rhead0(korb))*wnl(korb)
         enddo
         write(7,9395) etsum

         do korb=1,nval
            write(7,9393) nlm(korb),en(korb),
     $           (en(korb)-esave(korb))*1000,
     $           par(korb),qar(korb)
         enddo
      endif

c     ------------------------------------------------------------------------
c     ------------------------------------------------------------------------
c     write out some stuff
c     ------------------------------------------------------------------------
c      do nnn=1,nval
c         do i=1,np
c            rnlsave(i,nnn)=rnl(i,nnn)
c         enddo

c         ic1=nlm(nnn)/100
c         ll=(nlm(nnn)-ic1*100)/10
c         rc=rcall(nnn)

c         call fitwv      

c         do i=1,np
c            rnl(i,nnn)=rnlsave(i,nnn)
c         enddo

c      enddo

c      do korb=1,nval
c         do i=1,np
c            if (r(i).gt.15.0) goto 911
c            write(700+korb,*) r(i),rvcore(i,korb)
c            write(800+korb,*) r(i),rnl(i,korb),rnlo(i,korb)
c            write(810+korb,*) r(i),rnl(i,korb)-rnlo(i,korb)
c            write(820+korb,*) r(i),rnl(i,korb)**2-rnlo(i,korb)**2
c            write(750+korb,*) r(i),rvcore(i,korb),
c     $           rvcoresave(i,korb)
c            write(780+korb,*) r(i),(rvcore(i,korb)-
c     $           rvcoresave(i,korb))/r(i)
c            rdelt(i)=rdelt(i)+rnl(i,korb)**2-rnlo(i,korb)**2
c            rdelt2(i)=rdelt(i)*wnl(korb)
c         enddo
c 911     continue
c      enddo
c      do i=1,np
c         if (r(i).gt.15.0) goto 916
c         write(860,*) r(i),rdelt(i)
c 916     continue
c      enddo

      close(7)
      return

 9292 format(5x,'nlm',4x,'#it',2x
     $     '#q',4x,'#p',9x,'pp',13x,'qq',20x,'delta e(mRy)')
 9393 format(1x,'QP',i7,2f16.9,4x,2f16.8)
 9394 format(1x,'QP change in norm',i7,f16.9)
 9335 format(1x,'QP',4i5,2f14.8,4x,3d14.4)
 9395 format(1x,'QP Total #e transferred: ', f16.9)

 9000 format(1x,'--------------------Pseudizing state: ',
     $       '|',i3,'>',3x,'--------------------')
 9001 format(1x,'# basis functions        : ',i4)
 9002 format(1x,'rc                       : ',f8.4)
 9003 format(1x,'qc                       : ',f8.4)

 9011 format(1x,'Total Convergence Error: ',f10.6,1x,
     $       'mRy',3x,f10.6,1x,'meV') 
            
      close(unit=7)
      
      return      
      end


      subroutine setgam(iqpopt,ind,fgam)
c     *************************************************************************
c     smooth the hartree-fock potential
c     *************************************************************************

      implicit double precision(a-h,o-z)
      
#include "fortdim.h"
      common /grid/ h,r1,z,r(npdm),np      
      dimension fgam(npdm)

      if (iqpopt.eq.1) then
         do j=1,ind
            fgam(j)=r(j)**4 * (1-2.0/(3.0*r(ind)**2)*r(j)**2)
         enddo
         do j=ind+1,np
            fgam(j)=(r(ind)**4)/3.0
         enddo
      endif
      
      if (iqpopt.eq.2) then
         do j=1,ind
            fgam(j)=r(j)*(1.0-r(j)/(2.0*r(ind)))
         enddo
         do j=ind+1,np
            fgam(j)=r(ind)/2.0
         enddo
      endif
      
      if (iqpopt.eq.3) then
         do j=1,ind
            fgam(j)=r(j)
         enddo
         do j=ind+1,np
            fgam(j)=r(ind)
         enddo
      endif
      
      if (iqpopt.eq.4) then
         do j=1,ind
            fgam(j)=(1.0/r(j) + r(j)/r(ind)**2)
         enddo
         do j=ind+1,np
            fgam(j)=(2.0/r(ind))
         enddo
      endif

      return
      end

      subroutine applygamma(korb,ind,qq,pp,alp,
     $     zion,fgam,rvcoresave,rvcore)
c     *************************************************************************
c     smooth the hartree-fock potential
c     *************************************************************************

      implicit double precision(a-h,o-z)
      
#include "fortdim.h"
      common /grid/ h,r1,z,r(npdm),np      
      dimension fgam(npdm),zion(npdm),rgam(npdm)
      dimension rvcoresave(npdm,n0),rvcore(npdm,n0)

      do j=1,ind
         rgam(j)=qq+pp*fgam(j)
         rvcore(j,korb)=rvcoresave(j,korb)+rgam(j)*r(j)
c         rvcore(j,korb)=-(rvcoresave(j,korb)+rgam(j)*r(j))/2.d0
      enddo
            
      do j=ind+1,np-100
         rgam(j)=qq+pp*fgam(j)
         ggg=exp(-alp*(r(j)-r(ind))**2)
         rvv=ggg*(rgam(j)+rvcoresave(j,korb)/r(j)
     $        +zion(j)/r(j)*2.0)
     $        -zion(j)/r(j)*2.0
         rvcore(j,korb)=rvv*r(j)
c         rvcore(j,korb)=-rvv*r(j)/2.d0
      enddo
            
      do j=np-100+1,np
         rvcore(j,korb)=-zion(j)*2.0
c         rvcore(j,korb)=zion(j)
      enddo
c      write(7,*) 
      return
      end
            
      subroutine computep(korb,ilocalr,qq,fgam,rnlo,rnl,pp)
c     *************************************************************************
c     smooth the hartree-fock potential
c     *************************************************************************

      implicit double precision(a-h,o-z)
      
#include "fortdim.h"
      common /grid/ h,r1,z,r(npdm),np      
      dimension rnl(npdm,n0),rnlo(npdm,n0),ilocalr(n0)
      dimension wavea(npdm),waveb(npdm),fgam(npdm)

      do k=1,np
         wavea(k)=0.0
         waveb(k)=0.0
      enddo

      do k=1,ilocalr(korb)
         wavea(k)=rnlo(k,korb)*rnl(k,korb)*qq
      enddo
      pow1 = 2.0
      call radin(r,wavea,0,ilocalr(korb),h,pow1)
      
      do k=1,ilocalr(korb)
         waveb(k)=rnlo(k,korb)*rnl(k,korb)*fgam(k)
      enddo
      pow2 = 2.0
      call radin(r,waveb,0,ilocalr(korb),h,pow2)

      pp=-pow1/pow2
c      if (pp.gt.1) pp=1
c      if (pp.lt.-1) pp=-1
      return
      end

      subroutine headnorm(ilocalr,rnl,rhead)
c     *************************************************************************
c     smooth the hartree-fock potential
c     *************************************************************************

      implicit double precision(a-h,o-z)
      
#include "fortdim.h"
      common /grid/ h,r1,z,r(npdm),np      
      common /aorb/ ncore,nval,norb,nlm(n0),
     $     no(n0),lo(n0),nmax(n0),maxim
      common /adat/ wnl(n0),en(n0),so(n0),xion
      common /aval/ rcall(n0),rvap(n0),rnorm(n0),ibd(n0),etot
      dimension ilocalr(n0),rnl(npdm,n0)
      dimension rhead(n0),wavea(npdm)

      do korb=1,nval
         do k = 1,ilocalr(korb)
            wavea(k) = rnl(k,korb) * rnl(k,korb)
         enddo
         pow = 2 * lo(korb) + 2
         call radin(r,wavea,0,ilocalr(korb),h,pow)
         rhead(korb)=pow
      enddo

      return
      end
