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
      subroutine average(ixc,exccut,iexit)
      
c     *************************************************************************
c     This code reads in the relativ. atomic information from fort.22. and 
c     processes it for the scalar-relativistic pseudopotential procedure.
c     The file *.ae is generated.
c     *************************************************************************

      implicit double precision (a-h,o-z)
      
      external val
      
#include "fortdim.h"

      common /grid/ h,r1,z,r(npdm),np
      common /aorb/ ncore,nval,norb,nlm(n0),
     $     no(n0),lo(n0),nmax(n0),maxim
      common /adat/ wnl(n0),en(n0),so(n0),xion
      common /aval/ rcall(n0),rvap(n0),rnorm(n0),ibd(n0),etot

      common /wfn/ rnl(npdm,n0)
      common /wfnrel/ rnla(npdm,n0),rnlb(npdm,n0)

      common /rscore/ rscore(npdm),rdd(npdm),rddd(npdm),rscoretot(npdm),
     $     rsval(npdm)
      common /totpot/ rvcore(npdm,n0),rvps(npdm,n0),rvcoul(npdm)

      dimension rvxct(npdm),rexct(npdm),v(npdm)
      dimension ibp(n0),indall(n0),rsatom(npdm),rv(npdm)
      dimension fe_old(n0),fe(n0), rnlt(npdm,n0)
      dimension rvscr(npdm),rvh(npdm),phom(npdm)
      dimension xah(npdm),xdh(npdm),xch(npdm),xlh(npdm)
      dimension rhowf(npdm),temp1(npdm),temp2(npdm)
      dimension ratio(n0),tempx(npdm)

      dimension enrel(n0)
      common /filenames/ file_log

      character*80 file_log
      character*1 xc(4)
      
      xc(1)='s'
      xc(2)='p'
      xc(3)='d'
      xc(4)='f'

      open(unit=7,file=file_log,form='formatted',access='append')

      write(7,*) "-----------------------------------------------"
      write(7,*) "Starting j-averaging procedure                 "
      write(7,*) "-----------------------------------------------"

      zv=xion
      do iorb=ncore+1,norb
         zv=zv+wnl(iorb)
      enddo
      write (7,9000) z
      write (7,9010) zv
      write (7,9015)
 9015 format(5x,"Ideal tail norms: ...")

      norb_rel=0
      do k=1,ncore
         do j=1,2
            if (lo(k).gt.0.or.j.eq.1) then
               norb_rel=norb_rel+1
               enrel(norb_rel)=en(norb_rel)
            endif
         enddo
      enddo

      ncore_rel=norb_rel
      
      do k=ncore+1,norb
         do j=1,2
            if (lo(k).gt.0.or.j.eq.1) then
               norb_rel=norb_rel+1
               enrel(norb_rel)=en(norb_rel)
            endif
         enddo
      enddo
      nval_rel=norb_rel-ncore_rel

      do j=1,np
         if (rsval(j).gt.0.0) maxim=j
         
      enddo

      kk=ncore_rel
      do i=ncore+1,norb
         ic=i-ncore
         kk=kk+1
         
         if (lo(i).ne.0) kk=kk+1

         do j=1,maxim
            phom(j)=rnla(j,kk)**2
         enddo
         pow=2.d0*lo(i)+2
         call radin(r,phom,0,maxim,h,pow)
      enddo
         

      kk=ncore_rel
      do i=ncore+1,norb
         ic=i-ncore
         kk=kk+1
         fe_old(ic)=0.0
         
         xt=4.d0*float(lo(i))+2.d0
         xd=(2.d0*float(lo(i))+2.d0)/xt
         xu=(2.d0*float(lo(i)))/xt

         do j=1,maxim
            rhowf(j)=rnla(j,kk)**2*xu
         enddo
         
         if (lo(i).ne.0) kk=kk+1

         do j=1,maxim
            rhowf(j)=rhowf(j)+rnla(j,kk)**2*xd
         enddo
         if (ibd(ic).eq.0) goto 934

         pow=2*lo(i)+2
         call radin(r,rhowf,0,maxim,h,pow)
         totalint=pow
         do j = 1,np
            if (rcall(ic).gt.r(j)) indall(ic) = j
         enddo
         ibp(ic) = indall(ic) - 5
         
         pow=2*lo(i)+2
         call radin(r,rhowf,0,indall(ic),h,pow)
         rnorm(ic) = 1- pow
         write (7,9050) no(i),xc(lo(i)+1),rnorm(ic)
         pow=2*lo(i)+2
         call radin(r,rhowf,0,ibp(ic),h,pow)
         powall = 1- pow
         
         rcint=pow
         partint=abs(totalint-rcint)
         fe_old(ic)=partint
 934     continue
      enddo
      
      kk=ncore_rel
      do i=ncore+1,norb
         ic=i-ncore
         kk=kk+1
         do j=1,maxim
            rnlt(j,ic)=rnla(j,kk)**2*xu
         enddo
         if (lo(i).ne.0) kk=kk+1
         do j=1,maxim
            rnlt(j,ic)=rnlt(j,ic)+rnla(j,kk)**2*xd
         enddo
         do j=1,maxim
            rnlt(j,ic)=sqrt(rnlt(j,ic))
         enddo

         if (ibd(ic).eq.0) goto 935

         do j=1,ibp(ic)
            rl=r(j)**(lo(i)+1)
            chiro=rnlt(ibp(ic),ic)/r(ibp(ic))**(lo(i)+1)
            oneminus=(1-(r(j)/r(ibp(ic)))**4)
            temp1(j)=(rl*chiro)**2
            temp2(j)=(rl*oneminus)**2
            tempx(j)=2*rl*rl*chiro*oneminus
         enddo
         pow=2*lo(i)+2
         call radin(r,temp1,0,ibp(ic),h,pow)
         temp_int1=pow
         pow=2*lo(i)+2
         call radin(r,temp2,0,ibp(ic),h,pow)
         temp_int2=pow
         pow=2*lo(i)+2
         call radin(r,tempx,0,ibp(ic),h,pow)
         temp_intx=pow
         const1=(-temp_intx+sqrt(temp_intx**2-4*temp_int2
     $        *(temp_int1-1+fe_old(ic))))/2.d0/temp_int2
         
         const2=(-temp_intx-sqrt(temp_intx**2-4*temp_int2
     $        *(temp_int1-1+fe_old(ic))))/2.d0/temp_int2
         c=const1
         if(const1.ge.0.0) c=const1
         if(const2.ge.0.0) c=const2
         do j=1,ibp(ic)
            rl=r(j)**(lo(i)+1)
            chiro=rnlt(ibp(ic),ic)/r(ibp(ic))**(lo(i)+1)
            oneminus=(1-(r(j)/r(ibp(ic)))**4)
            rnlt(j,ic)=rl*(chiro+c*oneminus)
         enddo
 935     continue
      enddo
      
      c=1.d0
      do j=1, np
         rsval(j)=0
      enddo
      
      kk=ncore_rel
      do i=ncore+1,norb
         ic=i-ncore
         if (ibd(ic).ne.0) then
            do j=1, np
               rhowf(j)=0
            enddo
            
            do j=1,np
               rsval(j)=rsval(j)+rnlt(j,ic)**2*wnl(i)
            enddo
         endif
      enddo
      write (7,*)
      
c     eeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee
c     jjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjj
c     wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww

      do ijk=1,250
         do j=1,np
            rsatom(j)=rsval(j)+rscore(j)
            rvxct(j) = 0.0
            rexct(j) = 0.0
            rvh(j) = 0.0
         enddo

         irel=0
         call excorr (np,ixc,exccut,r,rsatom,rvxct,rexct,irel)
         call hrtree (np,h,r,rsatom,rvh)
         do j=1,np
            rvscr(j)=rvh(j)+rvxct(j)
            rv(j)=-2*z+rvscr(j)
         enddo

         kk=ncore_rel
         do i=ncore+1,norb
            ic=i-ncore
            kk=kk+1
            if (ibd(ic).eq.0) goto 933

            xt=4.d0*float(lo(i))+2.d0
            xd=(2.d0*float(lo(i))+2.d0)/xt
            xu=(2.d0*float(lo(i)))/xt
            
            x2=(lo(i)+0.5)**2
            h2=h*h/12.d0

            ei=enrel(kk)*xu
            if (lo(i).ne.0) kk=kk+1
            ei=ei+enrel(kk)*xd

            do j = 1, np
               v(j) = ((rv(j) - ei * r(j)) * r(j) + x2) * h2
            enddo
            do j = 1,np
               if (rsval(j).gt.0.0) maxim = j
            enddo
            imax=maxim
            imat=ibp(ic)+1
            imp1=imat+1
            xah(imat-1) = 1 - v(imat-1)
            xdh(imat-1) = -(2.d0 + 10.d0 * v(imat-1))
            xah(imat) = 1 - v(imat)
            xdh(imat) = -(2.d0 + 10.d0 * v(imat))
            xlh(imat) = xdh(imat) 
            xch(imat) = -(xah(imat-1)) * rnlt(imat-1,ic)/sqrt(r(imat-1))
            do j = imp1,imax+1
               xah(j) = 1.d0 - v(j)
               xch(j) = -xch(j-1) * xah(j-1)/xlh(j-1)
               xdh(j) = -(2.d0 + 10.d0 * v(j))
               xlh(j) = xdh(j) - xah(j) * xah(j-1)/xlh(j-1)
            enddo
            rootv = h * sqrt(v(imax)/h2)
            ap = exp(-rootv)
            rnlt(imax,ic) = xch(imax)/(xlh(imax) + xah(imax+1) * ap)
            do j = imax-1,imat,-1
               rnlt(j,ic) = (xch(j) - xah(j+1) * rnlt(j+1,ic))/xlh(j)
            enddo
 933        continue
         enddo
         
         kk=ncore_rel
         do i=ncore+1,norb
            ic=i-ncore
            if (ibd(ic).eq.0) goto 937
            fe(ic)=0
            do j=1,ibp(ic)-1
               rhowf(j)=rnlt(j,ic)*rnlt(j,ic)
            enddo
            do j=ibp(ic)+1,np
               rnlt(j,ic)=rnlt(j,ic)*sqrt(r(j))
            enddo
            do j=ibp(ic),np
               rhowf(j)=rnlt(j,ic)*rnlt(j,ic)
            enddo

            pow=2*lo(i)+2
            call radin(r,rhowf,0,np,h,pow)
            totalint=pow
            pow=2*lo(i)+2
            call radin(r,rhowf,0,ibp(ic),h,pow)
            rcint=pow
            partint=abs(totalint-rcint)

            fe(ic)=partint
            ratio(ic)=fe_old(ic)/fe(ic)
            rtsqrt=sqrt(ratio(ic))

            do j=ibp(ic),np
               rnlt(j,ic)=rnlt(j,ic)*rtsqrt
            enddo
            do j=1,np
               rhowf(j)=rnlt(j,ic)*rnlt(j,ic)
            enddo

            do j=1,ibp(ic)
               rl=r(j)**(lo(i)+1)
               chiro=rnlt(ibp(ic),ic)/r(ibp(ic))**(lo(i)+1)
               oneminus=(1-(r(j)/r(ibp(ic)))**4)
               temp1(j)=(rl*chiro)**2
               temp2(j)=(rl*oneminus)**2
               tempx(j)=2*rl*rl*chiro*oneminus
            enddo

            pow=2*lo(i)+2
            call radin(r,temp1,0,ibp(ic),h,pow)
            temp_int1=pow
            pow=2*lo(i)+2
            call radin(r,temp2,0,ibp(ic),h,pow)
            temp_int2=pow
            pow=2*lo(i)+2
            call radin(r,tempx,0,ibp(ic),h,pow)
            temp_intx=pow

            const1=(-temp_intx+sqrt(temp_intx**2-4*temp_int2
     $           *(temp_int1-1+fe_old(ic))))/2/temp_int2
            const2=(-temp_intx-sqrt(temp_intx**2-4*temp_int2
     $           *(temp_int1-1+fe_old(ic))))/2/temp_int2
            c=const1
            if(const1.ge.0.0) c=const1
            if(const2.ge.0.0) c=const2

            do j=1,ibp(ic)
               rl=r(j)**(lo(i)+1)
               chiro=rnlt(ibp(ic),ic)/r(ibp(ic))**(lo(i)+1)
               oneminus=(1.d0-(r(j)/r(ibp(ic)))**4)
               rnlt(j,ic)=rl*(chiro+c*oneminus)
            enddo
 937        continue
         enddo

         rttot=0.0
         do i=1,nval
            rttot=rttot+abs((ratio(i)-1.d0))
         enddo
         if(abs(rttot).lt.(1e-12)) then 
            goto 998
         endif

c     EJW:  now we exit when each norm error is < 1e-13
c         iexit=0
c         do i=1,nval
c            if (abs(ratio(i)-1.d0).lt.1d-13) iexit=iexit+1
c         enddo
c         if (iexit.eq.nval) goto 998
c     EJW

         do j=1, np
            rsval(j)=0
         enddo
         do i=ncore+1,norb
            ic=i-ncore
            if (ibd(ic).ne.0) then
               do j=1,np
                  rsval(j)=rsval(j)+rnlt(j,ic)**2*wnl(i)
               enddo
            endif
         enddo
      enddo
      
c     *************************************************************************
c     done with iterative solution and 
c     *************************************************************************

 998  write(7,9060) ijk

      do i=ncore+1,norb
         ic=i-ncore
         if (ibd(ic).ne.0) then
            write (7,9061) no(i),xc(lo(i)+1),abs(ratio(ic)-1.d0)
         else
            write (7,9077) no(i),xc(lo(i)+1)
         endif
      enddo

      write(7,*)
      write (7,9016)
 9016 format(5x,"Final averaged tail norms and eigenvalues: ...")

      kk=ncore_rel
      do i=ncore+1,norb
         ic=i-ncore
         kk=kk+1

         xt=4.d0*float(lo(i))+2.d0
         xd=(2.d0*float(lo(i))+2.d0)/xt
         xu=(2.d0*float(lo(i)))/xt
         
         en(ic)=enrel(kk)*xu
         if (lo(i).ne.0) kk=kk+1
         en(ic)=en(ic)+enrel(kk)*xd
         
         do j=1,np
            rnl(j,ic)=rnlt(j,ic)
            rhowf(j)=rnl(j,ic)**2
         enddo
         if (ibd(ic).ne.0) then
            pow=2*lo(i)+2
            call radin(r,rhowf,0,indall(ic),h,pow)
            rnorm(ic) = 1- pow
            write (7,9050) no(i),xc(lo(i)+1),rnorm(ic)

            pow=2*lo(i)+2
            call radin(r,rhowf,0,maxim,h,pow)
         else
c            write(7,*) 'NB2'
         endif
         
      enddo
      
      do i = 1,np
         rsatom(i) = rsval(i) + rscore(i)
         if (rsatom(i).gt.0.0) maxim = i
      enddo

      pow = 2.0
      call radin(r,rsatom,0,maxim,h,pow)
      pow = 2.0
      call radin(r,rscore,0,maxim,h,pow)
      irel=0
      call excorr (maxim,ixc,exccut,r,rsatom,rvxct,rexct,irel)
      call hrtree (maxim,h,r,rsatom,rvh)
      do  i = maxim+1,np
        rvh(i)  = z + z - xion - xion
      enddo
      xm = (z+z-xion-xion)/rvh(maxim)
      dif1 = abs(xm-1.0)

      if (dif1.gt.tol2) then
        write (7,*) 'Use more grid points. Potentials inaccurately '
        write (7,*) 'integrated from charge density. ',tol2
        write (7,*) dif1,(z+z-xion-xion),rvh(maxim),maxim 
      endif
      sum=0.0
      sum2=0.0
      do  i = 1,np
         rvh(i) = rvh(i) * xm
         rvcoul(i) = rvh(i) + rvxct(i)
         f1 = -z-z+rvcoul(i)
         do  j = 1,nval
            rvcore(i,j) = -2 * z
            rvps(i,j) = f1
         enddo
      enddo
      tch=2.d0

      do i = ncore+1,norb
         ic=i-ncore
         write(7,9070) no(i),xc(lo(i)+1),en(ic)
         
         do k=np,1,-1
            if (abs(rnl(k,ic)).gt.1e-10) then
               if (rnl(k,ic).lt.0.0) then
                  do j = 1,np
                     rnl(j,ic) = -rnl(j,ic)
                  enddo
               endif
               goto 770
            endif
         enddo
 770     continue

         do k=1,np
            if (rnl(k,ic).gt.1e-18) nmax(ic)=k
         enddo

      enddo

      
      
 9000 format(5x,"Z atom                            :",f10.6)
 9010 format(5x,"Z valence                         :",f10.6)
 9020 format(5x,"Total valence charge              :",f10.6)

 9030 format(1x,i1,a1,2x,"- eigenvalue                      :",f10.6)
 9040 format(1x,i1,a1,2x,"+ eigenvalue                      :",f10.6)

 9050 format(1x,i1,a1,2x,"Norm rc->oo                       :",f10.6)

 9060 format(5x,"Converged in ", i3,1x,"iterations")
 9061 format(5x,"Breakpoint norm error    :",1x,i1,a1,2x,e10.3)  
 9077 format(5x,"Breakpoint norm error    :",1x,i1,a1,2x, "-Not bound-")  

 9070 format(1x,i1,a1,2x,"Avg. eigenvalue                   :",f10.6)

c     *************************************************************************

      close(7)
      iexit=0
      return
      
 911  continue
      close(7)
      iexit=1
      return

      end
