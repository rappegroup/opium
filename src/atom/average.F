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
      subroutine average(ixc,exccut,iexit)
      
c     *************************************************************************
c     This code computes the j-averaged wfns from the all-electron ones
c     This is procedure comes from I. Grinbeg, N. J. Ramer, and A. M. Rappe, 
c     PRB 62, 2311 (1999). Referred to as "GRR" in this routine.

c     This routine was originally written by I. Grinberg and cleaned up by EJW
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

      common /iavgtype/ iavgtype

      dimension rvxct(npdm),rexct(npdm),v(npdm)
      dimension rvxt(npdm),rext(npdm)
      dimension rvct(npdm),rect(npdm)
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

      irelxc=0
      zv=xion
      do iorb=ncore+1,norb
         zv=zv+wnl(iorb)
      enddo
      write (7,9000) z
      write (7,9010) zv
      write (7,9015)
      write(7,*)

c     Keep a copy of rel eigenvalues, eigenvalues come into 
c     this routine in "en" and are copied into enrel.

c     "lo" is an array of the "l" values for each orbital

c     first get the core eigenvalues
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

c     now get the valence eigenvalues
      do k=ncore+1,norb
         do j=1,2
            if (lo(k).gt.0.or.j.eq.1) then
               norb_rel=norb_rel+1
               enrel(norb_rel)=en(norb_rel)
            endif
         enddo
      enddo
      nval_rel=norb_rel-ncore_rel

c     set maxim (virtual infinity for integrals) 
c     to be based on when charge density vanished
      do j=1,np
         if (rsval(j).gt.0.0) maxim=j
      enddo

      minrc=1e4
      kk=ncore_rel
      do i=ncore+1,norb
         ic=i-ncore
         kk=kk+1

         do j = 1,np
            if (rcall(ic).gt.r(j)) indall(ic) = j
         enddo
         
         minrc=min(indall(ic),minrc)
         
c     if we are NOT an 's' state, then grab next wfn         
         if (lo(i).ne.0) kk=kk+1

         do j=1,maxim
            phom(j)=rnla(j,kk)**2
         enddo
         pow=2.d0*lo(i)+2
         call radin(r,phom,0,maxim,h,pow)
      enddo
         
c     Step #1: get a trial j-average of the wavefns and set the
c     norm from 0->rc (=indall) and 0->"breakpoint" = ibp = MIN(over all rc) - 5 grid points
c     This section gets equations 3 and 4 from "GRR".

      kk=ncore_rel
      do i=ncore+1,norb
         ic=i-ncore
         
         kk=kk+1
         fe_old(ic)=0.0

c     get the j weights
         xt=4.d0*float(lo(i))+2.d0
         xu=(2.d0*float(lo(i))+2.d0)/xt
         xd=(2.d0*float(lo(i)))/xt

c     get j+ contribution
         do j=1,maxim
            rhowf(j)=rnla(j,kk)**2*xd
         enddo

c     if this an unbound state skip this next section

c     if we are NOT an 's' state, then grab next wfn
         if (lo(i).ne.0) kk=kk+1

c     add in the j- component
         do j=1,maxim
            rhowf(j)=rhowf(j)+rnla(j,kk)**2*xu
         enddo

         if (ibd(ic).eq.0) goto 934

c     compute the integral from 0 -> oo
         pow=2*lo(i)+2
         call radin(r,rhowf,0,maxim,h,pow)
         totalint=pow
         do j = 1,np
            if (rcall(ic).gt.r(j)) indall(ic) = j
         enddo

c     OLD way the can cause problems if rc's are not all the same
         if (iavgtype.eq.0) then
            ibp(ic) = indall(ic) - 5         
         else
c     NEW way, use the min of all rc's
            ibp(ic) = minrc-5
         endif

c     compute the integral from 0 -> rc 
         pow=2*lo(i)+2
         call radin(r,rhowf,0,indall(ic),h,pow)
         rnorm(ic) = 1- pow
         write (7,9050) no(i),xc(lo(i)+1),rnorm(ic),r(ibp(ic))

c     compute the integral from 0 -> "breakpoint"
         pow=2*lo(i)+2
         call radin(r,rhowf,0,ibp(ic),h,pow)
         powall = 1- pow
         
         rcint=pow
         partint=abs(totalint-rcint)

c     now we store the partial norm from breakpoint -> oo for each 'l'
         fe_old(ic)=partint
         
 934     continue
      enddo

c     Step #2: find the initial  "c_nl" from equation 5 in "GRR"
      kk=ncore_rel
      do i=ncore+1,norb
         ic=i-ncore

c     get j-averaged wfns again

         xt=4.d0*float(lo(i))+2.d0
         xu=(2.d0*float(lo(i))+2.d0)/xt
         xd=(2.d0*float(lo(i)))/xt

         kk=kk+1
         do j=1,maxim
            rnlt(j,ic)=rnla(j,kk)**2*xd
         enddo
         if (lo(i).ne.0) kk=kk+1
         do j=1,maxim
            rnlt(j,ic)=rnlt(j,ic)+rnla(j,kk)**2*xu
         enddo
         do j=1,maxim
            rnlt(j,ic)=sqrt(rnlt(j,ic))
         enddo

c     skip ahead if unbound 
         if (ibd(ic).eq.0) goto 935

c     find c_nl, using equation 5, solve |phi_nl(r)|^2 for c_nl
         do j=1,ibp(ic)

            rl=r(j)**(lo(i)+1)
c     r^(l+1)

            chiro=rnlt(ibp(ic),ic)/r(ibp(ic))**(lo(i)+1)
c     phi_nl(r_breakpoint)/(r_breakpoint^(l+1))

            oneminus=(1-(r(j)/r(ibp(ic)))**4)
c     (1-r/r_breakpoint)^4

c     multiply by r^(l+1)
            temp1(j)=(rl*chiro)**2
            temp2(j)=(rl*oneminus)**2

c     |phi_nl(r)|^2
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

c     find roots
         const1=(-temp_intx+sqrt(temp_intx**2-4*temp_int2
     $        *(temp_int1-1+fe_old(ic))))/2.d0/temp_int2
         
         const2=(-temp_intx-sqrt(temp_intx**2-4*temp_int2
     $        *(temp_int1-1+fe_old(ic))))/2.d0/temp_int2
         c=const1
         if(const1.ge.0.0) c=const1
         if(const2.ge.0.0) c=const2

c     now apply c_nl and find trial phi_nl(r)
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

c     collect new rho(r) (rsval)
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

c     Section #3:  Find set of {c_nl}'s s.t. norm from 0->r_breakpoint
c     is the same as it is for the initial j-averaged wavefunctions


      do ijk=1,250
         do j=1,np
            rsatom(j)=rsval(j)+rscore(j)
            rvxct(j) = 0.0
            rexct(j) = 0.0
            rvh(j) = 0.0
         enddo

c     3a: Find new xc and hartree potentials since rho has changed
         call excorr (np,ixc,exccut,rsatom,rvxct,rexct,
     $        rvxt,rext,rvct,rect,irelxc)
         call hrtree (np,h,r,rsatom,rvh)

c     rvscr is the screening potential
c     rv is the effective potential
         do j=1,np
            rvscr(j)=rvh(j)+rvxct(j)
            rv(j)=-2*z+rvscr(j)
         enddo

c     3b: Do inward integration to find new wfns using the new effective
c     potential and the reference eigenvalues
         kk=ncore_rel
         do i=ncore+1,norb
            ic=i-ncore
            kk=kk+1

c     Skip ahead if unbound
            if (ibd(ic).eq.0) then
               kk=kk+1
               if (lo(i).ne.0) kk=kk+1
               goto 933
            endif

            xt=4.d0*float(lo(i))+2.d0
            xu=(2.d0*float(lo(i))+2.d0)/xt
            xd=(2.d0*float(lo(i)))/xt

c     EJW: From here until 933, this is just code from schsl.F and should 
c     be a subroutine call NOT a duplicate of the code.  
            x2=(lo(i)+0.5)**2
            h2=h*h/12.d0

            ei=enrel(kk)*xd
            if (lo(i).ne.0) kk=kk+1
            ei=ei+enrel(kk)*xu
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

c     3c: piece together new rho and compute partial norms for each 'l'
         kk=ncore_rel
         do i=ncore+1,norb
            ic=i-ncore

c     Skip ahead if not bound
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


c     Now we have the new partial norm and original partial norm for this 'l'
c     find the ratio of the two:
            fe(ic)=partint
            ratio(ic)=fe_old(ic)/fe(ic)
            rtsqrt=sqrt(ratio(ic))

            do j=ibp(ic),np
               rnlt(j,ic)=rnlt(j,ic)*rtsqrt
            enddo
            do j=1,np
               rhowf(j)=rnlt(j,ic)*rnlt(j,ic)
            enddo

c     using the new wfns, find a new c_nl, this code should be a subroutine call 
c     since it is copied from above.

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

c     perhaps we should check this BEFORE computing the new trial wfns?
c     kind of wasteful but ok (since we are self-consistent at this point).

         rttot=0.0
         do i=1,nval
            rttot=rttot+abs((ratio(i)-1.d0))
         enddo
         if(abs(rttot).lt.(1e-12)) then 
            goto 998
         endif

         do j=1, np
            rsval(j)=0
         enddo
         do i=ncore+1,norb
            ic=i-ncore

c     compute new rho for bound states
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

c     Section #4: print out the results since we are done
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

c     recompute the partial norms one last time
      kk=ncore_rel
      do i=ncore+1,norb
         ic=i-ncore
         kk=kk+1

         xt=4.d0*float(lo(i))+2.d0
         xu=(2.d0*float(lo(i))+2.d0)/xt
         xd=(2.d0*float(lo(i)))/xt
         
         en(ic)=enrel(kk)*xd
         if (lo(i).ne.0) kk=kk+1
         en(ic)=en(ic)+enrel(kk)*xu
         
         do j=1,np
            rnl(j,ic)=rnlt(j,ic)
            rhowf(j)=rnl(j,ic)**2
         enddo
         if (ibd(ic).ne.0) then
            pow=2*lo(i)+2
            call radin(r,rhowf,0,indall(ic),h,pow)
            rnorm(ic) = 1- pow
            write (7,9058) no(i),xc(lo(i)+1),rnorm(ic)

            pow=2*lo(i)+2
            call radin(r,rhowf,0,maxim,h,pow)
         endif
         
      enddo
      
      do i = 1,np
         rsatom(i) = rsval(i) + rscore(i)
         if (rsatom(i).gt.0.0) maxim = i
      enddo
c     why recompute maxim again??

c     compute new v_xc and v_h
      pow = 2.0
      call radin(r,rsatom,0,maxim,h,pow)
      pow = 2.0
      call radin(r,rscore,0,maxim,h,pow)
      call excorr (maxim,ixc,exccut,rsatom,rvxct,rexct,
     $     rvxt,rext,rvct,rect,irelxc)
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

c     print final average eigenvalues
      do i = ncore+1,norb
         ic=i-ncore
         write(7,9070) no(i),xc(lo(i)+1),en(ic)

c     flip wfn if approaches 0 as r->oo from the negative side
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

c     This is not used later I don't think
         do k=1,np
            if (rnl(k,ic).gt.1e-18) nmax(ic)=k
         enddo

      enddo      
      
 9000 format(5x,"Z atom                            :",f10.6)
 9010 format(5x,"Z valence                         :",f10.6)
 9020 format(5x,"Total valence charge              :",f10.6)
 9015 format(5x,"Ideal tail norms and bp radius   :")

 9030 format(1x,i1,a1,2x,"- eigenvalue                      :",f10.6)
 9040 format(1x,i1,a1,2x,"+ eigenvalue                      :",f10.6)

 9050 format(1x,i1,a1,2x,"Norm rc->oo: ",f10.6,4x,"bp radius: ",f10.6)
 9058 format(1x,i1,a1,2x,"Norm rc->oo                       : ",f10.6)

 9060 format(5x,"Converged in ", i3,1x,"iterations")
 9061 format(5x,"Breakpoint norm error    :",1x,i1,a1,2x,e10.3)  
 9077 format(5x,"Breakpoint norm error    :",1x,i1,a1,2x, "-Not bound-")  

 9070 format(1x,i1,a1,2x,"Avg. eigenvalue                   : ",f10.6)

c     *************************************************************************

      close(7)
      iexit=0
      return
      
 911  continue
      close(7)
      iexit=1
      return

      end
