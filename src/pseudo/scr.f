      subroutine scr

c     *************************************************************************
c     This code finds an optimized pseudowavefunction following the paper
c       of Rappe et al. (PRB 41,1227 (1990)).
c     Modified from Andrew Rappe's subroutine scr to use the CONMAX
c       minimization routines rather than imsl's DNCONF.
c     The logarithmic derivative matching has been removed on Andrew's
c       advice.
c                                                           06/22/98. N. Hill
c     *************************************************************************

      implicit double precision (a-h,o-z)
      
#include "PARAMOPT"

      parameter(maxitn=1000)

      logical lbes
      common /roots/ xroot(numfn0)
      common /angm/ ll
      common /nmax/ nmax(n0),maxim
      common /numfn/ numfn
      common /wavrc/ wavrc, slope, curvrc, indrc
      common /nnn/ nnn
      common /grid/ h,r1,z,r(npdm),np
      common /atomic/ xion,rnl(npdm,n0),nlm(n0),wnl(n0),en(n0),norb
      common /atom3/ rsvale(npdm)
      common /bir/ bir(numfn0,npdm)
      common /frrv/ fr(npdm), rv(npdm)
      common /transum/ transumry

      common /totpot/ rvcore(npdm,n0),rvps(npdm,n0),rvcoul(npdm)
      
      common /scrconmax/ enc,tolc,lim,nstp,isw1,isw2,isw3
      
      common /convrpt/ wsumt(nvale0), lnode(nvale0), lghost(nvale0)
      
c     internal common block 
      common /xlog/ ffunc,sumo,sumt,sumn
      
      
c     *************************************************************************
c     local variables
c     *************************************************************************

      dimension fr2(npdm),xguess(numfn0),v0(npdm)
      dimension iwork(7*7 + 7*numfn0 +3)
      dimension work(2*numfn0**2+4*7*numfn0+11*7+27*numfn0+13)
      dimension error(10),fun(1),pttbl(1,1),icntyp(7),confun(7,numfn0+1)

c     *************************************************************************
c     Step 1:  Put the Bessel functions in real space into bir.
c     For a starting guess start with only the first Bessel function,
c     and make the wavefunction continuous.
c     *************************************************************************

      do i = 1,np
        rv(i) = -z-z+rvcoul(i)
        if (i.gt.maxim) rv(i) = -xion-xion
        if (nnn.eq.1) rsvale(i)=0.0
      enddo

      do j = 1,numfn
        do i = 1,np
          x = xroot(j) * r(i)
          bir(j,i) = besfn(x,ll)
        enddo
      enddo
 
      do i = 1,numfn0
        xguess(i) = 0.0
      enddo
      
      call guess(xguess)

c      do i=1,numfn
c        write(7,9011) xguess(i),xroot(i)
c      enddo
c 9011 format(1x,'Coeff guesses: ',2f20.10)
 
c     *************************************************************************
c     Start the CONMAX section
c     For info. on the CONMAX routines, see the comments in conmax.f or look at
c       http://www.netlib.org/opt/conmax.f
c     numgr is the number of constraints, which will always be equal to 7 for
c       Andrew's method (CONMAX thinks in a very strange way - the main
c       condition is treated as a constraint, and equality constraints
c       are imposed by setting + and - the expression .LE. 0.)
c     *************************************************************************

      nparm = numfn
      numgr = 7
      itlim = maxitn
      ifun = 1
      iptb = 1
      indm = 1
      liwrk = (7*numgr + 7*numfn0 +3)
      lwrk = 2*numfn0**2 + 4*NUMGR*numfn0+ 11*NUMGR + 27*numfn0+ 13
      
camr  This helps algorithm stop more promptly.
c      ioptn = 311
c      work(1) = 1.0e-5
c      iwork(1) = 5
c      work(2) = 1.0e-5
c      iwork(2) = 3

c -- gjt: now these values are set by parameters in the input file
      ioptn = 1
      if (isw1.gt.0) ioptn=ioptn+10
      if (isw2.gt.0) ioptn=ioptn+100
      if (isw3.gt.0) ioptn=ioptn+200
      work(1)  = enc
      iwork(1) = lim
      work(2)  = tolc
      iwork(2) = nstp
c      write(7,*) 'ioptn',ioptn,work(1),iwork(1),work(2),iwork(2)

cEJW added this for internal debugging (lets you load bessel functions)
      lbes = .false.
      inquire(file='BESSEL',exist=lbes)
      if (lbes) then
         open(unit=111,file='BESSEL',form='formatted')
         do i=1,nparm
            read(111,*,end=911) ijunk,xjunk,xguess(i)
         enddo
         goto 912

 911     stop'Bessel file broken!'

 912     continue
      else

         call newmin(numfn,xguess,sumt)
c         call conmax(ioptn,nparm,numgr,itlim,fun,ifun,pttbl,iptb,
c     $        indm,iwork,liwrk,work,lwrk,iter,xguess,error)
      endif
      call fnset(nparm,numgr,pttbl,iptb,indm,xguess,7,1,icntyp,confun)
     
c      write(7,800)iter,liwrk,lwrk
  800 FORMAT(/26H *****AFTER CONMAX ITER IS,I4,10H  LIWRK IS,I5,
     *9H  LWRK IS,I6)

      write(7,*) 'Bessel wavevectors(sqrt(Ry)) and coefficients'
      do i=1,nparm
        write(7,9012) i,xroot(i),xguess(i)
      enddo
 9012 format(1x,i3,2f20.10)

      write(7,*)
      write(7,*) 'Error in constraints'
      do i=1,numgr
        write(7,9013) i,error(i)
      enddo
      write(7,9013) numgr+1,error(numgr+1)
      write(7,9013) numgr+2,error(numgr+2)
      write(7,9013) numgr+3,error(numgr+3)

 9013 format(1x,i3,e20.5)

c     *************************************************************************
c     end CONMAX section.
c     *************************************************************************

      write(7,*)
      write(7,9102) sumt
      write(7,9132) wnl(nnn)
      wsumt(nnn) = sumt*wnl(nnn)
      write(7,9133) wsumt(nnn)*1000.0
      write(7,9134) wsumt(nnn)*13.6058*1000.0
      transumry=transumry+wsumt(nnn)

 9100 format(1x,'Minimized Function: ',f20.10)
 9101 format(1x,'Log deriv term    : ',f20.10)

 9102 format(1x,'Convergence term (Ry)    : ',f20.10)
 9132 format(1x,'Occupation               : ',f20.10)
 9133 format(1x,'Convergence error (mRy)  : ',f20.10)
 9134 format(1x,'Convergence error (meV)  : ',f20.10)
 9103 format(1x,'Norm conserv term : ',f20.10)

c     *************************************************************************
c     Invert the wavefunction 
c     *************************************************************************
      
      call optinvert(fr,r,h,ll,en(nnn),indrc,np,maxim,xroot,xguess,
     $     numfn,xion,v0)

c     *************************************************************************
c     Step 5: Now put the potential into rvps and rnl and rsvale.
c     *************************************************************************

      inode = 0
      do i = 1,np
c        if (fr(i).lt.0.0.and.i.le.indrc) then
c          inode = 1
c        endif
        fr2(i) = fr(i) * fr(i) * r(i) * r(i)
      enddo
c      lnode(nnn) = 0
      if (inode.eq.1) then
c       lnode(nnn) = 1
        write (7,*) 'This wavefunction has a node.  Stop.  nnn = ',nnn
        stop
      endif
      
      tch = float(ll+ll+2)
      call radin(r,fr2,0,maxim,h,tch)
 
      tch = sqrt(tch)
      do i = 1,np
        fr(i) = fr(i)/tch
      enddo

c      do i = 1,np
c        fr2(i) = fr(i) * fr(i) * r(i) * r(i)
c      enddo
c      tch = float(ll+ll+2)
c      call radin(r,fr2,0,maxim,h,tch)
 
      do i = 1,np
        rnl(i,nnn) = fr(i) * r(i)
        rsvale(i) = rsvale(i) + wnl(nnn) * rnl(i,nnn)**2
        if (i.gt.indrc+50) rvps(i,nnn) = -z-z+rvcoul(i)
        if (i.gt.maxim) rvps(i,nnn) = -xion-xion
      enddo
            
c      tch = float(ll+ll+2)
c      call radin(r,rsvale,0,maxim,h,tch)
c      write(7,9202) tch,wnl(nnn)
c 9202 format(1x,'Total charge so far is: ',2f20.10)

      return
      end
      
      
c     #########################################################################
      
      
      subroutine fnset(
     $          NPARM,NUMGR,PTTBL,IPTB,INDM,x,IPT,INDFN,ICNTYP,CONFUN)
     
c     *************************************************************************
c     Subroutine to calculate the primary and secondary minimization
c       constraints to send to CONMAX package.
c     Modified from Andrew Rappe's subroutine fcn2 in file scr.f
c                                                             06/22/98. N. Hill
c     *************************************************************************

      implicit double precision(a-h,o-z)
      
#include "PARAMOPT"


      common /bir/ bir(numfn0,npdm)
      common /numfn/ numfn
      common /wavrc/ wavrc, slope, curvrc, indrc
      common /nnn/ nnn
      common /grid/ h,r1,z,r(npdm),np
      common /atomic/ xion,rnl(npdm,n0),nlm(n0),wnl(n0),en(n0),norb
      common /angm/ ll
      common /frrv/ fr(npdm), rv(npdm)
      common /roots/ xroot(numfn0)
      common /a/ a(numfn0)
      common /d/ d(numfn0,numfn0)
      common /e/ e(numfn0)
      common /f/ fint
      common /g/ gint

c     internal common block      
      common /xlog/ ffunc,sumo,sumt,sumn

c     *************************************************************************
c     local variables
c     *************************************************************************

      dimension x(nparm), pttbl(iptb,indm)
      dimension icntyp(numgr), confun(numgr,nparm+1)
      dimension fr2(npdm)
      
c     *************************************************************************
c     Find the function to be minimized, including kinetic, and
c       norm-conservation terms. (Note that the logarithmic derivative
c       terms have been removed since Andrew advised against using them!)
c     x(j) are the Bessel functions alpha_i. In contrast to the original
c       paper, the pseudo-wavefunction is not split into F and C parts.
c     *************************************************************************

c     *************************************************************************
c     First calculate the pseudo-wavefunction (sum) and its Laplacian (sum2)
c       at rc, qi. (xroot is qi - the root of the transcendental equation).
c     *************************************************************************

      sum = 0.0
      sum2 = 0.0
      do j = 1,nparm
        sum = sum + x(j) * bir(j,indrc)
        sum2 = sum2 + x(j) * bir(j,indrc) * xroot(j) * abs(xroot(j))
c        write(7,*) j, x(j), bir(j,indrc), xroot(j)
      enddo
 
c -- gjt: this is debugging stuff ---
c      write(7,*) 'sum, sum2', sum, sum2    
c      stop

c     *************************************************************************
c     The first two constraints are the continuity of the wavefunction and
c       its derivatives. This is achieved by requiring that the pseudo-
c       and all-electron terms are equal at rc. For the weird formatting
c       of the CONMAX routine, this is achieved by setting the difference
c       and minus the difference each .LE. 0.0
c     *************************************************************************

      g2 = sum - wavrc
      g3 = sum2 + curvrc
      if (ipt.eq.1) then
        confun(ipt,1)=g2
        icntyp(ipt)=-2
        return
      end if
      if (ipt.eq.2) then
        confun(ipt,1)=-g2
        icntyp(ipt)=-2
        return
      end if
      if (ipt.eq.3) then
        confun(ipt,1)=g3
        icntyp(ipt)=-2
        return
      end if
      if (ipt.eq.4) then
        confun(ipt,1)=-g3
        icntyp(ipt)=-2
        return
      end if

c -- gjt: this is debugging stuff ---
c      write(7,*) 'g2, g3', g2, g3

c     *************************************************************************
c     Now construct the trial pseudo-wavefunction fr over the whole grid.
c     *************************************************************************

c     Up to rc this is the sum over Bessel functions:

      f = 0.0
      do i = 1,indrc
        fr(i) = 0.0
        do j = 1,numfn
          fr(i) = fr(i) + x(j) * bir(j,i)
        enddo
        fr2(i) = fr(i) * fr(i) * r(i) * r(i)
      enddo
      
c     Beyond rc this is the all-electron wavefunction:

      do i = indrc+1,np
        fr(i) = rnl(i,nnn)
        fr2(i) = fr(i) * fr(i) * r(i) * r(i)
      enddo
      
      xn9 = float(ll+ll+2)
      call radin(r,fr2,0,np,h,xn9)
      
c -- gjt: this is debugging stuff ---
c      do i=1, 20
c        write(7,*) 'fr2(',i,')=',fr2(i)
c      enddo
c      write(7,*) 'xn9=',xn9
      
c     *************************************************************************
c     xn9 is the norm of fr - this creates an additional constraint
c       (xn9 must equal 1!)
c     *************************************************************************

      g1 = xn9 - 1.0
      sumn = g1**2
      if (ipt.eq.5) then
        confun(ipt,1)=g1
        icntyp(ipt)=-2
        return
      end if
      if (ipt.eq.6) then
        confun(ipt,1)=-g1
        icntyp(ipt)=-2
        return
      end if

c     *************************************************************************
c     Here we find the part kinetic energy.
c     *************************************************************************

      sumt = 0.0
      do i = 1,numfn
        sum = 0.0
        do j = 1,numfn
          sum = sum + x(j) * d(i,j)
        enddo
        
c       Here xroot*abs(xroot) again corrects for xroot(1)<0 exception.
        fact = x(i) * xroot(i)*abs(xroot(i)) * a(i) - 2.0 * e(i) - sum
     
        sumt = sumt + fact * x(i)
      enddo
 
      sumt = sumt - fint - gint
      sumo = f
      f = f + sumt 

      
      if (f.lt.0.0) then
        write (7,*) 'less than zero',f
        write (7,*) fint,gint,(x(i),i=1,numfn),(xroot(i),i=1,numfn),
     $              (e(i),i=1,numfn)
        stop
      endif
      
      ffunc = f

c -- gjt: this is debugging stuff ---
c      write(7,*) 'stuff=',sumt,wke,fint,gint,ffunc

      if (ipt.eq.7) then
        confun(ipt,1)=ffunc
        icntyp(ipt)=1
      end if

      return
      end
