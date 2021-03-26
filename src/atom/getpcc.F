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
      subroutine getpcc(iexit)
      implicit double precision(a-h,o-z)
      
#include "fortdim.h"
      
      parameter (n=100)

c -------------------------------------------------------------------------
c     External (shared between C and Fortran) common blocks
c -------------------------------------------------------------------------
      common /grid/ h,r1,z,r(npdm),np
      common /rpcc/ rpcc,rpccz
      common /rscore/ rscore(npdm),rdd(npdm),rddd(npdm),rscoretot(npdm),
     $     rsval(npdm)
      common /filenames/ file_log
      common /ipccmeth/ ipccmeth
c -------------------------------------------------------------------------

      dimension rdensity(npdm)

      character*80 file_log

      open(unit=7,file=file_log,form='formatted',access='append')

      iexit=0
      pi=acos(-1.0)                                                   
      pi4=4.0*pi                                                      

      rsc2=0.0
      call radin(r,rscore,0,np,h,rsc2)

      irc=nint(log(rpcc/r(1))/h) +1

      do i=1,np
         rscoretot(i)=rscore(i)
         rdensity(i)=rscore(i)/(r(i)*r(i))
      enddo

      if (ipccmeth.eq.0) then
         write(7,*) 'Using LFC pcc scheme'
         call LFC(rdensity,rscore)
      else
         write(7,*) 'Using Fuchs pcc scheme'
         call fhipcc(rdensity,rscore,iexit)
         if (iexit.ne.0) goto 911
      endif

      do i=1,np-10
         x3old = 0.0
         diff3o = 1.0e-8
         do m = 3,10
            x2 = val2(rscore,r,np,r(i),m)
            diff2 = abs(x2-x2old)
            if (diff2.lt.diff2o) rdd(i) = x2
            x2old = x2
            diff2o = diff2
         enddo
         
         x3old = 0.0
         diff3o = 1.0e-8
         do m = 3,10
            x3 = val3(rscore,r,np,r(i),m)
            diff3 = abs(x3-x3old)
            if (diff3.lt.diff3o) rddd(i) = x3
            x3old = x3
            diff3o = diff3
         enddo
      enddo

      do i=1,np
         rscore(i)=rscore(i)*r(i)*r(i)
      enddo

      write(7,9010) rsc2

      write(7,9009) rpcc
      rsc=0.0
      call radin(r,rscore,0,np,h,rsc)
      write(7,9011) rsc
      rpccz = rsc

 9009 format(1x,'partial core radius : ', f10.4)
 9010 format(1x,'total core charge   : ', f20.10)
 9011 format(1x,'partial core charge : ', f20.10)

 911  continue

      close(7)
      return
      end

      subroutine LFC(rcore,rpcore)
      implicit double precision(a-h,o-z)
      
#include "fortdim.h"
      
      parameter (pi=3.141592653589793239)
      parameter (n=100)

c -------------------------------------------------------------------------
c     External (shared between C and Fortran) common blocks
c -------------------------------------------------------------------------
      common /grid/ h,r1,z,r(npdm),np
      common /rpcc/ rpcc,rpccz
      common /filenames/ file_log
c -------------------------------------------------------------------------

      dimension rleft(n),right(n),b(n)
      dimension rcore(npdm),rpcore(npdm)

      character*80 file_log

      ircpoint=nint(log(rpcc/r(1))/h) +1

c  get the charge density and the gradient of the charge density 
      rcgrad=(rcore(ircpoint+1)-rcore(ircpoint-1))
     $     /(r(ircpoint+1)-r(ircpoint-1)) 
      
      rconst=1.0/r(ircpoint)+rcgrad/rcore(ircpoint)

c--------------------------------------------------------------------------------------
c the following is solving for A and B of the SGL'moved because of error.
c
c--------------------------------------------------------------------------------------            
      if (rconst.gt.0) then
         
         right(1)=pi/(2.0*r(ircpoint))
         rleft(1)=0.0
         b(1)=(right(1)+rleft(1))/2.0
         
         do j=1,99
            
            eps=1/TAN(b(j)*r(ircpoint))-rconst/b(j)

            if  (eps.gt.0) then
               right(j+1)=right(j)
               rleft(j+1)=b(j)
               b(j+1)=(right(j+1)+rleft(j+1))/2.0

            else 
               right(j+1)=b(j)
               rleft(j+1)=rleft(j)
               b(j+1)=(right(j+1)+rleft(j+1))/2.0  

            endif
            if (AbS(eps).lt.1.0e-14) goto 912    
         enddo        
      endif
      
      if (rconst.lt.0) then
           
         right(1)=pi/r(ircpoint)
         rleft(1)=pi/(2.0*r(ircpoint))
         b(1)=(right(1)+rleft(1))/2.0
         
         do j=1,99
            
            eps=1/TAN(b(j)*r(ircpoint))-rconst/b(j) 
            if  (eps.gt.0) then
               right(j+1)=right(j)
              rleft(j+1)=b(j)
               b(j+1)=(right(j+1)+rleft(j+1))/2.0
            else 
               right(j+1)=b(j)
               rleft(j+1)=rleft(j)
               b(j+1)=(right(j+1)+rleft(j+1))/2.0  
            endif
            if (ABS(eps).lt.1.0e-14) goto 912    
         enddo        
      endif
      
      if (rconst.eq.0) then
         j=1
         b(j)=Pi/(2*r(ircpoint))
       endif
                
 912   continue
      sb=b(j)
      sa=r(ircpoint)*rcore(ircpoint)/SIN(sb*r(ircpoint))
      
      testr=sa*SIN(sb*r(ircpoint))/r(ircpoint)
      testl=rcore(ircpoint)
      dif=test-test1
      
      do  i=1,ircpoint 
         rpcore(i)=(sa*SIN(sb*r(i))/r(i))
      enddo
      do i=ircpoint+1,np
         rpcore(i)=rcore(i)
      enddo
      
      return
      end

      subroutine fhipcc(rcore,rpcore,iexit)
      implicit double precision(a-h,o-z)
      
#include "fortdim.h"
      
      parameter (pi=3.141592653589793239)
      parameter (n=100)

c -------------------------------------------------------------------------
c     External (shared between C and Fortran) common blocks
c -------------------------------------------------------------------------
      common /grid/ h,r1,z,r(npdm),np
      common /rpcc/ rpcc,rpccz
c -------------------------------------------------------------------------
      dimension rd1(npdm),rd2(npdm)
      dimension rcore(npdm),rpcore(npdm)
      dimension b(4),a(4,4)
      dimension ipvt(4)

c     This pcc method is used in the fhi98pp psp code

c     get the rpcc grid point
      irc=nint(log(rpcc/r(1))/h) +1

c     get the 1st and 2nd deriv. of the AE core density
      do i=1,irc+5
         x2old = 0.0
         diff2o = 1.0e-8
         do m = 3,10
            x2 = val2(rcore,r,np,r(i),m)
            diff2 = abs(x2-x2old)
            if (diff2.lt.diff2o) rd1(i) = x2
            x2old = x2
            diff2o = diff2
         enddo
         
         x3old = 0.0
         diff3o = 1.0e-8
         do m = 3,10
            x3 = val3(rcore,r,np,r(i),m)
            diff3 = abs(x3-x3old)
            if (diff3.lt.diff3o) rd2(i) = x3
            x3old = x3
            diff3o = diff3
         enddo
      enddo

c     Start the fhi98pp method
      i=irc

c     rpcc2 is the ACTUAL partial core radius (r of the closest g.p.)
      rpcc2=r(irc)

c     set the 1st coeff equal to the value of rho at rpcc2
      c0=rcore(irc)

c     loop until rho
      do it=1,30
         c0=1.8*rpden(r(i),c0,c3,c4,c5,c6)
         
         c1=0.d0
         c2=0.d0

         c3=(rd2(irc+1)-rd2(irc-1))/(r(irc+1)-r(irc-1))/rpcc2

         a(1,1)=1.d0
         a(1,2)=4*rpcc2
         a(1,3)=10*rpcc2**2
         a(1,4)=20*rpcc2**3

         a(2,1)=3.d0
         a(2,2)=6*rpcc2
         a(2,3)=10*rpcc2**2
         a(2,4)=15*rpcc2**3

         a(3,1)=1.d0
         a(3,2)=rpcc2
         a(3,3)=rpcc2**2
         a(3,4)=rpcc2**3

         a(4,1)=3.d0
         a(4,2)=4*rpcc2
         a(4,3)=5*rpcc2**2
         a(4,4)=6*rpcc2**3

c     constraint vector
         b(1)=c3/6.d0
         b(2)=rd2(irc)/(2*rpcc2)
         b(3)=(rcore(irc)-c0)/rpcc2**3
         b(4)=rd1(irc)/rpcc2**2

c     solve A . C = B
         call dgesv(4,1,a,4,ipvt,b,4,info)
         if(info .ne. 0) then
            if(info .lt. 0) write(7,*) 'stop pcc: singular matrix'      
            if(info .gt. 0) write(7,*) 'stop pcc: bad input'      
            goto 411
         endif

         c3=b(1)
         c4=b(2)
         c5=b(3)
         c6=b(4)
         
         do i=irc,2,-1
            if(rpden(r(i),c0,c3,c4,c5,c6) .gt. c0) goto 911
         enddo
 911     continue
         if(i.eq.1) goto 50
      enddo
      write(7,*) 'stop: partial core correction error',it 
      goto 411
 50   continue

      write(7,*) 'pcc matching at i',irc,' r ',r(irc)
      write(7,*) 'pcc parameters:'
      write(7,*) 'c0 ',c0
      write(7,*) 'c3 ',c3
      write(7,*) 'c4 ',c4
      write(7,*) 'c5 ',c5
      write(7,*) 'c6 ',c6

      do  i=1,irc
         rpcore(i)=rpden(r(i),c0,c3,c4,c5,c6)
      enddo
      do i=irc+1,np
         rpcore(i)=rcore(i)
      enddo

      return
      close(7)
      
 411  continue
      iexit=1
      close(7)
      return
      end

      double precision function rpden(rv,c0,c3,c4,c5,c6)
      implicit double precision(a-h,o-z)
      
      rv3=rv*rv*rv
      rv4=rv*rv3
      rv5=rv*rv4
      rv6=rv*rv5

      rpden=c0 + c3*rv3 + c4*rv4 + c5*rv5 + c6*rv6

      return
      end

