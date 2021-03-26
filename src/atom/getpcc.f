      subroutine getpcc
      implicit double precision(a-h,o-z)
      
#include "PARAMHFS"
      
      parameter (n=100)

c -------------------------------------------------------------------------
c     External (shared between C and Fortran) common blocks
c -------------------------------------------------------------------------
      common /grid/ h,r1,z,r(npdm),np
      common /rpcc/ rpcc,rpccz
      common /rscore/ rscore(npdm),rdd(npdm),rddd(npdm),rscoretot(npdm)
      common /filenames/ file_log
c -------------------------------------------------------------------------

      dimension rdensity(npdm)

      character*80 file_log

      open(unit=7,file=file_log,form='formatted',access='append')

      pi=acos(-1.0)                                                   
      pi4=4.0*pi                                                      

      rsc2=0.0
      call radin(r,rscore,0,np,h,rsc2)

      irc=nint(log(rpcc/r(1))/h) +1

      do i=1,np
         rscoretot(i)=rscore(i)
         rdensity(i)=rscore(i)/(r(i)*r(i))
      enddo

      call LFC(rdensity,rscore)

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


c      do i=1,np
c         if (r(i).lt.3) then
c            write(33,*) r(i),rscore(i)/r(i)**2
c            write(34,*) r(i),rdd(i)
c            write(35,*) r(i),rddd(i)
c            write(36,*) r(i),rscoretot(i)            
c         endif
c         if (r(i).lt.20) imax=i
c      enddo


c      qmax=100
c      ip=5000
c      dql=qmax/ip
c      do j=1,ip                                                    
c         ql=dql*dble(j-1)                                             
c         do k=1,imax                                                
c            qr=r(k)*ql                                             
c            rtemp(k) = rscore(k) * besfn(qr,l)
c         enddo
c         pow=1
c         call radin(r,rtemp,0,imax,h,pow)
c         pl(j)=pow * pi4    
c         write(70,*) ql,pl(j)
c      enddo

      close(7)
            
 999  return
      end

      subroutine LFC(rcore,rpcore)
      implicit double precision(a-h,o-z)
      
#include "PARAMHFS"
      
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
            if (AbS(eps).lt.1.0e-14) goto 911    
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
            if (ABS(eps).lt.1.0e-14) goto 911    
         enddo        
      endif
      
      if (rconst.eq.0) then
         j=1
         b(j)=Pi/(2*r(ircpoint))
       endif
                
 911  continue
      
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
