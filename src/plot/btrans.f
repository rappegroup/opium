      subroutine btrans(zeff,name)
      implicit double precision (a-h,o-z)
      
#include "PARAMHFS"

      common /atomic/ xion,rnl(npdm,n0),nlm(n0),wnl(n0),en(n0),norb
      common /grid/ h,r1,z,r(npdm),np
      common /box/ iboxstart(10),iboxend(10),boxheight(10),numbox
      common /ilogder/ ilogder
      common /ibound/ ibd(n0)
      common /local/ iloc,idesign
      common /logarith/ rphas,elogmax,elogmin,dlwf(npl0,n0)
      common /nlpot2/ inl,indrc(n0),IB,ID,IG
      common /npm/ ncores,nvales,iskip,maxip
      common /rpcc/ rpcc,rpccz
      common /rscore/ rscore(npdm),rdd(npdm),rddd(npdm),rscoretot(npdm)
      common /rcall/ rcall(n0)
      common /filenames/ file_log
      common /results/ rnorm(30), etot, lghost(30)
      common /nmax/ nmax(n0),maxim
      common /totpot/ rvcore(npdm,n0),rvps(npdm,n0),rvcoul(npdm)
      common /valden/ rsval(npdm)
      common /nlcore/ rvloc(npdm)
      common /nll/ nll

      character*1 name(2)

      dimension pl(1000),rtemp(npdm)

      if (name(2).eq."\0") then
         open(unit=70,file=name(1)//'.plt_wq',form='formatted')
         open(unit=60,file=name(1)//'.plt_vq',form='formatted')
      else
         open(unit=70,file=name(1)//name(2)//'.plt_wq',form='formatted')
         open(unit=60,file=name(1)//name(2)//'.plt_vq',form='formatted')
      endif

      rewind(70)
      rewind(60)
      pi=acos(-1.0)                                                   
      pi4=4.0*pi                                                      

      do m=1,nvales
         l = nlm(m)/10 - 10 * (nlm(m)/100)
         qmax=10
         qp=500
         dql=qmax/qp
         do j=1,qp                                                    
            ql=dql*dble(j-1)                                             
            do k=1,maxim                                                
               qr=r(k)*ql                                             
               rtemp(k) = rnl(k,m)*r(k) * besfn(qr,l)
            enddo
            pow=dble(l+1)**2
            call radin(r,rtemp,0,maxim,h,pow)
            pl(j)=pow * pi4    
            write(70,*) ql,pl(j)
         enddo
         write(70,*) '@'
      enddo

      do m=1,nll
         l = nlm(m)/10 - 10 * (nlm(m)/100)
         qmax=100
         qp=1000
         dql=qmax/qp
         do j=1,qp                                                    
            ql=dql*dble(j-1)                                             
            do k=1,maxim                                                
               qr=r(k)*ql                                             
               rtemp(k) = (rvcore(k,m)*r(k)+zeff*2.0*r(k)) 
     $              * besfn(qr,l)
            enddo
            pow=dble(l+1)**2
            call radin(r,rtemp,0,maxim,h,pow)
            pl(j)=pow * pi4    
            write(60,*) ql,pl(j)
         enddo
         write(60,*) '@'
      enddo

      do j=1,qp                                                    
         ql=dql*dble(j-1)                                             
         do k=1,maxim                                                
            qr=r(k)*ql                                             
            rtemp(k) = (rvloc(k)*r(k)+zeff*2.0*r(k)) 
     $           * besfn(qr,0)
         enddo
         pow=1.0
         call radin(r,rtemp,0,maxim,h,pow)
         pl(j)=pow * pi4    
         write(60,*) ql,pl(j)
      enddo
      write(60,*) '@'

      close(60)
      close(70)
      return

      end
