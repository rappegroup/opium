      subroutine etotal(etot,ev,ek,ep)
      implicit double precision(a-h,o-z)

#include "param.h"

      common /reli/ norb,ncore,nval,no(norbp),lo(norbp)
      common /reld/ zcore,zval,zo(norbp),so(norbp)
      common /rgrid/ aa,bb,r(nrmax),rab(nrmax)
      common /nrgrid/ nr

      dimension etot(10),ev(norbp),ek(norbp),ep(norbp),
     $     ev1(norbp),ek1(norbp),ep1(norbp)
c     
c     etot(i)    i=1,10 contains various contributions to the total
c     .          energy.
c     .          (1)   sum of eigenvalues ev
c     .          (2)   sum of orbital kinetic energies ek
c     .          (3)   el-ion interaction from sum of orbital
c     .                potential energies ep
c     .          (4)   electrostatic el-el interaction  (from velect)
c     .          (5)   vxc (exchange-correlation) correction to sum
c     .                of eigenvalues                   (from velect)
c     .          (6)   3 * vc - 4 * ec
c     .                correction term for virial theorem
c     .                when correlation is included     (from velect)
c     .          (7)   exchange and correlation energy  (from velect)
c     .          (8)   kinetic energy from eigenvalues  (1,3,4,5)
c     .          (9)   potential energy
c     .          (10)  total energy
c     
      dimension il(5)
      character il*1
c
c1    format(/,1x,a10,30(/,1x,10e13.4))
      pi = 4*atan(1.D0)
c     
c     sum up eigenvalues ev, kinetic energies ek, and
c     el-ion interaction ep
c     
      etot(1) = 0.D0
      etot(2) = 0.D0
      etot(3) = 0.D0
      do 100 i=1,norb
         etot(1) = etot(1) + zo(i)*ev(i)
         etot(2) = etot(2) + zo(i)*ek(i)
         etot(3) = etot(3) + zo(i)*ep(i)
  100 continue
c     
c     compute interaction shell - (nucleus-core)
c     

c     kinetic energy
c     
      etot(8) = etot(1) - etot(3) - 2*etot(4) - etot(5)
c     
c     potential energy
c     
      etot(9) = etot(3) + etot(4) + etot(7) 
c     
c     total energy
c     
      etot(10) = etot(2) + etot(3) + etot(4) + etot(7) 
      etst     = etot(1) - etot(4) - etot(5) + etot(7) 
c     
c     printout
c     
      il(1) = 's'
      il(2) = 'p'
      il(3) = 'd'
      il(4) = 'f'
      il(5) = 'g'
  110 format(' ',a2,' output data',/,1x,14('-'))
      if (norb .eq. 0) goto 150
      write(nout,120) 
  120 format('    output data for orbitals',/,1x,27('-'),//,
     +      ' nl    s      occ',9x,'eigenvalue',4x,'kinetic energy',
     +      6x,'pot energy',/)
      do 140 i=1,norb
         write(nout,130) no(i),il(lo(i)+1),so(i),zo(i),ev(i),ek(i),ep(i)
  130    format(1x,i1,a1,f6.1,f10.4,3f17.8)
  140 continue

      write(nout,*)
      write(nout,*) 'Averaged Values'
      write(nout,*)
      write(nout,141)
  141 format(' nl',6x,'occ',8x,'eigenvalue',4x,'kinetic energy',
     +	     6x,'pot energy')
      write(nout,*)
      do i=1,norb
         ev1(i)=0.0
         ek1(i)=0.0
         ep1(i)=0.0
       	 if (il(lo(i)+1).ne.'s') then
            if (no(i).eq.no(i+1).and.il(lo(i)+1).eq.il(lo(i+1)+1)) then
         	if (il(lo(i)+1).eq.'p') then
		      ev1(i) = ev(i)/3 + ev(i+1)*2/3
		      ek1(i) = ek(i)/3 + ek(i+1)*2/3
		      ep1(i) = ep(i)/3 + ep(i+1)*2/3
		elseif (il(lo(i)+1).eq.'d') then
		      ev1(i) = ev(i)*0.4 + ev(i+1)*0.6
		      ek1(i) = ek(i)*0.4 + ek(i+1)*0.6
		      ep1(i) = ep(i)*0.4 + ep(i+1)*0.6
		else
		      ev1(i) = ev(i)*3/7 + ev(i+1)*4/7
		      ek1(i) = ek(i)*3/7 + ek(i+1)*4/7
		      ep1(i) = ep(i)*3/7 + ep(i+1)*4/7
		endif
	    endif
	 else
	    ev1(i) = ev(i)
	    ek1(i) = ek(i)
	    ep1(i) = ep(i)
	 endif
	 if (abs(ev1(i)).gt.1e-18) then
            if (il(lo(i)+1).eq.'s') then
		write(nout,142) no(i),il(lo(i)+1),zo(i),ev1(i),
     +         	    ek1(i),ep1(i)
            else
		write(nout,142) no(i),il(lo(i)+1),zo(i)+zo(i+1),ev1(i),
     +         	    ek1(i),ep1(i)
	    endif
	 endif
  142    format(1x,i1,a1,f10.4,3f17.8)
      enddo

  150 write(nout,160) (etot(i),i=1,10),etst
  160 format(//,' total energies',/,1x,14('-'),/,
     +     /,' sum of eigenvalues        :',f18.8,
     +     /,' kinetic energy from ek    =',f18.8,
     +     /,' el-ion interaction energy =',f18.8,
     +     /,' el-el  interaction energy =',f18.8,
     +     /,' vxc    correction         :',f18.8,
     +     /,' virial correction         :',f18.8,
     +     /,' exchange + corr energy    =',f18.8,
     +     /,' kinetic energy from ev    :',f18.8,
     +     /,' potential energy          :',f18.8,/,1x,45('-'),
     +     /,' total energy              =',2f18.8)
c      if (itype .eq. 'fc' .or. itype .eq. 'pt' .or. itype .eq. 'pm'
c     +     .or. itype .eq. 'ps' .or. zsh .ne. 0.D0) return
c      if (ispp .eq. 'r') return

      return     

c     virial theorem
c     
      vsum = 2*etot(2) + etot(9) + etot(6)
      write(nout,170) 2*etot(2),etot(9),etot(6),vsum
  170 format(//,' virial theorem',/,1x,14('-'),/,
     +     /,' kinetic energy  *  2      =',f18.8,
     +     /,' potential energy          =',f18.8,
     +     /,' virial correction         =',f18.8,/,1x,45('-'),
     +     /,' virial sum                =',f18.8)
      return
      end

