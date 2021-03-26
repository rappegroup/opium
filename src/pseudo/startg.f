c ****************************************************************************
c Copyright (C) 2002 Andrew M Rappe group, University of Pennsylvania
c This file is distributed under the terms of the GNU General Public
c License as described in the file 'License' in the current directory.
c ****************************************************************************

      subroutine startg
      
c     *************************************************************************
c     This code generates starting guesses for xlam and for qc.
c     *************************************************************************

      implicit real (a-h,o-z)
      
#include "PARAMOPT"

      common /nnn/ nnn
      common /cuts/ qc,rc
      common /numfn/ numfn
      common /scat/ escat(numen0),xscat(numen0),rscat(numen0),numen
      common /scatw/ wscat(numen0),wnorm,wke,wghost
      common /nonloc/ nonloc,ifnl,numloc
      common /nlall/ nlall(nvale0)
      common /atomic/ xion,rnl(npdm,n0),nlm(n0),wnl(n0),en(n0),norb

      common /lparam/ qcl(10),nbl(10),ifnll(10),
     $                numenl(10),wnorml(10),wkel(10),
     $                escatl(numen0,10),rscatl(numen0,10),
     $                wscatl(numen0,10) 
      
c gjt: I introduced the /lparam/ block in order to get parameters from C

      qc = qcl(nnn)
      numfn = nbl(nnn)
      
c gjt: since I removed the nonloc loop in pseudo.f I need to set nonloc here

      if (numfn.gt.numfn0) then
        write (7,*) 'Too many basis functions requested.  Recompile.'
        write (7,*) 'Execution terminates.'
        stop
      endif

      write(7,9000) nlm(nnn) 
      write(7,9001) numfn
      write(7,9002) rc
      write(7,9003) qc

 9000 format(1x,'--------------------Pseudizing state: ',
     $       '|',i3,'>',3x,'--------------------')
 9001 format(1x,'# basis functions        : ',i4)
 9002 format(1x,'rc                       : ',f8.4)
 9003 format(1x,'qc                       : ',f8.4)

      return
      end
