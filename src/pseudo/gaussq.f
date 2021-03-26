c
c Copyright (c) 1998-2004 The OPIUM Group
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
      subroutine gaussq
      
      implicit double precision (a-h,o-z)           
      
#include "PARAMOPT"

c -------------------------------------------------------------------------
c     Internal (Fortran only) common blocks
c -------------------------------------------------------------------------
      common /gauss/ xquad(nquad0),wquad(nquad0)
      common /quads/ qq(nquad0)
      common /cuts/ qc,rc
c -------------------------------------------------------------------------
      
      xquad (  1 ) =-0.99886640442007105018e00
      xquad (  2 ) =-0.99403196943209071258e00
      xquad (  3 ) =-0.98535408404800588230e00    
      xquad (  4 ) =-0.97286438510669207371e00    
      xquad (  5 ) =-0.95661095524280794299e00    
      xquad (  6 ) =-0.93665661894487793378e00    
      xquad (  7 ) =-0.91307855665579189308e00    
      xquad (  8 ) =-0.88596797952361304863e00    
      xquad (  9 ) =-0.85542976942994608461e00    
      xquad ( 10 ) =-0.82158207085933594835e00    
      xquad ( 11 ) =-0.78455583290039926390e00    
      xquad ( 12 ) =-0.74449430222606853826e00    
      xquad ( 13 ) =-0.70155246870682225108e00    
      xquad ( 14 ) =-0.65589646568543936078e00    
      xquad ( 15 ) =-0.60770292718495023918e00    
      xquad ( 16 ) =-0.55715830451465005431e00    
      xquad ( 17 ) =-0.50445814490746420165e00    
      xquad ( 18 ) =-0.44980633497403878914e00    
      xquad ( 19 ) =-0.39341431189756512739e00    
      xquad ( 20 ) =-0.33550024541943735683e00    
      xquad ( 21 ) =-0.27628819377953199032e00    
      xquad ( 22 ) =-0.21600723687604175684e00    
      xquad ( 23 ) =-0.15489058999814590207e00    
      xquad ( 24 ) =-0.93174701560086140854e-1    
      xquad ( 25 ) =-0.31098338327188876112e-1    
      xquad ( 26 ) = 0.31098338327188876112e-1    
      xquad ( 27 ) = 0.93174701560086140854e-1    
      xquad ( 28 ) = 0.15489058999814590207e00    
      xquad ( 29 ) = 0.21600723687604175684e00    
      xquad ( 30 ) = 0.27628819377953199032e00    
      xquad ( 31 ) = 0.33550024541943735683e00    
      xquad ( 32 ) = 0.39341431189756512739e00    
      xquad ( 33 ) = 0.44980633497403878914e00    
      xquad ( 34 ) = 0.50445814490746420165e00    
      xquad ( 35 ) = 0.55715830451465005431e00    
      xquad ( 36 ) = 0.60770292718495023918e00    
      xquad ( 37 ) = 0.65589646568543936078e00    
      xquad ( 38 ) = 0.70155246870682225108e00    
      xquad ( 39 ) = 0.74449430222606853826e00    
      xquad ( 40 ) = 0.78455583290039926390e00    
      xquad ( 41 ) = 0.82158207085933594835e00    
      xquad ( 42 ) = 0.85542976942994608461e00    
      xquad ( 43 ) = 0.88596797952361304863e00    
      xquad ( 44 ) = 0.91307855665579189308e00    
      xquad ( 45 ) = 0.93665661894487793378e00    
      xquad ( 46 ) = 0.95661095524280794299e00    
      xquad ( 47 ) = 0.97286438510669207371e00    
      xquad ( 48 ) = 0.98535408404800588230e00    
      xquad ( 49 ) = 0.99403196943209071258e00    
      xquad ( 50 ) = 0.99886640442007105018e00    
                                             
      wquad (  1 ) = 0.29086225531551409584e-2    
      wquad (  2 ) = 0.67597991957454015027e-2    
      wquad (  3 ) = 0.10590548383650969263e-1    
      wquad (  4 ) = 0.14380822761485574419e-1    
      wquad (  5 ) = 0.18115560713489390351e-1    
      wquad (  6 ) = 0.21780243170124792981e-1    
      wquad (  7 ) = 0.25360673570012390440e-1    
      wquad (  8 ) = 0.28842993580535198029e-1    
      wquad (  9 ) = 0.32213728223578016648e-1    
      wquad ( 10 ) = 0.35459835615146154160e-1    
      wquad ( 11 ) = 0.38568756612587675244e-1    
      wquad ( 12 ) = 0.41528463090147697422e-1    
      wquad ( 13 ) = 0.44327504338803275492e-1    
      wquad ( 14 ) = 0.46955051303948432965e-1    
      wquad ( 15 ) = 0.49400938449466314921e-1    
      wquad ( 16 ) = 0.51655703069581138489e-1    
      wquad ( 17 ) = 0.53710621888996246523e-1    
      wquad ( 18 ) = 0.55557744806212517623e-1    
      wquad ( 19 ) = 0.57189925647728383723e-1    
      wquad ( 20 ) = 0.58600849813222445835e-1    
      wquad ( 21 ) = 0.59785058704265457509e-1    
      wquad ( 22 ) = 0.60737970841770216031e-1    
      wquad ( 23 ) = 0.61455899590316663756e-1    
      wquad ( 24 ) = 0.61936067420683243384e-1    
      wquad ( 25 ) = 0.62176616655347262321e-1    
      wquad ( 26 ) = 0.62176616655347262321e-1    
      wquad ( 27 ) = 0.61936067420683243384e-1    
      wquad ( 28 ) = 0.61455899590316663756e-1    
      wquad ( 29 ) = 0.60737970841770216031e-1    
      wquad ( 30 ) = 0.59785058704265457509e-1    
      wquad ( 31 ) = 0.58600849813222445835e-1    
      wquad ( 32 ) = 0.57189925647728383723e-1    
      wquad ( 33 ) = 0.55557744806212517623e-1    
      wquad ( 34 ) = 0.53710621888996246523e-1    
      wquad ( 35 ) = 0.51655703069581138489e-1    
      wquad ( 36 ) = 0.49400938449466314921e-1    
      wquad ( 37 ) = 0.46955051303948432965e-1    
      wquad ( 38 ) = 0.44327504338803275492e-1    
      wquad ( 39 ) = 0.41528463090147697422e-1    
      wquad ( 40 ) = 0.38568756612587675244e-1    
      wquad ( 41 ) = 0.35459835615146154160e-1    
      wquad ( 42 ) = 0.32213728223578016648e-1    
      wquad ( 43 ) = 0.28842993580535198029e-1    
      wquad ( 44 ) = 0.25360673570012390440e-1    
      wquad ( 45 ) = 0.21780243170124792981e-1    
      wquad ( 46 ) = 0.18115560713489390351e-1    
      wquad ( 47 ) = 0.14380822761485574419e-1    
      wquad ( 48 ) = 0.10590548383650969263e-1    
      wquad ( 49 ) = 0.67597991957454015027e-2    
      wquad ( 50 ) = 0.29086225531551409584e-2    
                                             
      do i = 1,nquad0
        qq(i) = (xquad(i) + 1.0e0) * qc/2.0
      enddo
 
      return                                  
      end                                     