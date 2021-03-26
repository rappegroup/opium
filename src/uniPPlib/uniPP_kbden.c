/*
 */

/****************************************************************************
* Part of uniPP Library providing access to universal pseudopotential files.*
*                                                                           *
* - uniPP_kbden() calculates the KB denominators of the non-local part of   *
*   the pseudopotential. Returned units: [ 1 / Hartree ]                    *
*                                                                           *
*                                                                           *
* last touched: 13.12.00 gjt                                                *
****************************************************************************/


/****************************************************************************
*  'uniPPlib' - Universal PseudoPotential Library                           *
*  Copyright (C) 2000  University of California, Santa Barbara              *
*                                                                           *
*  This program is free software; you can redistribute it and/or modify     *
*  it under the terms of the GNU General Public License as published by     *
*  the Free Software Foundation; either version 2 of the License, or        *
*  (at your option) any later version.                                      *
*                                                                           *
*  This program is distributed in the hope that it will be useful,          *
*  but WITHOUT ANY WARRANTY; without even the implied warranty of           *
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            *
*  GNU General Public License for more details.                             *
*                                                                           *
*  You should have received a copy of the GNU General Public License        *
*  along with this program; if not, write to the Free Software              *
*  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA  *
****************************************************************************/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "xmcomplex.h"
#include "xmintegral.h"

#include "uniPPlib.h"

#ifdef HAVE_SIMPSON

void uniPP_kbden(uniPP *unipp, double ***kbden, double nproj){

  int i, l, l_loc, j, j_max, proj;
  double pi = acos(-1.);
  double r, u, dv;
  double lambda, factor;
  double *f = (double *)malloc(unipp->m_mesh * sizeof(double));
    
  
  l_loc = unipp->l_loc;
    
  for (l=0; l<unipp->l_max; l++){
    if (l!=l_loc){
        
  
      if (unipp->rel && l){
        j_max  = 2;
        factor = 4.*pi;
      }else{
        j_max  = 1;
        factor = 4.*pi*(2.*(double)l+1.);
      }
 
      for (j=0; j<j_max; j++){
        for (proj=0; proj<nproj; proj++){
          if (proj==0){
            for (i=0; i<unipp->m_mesh; i++){
              r  = unipp->r_m[i];
              u  = unipp->u_ps[l][j][i];
 
 
              dv = unipp->v_ps[l][j][i] - unipp->v_loc[i];
 
 
              f[i] = r*u*u*dv;
 
            }
 
 
            kbden[l][j][0] = 1.
                 / (simpson(unipp->m_mesh, f, log(unipp->a_mesh)) + .5*f[0]);
 
          }else if (proj==1){
            for (i=0; i<unipp->m_mesh; i++){
              r = unipp->r_m[i];
              f[i] *= r*r;
            }
            lambda = simpson(unipp->m_mesh, f, log(unipp->a_mesh)) + .5*f[0];
 
            for (i=0; i<unipp->m_mesh; i++){
              r = unipp->r_m[i];
              f[i] *= r*r;
            }
 
            kbden[l][j][1] = 1. /
                    ((simpson(unipp->m_mesh, f, log(unipp->a_mesh)) + .5*f[0])
                     - lambda*lambda*kbden[l][j][0]);
 
          }else{
            printf("uniPP_kbden: requested projector of order"
                   " higher than two -> exit(0)\n");
            exit(0);
          }
        }
        for (proj=0; proj<nproj; proj++){
          kbden[l][j][proj] *= factor;
          printf("l=%d, j=%d, proj=%d, kbden=%g\n", l, j, proj,
                                                    kbden[l][j][proj]);
        }
        
        
      }
    }
  }
        
  free(f);
}

#endif /* HAVE_SIMPSON */
