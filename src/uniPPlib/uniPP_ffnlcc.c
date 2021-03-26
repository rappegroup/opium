/*
 */

/****************************************************************************
* Part of uniPP Library providing access to universal pseudopotential files.*
*                                                                           *
* - uniPP_ffnlcc() calculates the form factor of the partial core density   *
*   neccessary for non-linear core corrections.                             *
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

void uniPP_ffnlcc(uniPP *unipp, double *ffnlcc, double *g, int n){

  int i, j;
  double pi = acos(-1.);
  double r, z;
  double *f = (double *)malloc(unipp->m_mesh * sizeof(double));
  
  for (j=0; j<n; j++){
    if (g[j] == 0.){
      for (i=0; i<unipp->m_mesh; i++){
        r = unipp->r_m[i];
        f[i] = r*r*r*unipp->n_pc[i];
      }
    }else{
      for (i=0; i<unipp->m_mesh; i++){
        r = unipp->r_m[i]; 
        z = r*g[j];
        f[i] = (sin(z)/z)*r*r*r*unipp->n_pc[i];
      }
    }
    ffnlcc[j] = 4.*pi*
                (simpson(unipp->m_mesh, f, log(unipp->a_mesh)) + .5*f[0]);
  }
  
  free(f);
}

#endif /* HAVE_SIMPSON */
