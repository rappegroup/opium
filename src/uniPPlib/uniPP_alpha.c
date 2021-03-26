/*
 */

/****************************************************************************
* Part of uniPP Library providing access to universal pseudopotential files.*
*                                                                           *
* - uniPP_alpha() calculates the alpha term of the non-Coulomb part of the  *
*   local part of the pseudopotential. Returned units: [Hartree*(a.u.)^3]   *
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

double uniPP_alpha(uniPP *unipp){
  
  int i;
  double pi = acos(-1.);
  double r, z_ion;
  double alpha;
  double *f = (double *)malloc(unipp->m_mesh * sizeof(double));
  
  z_ion = unipp->z_ion;
  
  for (i=0; i<unipp->m_mesh; i++){
    r = unipp->r_m[i]; 
    f[i] = r*r*(r*unipp->v_loc[i] + z_ion);
  }
  
  
  /* the last term in the returned value comes from the [0;r_0] interval */
  
  alpha = 4.*pi*(simpson(unipp->m_mesh, f, log(unipp->a_mesh)) + .5*f[0]);
  
  free(f);
  
  return alpha;
}

#endif /* HAVE_SIMPSON */
