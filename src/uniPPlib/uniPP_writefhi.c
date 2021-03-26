/*
 */

/****************************************************************************
* Part of uniPP Library providing access to universal pseudopotential files.*
*                                                                           *
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

#include "uniPPlib.h"


int uniPP_writefhi(uniPP *unipp, FILE *fp){

  int i, l;

  if (fp==NULL) return 1;
  
  /* line #1 */
  fprintf(fp, "%13.9f%5d\n", unipp->z_ion, unipp->l_max);
  
  /* Gaussian block not used */
  fprintf(fp, " 0.0000    0.0000    0.0000   0.0000\n");
  fprintf(fp, " 0.0000    .00e+00   .00e+00\n");
  fprintf(fp, " 0.0000    .00e+00   .00e+00\n");
  fprintf(fp, " 0.0000    .00e+00   .00e+00\n");
  fprintf(fp, " 0.0000    .00e+00   .00e+00\n");
  fprintf(fp, " 0.0000    .00e+00   .00e+00\n");
  fprintf(fp, " 0.0000    .00e+00   .00e+00\n");
  fprintf(fp, " 0.0000    .00e+00   .00e+00\n");
  fprintf(fp, " 0.0000    .00e+00   .00e+00\n");
  fprintf(fp, " 0.0000    .00e+00   .00e+00\n");
  
  /* grid, wavefunction, potential */
  for (l=0; l<unipp->l_max; l++){
    fprintf(fp, "%5d%25.15e\n", unipp->m_mesh, unipp->a_mesh);
    for (i=0; i<unipp->m_mesh; i++)
      fprintf(fp, "%5d%25.15e%25.15e%25.15e\n", i+1, unipp->r_m[i],
        unipp->u_ps[l][0][i], unipp->v_ps[l][0][i]);
  }
  
  if (unipp->nlcc)
    for (i=0; i<unipp->m_mesh; i++)
      fprintf(fp, "%1.15e\t%1.15e\t%1.15e\t%1.15e\n", 
	      /*	      unipp->r_m[i], unipp->n_pc[i],0.0,0.0);*/
	      unipp->r_m[i], unipp->n_pc[i],unipp->n_pc1[i],unipp->n_pc2[i]);

  /* Notice that there is an inconsistency between the fhi98pp documentation
     and the actual implementation in their code. It might be necessary to
     to replace multiply unipp->n_pc[i] by 4.*pi in the above output! */
        
  return 0;
}
