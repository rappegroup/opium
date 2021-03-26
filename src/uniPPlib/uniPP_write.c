/*
 */

/****************************************************************************
* Part of uniPP Library providing access to universal pseudopotential files.*
*                                                                           *
*                                                                           *
* last touched: 30.04.02 gjt                                                *
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

#include "xmcomplex.h"

#include "uniPPlib.h"



int uniPP_write(uniPP *unipp, FILE *fp){

  int i, l;

  if (fp==NULL) return 1;
  
  /* line #1 */
  fprintf(fp, "%s\n", unipp->name);
  
  /* line #2 */
  fprintf(fp, "%5.1f\t", unipp->z_ion);
  fprintf(fp, "%d\t", unipp->l_max);
  fprintf(fp, "%d\t", unipp->l_loc);
  fprintf(fp, "%d\t", unipp->rel);
  fprintf(fp, "%d\t", unipp->nlcc);
  fprintf(fp, "%d\t", unipp->m_mesh);
  fprintf(fp, "%1.15e\n", unipp->a_mesh);
  
  /* lines #3... */  
  for (i=0; i<unipp->m_mesh; i++){
    fprintf(fp, "%1.15e\t", unipp->r_m[i]);
    fprintf(fp, "%1.15e\t", unipp->v_loc[i]);
    for (l=0; l<unipp->l_max; l++){
      if (unipp->rel && l){
        fprintf(fp, "%1.15e\t", unipp->v_ps[l][0][i]);
        fprintf(fp, "%1.15e\t", unipp->u_ps[l][0][i]);
        fprintf(fp, "%1.15e\t", unipp->v_ps[l][1][i]);
        fprintf(fp, "%1.15e\t", unipp->u_ps[l][1][i]);
      }else{
        fprintf(fp, "%1.15e\t", unipp->v_ps[l][0][i]);
        fprintf(fp, "%1.15e\t", unipp->u_ps[l][0][i]);
      }
    }
    if (unipp->nlcc)
      fprintf(fp, "%1.15e\n", unipp->n_pc[i]);
    else fprintf(fp, "\n");
  }
  
  return 0;
}
