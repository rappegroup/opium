/*
 */

/****************************************************************************
* Part of uniPP Library providing access to universal pseudopotential files.*
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

#include "xmcomplex.h"

#include "uniPPlib.h"



int uniPP_read(uniPP *unipp, FILE *fp){

  int i, l;

  if (fp==NULL) return 1;
  
  /* line #1 */
  if (fscanf(fp, "%[^\n]", unipp->name)==EOF) {fclose(fp); return 1;}
  
  /* line #2 */  
  if (fscanf(fp, "%lg", &unipp->z_ion)==EOF) {fclose(fp); return 1;}
  if (fscanf(fp, "%d",  &unipp->l_max)==EOF) {fclose(fp); return 1;}
  if (fscanf(fp, "%d",  &unipp->l_loc)==EOF) {fclose(fp); return 1;}
  if (fscanf(fp, "%d",  &unipp->rel)==EOF) {fclose(fp); return 1;}
  if (fscanf(fp, "%d",  &unipp->nlcc)==EOF) {fclose(fp); return 1;}
  if (fscanf(fp, "%d",  &unipp->m_mesh)==EOF) {fclose(fp); return 1;}
  if (fscanf(fp, "%lg", &unipp->a_mesh)==EOF) {fclose(fp); return 1;}
  
  /* allocate some memory */
  unipp->r_m = (double *)malloc(unipp->m_mesh*sizeof(double));
  unipp->v_loc = (double *)malloc(unipp->m_mesh*sizeof(double));
  unipp->u_ps = (double ***)malloc(unipp->l_max*sizeof(double **));
  unipp->v_ps = (double ***)malloc(unipp->l_max*sizeof(double **));
  for (l=0; l<unipp->l_max; l++){
    if (unipp->rel && l){
      unipp->u_ps[l] = (double **)malloc(2*sizeof(double));
      unipp->v_ps[l] = (double **)malloc(2*sizeof(double));      
      unipp->u_ps[l][0] = (double *)malloc(unipp->m_mesh*sizeof(double));
      unipp->v_ps[l][0] = (double *)malloc(unipp->m_mesh*sizeof(double));
      unipp->u_ps[l][1] = (double *)malloc(unipp->m_mesh*sizeof(double));
      unipp->v_ps[l][1] = (double *)malloc(unipp->m_mesh*sizeof(double));
    }else{
      unipp->u_ps[l] = (double **)malloc(sizeof(double));
      unipp->v_ps[l] = (double **)malloc(sizeof(double));      
      unipp->u_ps[l][0] = (double *)malloc(unipp->m_mesh*sizeof(double));
      unipp->v_ps[l][0] = (double *)malloc(unipp->m_mesh*sizeof(double));
    }
  }
  if (unipp->nlcc)
    unipp->n_pc = (double *)malloc(unipp->m_mesh*sizeof(double));

  /* lines #3... */  
  for (i=0; i<unipp->m_mesh; i++){
    if (fscanf(fp, "%lg", &unipp->r_m[i])==EOF) {fclose(fp); return 1;}
    if (fscanf(fp, "%lg", &unipp->v_loc[i])==EOF) {fclose(fp); return 1;}
    for (l=0; l<unipp->l_max; l++){
      if (unipp->rel && l){
        if (fscanf(fp, "%lg", &unipp->v_ps[l][0][i])==EOF) 
          {fclose(fp); return 1;}
        if (fscanf(fp, "%lg", &unipp->u_ps[l][0][i])==EOF) 
          {fclose(fp); return 1;}
        if (fscanf(fp, "%lg", &unipp->v_ps[l][1][i])==EOF) 
          {fclose(fp); return 1;}
        if (fscanf(fp, "%lg", &unipp->u_ps[l][1][i])==EOF) 
          {fclose(fp); return 1;}
      }else{
        if (fscanf(fp, "%lg", &unipp->v_ps[l][0][i])==EOF) 
          {fclose(fp); return 1;}
        if (fscanf(fp, "%lg", &unipp->u_ps[l][0][i])==EOF) 
          {fclose(fp); return 1;}
      }
    }
    if (unipp->nlcc)
      if (fscanf(fp, "%lg", &unipp->n_pc[i])==EOF) {fclose(fp); return 1;}
  }
  
  return 0;
}
