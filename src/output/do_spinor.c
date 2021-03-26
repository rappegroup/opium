/*  This is code for the "spinor" psp format. The suffix for this
format was originally .upf, but I am changing this to .spinor so 
the QE UPF output can be introduced.  This format also uses the UniPPlib 
directory. I am not sure if this output even works anymore?? */
/*
 * Copyright (c) 1998-2012 The OPIUM Group
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 */
/*
 */

/****************************************************************************
*                                                                           *
* generate *.spinor output                                                     *
*****************************************************************************
INPUT:  param_t structure, *.nlp, *.loc
OUTPUT: *.spinor
****************************************************************************/

/* standard libraries */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "parameter.h"        /* defines structure: 'param_t' */
#include "uniPPlib.h"         /* library for spinor pseudopotential format */
#include "cdim.h"        /* fortran code parameters */
#include "do_spinor.h"           /* the module's own header */


int do_spinor(param_t *param, char *logfile){

  int i, l;           /* loop counter */
  uniPP unipp;        /* a uniPP object */
  char filename[180]; /* filename */
  FILE *fp;           /* file pointer */
  FILE *fp_log;
  int nnval;
  
  static double rscore[NPDM];
  static double nlcore[NPDM];
  static double rvcore[NVALE0+1][NPDM];
  static double rnl[N0][NPDM];

  fp_log = fopen(logfile, "a");
  fprintf(fp_log,"<<<do_spinor>>>\n");
  
  /* set nnval to the correct number of orbitals, considering frl case */
  if (!strcmp(param->reltype, "frl"))
    nnval = 2 * param->nll -1;
  else
    nnval = param->nll;

  /* read in the local potential from a binary file created by do_nl() */
  sprintf(filename, "%s.loc", param->name);
  fp = fopen(filename, "rb");
  fread(nlcore, sizeof(double), param->ngrid, fp);
  fclose(fp);

  sprintf(filename, "%s.psi_nl", param->name);
  fp = fopen(filename, "rb");
  for (i=0; i<param->nval; i++)
    fread(rnl[i], sizeof(double), param->ngrid, fp);
  fclose(fp);

  sprintf(filename, "%s.pot_ps", param->name);
  fp = fopen(filename, "rb");
  for (i=0; i<param->nval; i++) {
    fread(rvcore[i], sizeof(double), param->ngrid, fp);
  }	
  fclose(fp);

  if (param->rpcc > 0.){
    sprintf(filename, "%s.rho_pcore", param->name);
    fp = fopen(filename, "rb");
    fread(rscore, sizeof(double), param->ngrid, fp);
    fclose(fp);
  }

  /* set the uniPPlib structure */
  strncpy(unipp.name, param->name, 80);
  unipp.z_ion = param->z;
  for (i=0; i<param->norb - param->nval; i++)
    unipp.z_ion -= param->wnl[i];
  unipp.l_max = param->nll;
  if (param->nboxes > 0)
    unipp.l_loc = param->nll;
  else
    unipp.l_loc = param->local;
  if (!strcmp(param->reltype, "frl"))
    unipp.rel = 1;
  else
    unipp.rel = 0;
  if (param->rpcc > 0.)
    unipp.nlcc = 1;
  else 
    unipp.nlcc = 0;
  
  unipp.m_mesh = param->ngrid;
  unipp.a_mesh = exp(param->b);
  
  /* allocate some memory for the radial arrays */
  unipp.r_m = (double *)malloc(unipp.m_mesh*sizeof(double));
  unipp.v_loc = (double *)malloc(unipp.m_mesh*sizeof(double));
  unipp.u_ps = (double ***)malloc(unipp.l_max*sizeof(double **));
  unipp.v_ps = (double ***)malloc(unipp.l_max*sizeof(double **));
  for (l=0; l<unipp.l_max; l++){
    if (unipp.rel && l){
      unipp.u_ps[l] = (double **)malloc(2*sizeof(double));
      unipp.v_ps[l] = (double **)malloc(2*sizeof(double));
      unipp.u_ps[l][0] = (double *)malloc(unipp.m_mesh*sizeof(double));
      unipp.v_ps[l][0] = (double *)malloc(unipp.m_mesh*sizeof(double));
      unipp.u_ps[l][1] = (double *)malloc(unipp.m_mesh*sizeof(double));
      unipp.v_ps[l][1] = (double *)malloc(unipp.m_mesh*sizeof(double));
    }else{
      unipp.u_ps[l] = (double **)malloc(sizeof(double));
      unipp.v_ps[l] = (double **)malloc(sizeof(double));      
      unipp.u_ps[l][0] = (double *)malloc(unipp.m_mesh*sizeof(double));
      unipp.v_ps[l][0] = (double *)malloc(unipp.m_mesh*sizeof(double));
    }
  }
  if (unipp.nlcc)
    unipp.n_pc = (double *)malloc(unipp.m_mesh*sizeof(double));

  /* fill the radial arrays with information out of the previous pseudo call */
  unipp.r_m[0] = param->a * pow(param->z, -1./3.);
  for (i=0; i<param->ngrid; i++){
    unipp.r_m[i] = unipp.r_m[0] * exp(param->b * i);
    if (unipp.l_loc < 0 || unipp.l_loc >= unipp.l_max)
      unipp.v_loc[i] = nlcore[i]/(2.*unipp.r_m[i]);
    else
      unipp.v_loc[i] = rvcore[unipp.l_loc][i]/(2.*unipp.r_m[i]);
    for (l=0; l<unipp.l_max; l++){
      if (unipp.rel && l){
        unipp.v_ps[l][0][i] = rvcore[2*l][i]/(2.*unipp.r_m[i]);
        unipp.u_ps[l][0][i] = rnl[2*l][i];
        unipp.v_ps[l][1][i] = rvcore[2*l-1][i]/(2.*unipp.r_m[i]);
        unipp.u_ps[l][1][i] = rnl[2*l-1][i];
      }else{
        unipp.v_ps[l][0][i] = rvcore[l][i]/(2.*unipp.r_m[i]);
        unipp.u_ps[l][0][i] = rnl[l][i];
      }
    }
    if (unipp.nlcc)
      unipp.n_pc[i] = rscore[i] / (unipp.r_m[i] * unipp.r_m[i]);
  }
  
  /* open the file and call the method to write spinor format */  
  sprintf(filename, "%s.spinor", param->name);
  fp = fopen(filename, "w");
  uniPP_write(&unipp, fp);
  fclose(fp);  
  
  fprintf(fp_log,"   ================================================\n");
  fclose(fp_log);
  
  return 0;
}
