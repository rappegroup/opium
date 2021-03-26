/*
 * Copyright (c) 1998-2004 The OPIUM Group
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
 * $Id: do_fhi.c,v 1.11 2004/10/02 18:34:49 ewalter Exp $
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "parameter.h"
#include "uniPPlib.h"
#include "fortparam.h"        /* fortran code parameters */
#include "do_fhi.h"
#include "nlm.h"

int do_fhi(param_t *param, FILE *fp_param, char *logfile){

  int i, l, k;
  uniPP unipp;
  char filename[180];
  FILE *fp;
  FILE *fp_log;
  time_t t;
  int icount,lixc;
  char datestring[7];
  double rpccmax;     /* where the pcc goes to virtually zero */
  int ill[4];
  int ncore,c;
  int kk=0;
  
  static double rscore[NPDM],rdd[NPDM],rddd[NPDM];
  static double rvcore[NVALE0+1][NPDM],rvloc[NPDM];
  static double rnl[N0][NPDM];

  fp_log = fopen(logfile, "a");
  fprintf(fp_log,"<<<do_fhi>>>\n");

  ncore=param->norb - param->nval;
  
  /* read in the local potential from a binary file created by do_nl() */
  /*  sprintf(filename, "%s.loc", param->name);
  fp = fopen(filename, "rb");
  fread(nlcore, sizeof(double), param->ngrid, fp);
  fclose(fp); */
  
  /* read in the non-local potential from a binary file created by do_ps() */

  sprintf(filename, "%s.loc", param->name);
  fp = fopen(filename, "rb");
  fread(rvloc, sizeof(double), param->ngrid, fp);
  fclose(fp);

  sprintf(filename, "%s.pot_ps", param->name);
  fp = fopen(filename, "rb");
  for (i=0; i<param->nval; i++)
    fread(rvcore[i], sizeof(double), param->ngrid, fp);
  fclose(fp);

  sprintf(filename, "%s.psi_nl", param->name);
  fp = fopen(filename, "rb");
  for (i=0; i<param->nval; i++)
    fread(rnl[i], sizeof(double), param->ngrid, fp);
  fclose(fp);

  /* if NLCC then read in the PCC from a binary file created by do_ae() */
  if (param->rpcc > 0.){
    sprintf(filename, "%s.rho_core", param->name);
    fp = fopen(filename, "rb");
    fread(rscore, sizeof(double), param->ngrid, fp);
    fread(rdd, sizeof(double), param->ngrid, fp);
    fread(rddd, sizeof(double), param->ngrid, fp);
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
    unipp.l_loc = nlm_label(param->nlm[param->localind+ncore]).l;
  unipp.rel = 0;        /* this format does not supports fully relativistic */
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
  for (l=0; l<unipp.l_max+1; l++){
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
  if (unipp.nlcc) {
    unipp.n_pc = (double *)malloc(unipp.m_mesh*sizeof(double));
    unipp.n_pc1 = (double *)malloc(unipp.m_mesh*sizeof(double));
    unipp.n_pc2 = (double *)malloc(unipp.m_mesh*sizeof(double));
  }

  for (i=0;i<4;i++)
    ill[i]=0;
  ncore=param->norb - param->nval;

  /* fill the radial arrays with information out of the previous pseudo call */
  unipp.r_m[0] = param->a * pow(param->z, -1./3.);
  for (i=0; i<param->ngrid; i++)
    unipp.r_m[i] = unipp.r_m[0] * exp(param->b * i);


  /* EJW: now we set the unipp variables equal to the 'correct' v_ps and u_ps 
     correct means:  the first instance of each l (for single projector semicore) 
     and, of these, arrange them in increasing l-order */

  icount=0; 
  for (kk=0; kk<param->nll;kk++)
    for (k=0; k<param->nval; k++){
      if ((ill[nlm_label(param->nlm[k+ncore]).l]==0) && (nlm_label(param->nlm[k+ncore]).l == kk)) {
	ill[nlm_label(param->nlm[k+ncore]).l]++;
	for (i=0; i<param->ngrid; i++){
	  if (unipp.rel && k){
	    /* set different j-components */
	  }else{
	    unipp.v_ps[icount][0][i] = rvcore[k][i]/(2.*unipp.r_m[i]);
	    unipp.u_ps[icount][0][i] = rnl[k][i];
	    
	  }
	}
	icount++;
      }
    }

  /* New section to add DNL to abinit */
  if (param->nboxes > 0) {
    /*    fprintf(fp_log," fhi format does not support the use of augmentation operators\n");*/
    fprintf(fp_log," Making l+1 the local potential %d\n",icount);
    for (i=0; i<param->ngrid; i++){
      unipp.v_ps[icount][0][i] = rvloc[i]/(2.*unipp.r_m[i]);
      unipp.u_ps[icount][0][i] = 0.0;
    }
    unipp.l_max++;
    unipp.l_loc = param->nll;
  }
  /* End new section */  

  rpccmax=0.0;
  if (unipp.nlcc) {
    for (i=0;i<param->ngrid;i++) {
      unipp.n_pc[i] = rscore[i] / (unipp.r_m[i] * unipp.r_m[i]);
      unipp.n_pc1[i] = rdd[i];
      unipp.n_pc2[i] = rddd[i];
      if ((i > 1) && (rpccmax == 0.0)) {
	if ((unipp.n_pc[i]<1e-15)&&(unipp.n_pc[i]<unipp.n_pc[i-1])) {
	  rpccmax = unipp.r_m[i];
	}
      }
    }
  }

  if (param->ixc == 0) lixc=2;
  if (param->ixc == 1) lixc=7;
  if (param->ixc == 2) lixc=11;

  /* open the file and call the method to write cpi format */  
  sprintf(filename, "%s.cpi", param->name);
  fp = fopen(filename, "w");
  uniPP_writefhi(&unipp, fp);
  fclose(fp);
  
  /* open the file and call write the fhi format for use as abinit format 6 */
  sprintf(filename, "%s.fhi", param->name);
  fp = fopen(filename, "w");
  fprintf(fp, "OPIUM generated %s potential\n", param->symbol);
  time(&t);
  strftime(datestring, 7, "%y%m%d", localtime(&t));
  fprintf(fp, "%10.5f%10.5f  %s\t\t\tzatom,zion,pspdat\n",
    param->z, unipp.z_ion, datestring);
  fprintf(fp, "%5d%5d%5d%5d%10d%10.5f\tpspcod,pspxc,lmax,lloc,mmax,r2well\n",
    6, lixc, unipp.l_max-1, unipp.l_loc,
    unipp.m_mesh, 0.);
  fprintf(fp, "%10.5f%10.5f%10.5f\t\t\trchrg,fchrg,qchrg\n",
    (unipp.nlcc)?rpccmax:0., (unipp.nlcc)?1.:0., 0.);
  fprintf(fp, "5 --- reserved for future features\n");
  fprintf(fp, "6 --- reserved for future features\n");
  fprintf(fp, "7 --- Here follows the cpi file in the fhi98pp format -\n");
  
  uniPP_writefhi(&unipp, fp);

  /* append the parameter file to the end of the file */
  fprintf(fp, "\n");
  fprintf(fp, "############################################################\n");
  fprintf(fp, "#    Opium Parameter File                                  #\n");
  fprintf(fp, "############################################################\n");
  fprintf(fp, "\n");
  rewind(fp_param);
  while((c=fgetc(fp_param)) != EOF) fputc(c, fp);
  fprintf(fp, "############################################################\n");
  fclose(fp);
  
  fp_log = fopen(logfile, "a");  
  fprintf(fp_log,"   ================================================\n");
  fclose(fp_log);
  
  return 0;
}
