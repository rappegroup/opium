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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "parameter.h"
#include "uniPPlib.h"
#include "cdim.h"        /* fortran code parameters */
#include "common_blocks.h"
#include "do_cpmd.h"
#include "nlm.h"

void writeparam(param_t *param, FILE *fp, FILE *fp_param);
void nrelsproj(param_t *param, char *);
int do_cpmd(param_t *param, FILE *fp_param, char *logfile){

  int i, l, k, j;
  uniPP unipp;
  char filename[180];
  FILE *fp;
  FILE *fp_log;
  time_t t;
  int icount,lixc;
  char datestring[7];
  double rpccmax;     /* where the pcc goes to virtually zero */
  double zeff;
  int ill[4];
  int ncore,c;
  int kk=0;
  
  static double rscore[NPDM],rdd[NPDM],rddd[NPDM];
  static double rvcore[N0][NPDM],rvloc[NPDM];
  static double rnl[N0][NPDM];

  fp_log = fopen(logfile, "a");
  fprintf(fp_log,"<<<do_cpmd>>>\n");
  fclose(fp_log);

  ncore=param->norb - param->nval;

  nrelsproj(param,logfile);

  for (i=0; i<param->nll; i++) {
    sprintf(filename, "%s.pot.ps.l=%d", param->name,i);
    fp = fopen(filename, "rb");
    fread(rvcore[i], sizeof(double), param->ngrid, fp);
    fseek(fp,sizeof(double) ,param->ngrid);
    fclose(fp);
  }

  for (i=0; i<param->nll; i++) {
    sprintf(filename, "%s.psi.ps.l=%d", param->name,i);
    fp = fopen(filename, "rb");
    fread(rnl[i], sizeof(double), param->ngrid, fp);
    fclose(fp);
  }
  
  if (param->rpcc > 0.){
    fprintf(stderr, "CPMD format not working for pcc yet, sorry!\n");
    return 0;
  }

  if (param->rpcc > 0.){
    sprintf(filename, "%s.rho_pcore", param->name);
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
  unipp.rel = 0;        /* this format does not support fully relativistic */
  if (param->rpcc > 0.)
    unipp.nlcc = 1;
  else 
    unipp.nlcc = 0;
  
  unipp.m_mesh = param->ngrid;
  unipp.a_mesh = exp(param->b);

  /* allocate some memory for the radial arrays */
  unipp.r_m = (double *)malloc(unipp.m_mesh*sizeof(double));
  unipp.v_loc = (double *)malloc(unipp.m_mesh*sizeof(double));
  unipp.u_ps = (double ***)malloc((unipp.l_max+1)*sizeof(double **));
  unipp.v_ps = (double ***)malloc((unipp.l_max+1)*sizeof(double **));
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


  for (kk=0; kk<param->nll;kk++){
    /*    for (k=0; k<param->nval; k++){
	  if ((ill[nlm_label(param->nlm[k+ncore]).l]==0) && (nlm_label(param->nlm[k+ncore]).l == kk)) {
	  ill[nlm_label(param->nlm[k+ncore]).l]++;

	  if (unipp.rel && k){
	  }else{*/
    for (i=0; i<param->ngrid; i++){
      unipp.v_ps[kk][0][i] = rvcore[kk][i]/(2.*unipp.r_m[i]);
      unipp.u_ps[kk][0][i] = rnl[kk][i];
    }
  }
  /* New section to add DNL to cpmd */
  if (param->nboxes > 0) {
    sprintf(filename, "%s.loc", param->name);
    if (fp = fopen(filename, "rb")) {
      fread(rvloc, sizeof(double), param->ngrid, fp);
      fclose(fp);
    } else {
      fp_log = fopen(logfile, "a");
      fprintf(fp_log,"Looks like you never ran nl yet you have augmentation functions :( --EXIT!\n");
      printf("Looks like you never ran nl yet you have augmentation functions :( --EXIT!\n");
      fclose(fp_log);
      exit(1);
    }
    fp_log = fopen(logfile, "a");
    fprintf(fp_log," Making l+1 the local potential %d\n",kk);
    fclose(fp_log);

    for (i=0; i<param->ngrid; i++){
      unipp.v_ps[kk][0][i] = rvloc[i]/(2.*unipp.r_m[i]);
      unipp.u_ps[kk][0][i] = 0.0;
    }
    unipp.l_max++;

    /*    unipp.l_loc = param->nll;*/
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

  lixc=2;
  if (param->ixc == 0) lixc=2;
  if (param->ixc == 1) lixc=7;
  if (param->ixc == 2) lixc=11;

  sprintf(filename, "%s.cpmd", param->name);
  fp = fopen(filename, "w");
  zeff = param->z;
  for (i=0; i<param->norb-param->nval; i++)
    zeff -= param->wnl[i];

  fprintf(fp, "&ATOM\n");  
  fprintf(fp, " Z  =%5.0f\n",param->z);  
  fprintf(fp, " ZV =%5.0f\n",zeff);  
  fprintf(fp, "  XC = 1312        .666670\n");
  fprintf(fp, "&END\n");  
  fprintf(fp, "&INFO\n");  
  fprintf(fp, "OPIUM generated %s potential\n", param->symbol);
  fprintf(fp, "&END\n");  

  fprintf(fp, " &POTENTIAL\n");

  fprintf(fp, "%5d%25.15e\n", unipp.m_mesh, unipp.a_mesh);  
  for (i=0; i<unipp.m_mesh; i++) {
    fprintf(fp, "%25.15e", unipp.r_m[i]);
    for (l=0; l<unipp.l_max; l++){    
      fprintf(fp, "%25.15e", unipp.v_ps[l][0][i]);
    }
    fprintf(fp, "\n");
  }
  fprintf(fp, " &END\n");
  fprintf(fp, " &WAVEFUNCTION\n");
  fprintf(fp, "%5d%25.15e\n", unipp.m_mesh, unipp.a_mesh);  
  for (i=0; i<unipp.m_mesh; i++) {
    fprintf(fp, "%25.15e", unipp.r_m[i]);
      for (l=0; l<unipp.l_max; l++){    
	fprintf(fp, "%25.15e", unipp.u_ps[l][0][i]);
      }
      fprintf(fp, "\n");
  }
  fprintf(fp, " &END\n");

  writeparam(param, fp, fp_param);  

  fp_log = fopen(logfile, "a");  
  fprintf(fp_log,"   ================================================\n");
  fclose(fp_log);

  uniPP_free(&unipp);
  
  return 0;
}
