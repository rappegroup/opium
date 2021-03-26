/* 
 * $Id: do_fhi.c,v 1.6 2004/06/16 20:46:17 mbarnes Exp $
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

int do_fhi(param_t *param, char *logfile){

  int i, l, k;
  uniPP unipp;
  char filename[180];
  FILE *fp;
  FILE *fp_log;
  time_t t;
  int icount;
  char datestring[7];
  double rpccmax;     /* where the pcc goes to virtually zero */
  int ill[4];
  int ncore;
  
  static double rscore[NPDM],rdd[NPDM],rddd[NPDM];
  static double rvcore[NVALE0+1][NPDM];
  static double rnl[N0][NPDM];

  fp_log = fopen(logfile, "a");
  fprintf(fp_log,"<<<do_fhi>>>\n");
  
  if (param->nboxes > 0) {
    fprintf(fp_log," fhi format does not support the use of augmentation operators\n");
    return 1;
  }

  /* read in the local potential from a binary file created by do_nl() */
  /*  sprintf(filename, "%s.loc", param->name);
  fp = fopen(filename, "rb");
  fread(nlcore, sizeof(double), param->ngrid, fp);
  fclose(fp); */
  
  /* read in the non-local potential from a binary file created by do_ps() */

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
    unipp.l_loc = param->local;
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

  /* EJW this is not needed and breaks things
  if (unipp.l_loc < 0 || unipp.l_loc >= unipp.l_max)
    unipp.v_loc[i] = nlcore[i]/(2.*unipp.r_m[i]);
  else
    unipp.v_loc[i] = rvcore[unipp.l_loc][i]/(2.*unipp.r_m[i]);
  */

  icount=0; 
  for (k=0; k<param->nval; k++){
    if (ill[nlm_label(param->nlm[k+ncore]).l]==0) {
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
    
    /* compute the first and second deriv. of the pcc */
    /* this should probably be its own fuction -- EJW */
    /*    
    for (i=1;i<param->ngrid-1;i++) {
      unipp.n_pc1[i] = (unipp.n_pc[i+1] - unipp.n_pc[i-1]);
      unipp.n_pc1[i] /= param->b * unipp.r_m[i];
    }
    unipp.n_pc1[0] = 2.0*unipp.n_pc1[1] - unipp.n_pc1[2];
    unipp.n_pc1[param->ngrid-1] = 2.0*unipp.n_pc1[param->ngrid-2] - unipp.n_pc1[param->ngrid-3];
    
    for (i=1;i<param->ngrid-1;i++) {
     unipp.n_pc2[i] = (unipp.n_pc1[i+1] - unipp.n_pc1[i-1]);
     unipp.n_pc2[i] /= param->b * unipp.r_m[i];
    }
    unipp.n_pc2[0] = 2.0*unipp.n_pc2[1] - unipp.n_pc2[2];
    unipp.n_pc2[param->ngrid-1] = 2.0*unipp.n_pc2[param->ngrid-2] - unipp.n_pc2[param->ngrid-3];
    */
  }

  

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
    6, (!strcmp(param->xcparam, "lda")?2:11), unipp.l_max-1, unipp.l_loc,
    unipp.m_mesh, 0.);
  fprintf(fp, "%10.5f%10.5f%10.5f\t\t\trchrg,fchrg,qchrg\n",
    (unipp.nlcc)?rpccmax:0., (unipp.nlcc)?1.:0., 0.);
  fprintf(fp, "5 --- reserved for future features\n");
  fprintf(fp, "6 --- reserved for future features\n");
  fprintf(fp, "7 --- Here follows the cpi file in the fhi98pp format -\n");
  
  uniPP_writefhi(&unipp, fp);
  fclose(fp);
  
  
  fprintf(fp_log,"   ================================================\n");
  fclose(fp_log);
  
  return 0;
}
