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

/* standard libraries */
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "parameter.h"        /* defines structure: 'param_t' */
#include "nlm.h"              /* nlm_label call */
#include "cdim.h"        /* fortran code parameters */
#include "do_fc.h"            /* the module's own header */
#include "common_blocks.h"    /* fortran common blocks */
#include "energy.h"           /* this is for saving the energies */

/* report feature */

static char report[8000];
#define streq(a,b) (!strcasecmp(a,b))
void nrelorbnl(param_t *param, int); 
char * write_reportfc(param_t *param, char *rp,int,double temp_eigen[], double temp_norm[]);
void dftsolve_(double  *, int * ,double *, int *, int *, int *, int *, int *, int *);
void startae(param_t *param, int); 
void relorbae(param_t *param, int, char *); 
void nrelorbae(param_t *param, int, char *); 

int do_fc(param_t *param, char *logfile){

  int i;  
  FILE *fp_log;
  double zeff;
  int config=-1;
  int ipsp=1;
  int ifc=1;
  int iexit=0;
  int iprint=1;
  int irel=1;
  int irelxc=1;
  int ncore;
  FILE *fp;
  char filename[160];
  char *rp=report;
  double temp_eigen[10];
  double temp_norm[10];

  /* set the log file */
  sprintf(filenames_.file_log, "%s", logfile);

  fp_log = fopen(logfile, "a");

  fprintf(fp_log,"\n ======================================================================== \n");
  fprintf(fp_log," Begin FC calculation\n");
  fprintf(fp_log," ======================================================================== \n");
  
  ipsp = 0;   

  ncore = param->norb-param->nval;

  fprintf(fp_log," Performing non-relativistic FC calculation...  \n" );
  fclose(fp_log);

  if (!strcmp(param->reltype, "nrl")){
    /* set up the nlm's ; config=-1 is the reference state */
    nrelorbae(param,config,logfile);
  }else{
    relorbae(param,config,logfile);
  }

  /* starting guess for the AE potential */
  startae(param, aorb_.norb);

  sprintf(filename, "%s.psi_ae_core", param->name);
  fp = fopen(filename, "rb");
  for (i=0; i<ncore; i++) 
    fread(wfn_.rnl[i], sizeof(double), param->ngrid, fp);
  fclose(fp);

  sprintf(filename, "%s.eig_ae_core", param->name);
  fp = fopen(filename, "rb");

   for (i=0; i<ncore; i++) {
    fread(&adat_.en[i], sizeof(double),1, fp);
    fread(&adat_.wnl[i], sizeof(double),1, fp);
    fread(&aorb_.nlm[i], sizeof(int),1, fp);
    fread(&aorb_.nmax[i], sizeof(int),1, fp);
    fread(&aorb_.maxim, sizeof(int),1, fp);
    fread(&adat_.xion, sizeof(double),1, fp);
    fread(&aval_.ibd[i], sizeof(int),1, fp);
    fread(&etrial_.etrial[i], sizeof(double),1, fp);

    printf(" en: %d %lg \n", i,adat_.en[i]);
  }
  fclose(fp);
  
  zeff=adat_.xion;
  for (i=0; i<param->nval; i++) {
    if (aval_.ibd[i]==0) {
      fp_log = fopen(logfile, "a");
      fprintf(fp_log," !NOTE! State: |%3d> marked as unbound, using reference eigenvalue of %6.3f \n",
	      aorb_.nlm[i],etrial_.etrial[i]);
      fclose(fp_log);
    }
    zeff +=adat_.wnl[i];
  }


  if (!strcmp(param->reltype, "nrl")){
    irel=0;
    irelxc=0;
  }else{
    irel=1;
    if (streq(param->relxc,"rxc")){ 
      irelxc=1;
    }else{
      irelxc=0;
    }
  }

  dftsolve_(&param->z,&param->ixc,&param->exccut,&ipsp,&ifc,&iexit,&irel,&irelxc,&iprint);
  
  rp = write_reportfc(param,rp,config,temp_eigen,temp_norm);

  enfc[0] = aval_.etot;

  fp_log = fopen(logfile, "a");

  fprintf(fp_log,"\n ======================================================================== \n");
  fprintf(fp_log," End FC calculation\n");
  fprintf(fp_log," ======================================================================== \n");
  fclose(fp_log);
  return 0;
}

/* report section */

void do_fc_report(FILE *fp){
  fprintf(fp, "%s", report);
}

char * write_reportfc(param_t *param, char *rp, int config, double temp_eigen[], double temp_norm[]) {
  
  int i,j;
  double temp_eigen_tot, temp_norm_tot;
  int ncore;

  ncore=param->norb-param->nval;

  if (config<0) {
    rp+=sprintf(rp, 
		"\n FC:Orbital    Filling       Eigenvalues[Ry]          Norm       Ghost\n"
		"    ------------------------------------------------------------------\n");
    for (i=ncore; i<param->norb; i++)
      rp+=sprintf(rp, "\t%3d\t%6.3f\t%19.10f\t%14.10f\t%6s\n",
		  aorb_.nlm[i], adat_.wnl[i], adat_.en[i],
		  aval_.rnorm[i],(local_.nlghost[i]>0)?"yes":((local_.nlghost[i]<0)?"?":"no"));
    
  }else{
    rp+=sprintf(rp, 
		"\n FC:Orbital    Filling       Eigenvalues[Ry]          Norm       Ghost\n"
		"    ------------------------------------------------------------------\n");
    for (i=ncore; i<param->norb; i++)
      rp+=sprintf(rp, "\t%3d\t%6.3f\t%19.10f\t%14.10f\t%6s\n",
		  aorb_.nlm[i], adat_.wnl[i], adat_.en[i],
		  aval_.rnorm[i],(local_.nlghost[i]>0)?"yes":((local_.nlghost[i]<0)?"?":"no"));
    
  }
  rp+=sprintf(rp, "\n      E_tot = %19.10f Ry\n",
	      aval_.etot);

  if (config>=0) {
    rp+=sprintf(rp, 
		"\n AE-FC:Orbital Filling       Eigenvalues[mRy]         Norm[1e-3] \n"
		" AE-FC- --------------------------------------------------------------\n");
    temp_eigen_tot = temp_norm_tot = 0.;

    for (i=ncore; i<param->norb; i++) {
      j=i-ncore;
      rp+=sprintf(rp, " AE-FC- %3d\t%6.3f\t%19.10f\t%14.10f\t\n",
		  aorb_.nlm[i], adat_.wnl[i], 1000.0*(temp_eigen[j]-adat_.en[i]),
		  1000.0*(temp_norm[j]-aval_.rnorm[i]));
      
      temp_eigen_tot += fabs(temp_eigen[j]-adat_.en[i]);
      temp_norm_tot += fabs(temp_norm[j]-aval_.rnorm[i]);
    }
    rp+=sprintf(rp, " AE-FC-  total error =\t%19.10f\t%14.10f\n\n",
		1000*temp_eigen_tot, 1000*temp_norm_tot);

    rp+=sprintf(rp," =====================================================================\n");
  }    

  return rp;
}


