/*
 * Copyright (c) 1998-2005 The OPIUM Group
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
 * $Id: do_nl.c,v 1.10 2004/10/02 18:34:48 ewalter Exp $
 */

/* standard libraries */
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "parameter.h"        /* defines structure: 'param_t' */
#include "nlm.h"              /* nlm_label call */
#include "cdim.h"        /* fortran code parameters */
#include "do_nl.h"            /* the module's own header */
#include "common_blocks.h"    /* fortran common blocks */
#include "energy.h"           /* this is for saving the energies */

/* report feature */

static char report[800];
void nrelorbnl(param_t *param, int); 
char * write_reportnl(param_t *param, char *rp,int,double temp_eigen[], double temp_norm[]);
void scpot_(double  *, int * ,int *, int *, int *);
void readPS(param_t *param);
void writeNL(param_t *param);

int do_nl(param_t *param, char *logfile){

  int i;  
  FILE *fp_log;
  double zeff;
  int config=-1;
  int ipsp=1;
  int ifc=0;
  int iexit=0;
  char *rp=report;
  double temp_eigen[10];
  double temp_norm[10];


  /* set the log file */
  sprintf(filenames_.file_log, "%s", logfile);

  fp_log = fopen(logfile, "a");

  fprintf(fp_log,"\n ======================================================================== \n");
  fprintf(fp_log," Begin NL calculation\n");
  fprintf(fp_log," ======================================================================== \n");
  fclose(fp_log);
  
  atomic_.norb = param->nval;
  npm_.ncores = 0; 
  ipsp = 1;   
  nlpot2_.inl = 1;  

  readPS(param);
  nrelorbnl(param,config);

  zeff=atomic_.xion;
  for (i=0; i<param->nll; i++) {
    if (ibound_.ibd[i]==0) {
      fp_log = fopen(logfile, "a");
      fprintf(fp_log," !NOTE! State: |%3d> marked as unbound, using reference eigenvalue of %6.3f \n",
	      atomic_.nlm[i],ensave_.ensave[i]);
      fclose(fp_log);
    }
    zeff +=atomic_.wnl[i];
  }

  scpot_(&zeff,&param->ixc,&ipsp,&ifc, &iexit); 
  if (iexit) {
    printf("Terminal error in: scpot <-- do_nl\n EXITING OPIUM \n");
    exit(1);
  }


  writeNL(param);

  rp = write_reportnl(param,rp,config,temp_eigen,temp_norm);

  ennl[0] = results_.etot;

  fp_log = fopen(logfile, "a");

  fprintf(fp_log,"\n ======================================================================== \n");
  fprintf(fp_log," End NL calculation\n");
  fprintf(fp_log," ======================================================================== \n");
  fclose(fp_log);
  return 0;
}

/* report section */

void do_nl_report(FILE *fp){
  fprintf(fp, "%s", report);
}

void nrelorbnl(param_t *param, int config) {

  int i,ii,j,k,l;
  int nrcore;

  nrcore=param->norb - param->nval;
  atomic_.norb = param->nval;
  atomic_.xion = param->z;

  for (i=0; i<nrcore; i++)
    atomic_.xion -= param->wnl[i];

  /* map the AE |nlm> labels onto the NL |nlm> labels */
  ii = 0;
  for (i=0; i<param->nval; i++){
    if (config<0) {
      atomic_.nlm[ii] = param->nlm[nrcore + i];
      atomic_.wnl[ii] = param->wnl[nrcore + i];
    }else{
      atomic_.nlm[ii] = param->nlm_conf[config][i];
      atomic_.wnl[ii] = param->wnl_conf[config][i];
    }      
    /*    ibound_.ibd[ii] = param->ibound[i+nrcore];*/
    atomic_.xion -= atomic_.wnl[ii];
    k = 0;
    l = nlm_label(atomic_.nlm[ii]).l;
    for (j=0; j<ii; j++)
      if (nlm_label(atomic_.nlm[j]).l == l) k++;
    atomic_.nlm[ii] = (l+k+1)*100 + l*10
      + nlm_label(atomic_.nlm[ii]).m;
    ++ii;
  }
}

char * write_reportnl(param_t *param, char *rp, int config, double temp_eigen[], double temp_norm[]) {
  
  int i;
  double temp_eigen_tot, temp_norm_tot;
  int igh=0;

  if (config<0) {
    rp+=sprintf(rp, 
		"\n NL:Orbital    Filling       Eigenvalues[Ry]          Norm       Ghost\n"
		"    ------------------------------------------------------------------\n");
    for (i=0; i<param->nval; i++){
      if (ibound_.ibd[i]==1) {
	rp+=sprintf(rp, "\t%3d\t%6.3f\t%19.10f\t%14.10f\t%6s\n",
		    atomic_.nlm[i], atomic_.wnl[i], atomic_.en[i],
		    results_.rnorm[i],(results_.lghost[i]>0)?"yes":((results_.lghost[i]<0)?"?":"no"));
      }else{
	rp+=sprintf(rp, "\t%3d\t%6.3f\t%19.10f\t\t\t%6s\n",
		    atomic_.nlm[i], atomic_.wnl[i], atomic_.en[i],
		    (results_.lghost[i]>0)?"yes":((results_.lghost[i]<0)?"?":"no"));
	
      }
      igh+=results_.lghost[i];
      
    }
    if (igh > 0) {
      rp+=sprintf(rp, "\n");
      rp+=sprintf(rp, "  !!ERROR!! Ghosts are present in pseudopotential \n");
      rp+=sprintf(rp, "  !!ERROR!! See log file for more information \n");
    }else{
      rp+=sprintf(rp, "\n\t ========== No ghosts in potential!!========== \n");
    }

  }else{
    rp+=sprintf(rp, 
		"\n NL:Orbital    Filling       Eigenvalues[Ry]          Norm       Ghost\n"
		"    ------------------------------------------------------------------\n");
    for (i=0; i<param->nval; i++)
      if (ibound_.ibd[i]==1) {
	rp+=sprintf(rp, "\t%3d\t%6.3f\t%19.10f\t%14.10f\t%6s\n",
		    atomic_.nlm[i], atomic_.wnl[i], atomic_.en[i],
		    results_.rnorm[i],(results_.lghost[i]>0)?"yes":((results_.lghost[i]<0)?"?":"no"));
      }else{
	rp+=sprintf(rp, "\t%3d\t%6.3f\t%19.10f\t\t\t%6s\n",
		    atomic_.nlm[i], atomic_.wnl[i], atomic_.en[i],
		    (results_.lghost[i]>0)?"yes":((results_.lghost[i]<0)?"?":"no"));
	
      }
  }
  rp+=sprintf(rp, "\n      E_tot = %19.10f Ry\n",
	      results_.etot);

  if (config>=0) {
    rp+=sprintf(rp, 
		"\n AE-NL:Orbital Filling       Eigenvalues[mRy]         Norm[1e-3] \n"
		" AE-NL- --------------------------------------------------------------\n");
    temp_eigen_tot = temp_norm_tot = 0.;

    for (i=0; i<param->nval; i++) {
      if (ibound_.ibd[i]==1) {
	
	rp+=sprintf(rp, " AE-NL- %3d\t%6.3f\t%19.10f\t%14.10f\t\n",
		    atomic_.nlm[i], atomic_.wnl[i], 1000.0*(temp_eigen[i]-atomic_.en[i]),
		    1000.0*(temp_norm[i]-results_.rnorm[i]));

	temp_eigen_tot += fabs(temp_eigen[i]-atomic_.en[i]);
	temp_norm_tot += fabs(temp_norm[i]-results_.rnorm[i]);
      }
    }
    
    rp+=sprintf(rp, " AE-NL-  total error =\t%19.10f\t%14.10f\n\n",
		1000*temp_eigen_tot, 1000*temp_norm_tot);
    
    rp+=sprintf(rp," =====================================================================\n");
  }    
  
  return rp;
}

void readPS(param_t *param) {
  int i,j,junk;
  FILE *fp;
  char filename[160];
  double junk2;
  
  sprintf(filename, "%s.psi_ps", param->name);
  fp = fopen(filename, "rb");
  for (i=0; i<param->nval; i++) {
    fread(atomic_.rnl[i], sizeof(double), param->ngrid, fp);
  }
  fclose(fp);

  sprintf(filename, "%s.pot_ps", param->name);
  fp = fopen(filename, "rb");
  for (i=0; i<param->nval; i++) {
    fread(totpot_.rvcore[i], sizeof(double), param->ngrid, fp);
  }	
  fread(totpot_.rvcoul, sizeof(double), param->ngrid, fp);
  fclose(fp);

  for (i=0; i<param->nval; i++) 
    for (j=0; j<param->ngrid; j++) 
      totpot_.rvps[i][j]=totpot_.rvcore[i][j]+totpot_.rvcoul[j] ;

  sprintf(filename, "%s.eig_ae", param->name);
  fp = fopen(filename, "rb");
  for (i=0; i<param->nval; i++){
    fread(&atomic_.en[i], sizeof(double),1, fp);
    fread(&junk2, sizeof(double),1, fp);
    fread(&junk, sizeof(int),1, fp);
    fread(&nmax_.nmax[i], sizeof(int),1, fp);
    fread(&nmax_.maxim, sizeof(int),1, fp);
    fread(&atomic_.xion, sizeof(double),1, fp);
    fread(&ibound_.ibd[i], sizeof(int),1, fp);
    fread(&ensave_.ensave[i], sizeof(double),1, fp);
  }
  fclose(fp);

  if (param->rpcc>1e-12){
    sprintf(filename, "%s.rho_pcore", param->name);
    fp = fopen(filename, "rb");
    fread(rscore_.rscore, sizeof(double), param->ngrid, fp);
    fread(rscore_.rdd, sizeof(double), param->ngrid, fp);
    fread(rscore_.rddd, sizeof(double), param->ngrid, fp);
    fclose(fp);

    sprintf(filename, "%s.rho_fcore", param->name);
    fp = fopen(filename, "wb");
    fread(rscore_.rscoretot, sizeof(double), param->ngrid, fp);
    fclose(fp);
  }else{
    sprintf(filename, "%s.rho_fcore", param->name);
    fp = fopen(filename, "wb");
    fread(rscore_.rscore, sizeof(double), param->ngrid, fp);
    fclose(fp);
  }
}

void writeNL(param_t *param) {

  int i,j;
  int iset=0;
  FILE *fp;
  char filename[160];

  /* now dump the local potential into a binary file */
  sprintf(filename, "%s.loc", param->name);
  fp = fopen(filename, "wb");
  fwrite(nlcore_.rvloc, sizeof(double), param->ngrid, fp);
  fclose(fp);
  
  sprintf(filename, "%s.psi_nl", param->name);
  fp = fopen(filename, "wb");
  for (i=0; i<param->nval; i++)
    fwrite(atomic_.rnl[i], sizeof(double), param->ngrid, fp);
  fclose(fp);

  sprintf(filename, "%s.pcc_plt", param->name);
  fp = fopen(filename, "a");
  for (j=0;j<param->ngrid;j++)
    fprintf(fp,"%lg %lg \n",grid_.r[j],valden_.rsval[j]);
  fprintf(fp,"@ \n");
  fclose(fp);

  sprintf(filename, "%s.nl_plt", param->name);
  fp = fopen(filename, "w");
  for (i=0; i<param->nval;i++) {
    if (ibound_.ibd[i]==1) {
      for (j=0;j<param->ngrid;j++)
	if ((grid_.r[j]<param->rc[i])||(iset)) {
	  fprintf(fp,"%lg %lg 0.0\n",grid_.r[j],atomic_.rnl[i][j]);
	}else{
	  fprintf(fp,"%lg %lg 1e-8\n",grid_.r[j],atomic_.rnl[i][j]);
	  iset=1;
	}
    } else {
      j=0;
      while(grid_.r[j] < param->rc[i]) {
	fprintf(fp,"%lg %lg 0.0\n",grid_.r[j],atomic_.rnl[i][j]);
	j++;
      }
    }
    fprintf(fp,"@ \n");
    iset=0;
  }
  fclose(fp);

  sprintf(filename, "%s.vi_plt", param->name);
  fp = fopen(filename, "a");
  for (j=0;j<param->ngrid;j++){
    fprintf(fp,"%lg %lg 0.0\n",grid_.r[j],nlcore_.rvloc[j]/grid_.r[j]);
  }
  fprintf(fp,"@ \n");
  fclose(fp);
  
  sprintf(filename, "%s.vs_plt", param->name);
  fp = fopen(filename, "a");
  for (j=0;j<param->ngrid;j++){
    fprintf(fp,"%lg %lg 0.0\n",grid_.r[j],(nlcore_.rvloc[j]+totpot_.rvcoul[j])/grid_.r[j]);
  }
  fprintf(fp,"@ \n");
  fclose(fp);
  
}
