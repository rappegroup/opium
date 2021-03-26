/* 
 * $Id: do_nl.c,v 1.6 2004/06/16 20:46:17 mbarnes Exp $
 */

/* standard libraries */
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "parameter.h"        /* defines structure: 'param_t' */
#include "nlm.h"              /* nlm_label call */
#include "fortparam.h"        /* fortran code parameters */
#include "do_nl.h"            /* the module's own header */
#include "common_blocks.h"    /* fortran common blocks */
#include "energy.h"           /* this is for saving the energies */

/* report feature */

static char report[800];
void nrelorbnl(param_t *param, int); 
char * write_reportnl(param_t *param, char *rp,int,double temp_eigen[], double temp_norm[]);
void scpot_(double  *, int * ,int *);
void readPS(param_t *param);
void writeNL(param_t *param);

int do_nl(param_t *param, char *logfile){

  int i;  
  FILE *fp_log;
  double zeff;
  int config=-1;
  int ipsp=1;
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
    ibound_.ibd[i]=param->ibound[param->norb-param->nval + i];
    zeff +=atomic_.wnl[i];
  }

  scpot_(&zeff,&param->ixc,&ipsp); 

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
    ibound_.ibd[i]=ibound_.ibd[i+nrcore];
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

  if (config<0) {
    rp+=sprintf(rp, 
		"\n NL:Orbital    Filling       Eigenvalues[Ry]          Norm       Ghost\n"
		"    ------------------------------------------------------------------\n");
    for (i=0; i<param->nval; i++)
      rp+=sprintf(rp, "\t%3d\t%6.3f\t%19.10f\t%14.10f\t%6s\n",
		  atomic_.nlm[i], atomic_.wnl[i], atomic_.en[i],
		  results_.rnorm[i], (results_.lghost[i]==0)?"no":"yes");
  }else{
    rp+=sprintf(rp, 
		"\n NL:Orbital    Filling       Eigenvalues[Ry]          Norm       Ghost\n"
		"    ------------------------------------------------------------------\n");
    for (i=0; i<param->nval; i++)
      rp+=sprintf(rp, "\t%3d\t%6.3f\t%19.10f\t%14.10f\t%6s\n",
		  atomic_.nlm[i], atomic_.wnl[i], atomic_.en[i],
		  results_.rnorm[i],(results_.lghost[i]==0)?"no":"yes");
    
  }
  rp+=sprintf(rp, "\n      E_tot = %19.10f Ry\n",
	      results_.etot);

  if (config>=0) {
    rp+=sprintf(rp, 
		"\n AE-NL:Orbital Filling       Eigenvalues[mRy]         Norm[1e-3] \n"
		" AE-NL- --------------------------------------------------------------\n");
    temp_eigen_tot = temp_norm_tot = 0.;

    for (i=0; i<param->nval; i++) {
      rp+=sprintf(rp, " AE-NL- %3d\t%6.3f\t%19.10f\t%14.10f\t\n",
		  atomic_.nlm[i], atomic_.wnl[i], 1000.0*(temp_eigen[i]-atomic_.en[i]),
		  1000.0*(temp_norm[i]-results_.rnorm[i]));
      
      temp_eigen_tot += fabs(temp_eigen[i]-atomic_.en[i]);
      temp_norm_tot += fabs(temp_norm[i]-results_.rnorm[i]);
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
  }
  fclose(fp);

  sprintf(filename, "%s.rho_core", param->name);
  fp = fopen(filename, "rb");
  fread(rscore_.rscore, sizeof(double), param->ngrid, fp);
  if (param->rpcc>1e-12){
    fread(rscore_.rscoretot, sizeof(double), param->ngrid, fp);
  }
  fclose(fp);
}

void writeNL(param_t *param) {

  int i,j;
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

  sprintf(filename, "%s.plt_pcc", param->name);
  fp = fopen(filename, "a");
  for (j=0;j<param->ngrid;j++)
    fprintf(fp,"%lg %lg \n",grid_.r[j],valden_.rsval[j]);
  fprintf(fp,"@ \n");
  fclose(fp);

  sprintf(filename, "%s.plt_nl", param->name);
  fp = fopen(filename, "w");
  for (i=0; i<param->nval;i++) {
    if (ibound_.ibd[i]==1) {
      for (j=0;j<param->ngrid;j++)
	fprintf(fp,"%lg %lg \n",grid_.r[j],atomic_.rnl[i][j]);
    } else {
      j=0;
      while(grid_.r[j] < param->rc[i]) {
	fprintf(fp,"%lg %lg \n",grid_.r[j],atomic_.rnl[i][j]);
	j++;
	
      }
    }
    fprintf(fp,"@ \n");
  }

  sprintf(filename, "%s.plt_ips", param->name);
  fp = fopen(filename, "a");
  for (j=0;j<param->ngrid;j++){
    fprintf(fp,"%lg %lg \n",grid_.r[j],nlcore_.rvloc[j]/grid_.r[j]);
  }
  fprintf(fp,"@ \n");
  fclose(fp);
  
  sprintf(filename, "%s.plt_sps", param->name);
  fp = fopen(filename, "a");
  for (j=0;j<param->ngrid;j++){
    fprintf(fp,"%lg %lg \n",grid_.r[j],(nlcore_.rvloc[j]+totpot_.rvcoul[j])/grid_.r[j]);
  }
  fprintf(fp,"@ \n");
  fclose(fp);
  
}
