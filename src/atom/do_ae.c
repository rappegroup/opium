/* 
 * $Id: do_ae.c,v 1.6 2004/07/06 15:29:50 ewalter Exp $
 */

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>

#include "parameter.h"
#include "nlm.h"
#include "fortparam.h"
#include "do_ae.h"
#include "common_blocks.h"    /* fortran common blocks */
#include "energy.h"

/* report feature */

static char report[800];
void startae(param_t *param); 
void relorb(param_t *param, int); 
void nrelorbae(param_t *param, int); 
char * write_reportae(param_t *param, char *rp,int,double temp_eigen[], double temp_norm[]);
void writeAE(param_t *param);

/* fortran prototypes Should this be somewher else?????? */
void scpot_(double  *, int * ,int *);
void atm_(double *, int *);
void interp_(int *, int *, int *);
void getpcc_(void);

int do_ae(param_t *param, char *logfile){
  
  int i;
  int ncore;
  FILE *fp_log;
  char *rp = report;
  int ifrl;
  int config=-1;
  int ipsp=0;
  double temp_eigen[10];
  double temp_norm[10];


  /* clear report buffer */
  *rp = 0;

  /* set log file name */
  sprintf(filenames_.file_log, "%s", logfile);

  ncore = param->norb-param->nval;
  npm_.ncores=ncore;

  fp_log=fopen(logfile,"a");

  fprintf(fp_log,"\n ======================================================================== \n");
  fprintf(fp_log," Begin AE calculation\n");
  fprintf(fp_log," ======================================================================== \n");

  if (!strcmp(param->reltype, "nrl")){
    
    /* do the NRL type AE solve */

    fprintf(fp_log," Performing non-relativistic AE calculation...  \n" );
    fclose(fp_log);

    /* set up the nlm's ; config=-1 is the reference state */
    nrelorbae(param,config);
    
    /* starting guess for the AE potential */
    startae(param);

    /* find the SCF solution */
    scpot_(&param->z,&param->ixc,&ipsp);

    for (i=0; i<param->norb; i++){
      param->ibound[i]=ibound_.ibd[i];
    }

  } else {

    fprintf(fp_log," Performing relativistic AE calculation...  \n" );
    fclose(fp_log);

    /* turn the l-orbitals into j-orbitals */
    relorb(param,config);

    /* find the SCF solution */
    atm_(&param->z,&param->ixc);

    ifrl = (!strcmp(param->reltype,"frl"))?1:0;
    param->aereltransf = 0;
    nrelorbae(param,config);
    interp_(&ifrl, &param->aereltransf,&param->ixc);
    
  }

  enae[0] = results_.etot;

  if (param->rpcc>1e-12) {
    fp_log = fopen(logfile, "a");
    fprintf(fp_log, "   ================================================\n");
    fprintf(fp_log,"<<<pseudizing core>>>\n");
    fclose(fp_log);
    getpcc_();
    fp_log = fopen(logfile, "a");
    fprintf(fp_log, "   ================================================\n");
    fclose(fp_log);
  }
  /* record data in report */
  rp = write_reportae(param,rp,config,temp_eigen,temp_norm);
  writeAE(param);
  
  fp_log = fopen(logfile, "a");
  fprintf(fp_log,"\n ======================================================================== \n");
  fprintf(fp_log," End AE calculation\n");
  fprintf(fp_log," ======================================================================== \n");

  fclose(fp_log);

  return 0;
}

/* report section */

void do_ae_report(FILE *fp){
  fprintf(fp, "%s", report);
}

void startae(param_t *param) {
  
  double y=0.0;
  int i,j;
  
  if ((param->z - atomic_.xion - 1.0) > 0.0) 
    y=pow((param->z - atomic_.xion - 1.0),0.4);
  
  for (i=0; i<param->ngrid; i++) {
    double t=grid_.r[i]/0.675;
    if (t>100) 
      t=100.0;
    t=2.0*(param->z - atomic_.xion - 1.0) * (1.0 - 1.0/(0.675*y*(exp(t)-1.0)+1.0))-2*param->z;
    for (j=0; j<atomic_.norb; j++){
      totpot_.rvps[j][i] = t;
      totpot_.rvcore[j][i] = -2.0* param->z; 
    }
  }
}

void nrelorbae(param_t *param, int config) {
  
  int i;
  int nrcore;

  nrcore=param->norb - param->nval;
  
  atomic_.norb = param->norb;
  atomic_.xion = param->z;

  for (i=0; i<nrcore; i++){
    atomic_.nlm[i] = param->nlm[i];
    atomic_.wnl[i] = param->wnl[i];
    atomic_.en[i] = param->en[i];    
    atomic_.xion -= param->wnl[i];
    ensave_.ensave[i]=param->ensave[i];
    ibound_.ibd[i]=param->ibound[i];
  }
  for (i=nrcore; i<param->norb; i++){
    if (config<0) {  
      atomic_.nlm[i] = param->nlm[i];
      atomic_.wnl[i] = param->wnl[i];
      atomic_.en[i] = param->en[i];    
      ensave_.ensave[i]=param->ensave[i];
    }else{
      atomic_.nlm[i] = param->nlm_conf[config][i-nrcore];
      atomic_.wnl[i] = param->wnl_conf[config][i-nrcore];
      atomic_.en[i] = param->en_conf[config][i-nrcore];
      ensave_.ensave[i]=param->ensave_conf[config][i-nrcore];
    }
    ibound_.ibd[i]=param->ibound[i];
    atomic_.xion -= atomic_.wnl[i];
  }
}


void relorb(param_t *param, int config) {
  
  int i,j;
  double sc;
  int nrcore;
  
  nrcore=param->norb - param->nval;
  atomic_.norb = param->norb;

  sc=-0.5;
  reli_.norb=0;
  reld_.zval=0;
  reld_.zcore=0;

  for (i=0; i<nrcore; i++){
    for (j=0; j<2; j++){
      reli_.no[reli_.norb] = nlm_label(param->nlm[i]).n;
      reli_.lo[reli_.norb] = nlm_label(param->nlm[i]).l;
      reld_.so[reli_.norb] = sc;
      reld_.zo[reli_.norb] = 2.0*(nlm_label(param->nlm[i]).l+sc)+1.0;
      reld_.zcore += reld_.zo[reli_.norb];
      if (fabs(reld_.zo[reli_.norb])>0.1) reli_.norb++;
      sc=-1.0*sc;
    }
  }
  sc=-0.5;
  reli_.ncore = reli_.norb;

  for (i=nrcore; i<param->norb; i++){
    for (j=0; j<2; j++){
      if (config <0) {
	reli_.no[reli_.norb] = nlm_label(param->nlm[i]).n;
	reli_.lo[reli_.norb] = nlm_label(param->nlm[i]).l;
	reld_.so[reli_.norb] = sc;
	reld_.zo[reli_.norb] = param->wnl[i]*(2*(reli_.lo[reli_.norb]+reld_.so[reli_.norb])+1)/(4*reli_.lo[reli_.norb]+2);
      }else{
	reli_.no[reli_.norb] = nlm_label(param->nlm_conf[config][i-nrcore]).n;
	reli_.lo[reli_.norb] = nlm_label(param->nlm_conf[config][i-nrcore]).l;
	reld_.so[reli_.norb] = sc;
	reld_.zo[reli_.norb] = param->wnl_conf[config][i-nrcore]*(2*(reli_.lo[reli_.norb]
					 +reld_.so[reli_.norb])+1)/(4*reli_.lo[reli_.norb]+2);
      }
      reld_.zval += reld_.zo[reli_.norb];
      if (reli_.lo[reli_.norb]+reld_.so[reli_.norb]>0.0) reli_.norb++;
      sc=-1.0*sc;
    }
  }

  reli_.nval =  reli_.norb - reli_.ncore;
  
  for (i=0; i< reli_.norb; i++) {
    for (j=0; j< reli_.norb; j++) {
      if (i==j) continue;
      if (reli_.no[i] != reli_.no[j]) continue;
      if (reli_.lo[i] != reli_.lo[j]) continue;
      if (fabs(reld_.so[i] - reld_.so[j])>0.001) continue;
      printf("%d %d two orbitals the same! -- abort! \n",i,j);
      exit(1);
    }
  }
}

char * write_reportae(param_t *param, char *rp,int config, double temp_eigen[], double temp_norm[]) {
  
  int i,j,l,ncore;
  double e_rel,n_rel;

  ncore=param->norb - param->nval;

  rp += sprintf(rp, 
		" AE:Orbital    Filling       Eigenvalues[Ry]          Norm\n"
		"    ----------------------------------------------------------\n");

  if (!strcmp(param->reltype, "nrl")){
    for (i=0; i<param->norb; i++)
      if (i>ncore-1) {
	rp += sprintf(rp, "\t%3d\t%6.3f\t%19.10f\t%14.10f\n",
		      atomic_.nlm[i], atomic_.wnl[i], atomic_.en[i],results_.rnorm[i]);
	temp_eigen[i-ncore]=atomic_.en[i];
	temp_norm[i-ncore]=results_.rnorm[i];
      }else{
	rp += sprintf(rp, "\t%3d\t%6.3f\t%19.10f\n",
		      atomic_.nlm[i], atomic_.wnl[i], atomic_.en[i]);
      }

  }else if (!strcmp(param->reltype, "srl")){
    j=0;

    if (config<0) {
      for (i=0; i<param->nval; i++)
	rp += sprintf(rp, "\t%3d\t%6.3f\t%19.10f\t%14.10f\n",
		      atomic_.nlm[i], atomic_.wnl[i], atomic_.en[i],
		      frlwf_.rnor[i]);
    }else{
      for (i=0; i<param->nval; i++){
	if ((l=nlm_label(param->nlm_conf[config][i]).l)){
	  /* need to average the eigenvalue */
	  e_rel = (l*atmwave_.eev[j] + (l+1)*atmwave_.eev[j+1])
	    / (double)(2.*l+1.);
	  n_rel = (l*rcrel_.relnorm[j] + (l+1)*rcrel_.relnorm[j+1])
	    / (double)(2.*l+1.);
	  ++j;
	}else{
	  e_rel = atmwave_.eev[j];
	  n_rel = rcrel_.relnorm[j];
	}
	temp_eigen[i]=e_rel;
	temp_norm[i]=n_rel;

	rp+=sprintf(rp, "\t%3d\t%6.3f\t%19.10f\t%14.10f\n",
		param->nlm_conf[config][i], param->wnl_conf[config][i], e_rel, n_rel);
	j++;
      }
    }
  }
  
  rp += sprintf(rp, "\n      E_tot = %19.10f Ry\n", results_.etot);

  return rp;
  
}

void writeAE(param_t *param) {

  int i,j,ncore;
  FILE *fp;
  char filename[160];

  /* New plot section */

  sprintf(filename, "%s.plt_pcc", param->name);
  fp = fopen(filename, "w");

  if (param->rpcc>1e-12){
    for (j=0;j<param->ngrid;j++)
      fprintf(fp,"%lg %lg \n",grid_.r[j],rscore_.rscoretot[j]);
    fprintf(fp,"@ \n");
  }else{  
    for (j=0;j<param->ngrid;j++)
      fprintf(fp,"%lg %lg \n",grid_.r[j],rscore_.rscore[j]);
    fprintf(fp,"@ \n");
  }
  if (param->rpcc>1e-12){
    
    for (j=0;j<param->ngrid;j++)
      fprintf(fp,"%lg %lg \n",grid_.r[j],rscore_.rscore[j]);
    fprintf(fp,"@ \n");
    
    /*    for (j=0;j<param->ngrid;j++)
	  fprintf(fp,"%lg %lg \n",grid_.r[j],rscore_.rdd[j]*grid_.r[j]*grid_.r[j]);
	  fprintf(fp,"@ \n");
	  
	  for (j=0;j<param->ngrid;j++)
	  fprintf(fp,"%lg %lg \n",grid_.r[j],rscore_.rddd[j]*grid_.r[j]*grid_.r[j]);
	  fprintf(fp,"@ \n"); */

  }
  fclose(fp);

  sprintf(filename, "%s.plt_ae", param->name);
  fp = fopen(filename, "w");

  ncore = param->norb-param->nval;  

  if (!strcmp(param->reltype, "nrl")){
    for (i=0; i<param->nval;i++) {

      if (ibound_.ibd[i+ncore]==1) {
	for (j=0;j<param->ngrid;j++)
	   fprintf(fp,"%lg %lg \n",grid_.r[j],atomic_.rnl[i+ncore][j]);
      } else {
	j=0;
	while(grid_.r[j] < param->rc[i]) {
	  fprintf(fp,"%lg %lg \n",grid_.r[j],atomic_.rnl[i+ncore][j]);
	  j++;

	}
      }
      fprintf(fp,"@ \n");
    }
  } else {
    for (i=0; i<param->nval; i++) {
      for (j=0;j<param->ngrid;j++)
	if (param->rc[i] < grid_.r[j]) fprintf(fp,"%lg %lg \n",grid_.r[j],atomic_.rnl[i][j]);
      fprintf(fp,"@ \n");
    }

    for (i=0; i<atmwave_.numorb;i++) {
      for (j=0;j<atmwave_.nnr;j++) {
	fprintf(fp,"%lg %lg \n",atmwave_.rr[j],atmwave_.aar[i][j]);
      }
      fprintf(fp,"@ \n");
    }
    /*    for (i=0; i<atmwave_.numorb;i++) {
	  for (j=0;j<atmwave_.nnr;j++) {
	  fprintf(fp,"%lg %lg \n",atmwave_.rr[j],atmwave_.bbr[i][j]);
	  }
	  fprintf(fp,"@ \n");
	  } */
  }
  fclose(fp);

  /* now dump the all-electron wave function into a binary file */
    
  sprintf(filename, "%s.psi_ae", param->name);
  fp = fopen(filename, "wb");
  if (!strcmp(param->reltype, "nrl")){
    ncore = param->norb-param->nval;
    for (i=0; i<param->nval; i++) 
      fwrite(atomic_.rnl[ncore+i], sizeof(double), param->ngrid, fp);
  }else{
    ncore = 0;
    for (i=0; i<param->nval; i++)
      fwrite(atomic_.rnl[ncore+i], sizeof(double), param->ngrid, fp);
    for (i=0; i<param->nval; i++){
      fwrite(relwf_.wfu[ncore+i], sizeof(double), param->ngrid, fp);
      fwrite(relwf_.wfd[ncore+i], sizeof(double), param->ngrid, fp);
    }
  }
  fclose(fp);
  
  /* now dump the all-electron potential into a binary file */
  
  sprintf(filename, "%s.pot_ae", param->name);
  fp = fopen(filename, "wb");
  fwrite(totpot_.rvcore[0], sizeof(double), param->ngrid, fp);
  fwrite(totpot_.rvcoul, sizeof(double), param->ngrid, fp);
  fclose(fp);

  sprintf(filename, "%s.eng", param->name);
  fp = fopen(filename, "wb");
  fwrite(&results_.etot, sizeof(double),1, fp);
  fclose(fp);

  sprintf(filename, "%s.eig_ae", param->name);
  fp = fopen(filename, "wb");

  /*this is a hack until interp is changed */
  if (!strcmp(param->reltype, "nrl")){
    ncore = param->norb-param->nval;
  }else{
    ncore = 0;
  }
  for (i=0; i<param->nval; i++) {
    fwrite(&atomic_.en[ncore+i], sizeof(double),1, fp);
    fwrite(&atomic_.wnl[ncore+i], sizeof(double),1, fp);
    fwrite(&atomic_.nlm[ncore+i], sizeof(int),1, fp);
    fwrite(&nmax_.nmax[ncore+i], sizeof(int),1, fp);
    fwrite(&nmax_.maxim, sizeof(int),1, fp);
    fwrite(&atomic_.xion, sizeof(double),1, fp);
  }
  fclose(fp);
  
  /* now dump the PCC into a binary file if NLCC requested */
  
  sprintf(filename, "%s.rho_core", param->name);
  fp = fopen(filename, "wb");
  fwrite(rscore_.rscore, sizeof(double), param->ngrid, fp);
  if (param->rpcc>1e-12){
    fwrite(rscore_.rdd, sizeof(double), param->ngrid, fp);
    fwrite(rscore_.rddd, sizeof(double), param->ngrid, fp);
    fwrite(rscore_.rscoretot, sizeof(double), param->ngrid, fp);
  }
  fclose(fp);

}
