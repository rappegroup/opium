/*
 * $Id: do_ps.c,v 1.6 2004/06/16 20:46:17 mbarnes Exp $
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "parameter.h"
#include "fortparam.h"        /* fortran code parameters */
#include "do_ps.h"
#include "common_blocks.h"
#include "nlm.h"

/* fortran prototypes Should this be somewher else?????? */
void optim_(int *);
void kerker_(int *);

/* report feature */

static char report[800];

void readAE(param_t *param);
char * write_reportps(param_t *param , char *rp);
void writePS(param_t *param);

int do_ps(param_t *param, char *logfile){

  int i;  
  FILE *fp_log;
  char *rp=report;

  fp_log = fopen(logfile, "a");

  fprintf(fp_log,"\n ======================================================================== \n");
  fprintf(fp_log," Begin PS construction\n");
  fprintf(fp_log," ======================================================================== \n");
  fclose(fp_log);  

  /* set the log file */
  sprintf(filenames_.file_log, "%s", logfile);

  for (i=0; i<param->nval; i++) {
    ibound_.ibd[i]=param->ibound[param->norb-param->nval + i];
    ensave_.ensave[i]=param->ensave[param->norb-param->nval + i];
  }

  atomic_.norb=param->nll;
  np_.nvales=param->nll;

  readAE(param);

  if (param->psmeth == 'o') {
    if (param->optmeth == 'c') {
      opt_.meth=0;
    } else { 
      opt_.meth=1;
    }
    optim_(&param->ixc);
  } else if(param->psmeth == 'k') {
    kerker_(&param->ixc);
  } else {
    fp_log = fopen(logfile, "a");
    fprintf(fp_log, "   PS constuction method not known (must be optimized or kerker) \n");
    fclose(fp_log);
    exit(1);
  }

  rp=write_reportps(param,rp);
  writePS(param);

  
  fp_log = fopen(logfile, "a");

  fprintf(fp_log,"\n ======================================================================== \n");
  fprintf(fp_log," End PS construction\n");
  fprintf(fp_log," ======================================================================== \n");


  fclose(fp_log);
  return 0;
}


/* report section */

void do_ps_report(FILE *fp){
  fprintf(fp, "%s", report);
}

void readAE(param_t *param) {
  int i,j,ic;
  FILE *fp;
  char filename[160];
  
  ic=0;
  sprintf(filename, "%s.psi_ae", param->name);
  fp = fopen(filename, "rb");
  for (i=0; i<param->nval; i++) {
    if (i == param->ipot[ic]) {
      fread(atomic_.rnl[ic], sizeof(double), param->ngrid, fp);
      ic++;
    }else{
      fseek(fp,sizeof(double)*param->ngrid,1);
    }
  }
  fclose(fp);

  ic=0;
  sprintf(filename, "%s.pot_ae", param->name);
  fp = fopen(filename, "rb");
  for (i=0; i<param->nval; i++) {
    if (i == param->ipot[ic]) {
      fread(totpot_.rvcore[ic], sizeof(double), param->ngrid, fp);
      ic++;
    }else{
      fseek(fp,sizeof(double)*param->ngrid,1);
    }
    fread(totpot_.rvcoul, sizeof(double), param->ngrid, fp);
  }
  fclose(fp);

  for (i=0; i<param->nll; i++) 
    for (j=0; j<param->ngrid; j++) 
      totpot_.rvps[i][j]=totpot_.rvcore[i][j]+totpot_.rvcoul[j] ;

  ic=0;
  sprintf(filename, "%s.eig_ae", param->name);
  fp = fopen(filename, "rb");
  for (i=0; i<param->nval; i++){
    if (i == param->ipot[ic]) {      
      fread(&atomic_.en[ic], sizeof(double),1, fp);
      fread(&atomic_.wnl[ic], sizeof(double),1, fp);
      fread(&atomic_.nlm[ic], sizeof(int),1, fp);
      fread(&nmax_.nmax[ic], sizeof(int),1, fp);
      fread(&nmax_.maxim, sizeof(int),1, fp);
      fread(&atomic_.xion, sizeof(double),1, fp);
      ic++;
    }else{
      fseek(fp,3*(sizeof(double)+sizeof(int)),1);
    }
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

char * write_reportps(param_t *param , char *rp) {
  
  int i;
  double tot_conv_error = 0.;

  if (param->psmeth=='o') {
  
    rp+=sprintf(rp, 
		"    Orbital  Conv. error: [mRy/e]            [mRy]             [meV]     Ghost\n"
		"    --------------------------------------------------------------------------\n");
    for (i=0; i<param->nval; i++){
      if (convrpt_.sumc[i] > 1e-12) {
	rp+=sprintf(rp, "\t%3d        %16.10f  %16.10f  %16.10f\t%6s\n",
		    atomic_.nlm[i], convrpt_.sumc[i]*1000., convrpt_.wsumt[i]*1000., 
		    convrpt_.wsumt[i]*13.6058*1000.,
		    (convrpt_.lghost[i]==0)?"no":"yes");
	tot_conv_error += convrpt_.wsumt[i]*1000.;
      } else {
	/*	rp+=sprintf(rp, "\t%3d  (unbound)   --------          --------          --------       -- \n",
		atomic_.nlm[i]); */
      }
    }
    rp+=sprintf(rp, 
		"\n                  Tot. error =     %18.10f  %16.10f\n",
		tot_conv_error, tot_conv_error*13.6058);
  } else {

    rp+=sprintf(rp, 
		"    Orbital       Ghost\n"
		"    --------------------------------------------------------------------------\n");
    for (i=0; i<param->nval; i++){
      rp+=sprintf(rp, "\t%3d\t%6s\n",
		  atomic_.nlm[i],(convrpt_.lghost[i]==0)?"no":"yes");
    }
  }
  
  return rp;
}

void writePS(param_t *param) {

  int i,j,ncore; 
  FILE *fp;
  char filename[160];
  double zeff;

  /* This is calculated in too many places */
  ncore=param->norb-param->nval;
  zeff=atomic_.xion;
  for (i=0; i<param->nll; i++) {
    zeff +=atomic_.wnl[i];
  }

  sprintf(filename, "%s.psi_ps", param->name);
  fp = fopen(filename, "wb");
  for (i=0; i<param->nval; i++) 
    for (j=0; j<param->nll; j++) 
      if (nlm_label(atomic_.nlm[j]).l == nlm_label(param->nlm[i+ncore]).l) {
	fwrite(atomic_.rnl[j], sizeof(double), param->ngrid, fp);
      }
  fclose(fp);  

  sprintf(filename, "%s.plt_ips", param->name);
  fp = fopen(filename, "w");

  for (i=0; i<param->nll;i++) {
    for (j=0;j<param->ngrid;j++) {
      fprintf(fp,"%lg %lg \n",grid_.r[j],totpot_.rvcore[i][j]/grid_.r[j]);
    }
    fprintf(fp,"@ \n");
  }
  /* add the AE ionic potential */
  for (j=0;j<param->ngrid;j++) {
    fprintf(fp,"%lg %lg \n",grid_.r[j],-2.0*zeff/grid_.r[j]);
  }
  fprintf(fp,"@ \n");
  fclose(fp);

  sprintf(filename, "%s.plt_sps", param->name);
  fp = fopen(filename, "w");

  for (i=0; i<param->nll;i++) {
    for (j=0;j<param->ngrid;j++) {
      fprintf(fp,"%lg %lg \n",grid_.r[j],(totpot_.rvcore[i][j]+totpot_.rvcoul[j])/grid_.r[j]);
    }
    fprintf(fp,"@ \n");
  }
  fclose(fp);

  sprintf(filename, "%s.pot_ps", param->name);
  fp = fopen(filename, "wb");
  for (i=0; i<param->nval; i++) 
    for (j=0; j<param->nll; j++) 
      if (nlm_label(atomic_.nlm[j]).l == nlm_label(param->nlm[i+ncore]).l) {
	fwrite(totpot_.rvcore[j], sizeof(double), param->ngrid, fp);

      }
  fwrite(totpot_.rvcoul, sizeof(double), param->ngrid, fp);
  fclose(fp);
  
}
