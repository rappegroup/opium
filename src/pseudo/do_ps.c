/*
 * Copyright (c) 1998-2008 The OPIUM Group
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

#include "parameter.h"
#include "cdim.h"        /* fortran code parameters */
#include "do_ps.h"
#include "common_blocks.h"
#include "nlm.h"

#define streq(a,b) (!strcasecmp(a,b))

/* fortran prototypes Should this be somewhere else?????? */
void optim_(int *, int *, double *);
void kerker_(int *, int *, double *);
void tmsub_(int *, int *, double *);
/* report feature */

static char report[8000];

void readAE(param_t *param);
char * write_reportps(param_t *param , char *rp);
void writePS(param_t *param);
void nrelsproj(param_t *param, char *);
void relsproj(param_t *param, char *);
void hfsmooth_(int * ,double *, int *);

int do_ps(param_t *param, char *logfile){

  int i;  
  FILE *fp_log;
  char *rp=report;
  int irel;
  char filename[80];
  int qpopt;
  double rlocalr;

  fp_log = fopen(logfile, "a");

  fprintf(fp_log,"\n ======================================================================== \n");
  fprintf(fp_log," Begin PS construction\n");
  fprintf(fp_log," ======================================================================== \n");
  fclose(fp_log);  

  /* set the log file */
  sprintf(filenames_.file_log, "%s", logfile);

  irel= (!strcmp(param->reltype, "nrl")) ? 0:1;

  if ((irel==0)||(param->ixc >= 0)) {
    nrelsproj(param,logfile);
  } else {
    relsproj(param,logfile);
  }
  readAE(param);

  aorb_.nval=param->nll;
  aorb_.norb=param->nll;
  aorb_.ncore=0;

  if (param->psmeth == 'o') {
    if (param->optmeth == 'c') {
      psdat_.meth=0;
    } else { 
      psdat_.meth=1;
    }
    optim_(&param->ixc, &irel, &param->exccut);
  } else if(param->psmeth == 'k') {
    kerker_(&param->ixc,&irel, &param->exccut);
  } else if(param->psmeth == 't') {
    tmsub_(&param->ixc, &irel, &param->exccut);
  } else {
    fp_log = fopen(logfile, "a");
    fprintf(fp_log, "   PS constuction method not known (must be (o)ptimized, (t)m, or (k)erker) \n");
    fclose(fp_log);
    exit(1);
  }

  param->nll=aorb_.nval;
  psdat_.nll = param->nll;    

  if ( (streq(param->xcparam,"hf"))&&(param->qpopt > 0)) 
    hfsmooth_(&param->qpopt,param->rlocalr,&param->ixc);

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


void nrelsproj(param_t *param, char *logfile) {
  
  int i,j,ic,jc,ii,jj,ncore;
  FILE *fp_log;
  int lc[7];

  ncore=param->norb - param->nval;
  param->nll=0;

  for (i=0;i<7;i++) {
    param->npot[i]=0;
    lc[i]=0;
  }

  for (i=0;i<param->nval;i++) {
    param->lpot[i]=nlm_label(param->nlm[i+ncore]).l;
  }

  for (i=0;i<param->nval;i++) { 
    for (j=0;j<4;j++) {
      if (param->lpot[i]==j) {
	param->npot[i]=lc[j];
	lc[j]++;
      }
    }
  }

  j=0;
  for (i=0;i<param->nval;i++) { 
    if (param->npot[i]==0) {
      aval_.rcall[j] = param->rc[i];
      optparam_.qcl[j] = param->qc[i];
      optparam_.nbl[j] = param->nb[i];
      param->nll++;
      j++;
    }

  }

}

void relsproj(param_t *param, char *logfile) {
  
  int i,j,ic,jc,ii,jj,ncore;
  FILE *fp_log;
  int lc[4];

  ncore=aorb_.norb - aorb_.nval;
  param->nll=0;

  for (i=0;i<4;i++) {
    param->npot[i]=0;
    lc[i]=0;
  }

  for (i=0;i<aorb_.nval;i++) {
    param->lpot[i]=aorb_.lo[i+aorb_.ncore];
  }

  for (i=0;i<aorb_.nval;i++) { 
    for (j=0;j<4;j++) {
      if (param->lpot[i]==j) {
	param->npot[i]=lc[j];
	lc[j]++;
      }
    }
  }

  for (i=0;i<aorb_.nval;i++) { 
    if (param->npot[i]==0) param->nll++;
    if ((param->npot[i]==1)&&(param->lpot[i]>0)){
      param->nll++;
      param->npot[i]=0;
    }
  }
}


void readAE(param_t *param) {
  int i,j,ic;
  FILE *fp;
  char filename[160];

  sprintf(filename, "%s.psi_ae", param->name);
  fp = fopen(filename, "rb");
  ic=0;
  for (i=0; i<aorb_.nval; i++) {
    if (param->npot[i]==0) {
      fread(wfn_.rnl[ic], sizeof(double), param->ngrid, fp);
      ic++;
    }else{
      fseek(fp,sizeof(double)*param->ngrid,1);
    }
  }
  fclose(fp);

  sprintf(filename, "%s.pot_ae", param->name);
  fp = fopen(filename, "rb");
  ic=0;
  for (i=0; i<aorb_.nval; i++) {
    if (param->npot[i]==0) {
      fread(totpot_.rvcore[ic], sizeof(double), param->ngrid, fp);
      fread(totpot_.rvps[ic], sizeof(double), param->ngrid, fp);
      ic++;
    }else{
      fseek(fp,sizeof(double)*param->ngrid,2);
    }
  }
  fclose(fp);

  /*  for (j=0; j<param->ngrid; j++) 
      totpot_.rvcoul[j] = totpot_.rvps[0][j]-totpot_.rvcore[0][j];*/
  
  sprintf(filename, "%s.eig_ae", param->name);
  fp = fopen(filename, "rb");
  ic=0;
  for (i=0; i<aorb_.nval; i++){
    if (param->npot[i]==0) {
      fread(&adat_.en[ic], sizeof(double),1, fp);
      fread(&adat_.wnl[ic], sizeof(double),1, fp);
      fread(&aorb_.nlm[ic], sizeof(int),1, fp);
      fread(&aorb_.no[ic], sizeof(int),1, fp);
      fread(&aorb_.lo[ic], sizeof(int),1, fp);
      fread(&adat_.so[ic], sizeof(double),1, fp);
      fread(&aorb_.nmax[ic], sizeof(int),1, fp);
      fread(&aorb_.maxim, sizeof(int),1, fp);
      fread(&adat_.xion, sizeof(double),1, fp);
      fread(&aval_.ibd[ic], sizeof(int),1, fp);
      ic++;
    }else{
      fseek(fp,sizeof(double)*4+sizeof(int)*6,1);
    }
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

char * write_reportps(param_t *param , char *rp) {
  
  int i,igh;
  double tot_conv_error = 0.;
  double tot_conv_p = 0.;

  if (param->ixc<0) 
    rp+=sprintf(rp, "\n\n NOTICE!! :Ghost testing not done for HF psps yet, sorry :( \n");

  if (param->psmeth=='o') {
    
    rp+=sprintf(rp,
		"    ====================Optimized pseudopotential method====================\n\n");   
    rp+=sprintf(rp, 
		"                       Pseudopotential convergence error                      \n");
    rp+=sprintf(rp, 
		"    Orbital      [mRy/e]       [meV/e]         [mRy]        [meV]        Ghost\n"
		"    --------------------------------------------------------------------------\n");
    for (i=0; i<param->nll; i++){
      if (param->ixc<0) {	
	rp+=sprintf(rp, "\t%3d  %12.6f  %12.6f  %12.6f  %12.6f\t%6s\n",
		    aorb_.nlm[i], psout_.sumc[i]*1000.,psout_.sumc[i]*1000*13.6057,
                    psout_.wsumt[i]*1000., 
		    psout_.wsumt[i]*13.6057*1000.,"???");
	tot_conv_error += psout_.wsumt[i]*1000.;

      } else if (aval_.ibd[i]==1) {
	rp+=sprintf(rp, "\t%3d  %12.6f  %12.6f  %12.6f  %12.6f\t%6s\n",
		    aorb_.nlm[i], psout_.sumc[i]*1000.,psout_.sumc[i]*1000*13.6057,
                    psout_.wsumt[i]*1000., 
		    psout_.wsumt[i]*13.6057*1000.,
		    (psout_.npsghost[i]>0)?"yes":((psout_.npsghost[i]<0)?"?":"no"));
	tot_conv_error += psout_.wsumt[i]*1000.;
      } else {
	rp+=sprintf(rp, "\t%3d     (unbound)      --------      --------      --------     %6s\n",
		    aorb_.nlm[i],(psout_.npsghost[i]>0)?"yes":((psout_.npsghost[i]<0)?"?":"no"));
      }
    }
    rp+=sprintf(rp, 
		"\n                  Tot. error =           %12.6f  %12.6f\n",
		tot_conv_error, tot_conv_error*13.6057);
  } else {
    
    if (param->psmeth=='t') {
      rp+=sprintf(rp,
		  "    ==========Troullier-Martins Pseudodoptential method==========\n\n");   
    }
    
    if (param->psmeth=='k') {
      rp+=sprintf(rp,
		  "    ===============Kerker Pseudodoptential method================\n\n");   
    }
    
    rp+=sprintf(rp, 
		"    Orbital       Ghost\n"
		"    -----------------------------------------------------------------\n");
    for (i=0; i<param->nll; i++){
      if (param->ixc<0) {	
	rp+=sprintf(rp, "\t%3d\t%6s\n",
		    aorb_.nlm[i],"???");
      } else {
	rp+=sprintf(rp, "\t%3d\t%6s\n",
		    aorb_.nlm[i],(psout_.npsghost[i]==0)?"no":"yes");
      }
    }    
  }
  rp+=sprintf(rp, "\n");
  igh=0;
  for (i=0; i<param->nll; i++){
    if (psout_.npsghost[i]==0) {
      /*rp+=sprintf(rp, "\t%3d ok as the local potential \n",aorb_.nlm[i]);*/
    } else {
      rp+=sprintf(rp, "\t%3d SHOULD NOT be used as the local potential \n",aorb_.nlm[i]);
      igh+=1;
    }
  }
  if (igh==param->nll) {
    rp+=sprintf(rp, "\t !!ERROR!! There are no choices for local potential\n");
  }
  
  return rp;
}

  
void writePS(param_t *param) {

  int i,ii,j,ncore; 
  int iset=0;
  FILE *fp;
  char filename[160];
  double zeff;

  ncore=param->norb-param->nval;
  zeff=adat_.xion;
  for (i=0; i<param->nll; i++) {
    zeff +=adat_.wnl[i];
  }

  sprintf(filename, "%s.vi_plt", param->name);
  fp = fopen(filename, "w");

  for (i=0; i<param->nll;i++) {
    for (j=0;j<param->ngrid;j++) {
      if ((grid_.r[j]<param->rc[i])||(iset)) {
	fprintf(fp,"%lg %lg 0.0\n",grid_.r[j],totpot_.rvcore[i][j]/grid_.r[j]);
      }else{
	fprintf(fp,"%lg %lg 1e-8\n",grid_.r[j],totpot_.rvcore[i][j]/grid_.r[j]);
	iset=1;
      }
    }
    fprintf(fp,"@ \n");
    iset=0;
  }
  /* add the AE ionic potential */
  for (j=0;j<param->ngrid;j++) {
    fprintf(fp,"%lg %lg 0.0\n",grid_.r[j],-2.0*zeff/grid_.r[j]);
  }
  fprintf(fp,"@ \n");
  fclose(fp);

  sprintf(filename, "%s.vs_plt", param->name);
  fp = fopen(filename, "w");

  for (i=0; i<param->nll;i++) {
    for (j=0;j<param->ngrid;j++) {
      if ((grid_.r[j]<param->rc[i])||(iset)) {
	fprintf(fp,"%lg %lg 0.0\n",grid_.r[j],(totpot_.rvcore[i][j]+totpot_.rvcoul[j])/grid_.r[j]);
      }else{
	fprintf(fp,"%lg %lg 1e-8\n",grid_.r[j],(totpot_.rvcore[i][j]+totpot_.rvcoul[j])/grid_.r[j]);
	iset=1;
      }
    }
    fprintf(fp,"@ \n");
    iset=0;
  }
  fclose(fp);


  for (j=0; j<param->nll; j++) {
    sprintf(filename, "%s.psi.ps.l=%d", param->name,nlm_label(param->nlm[j+ncore]).l);
    fp = fopen(filename, "wb");
    fwrite(wfn_.rnl[j], sizeof(double), param->ngrid, fp);
    fclose(fp);
  }      

  for (j=0; j<param->nll; j++) {
    sprintf(filename, "%s.pot.ps.l=%d", param->name,nlm_label(param->nlm[j+ncore]).l);
    fp = fopen(filename, "wb");
    fwrite(totpot_.rvcore[j], sizeof(double), param->ngrid, fp);
    fwrite(totpot_.rvps[j], sizeof(double), param->ngrid, fp);
    fclose(fp);
    /*    printf(" at dops: j psi pot %d %lg %lg \n", j, wfn_.rnl[j][0],totpot_.rvcore[j][0]);*/
  }

}

