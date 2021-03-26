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
 * $Id: do_tc.c,v 1.11 2004/10/10 22:10:06 ewalter Exp $
 */

/* standard libraries */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "parameter.h"        /* defines structure: 'param_t' */
#include "nlm.h"              /* nlm_label call */
#include "cdim.h"        /* fortran code parameters */
#include "do_tc.h"            /* the module's own header */
#include "common_blocks.h"    /* fortran common blocks */
#include "energy.h"           /* this is for saving the energies */

static char report[8000];
void startae(param_t *param);
void relorb(param_t *param, int); 
void nrelorbnl(param_t *param, int); 
void nrelorbae(param_t *param, int); 
char * write_reportae(param_t *param, char *rp,int, double temp_eigen[], double temp_norm[]);
char * write_reportnl(param_t *param, char *rp,int, double temp_eigen[], double temp_norm[]);
char * write_reportfc(param_t *param, char *rp,int, double temp_eigen[], double temp_norm[]);
void readPS(param_t *param);
void readAE(param_t *param);
void interp_(int *, int *, int *);
void scpot_(double  *, int * ,double *, int *, int *, int *);
void atm_(double *, int *, int *);

int do_tc(param_t *param, char *logfile, int job, int doifc){

  int i, j,k, kk; 
  FILE *fp_log;
  FILE *fp;
  static char filename[160];
  char *rp=report;
  int ncore,ifrl; 
  double zeff;
  double e,dele;
  int ifc=0;
  int iexit=0;
  int config,con1,con2;
  int ipsp,ncoreAE;
  double temp_eigen[10];
  double temp_norm[10];
  int npot[10];
  int ntpot;
  int ill[10];

  sprintf(filename, "%s.logd_plt", param->name);
  fp=fopen(filename,"w");
  fclose(fp);

  ncore = param->norb-param->nval;

  for (i=0; i<10; i++){
    npot[i]=0;
  }

  ntpot=0;

  for (i=0; i<param->nval; i++){
    npot[nlm_label(param->nlm[i+ncore]).l]++;
    if (npot[nlm_label(param->nlm[i+ncore]).l] == 1) {
      ntpot++;
    }
  }
  param->nll=ntpot;
  
  fp_log = fopen(logfile, "a");
  fprintf(fp_log,"<<<do_tc>>>\n");
  fclose(fp_log);

  /* set the log file */
  sprintf(filenames_.file_log, "%s", logfile);

  /* now loop over all the test configurations and perform AE and NL atom runs*/

  if (job==-67) {
    con1=0;
    con2=param->nconfigs;
  } else {
    con1=job-1;
    con2=job;
  }
  
  for (config=con1;config<con2;config++) {
    
    fp_log = fopen(logfile, "a");
    fprintf(fp_log,"\n\n ===============Configuration %d AE Calc=============== \n",config+1);
    fclose(fp_log);
    
    npm_.ncores = ncore;
    atomic_.norb = param->norb;
    ipsp = 0;   
    nlpot2_.inl = 0;  

    ilogder_.ilogder = 0;
    if (job != -67) ilogder_.ilogder = 1;
    
    logarith_.rphas = param->rphas;
    logarith_.elogmin = param->emin;
    logarith_.elogmax = param->emax;

    if (!strcmp(param->reltype, "nrl")){

      nrelorbae(param,config);      
      startae(param);

      for (i=0; i<param->norb; i++) {
	if (ibound_.ibd[i]==0) {
	  fp_log = fopen(logfile, "a");
	  fprintf(fp_log," !NOTE! State: |%3d> marked as unbound, using reference eigenvalue of %6.3f \n",
		  atomic_.nlm[i],ensave_.ensave[i]);
	  fclose(fp_log);
	}
      }

      scpot_(&param->z,&param->ixc,&param->exccut,&ipsp,&ifc, &iexit);
      if (iexit) {
	printf("Terminal error in: scpot <-- config #%d AE <-- do_tc\n EXITING OPIUM \n",config+1);
	exit(1);
      }
      ncoreAE=ncore;
      
    }else{
      
      if (job != -67) {
	fprintf(stderr," No Log derivs for scalar relativistic atoms yet :( \n");
	exit(1);
      }
      relorb(param,config);
      atm_(&param->z,&param->ixc,&iexit);
      if (iexit) {
	printf("Terminal error in: atm <-- config #%d AE <-- do_tc\n EXITING OPIUM \n",config+1);
	exit(1);
      }


      /*      ifrl = (!strcmp(param->reltype,"frl"))?1:0;
	      param->aereltransf = 0;
	      nrelorbae(param,config);
	      interp_(&ifrl, &param->aereltransf,&param->ixc);*/

      ncoreAE=0;
    }

    /* now dump the Log derivs from the AE calc */
    if (job != -67){
      for (i=0;i<4;i++)
	ill[i]=0;

      sprintf(filename, "%s.logd.%d", param->name,config+1);
      fp = fopen(filename, "wb");
      for (kk=0; kk<param->nll;kk++)
	for (k=0; k<param->nval; k++)
	  if ((ill[nlm_label(param->nlm[k+ncore]).l]==0) && (nlm_label(param->nlm[k+ncore]).l == kk)) {
	    ill[nlm_label(param->nlm[k+ncore]).l]++;
	    fwrite(logarith_.dlwf[k], sizeof(double), NPL0, fp);
	  }
      fclose(fp);
    
      sprintf(filename, "%s.logd_plt", param->name);
      fp = fopen(filename, "a");

      for (i=0;i<4;i++)
	ill[i]=0;
      
      dele = (param->emax - param->emin)/(NPL0-1);
      for (kk=0; kk<param->nll;kk++)
	for (k=0; k<param->nval; k++)
	  if ((ill[nlm_label(param->nlm[k+ncore]).l]==0) && (nlm_label(param->nlm[k+ncore]).l == kk)) {
	    ill[nlm_label(param->nlm[k+ncore]).l]++;
	    
	    for (j=0; j < NPL0; j++) {
	      e=param->emin + dele * j; 
	      fprintf(fp,"%lg %lg \n",e,logarith_.dlwf[k][j]);
	    }
	    fprintf(fp,"@ \n");
	  }
      fclose(fp);

      for (i=0;i<4;i++)
	ill[i]=0;
      
      sprintf(filename, "%s.logdeAE", param->name);
      fp = fopen(filename, "wb");
      for (kk=0; kk<param->nll;kk++)
	for (k=0; k<param->nval; k++)
	  if ((ill[nlm_label(param->nlm[k+ncore]).l]==0) && (nlm_label(param->nlm[k+ncore]).l == kk)) {
	    ill[nlm_label(param->nlm[k+ncore]).l]++;
	    /* Here we use ncoreAE */
	    fwrite(&atomic_.en[k+ncoreAE], sizeof(double), 1, fp);
	  }
      fclose(fp);
    }

    rp=write_reportae(param,rp,config,temp_eigen,temp_norm);
    enae[config+1] = results_.etot;

    /*    fp_log = fopen(logfile, "a");
    fprintf(fp_log,"\n\n ===============Configuration %d FC: Calc ===============\n",config+1);
    fclose(fp_log);*/

    /* FC */
    if (doifc != 0) {

      nlpot2_.inl = 0;  
      ipsp=0;
      
      ncore = param->norb-param->nval;
      
      sprintf(filename, "%s.psi_ae_core", param->name);
      fp = fopen(filename, "rb");
      for (i=0; i<ncore; i++) 
	fread(atomic_.rnl[i], sizeof(double), param->ngrid, fp);
      fclose(fp);
      
      sprintf(filename, "%s.eig_ae_core", param->name);
      fp = fopen(filename, "rb");
      
      for (i=0; i<ncore; i++) {
	fread(&atomic_.en[i], sizeof(double),1, fp);
	fread(&atomic_.wnl[i], sizeof(double),1, fp);
	fread(&atomic_.nlm[i], sizeof(int),1, fp);
	fread(&nmax_.nmax[i], sizeof(int),1, fp);
	fread(&nmax_.maxim, sizeof(int),1, fp);
	fread(&atomic_.xion, sizeof(double),1, fp);
	fread(&ibound_.ibd[i], sizeof(int),1, fp);
	fread(&ensave_.ensave[i], sizeof(double),1, fp);
      }
      fclose(fp);

      scpot_(&zeff,&param->ixc,&param->exccut,&ipsp,&ifc, &iexit); 
    rp = write_reportfc(param,rp,config,temp_eigen,temp_norm);
    enfc[config+1] = results_.etot;
    }

    /* NL */

    fp_log = fopen(logfile, "a");
    fprintf(fp_log,"\n\n ===============Configuration %d NL: Calc ===============\n",config+1);
    fclose(fp_log);

    zeff=atomic_.xion;
    for (i=0; i<param->nval; i++){
      if (ibound_.ibd[i]==0) {
	fp_log = fopen(logfile, "a");
	fprintf(fp_log," !NOTE! State: |%3d> marked as unbound, using reference eigenvalue of %6.3f \n",
		atomic_.nlm[i],ensave_.ensave[i]);
	fclose(fp_log);
      }
      zeff +=atomic_.wnl[i];
    }

    ilogder_.ilogder = 0;
    if (job != -67) ilogder_.ilogder = 1;
    logarith_.rphas = param->rphas;
    logarith_.elogmin = param->emin;
    logarith_.elogmax = param->emax;

    readPS(param);
    nrelorbnl(param,config);

    atomic_.norb = param->nval;
    npm_.ncores = 0; 
    nlpot2_.inl = 1;  
    ifc=0;
    ipsp=1;
    scpot_(&zeff,&param->ixc,&param->exccut,&ipsp,&ifc, &iexit); 

    sprintf(filename, "%s.psi_last", param->name);
    fp = fopen(filename, "wb");
    for (i=0; i<param->nval; i++){
      fwrite(&atomic_.nlm[i], sizeof(int), 1, fp);
      fwrite(&atomic_.wnl[i], sizeof(double), 1, fp);
      fwrite(&atomic_.en[i], sizeof(double), 1, fp);
      fwrite(atomic_.rnl[i], sizeof(double), param->ngrid, fp);
    }
    fclose(fp);

    if (iexit) {
      printf("Terminal error in: scpot <-- config #%d NL <-- do_tc\n EXITING OPIUM \n",config+1);
      exit(1);
    }

    rp = write_reportnl(param,rp,config,temp_eigen,temp_norm);
    ennl[config+1] = results_.etot;

  /* now dump the Log derivs from the NL calc */

    if (job != -67){
    for (i=0;i<4;i++)
      ill[i]=0;

      sprintf(filename, "%s.logd.%d", param->name,config+1);
      fp = fopen(filename, "ab");
      for (kk=0; kk<param->nll;kk++)
	for (k=0; k<param->nval; k++)
	  if ((ill[nlm_label(param->nlm[k+ncore]).l]==0) && (nlm_label(param->nlm[k+ncore]).l == kk)) {
	    ill[nlm_label(param->nlm[k+ncore]).l]++;
	    fwrite(logarith_.dlwf[k], sizeof(double), NPL0, fp);
	  }
      fclose(fp);
      
      sprintf(filename, "%s.logd_plt", param->name);
      fp = fopen(filename, "a");
      
      dele = (param->emax - param->emin)/(NPL0-1);

      for (i=0;i<4;i++)
	ill[i]=0;

      for (kk=0; kk<param->nll;kk++)
	for (k=0; k<param->nval; k++)
	  if ((ill[nlm_label(param->nlm[k+ncore]).l]==0) && (nlm_label(param->nlm[k+ncore]).l == kk)) {
	    ill[nlm_label(param->nlm[k+ncore]).l]++;
	    for (j=0; j < NPL0; j++) {
	      e=param->emin + dele * j; 
	      fprintf(fp,"%lg %lg \n",e,logarith_.dlwf[k][j]);
	    }
	    fprintf(fp,"@ \n");
	  }
      fclose(fp);

      for (i=0;i<4;i++)
	ill[i]=0;

      sprintf(filename, "%s.logdeNL", param->name);
      fp = fopen(filename, "wb");
      for (kk=0; kk<param->nll;kk++)
	for (k=0; k<param->nval; k++)
	  if ((ill[nlm_label(param->nlm[k+ncore]).l]==0) && (nlm_label(param->nlm[k+ncore]).l == kk)) {
	    ill[nlm_label(param->nlm[k+ncore]).l]++;
	    fwrite(&atomic_.en[k], sizeof(double), 1, fp);
	  }
      fclose(fp);
    }
  }

  rp+=sprintf(rp,
	      "\n  Comparison of total energy differences.           \n "
	      "  DD_ij = (E_i - E_j)_AE - (E_i-E_j)_NL     \n\n "
	      "AE-NL-   i   j         DD[mRy]        DD[meV] \n"
	      " AE-NL- ------------------------------------------\n");
  
  for(k=0;k<param->nconfigs+1;k++){
    for(kk=k+1;kk<param->nconfigs+1;kk++){
      sprintf(report+strlen(report), " AE-NL- %3d %3d   %14.6f %14.6f\n",
      	      k,kk,((enae[k] - enae[kk]) - (ennl[k] - ennl[kk]))*1000.,((enae[k] - enae[kk]) - (ennl[k] - ennl[kk]))*1000*13.6057);
    }    

  }    
  fp_log = fopen(logfile, "a");
  fprintf(fp_log,"   ------------------------------------------------\n");
  fclose(fp_log);


  fp_log = fopen(logfile, "a");
  fprintf(fp_log,"   ================================================\n");
  fclose(fp_log);

  return 0;
}

/* report section */

void do_tc_report(FILE *fp){
  fprintf(fp, "%s", report);
}

