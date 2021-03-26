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
 * $Id: do_ps.c,v 1.10 2004/10/02 18:34:49 ewalter Exp $
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "parameter.h"
#include "cdim.h"        /* fortran code parameters */
#include "do_ps.h"
#include "common_blocks.h"
#include "nlm.h"

/* fortran prototypes Should this be somewher else?????? */
void optim_(int *);
void kerker_(int *);
void tmsub_(int *);
/* report feature */

static char report[800];

void readAE(param_t *param);
char * write_reportps(param_t *param , char *rp);
void writePS(param_t *param);

int do_ps(param_t *param, char *logfile){

  int i;  
  FILE *fp_log;
  char *rp=report;
  char filename[80];

  fp_log = fopen(logfile, "a");

  fprintf(fp_log,"\n ======================================================================== \n");
  fprintf(fp_log," Begin PS construction\n");
  fprintf(fp_log," ======================================================================== \n");
  fclose(fp_log);  

  /* set the log file */
  sprintf(filenames_.file_log, "%s", logfile);

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
  } else if(param->psmeth == 't') {
    tmsub_(&param->ixc);
  } else {
    fp_log = fopen(logfile, "a");
    fprintf(fp_log, "   PS constuction method not known (must be optimized, tm, or kerker) \n");
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
  fread(totpot_.rvcore[0], sizeof(double), param->ngrid, fp);
  fread(totpot_.rvcoul, sizeof(double), param->ngrid, fp);
  fclose(fp);

  for (i=0; i<param->nval; i++) {
    if (i == param->ipot[ic]) {
      /*      fread(totpot_.rvcore[ic], sizeof(double), param->ngrid, fp);*/
      for (j=0; j<param->ngrid; j++) 
	totpot_.rvcore[ic][j]=totpot_.rvcore[0][j];
      ic++;

      /*    }else{
	    fseek(fp,sizeof(double)*param->ngrid,1);
	    }
	    fread(totpot_.rvcoul, sizeof(double), param->ngrid, fp);*/
    }
  }
    

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
      fread(&ibound_.ibd[ic], sizeof(int),1, fp);
      fread(&ensave_.ensave[ic], sizeof(double),1, fp);
      ic++;
    }else{
      fseek(fp,4*(sizeof(double)+sizeof(int)),1);
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

  if (param->psmeth=='o') {
    
    rp+=sprintf(rp,
		"    ====================Optimized pseudopotential method====================\n\n");   
    rp+=sprintf(rp, 
		"                       Pseudopotential convergence error                      \n");
    rp+=sprintf(rp, 
		"    Orbital      [mRy/e]       [meV/e]         [mRy]        [meV]        Ghost\n"
		"    --------------------------------------------------------------------------\n");
    for (i=0; i<param->nll; i++){
      if (ibound_.ibd[i]==1) {
	rp+=sprintf(rp, "\t%3d  %12.6f  %12.6f  %12.6f  %12.6f\t%6s\n",
		    atomic_.nlm[i], convrpt_.sumc[i]*1000.,convrpt_.sumc[i]*1000*13.6057,
                    convrpt_.wsumt[i]*1000., 
		    convrpt_.wsumt[i]*13.6057*1000.,
		    (convrpt_.lghost[i]>0)?"yes":((convrpt_.lghost[i]<0)?"?":"no"));
	tot_conv_error += convrpt_.wsumt[i]*1000.;
      } else {
	rp+=sprintf(rp, "\t%3d  (unbound)   --------          --------          --------   %6s\n",
		    atomic_.nlm[i],(convrpt_.lghost[i]>0)?"yes":((convrpt_.lghost[i]<0)?"?":"no"));
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
      rp+=sprintf(rp, "\t%3d\t%6s\n",
		  atomic_.nlm[i],(convrpt_.lghost[i]==0)?"no":"yes");
    }
  }    
    rp+=sprintf(rp, "\n");
    igh=0;
    for (i=0; i<param->nll; i++){
      if (convrpt_.lghost[i]==0) {
	/*rp+=sprintf(rp, "\t%3d ok as the local potential \n",atomic_.nlm[i]);*/
      } else {
	rp+=sprintf(rp, "\t%3d SHOULD NOT be used as the local potential \n",atomic_.nlm[i]);
	igh+=1;
      }
    }
    if (igh==param->nll) {
      rp+=sprintf(rp, "\t !!ERROR!! There are no choices for local potential\n");
    }

  return rp;
}

void writePS(param_t *param) {

  int i,j,ncore; 
  int iset=0;
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
