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

static char report[16000];
void startae(param_t *param, int);
void relorbae(param_t *param, int, char *); 
void nrelorbnl(param_t *param, int, char *); 
void nrelorbae(param_t *param, int, char *); 
char * write_reportae(param_t *param, char *rp,int, double temp_eigen[], double temp_norm[],int);
char * write_reportnl(param_t *param, char *rp,int, double temp_eigen[], double temp_norm[], int, int);
char * write_reportfc(param_t *param, char *rp,int, double temp_eigen[], double temp_norm[]);
void hfsolve_(double  *, int * ,double *, int * , int *, int *, int *, int *);
void dfsolve_(double  *, int * ,double *, int * , int *, int *, int *, int *);
void readPS(param_t *param);
void readAE(param_t *param);
void dftsolve_(double  *, int * ,double *, int *, int *, int *, int *, int *, int *);

#define streq(a,b) (!strcasecmp(a,b))

int do_tc(param_t *param, char *logfile, int job, int doifc, int donl){

  int i, j,k, kk; 
  FILE *fp_log;
  FILE *fp;
  static char filename[160];
  char *rp=report;
  int ncore,ifrl; 
  double zeff;
  int irel=0;
  int irelxc=0;
  int iminne=0;
  double e,dele,minne;
  int ifc=0;
  int iexit=0;
  int iprint=1;
  int config,con1,con2;
  int ipsp,ncoreAE;
  double temp_eigen[10];
  double temp_norm[10];
  double exccut_temp;
  int ntpot;
  int ill[10];
  int insl=1;
  char box[80];
  double boxval;

  sprintf(filename, "%s.logd_plt", param->name);
  fp=fopen(filename,"w");
  fclose(fp);

  if (param->exccut < 0) {
    exccut_temp=0.0;
  } else {
    exccut_temp=param->exccut;
  }

  /*  fp = fopen("boxfile", "r");
  while(fgets(box, 80, fp) != NULL)
    sscanf (box, "%lg", &boxval);
  fclose(fp);
  printf("boxval %lg \n",boxval);
  */
  ncore = param->norb-param->nval;

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
    
    aorb_.ncore = ncore;
    aorb_.norb = param->norb;
    ipsp = 0;   
    nlpot2_.inl = 0;  

    ilogder_.ilogder = 0;
    if (job != -67) {
      if (param->ixc<0) {
	fprintf(stderr," No Log derivs for HF/DF yet :( \n");
	exit(1);
      }
      ilogder_.ilogder = 1;
    }
    
    logarith_.rphas = param->rphas;
    logarith_.elogmin = param->emin;
    logarith_.elogmax = param->emax;

    if (!strcmp(param->reltype, "nrl")){

      nrelorbae(param,config,logfile);      
      startae(param,aorb_.norb);
      
      irel=0;
      iprint=1;
      irelxc=0;
      if (param->ixc < 0){ 
	hfsolve_(&param->z,&param->ixc,&exccut_temp,&ipsp,&ifc,&iexit,&irel,&iprint);
      } else {

	/*        for (i=0; i<param->ngrid+1; i++)
          for (j=0; j<aorb_.norb; j++)
	    if (aorb_.lo[j]  == 3)
	      totpot_.rvcore[j][i] += boxval*grid_.r[i];
	
        printf("box %lg \n", boxval);
	*/
	dftsolve_(&param->z,&param->ixc,&exccut_temp,&ipsp,&ifc,&iexit,&irel,&irelxc,&iprint);

	/*
        for (i=0; i<param->ngrid+1; i++)
          for (j=0; j<aorb_.norb; j++)
	    if (aorb_.lo[j]  == 3)
            totpot_.rvcore[j][i] -= boxval*grid_.r[i];
        printf("box %lg \n", boxval);
	*/
      }
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

      irel=1;
      if (streq(param->relxc,"rxc")){ 
	irelxc=1;
      }else{
	irelxc=0;
      }
      relorbae(param,config,logfile);
      startae(param,aorb_.norb);

      if (param->ixc < 0) {
	dfsolve_(&param->z,&param->ixc,&exccut_temp,&ipsp,&ifc,&iexit,&irel,&iprint);
      }else {
	/*
        for (i=0; i<param->ngrid+1; i++)
          for (j=0; j<aorb_.norb; j++)
	    if (aorb_.lo[j]  == 3)
            totpot_.rvcore[j][i] += boxval*grid_.r[i];
	*/
	dftsolve_(&param->z,&param->ixc,&exccut_temp,&ipsp,&ifc,&iexit,&irel,&irelxc,&iprint);
	
	/*
        for (i=0; i<param->ngrid+1; i++)
          for (j=0; j<aorb_.norb; j++)
	    if (aorb_.lo[j]  == 3)
            totpot_.rvcore[j][i] -= boxval*grid_.r[i];
	
	*/
      }

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

      sprintf(filename, "%s.logd.%d", param->name,config+1);
      fp = fopen(filename, "wb");
      for (k=0; k<param->nll; k++)
	fwrite(logarith_.dlwf[k], sizeof(double), NPL0, fp);
      fclose(fp);
    
      sprintf(filename, "%s.logd_plt", param->name);
      fp = fopen(filename, "a");

      dele = (param->emax - param->emin)/(NPL0-1);
      for (k=0; k<param->nll; k++) {
	for (j=0; j < NPL0; j++) {
	  e=param->emin + dele * j; 
	  fprintf(fp,"%lg %lg \n",e,logarith_.dlwf[k][j]);
	}
	fprintf(fp,"@ \n");
      }
      fclose(fp);

      sprintf(filename, "%s.logdeAE", param->name);
      fp = fopen(filename, "wb");
      for (k=0; k<param->nval; k++) {
	minne=1000000.0;
	for (j=0; j < NPL0; j++) {
	  e=param->emin + dele * j; 
	  if ((fabs(e-adat_.en[k+ncoreAE])) < fabs(minne-adat_.en[k+ncoreAE])) {
	    minne=e;
	    iminne=j;
	  }
	}
	fwrite(&adat_.en[k+ncoreAE], sizeof(double), 1, fp);
	fwrite(&logarith_.dlwf[k][iminne], sizeof(double), 1, fp);
      }
      fclose(fp);    
    }

    rp=write_reportae(param,rp,config,temp_eigen,temp_norm,aorb_.norb);
    enae[config+1] = aval_.etot;

    /* NL or SL */

    if (param->ixc <0) {
      if (donl != 0) {
	printf("!!NOTE!!: All Hartree-Fock tests are done with the semi-local form of the potential \n");
	printf("!!NOTE!!: Use the 'sl' command instead of the 'nl' command to remove this warning...\n");

	fp_log = fopen(logfile, "a");
	
	fprintf(fp_log,"!!NOTE!!: All Hartree-Fock tests are done with the semi-local form of the potential \n");
	fprintf(fp_log,"!!NOTE!!: Use the 'sl' command instead of the 'nl' command to remove this warning...\n");
	
	fclose(fp_log);
	donl=0;	
      }
    }

    nlpot2_.inl = donl;  

    fp_log = fopen(logfile, "a");
    if (donl==0) fprintf(fp_log,"\n\n ===============Configuration %d SL: Calc ===============\n",config+1);
    if (donl==1) fprintf(fp_log,"\n\n ===============Configuration %d NL: Calc ===============\n",config+1);
    fclose(fp_log);

    ilogder_.ilogder = 0;
    if (job != -67) ilogder_.ilogder = 1;
    logarith_.rphas = param->rphas;
    logarith_.elogmin = param->emin;
    logarith_.elogmax = param->emax;


    if (!strcmp(param->reltype, "frl")) { 
      relorbnl(param,config,logfile);
    }else{
      nrelorbnl(param,config,logfile);
    }

    readPS(param);


    zeff=adat_.xion;
    for (i=0; i<param->nval; i++){
      zeff +=adat_.wnl[i];
    }

    irelxc=0;
    ifc=0;
    ipsp=1;
    irel=0;
    if (param->ixc < 0) {
      hfsolve_(&zeff,&param->ixc,&exccut_temp,&ipsp,&ifc,&iexit,&irel,&iprint);
    }else {
      /*
        for (i=0; i<param->ngrid+1; i++)
          for (j=0; j<aorb_.norb; j++)
	    if (aorb_.lo[j]  == 3)
            totpot_.rvcore[j][i] += boxval*grid_.r[i];
      */
	dftsolve_(&param->z,&param->ixc,&exccut_temp,&ipsp,&ifc,&iexit,&irel,&irelxc,&iprint);
	/*
        for (i=0; i<param->ngrid+1; i++)
          for (j=0; j<aorb_.norb; j++)
	    if (aorb_.lo[j]  == 3)
            totpot_.rvcore[j][i] -= boxval*grid_.r[i];
	*/
    }

    sprintf(filename, "%s.psi_last", param->name);
    fp = fopen(filename, "wb");
    for (i=0; i<param->nval; i++){
      fwrite(&aorb_.nlm[i], sizeof(int), 1, fp);
      fwrite(&adat_.wnl[i], sizeof(double), 1, fp);
      fwrite(&adat_.en[i], sizeof(double), 1, fp);
      fwrite(wfn_.rnl[i], sizeof(double), param->ngrid, fp);
    }
    fclose(fp);

    if (iexit) {
      if (donl==0) printf("Terminal error in: scpot <-- config #%d SL <-- do_tc\n EXITING OPIUM \n",config+1);
      if (donl==1) printf("Terminal error in: scpot <-- config #%d NL <-- do_tc\n EXITING OPIUM \n",config+1);
      exit(1);
    }


    rp = write_reportnl(param,rp,config,temp_eigen,temp_norm,donl,aorb_.nval);
    ennl[config+1] = aval_.etot;

  /* now dump the Log derivs from the NL calc */

    if (job != -67){
      sprintf(filename, "%s.logd.%d", param->name,config+1);
      fp = fopen(filename, "ab");
      for (k=0; k<param->nll; k++)
	fwrite(logarith_.dlwf[k], sizeof(double), NPL0, fp);
      fclose(fp);
      
      sprintf(filename, "%s.logd_plt", param->name);
      fp = fopen(filename, "a");
      
      dele = (param->emax - param->emin)/(NPL0-1);

      for (k=0; k<param->nll; k++){
	for (j=0; j < NPL0; j++) {
	  e=param->emin + dele * j; 
	  fprintf(fp,"%lg %lg \n",e,logarith_.dlwf[k][j]);
	}
	fprintf(fp,"@ \n");
      }
      fclose(fp);
      
      sprintf(filename, "%s.logdeNL", param->name);
      fp = fopen(filename, "wb");
      for (k=0; k<param->nval; k++) {
	minne=1000000.0;
	for (j=0; j < NPL0; j++) {
	  e=param->emin + dele * j; 
	  if ((fabs(e-adat_.en[k])) < fabs(minne-adat_.en[k])) {
	    minne=e;
	    iminne=j;
	  }
	}
	fwrite(&adat_.en[k], sizeof(double), 1, fp);
	fwrite(&logarith_.dlwf[k][iminne], sizeof(double), 1, fp);
      }
      fclose(fp);    
    }

  }
  
  if (donl==0) {
    rp+=sprintf(rp,
		"\n  Comparison of total energy differences.           \n "
		"  DD_ij = (E_i - E_j)_AE - (E_i-E_j)_SL     \n\n "
		"AE-SL-   i   j         DD[mRy]        DD[meV] \n"
		" AE-SL- ------------------------------------------\n");
  }else{
    rp+=sprintf(rp,
		"\n  Comparison of total energy differences.           \n "
		"  DD_ij = (E_i - E_j)_AE - (E_i-E_j)_NL     \n\n "
		"AE-NL-   i   j         DD[mRy]        DD[meV] \n"
		" AE-NL- ------------------------------------------\n");
  }
  
  for(k=0;k<param->nconfigs+1;k++){
    for(kk=k+1;kk<param->nconfigs+1;kk++){
      if (donl==0) {
	sprintf(report+strlen(report), " AE-SL- %3d %3d   %14.6f %14.6f\n",
		k,kk,((enae[k] - enae[kk]) - (ennl[k] - ennl[kk]))*1000.,((enae[k] - enae[kk]) - (ennl[k] - ennl[kk]))*1000*13.6057);
      }else{
	sprintf(report+strlen(report), " AE-NL- %3d %3d   %14.6f %14.6f\n",
		k,kk,((enae[k] - enae[kk]) - (ennl[k] - ennl[kk]))*1000.,((enae[k] - enae[kk]) - (ennl[k] - ennl[kk]))*1000*13.6057);
      }
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

