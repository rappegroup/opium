/* 
 * $Id: do_tc.c,v 1.5 2004/06/16 20:46:17 mbarnes Exp $
 */

/* standard libraries */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "parameter.h"        /* defines structure: 'param_t' */
#include "nlm.h"              /* nlm_label call */
#include "fortparam.h"        /* fortran code parameters */
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
void readPS(param_t *param);

/* fortran prototypes Should this be somewher else?????? */
void scpot_(double  *, int * ,int *);
void atm_(double *, int *);


int do_tc(param_t *param, char *logfile, int job){

  int i, k, kk; 
  FILE *fp_log;
  FILE *fp;
  static char filename[160];
  char *rp=report;
  int ncore; 
  double zeff;
  double e,dele;
  int config,con1,con2;
  int ipsp;
  double temp_eigen[10];
  double temp_norm[10];
  int npot[10];
  int ntpot;

  sprintf(filename, "%s.plt_logd", param->name);
  fp=fopen(filename,"w");
  fclose(fp);

  for (i=0; i<10; i++){
    npot[i]=0;
  }

  ntpot=0;

  for (i=0; i<param->nval; i++){
    npot[nlm_label(atomic_.nlm[i]).l]++;
    if (npot[nlm_label(atomic_.nlm[i]).l] == 1) {
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

    ncore = param->norb-param->nval;
    npm_.ncores = param->norb-param->nval; 
    atomic_.norb = param->norb;
    ipsp = 0;   
    nlpot2_.inl = 0;  

    if (!strcmp(param->reltype, "nrl")){

      nrelorbae(param,config);      
      startae(param);

      ilogder_.ilogder = 0;
      if (param->ilogder != -67) ilogder_.ilogder = 1;
      
      logarith_.rphas = param->rphas;
      logarith_.elogmin = param->emin;
      logarith_.elogmax = param->emax;

      scpot_(&param->z,&param->ixc,&ipsp);

      /* now dump the Log derivs from the AE calc */
      if (param->ilogder != -67){
	sprintf(filename, "%s.logd.%d", param->name,config);
	fp = fopen(filename, "wb");
	for (i=0; i<param->nll; i++)
	  fwrite(logarith_.dlwf[i], sizeof(double), NPL0, fp);
	fclose(fp);

	sprintf(filename, "%s.plt_logd", param->name);
	fp = fopen(filename, "a");
	
	dele = (param->emax - param->emin)/(NPL0-1);
	
	for (i=0; i<param->nll;i++) {
	  for (k=0; k < NPL0; k++) {
	    e=param->emin + dele * k; 
	    fprintf(fp,"%lg %lg \n",e,logarith_.dlwf[i][k]);
	  }
	  fprintf(fp,"@ \n");
	}
	fclose(fp);
      }
    }else{

      relorb(param,config);
      atm_(&param->z,&param->ixc);
      
    }

  fp_log = fopen(logfile, "a");
  fprintf(fp_log,"\n\n ===============Configuration %d NL: Calc ===============\n",config+1);
  fclose(fp_log);


    rp=write_reportae(param,rp,config,temp_eigen,temp_norm);
    enae[config+1] = results_.etot;

    /* NL */

    atomic_.norb = param->nval;
    npm_.ncores = 0; 
    ipsp = 1;   
    nlpot2_.inl = 1;  

    readPS(param);
    nrelorbnl(param,config);

    zeff=atomic_.xion;
    for (i=0; i<param->nval; i++)
      zeff +=atomic_.wnl[i];

    ilogder_.ilogder = 0;
    if (param->ilogder != -67) ilogder_.ilogder = 1;
    logarith_.rphas = param->rphas;
    logarith_.elogmin = param->emin;
    logarith_.elogmax = param->emax;
    
    scpot_(&zeff,&param->ixc,&ipsp); 

    rp = write_reportnl(param,rp,config,temp_eigen,temp_norm);
    ennl[config+1] = results_.etot;

  /* now dump the Log derivs from the NL calc */
    if (param->ilogder != -67){
      sprintf(filename, "%s.logd.%d", param->name,config);
      fp = fopen(filename, "ab");
      for (i=0; i<param->nll; i++)
	fwrite(logarith_.dlwf[i], sizeof(double), NPL0, fp);
      fclose(fp);

      sprintf(filename, "%s.plt_logd", param->name);
      fp = fopen(filename, "a");
      
      dele = (param->emax - param->emin)/(NPL0-1);
      
      for (i=0; i<param->nll;i++) {
	for (k=0; k < NPL0; k++) {
	  e=param->emin + dele * k; 
	  fprintf(fp,"%lg %lg \n",e,logarith_.dlwf[i][k]);
	}
	fprintf(fp,"@ \n");
      }
      fclose(fp);
    }

  }

  rp+=sprintf(rp,
	      "\n  Comparison of total energy differences.           \n "
	      "  DD_ij = (E_i - E_j)_all-electron - (E_i - E_j)_pseudo     \n\n "
	      "AE-NL-  i   j          DD[mRy]        DD[meV] \n"
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


