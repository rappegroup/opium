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
/*                                                                          */
/*generate *.pwf output                                                     */
/****************************************************************************/
/*INPUT:  param_t structure, *.ps, *.loc                                    */
/*OUTPUT: *.pwf                                                             */
/****************************************************************************/

/* standard libraries */
#include <stdio.h>
#include <stdlib.h>

#include "parameter.h"        /* defines structure: 'param_t' */
#include "cdim.h"        /* fortran code parameters */
#include "do_pwf.h"           /* the module's own header */
#include "common_blocks.h"    /* fortran common blocks */

/* from atom/do_nl.c */
void nrelorbnl(param_t *param, int, char *);
void nrelsproj(param_t *param, char *);
/* creative passing of the 1st argument to follow */
void writepwf_(void *, char *);
void writeparam(param_t *param, FILE *fp, FILE *fp_param);

int do_pwf(param_t *param, FILE *fp_param, char *logfile){

  int i,j,ic;              /* loop counter */
  char filename[180]; /* filename */
  int c,config;              /* dummy character */
  FILE *fp;           /* file pointer */
  FILE *fp2;           /* file pointer */
  FILE *fp_log; 

  fp_log = fopen(logfile, "a");
  fprintf(fp_log,"<<<do_pwf>>>\n");
  fclose(fp_log);
  
  /* set the log file */
  sprintf(filenames_.file_log, "%s", logfile);

  nrelsproj(param,logfile);
  /* new routine to set the arrays for semicore states */

  /* set the pwf common block */
  /*  sprintf(filenames_.file_ps, "%s.ps", param->name);*/

  /* This section was put here to support the old pseudo code 
     It should be removed !!! - EJW */
  pwf_.zeff = param->z;
  for (i=0; i<param->norb-param->nval; i++)
    pwf_.zeff -= param->wnl[i];
  pwf_.rpcc = param->rpcc;  
  pwf_.igrid = param->ngrid;
  pwf_.spacing = param->b;
  pwf_.first = param->a;
  pwf_.nval = param->nll;
  if (param->nboxes > 0)
    pwf_.inl = 1;
  else
    pwf_.inl = 0;

  if (param->nll == 1)
    param->ist1 = 3;
  else
    param->ist1 = 2;

  pwf_.ist1 = param->ist1;
  pwf_.ilocal = param->local;
  pwf_.ilocalind = param->localind;

  param->ncut = param->ngrid;
  pwf_.ncut = grid_.np;

  aorb_.nval=param->nval;
  aorb_.ncore=param->norb-param->nval;
  psdat_.nll=param->nll;
  
  /* END bad code */
  
  config=-1;
  nrelorbnl(param,config,logfile);

  sprintf(filename, "%s.loc", param->name);
  if (fp = fopen(filename, "rb")) {
    fread(nlcore_.rvloc, sizeof(double), param->ngrid, fp);
    fclose(fp);
  } else {
    fp_log = fopen(logfile, "a");
    fprintf(fp_log,"Looks like you never ran nl yet you have augmentation functions :( --EXIT!\n");
    printf("Looks like you never ran nl yet you have augmentation functions :( --EXIT!\n");
    fclose(fp_log);
    exit(1);
    }

  for (i=0; i<param->nll; i++) {
    sprintf(filename, "%s.pot.ps.l=%d", param->name,i);
    fp = fopen(filename, "rb");
    fread(totpot_.rvcore[i], sizeof(double), param->ngrid, fp);
    fseek(fp,sizeof(double) ,param->ngrid);
    fclose(fp);
  }

  for (i=0; i<param->nll; i++) {
    sprintf(filename, "%s.psi.ps.l=%d", param->name,i);
    fp = fopen(filename, "rb");
    fread(wfn_.rnl[i], sizeof(double), param->ngrid, fp);
    fclose(fp);
    }

  if (param->rpcc > 0.){
    sprintf(filename, "%s.rho_pcore", param->name);
    fp = fopen(filename, "rb");
    fread(rscore_.rscore, sizeof(double), param->ngrid, fp);
    fclose(fp);
    }

  fp_log = fopen(logfile, "a");
  fprintf(fp_log,"<<calling: writepwf>>\n");
  fclose(fp_log);

  sprintf(filename, "%s.pwf", param->name);
  fp = fopen(filename, "w");

  writepwf_(fp, param->symbol);

  sprintf(filename, "%s.psi_last", param->name);
  fp2 = fopen(filename, "rb");
  for (i=0; i<param->nval; i++){
    fread(&aorb_.nlm[i], sizeof(int), 1, fp2);
    fread(&adat_.wnl[i], sizeof(double), 1, fp2);
    fread(&adat_.en[i], sizeof(double), 1, fp2);
    fread(wfn_.rnl[i], sizeof(double), param->ngrid, fp2);
  }
  fclose(fp2);

  for (i=0; i<param->ngrid; i++) {
    fprintf(fp,"%20.15f",grid_.r[i]);
    for (j=0; j<param->nval; j++){    
      fprintf(fp,"%20.15f",wfn_.rnl[j][i]/grid_.r[i]);
    }    
    fprintf(fp,"\n");
  }

  if (param->rpcc > 0.){
    for (i=0; i<param->ngrid; i++) {
      fprintf(fp,"%26.16E %26.16E \n",grid_.r[i],rscore_.rscore[i]);
    }
  }

  fp_log = fopen(logfile, "a");

  fprintf(fp_log,"\n The following configuration was used for the wavefunction output:\n");  

  for (j=0; j<param->nval; j++){    
    fprintf(fp,"%5d",aorb_.nlm[j]);
    /* fprintf(fp,"%15.8f",adat_.wnl[j]);*/
    fprintf(fp,"%15.8f\n",adat_.en[j]*13.6057);
    fprintf(fp_log,"%5d",aorb_.nlm[j]);
    fprintf(fp_log,"%15.8f",adat_.wnl[j]);
    fprintf(fp_log,"%15.8f Ry        %15.8f eV\n",adat_.en[j],adat_.en[j]*13.6057);
  }    

  fclose(fp_log);  


  writeparam(param, fp, fp_param);

  fp_log = fopen(logfile, "a");
  fprintf(fp_log, "   ================================================\n");
  fprintf(fp_log,"<<end:  do_pwf>>\n");
  fclose(fp_log);

  return 0;


}


void writeparam(param_t *param, FILE *fp, FILE *fp_param){

  int c;

  /* append the parameter file to the end of the pwf file */
  fprintf(fp, "\n");
  fprintf(fp, "############################################################\n");
  fprintf(fp, "#    Opium Parameter File                                  #\n");
  fprintf(fp, "############################################################\n");
  fprintf(fp, "#-----------------------------------------------------------\n");
  fprintf(fp, "#Command line   : %s \n",param->execstring);
  fprintf(fp, "#Opium version  : %s \n",param->version);
  fprintf(fp, "#Compile Host   : %s \n",param->chost);
  fprintf(fp, "#Compile OS     : %s \n",param->csys);
  fprintf(fp, "#Compile Date   : %s \n",param->cdate);
  fprintf(fp, "#-----------------------------------------------------------\n");
  fprintf(fp, "#Execution Host : %s \n",param->ehost);
  fprintf(fp, "#Execution Date : %s \n",param->edate);
  fprintf(fp, "#-----------------------------------------------------------\n");  
  fprintf(fp, "############################################################\n");
  rewind(fp_param);
  while((c=fgetc(fp_param)) != EOF) fputc(c, fp);
  fprintf(fp, "############################################################\n");
  fclose(fp);
  
}
