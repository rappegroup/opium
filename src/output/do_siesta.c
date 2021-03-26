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

/****************************************************************************
 * generate *.siesta output (Siesta compatible) *
 *****************************************************************************
 * INPUT:  param_t structure, *.ps, *.loc
 * OUTPUT: *.siesta
 ****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#ifndef M_PI 
#define M_PI 3.14159265358979323846
#endif

#define MAX(a, b)   (((a) > (b)) ? (a):(b))

#include "parameter.h"        /* defines structure: 'param_t' */
#include "cdim.h"        /* fortran code parameters */
#include "do_siesta.h"          /* the module's own header */
#include "nlm.h"
#include "common_blocks.h"

void interp2_(int *, double *,double *, double *);
void writeparam(param_t *param, FILE *fp, FILE *fp_param);
void nrelsproj(param_t *param, char *);
int do_siesta(param_t *param, FILE *fp_param, char *logfile){

  int i,j,k,ic,icount,ncore,ii;
  char filename[180];
  int nval2=0;
  FILE *fp_file;
  FILE *fp_out;
  FILE *fp_log;
  double junk=0.0;
  double zeff;	                /* effective Z */
  int lmax;                     /* max l value */   
  int lloc;
  int ngg=nrgrid_.nr-10;
/* fix */
  int kk;
  int c;
  char xao[5];
  char xstr[3];
  char relstr[4];
  char pccstr[5];
  
  static double dumm[NPDM];
  static double dumm2[NPDM];

  static double rho[NPDM];
  static double rho2[NPDM];

  static double rscore[NPDM];
  static double rscore2[NPDM];

  static double rvv[N0][NPDM];
  static double rvv2[N0][NPDM];

  static double rvloc[NPDM]; 
  static double rvloc2[NPDM];  

  xao[0]='s';
  xao[1]='p';
  xao[2]='d';
  xao[3]='f';
  xao[4]='g';

  ncore=param->norb-param->nval;

  nrelsproj(param,logfile);

  fp_log = fopen(logfile, "a");
  fprintf(fp_log,"<<<do_siesta -- psf>>>\n");  
  fclose(fp_log);

  if (param->nboxes > 0) {
    fp_log = fopen(logfile, "a");
    fprintf(fp_log,"!!ERROR!!: The Siesta format does not support the use of augmentation operators \n");
    fclose(fp_log);
    return 1;
  }

  fp_log = fopen(logfile, "a");
  fprintf(fp_log,"!!NOTE!!: the Siesta code redefines the local part of the psp, proceed with caution ...\n");
  /*printf(" Note: the Siesta code redefines the local part of the psp, proceed with caution ...\n");*/
  fclose(fp_log);

  zeff = param->z;
  for (i=0; i<param->norb - param->nval; i++)
    zeff -= param->wnl[i];

  /*  sprintf(filename, "%s.loc", param->name);
  fp_file = fopen(filename, "rb");
  fread(rvloc, sizeof(double), param->ngrid, fp_file);
  fclose(fp_file);

  interp2_(rvloc,rvloc2);*/
    
  for (i=0; i<param->nll; i++) {
    sprintf(filename, "%s.pot.ps.l=%d", param->name,i);
    
    fp_file = fopen(filename, "rb");
    fread(dumm, sizeof(double), param->ngrid, fp_file);    

    interp2_(&nrgrid_.nr,dumm,dumm2,rgrid_.r);

    for (j=0; j<ngg; j++){
      rvv[i][j]=dumm2[j];
    }

    fclose(fp_file);
  }

  icount=0; 
  for (k=0; k<param->nll;k++){
    for (i=0; i<ngg; i++){
      rvv2[k][i] = rvv[k][i];
    }
    icount++;
  }

  lmax=param->nll-1;
  lloc=nlm_label(param->nlm[param->localind+ncore]).l;
  
  
  if (param->rpcc > 0.){
      sprintf(filename, "%s.rho_pcore", param->name);
      fp_file = fopen(filename, "rb");
      fread(rscore, sizeof(double), param->ngrid, fp_file);
      interp2_(&nrgrid_.nr,rscore,rscore2,rgrid_.r);
      fclose(fp_file);
  }

  sprintf(filename, "%s.rho_nl", param->name);
  if (fp_file = fopen(filename, "rb")) {
    fread(rho, sizeof(double), param->ngrid, fp_file);
    fclose(fp_file);
  } else {
    sprintf(filename, "%s.rho_sl", param->name);
    if (fp_file = fopen(filename, "rb")) {
      fread(rho, sizeof(double), param->ngrid, fp_file);
      fclose(fp_file);
    } else {
      fp_log = fopen(logfile, "a");
      fprintf(fp_log,"No valence charge density avaliable :( --EXIT!\n");
      printf("No valence charge density avaliable :( --EXIT!\n");
      fclose(fp_log);
      exit(1);
    }
  }

  interp2_(&nrgrid_.nr,rho,rho2,rgrid_.r);

  /* section 1 : header */

  sprintf(filename, "%s.psf", param->name);
  fp_out = fopen(filename, "w");


  /* 1st line of siesta format: The atom symbol, icorr, irel, nicore */
 
  switch (param->ixc) {
    
  case 0 : sprintf(xstr,"%s","ca"); break;
  case 1 : sprintf(xstr,"%s","pw"); break;
  case 2 : sprintf(xstr,"%s","pb"); break;

  default :
    fp_log = fopen(logfile, "a");
    fprintf(fp_log,"!!WARNING!!: Can not determine XC string for Siesta format, using ca \n");
    fclose(fp_log);
    printf("!!WARNING!!: Can not determine XC string for Siesta format, using ca \n");

    sprintf(xstr,"%s","ca");
  }

  sprintf(relstr,"%s","nrl");

  if (param->rpcc>0) {
    sprintf(pccstr,"%s","pcec");
  } else {
    sprintf(pccstr,"%s","nc  ");
  }

  fprintf(fp_out,"%2s %3s %3s %4s \n",param->symbol,xstr,relstr,pccstr);
  
  /* 2nd line of siesta format: method : 6 strings each 10chars */
  fprintf(fp_out," OPIUM generated potential, local l is: %d \n",lloc); 

  /* 3rd line of siesta format: text : 70 chars */
  fprintf(fp_out," ");
  for (i=0; i<param->nval;i++){
    if (param->npot[i]==0) {
      ii=i+ncore;
      fprintf(fp_out,"%d%c%5.2f  r=%5.2f/",nlm_label(param->nlm[ii]).n,xao[nlm_label(param->nlm[ii]).l],MAX(param->wnl[ii],0),param->rc[i]);
    }
  }
  fprintf(fp_out,"\n");

  /* 4th line of siesta format: text : Ndown Nup ngrid a b zval */
  fprintf(fp_out," %3d%3d%5d%20.12E%20.12E%20.12E\n",lmax+1,nval2,ngg-1,rgrid_.a,rgrid_.b,zeff); 
  
  /* section 2: radial grid */
  fprintf(fp_out," Radial grid follows   \n"); 

  for (j=1; j<ngg; j=j+4) {
    fprintf(fp_out, "%#20.12E%#20.12E%#20.12E%#20.12E\n",rgrid_.r[j],rgrid_.r[j+1],rgrid_.r[j+2],rgrid_.r[j+3]);
  }

  /* section 3: Potentials */

  for (i=0; i<lmax+1; i++){
    fprintf(fp_out," Down Pseudopotential follows (l on next line) \n"); 
    fprintf(fp_out,"  %d \n",i);
      for (j=1; j<ngg; j=j+4) {
	  fprintf(fp_out, "%#20.12E%#20.12E%#20.12E%#20.12E\n",rvv2[i][j],rvv2[i][j+1],rvv2[i][j+2],rvv2[i][j+3]);
      }
  }

  /* section 4: Core Charge */

  fprintf(fp_out,"Core charge follows \n");
  if (param->rpcc > 0.) {
    for (j=1; j<ngg; j=j+4) {
      fprintf(fp_out, "%#20.12E%#20.12E%#20.12E%#20.12E\n",rscore2[j],rscore2[j+1],rscore2[j+2],rscore2[j+3]);
    } 
  } else {
    for (j=1; j<ngg; j=j+4) {
      fprintf(fp_out, "%#20.12E%#20.12E%#20.12E%#20.12E\n",junk,junk,junk,junk);
    }	
  }

  /* section 5: Valence Charge */
  fprintf(fp_out,"Valence charge follows \n");
  for (j=1; j<ngg; j=j+4) {
    fprintf(fp_out, "%#20.12E%#20.12E%#20.12E%#20.12E\n",rho2[j],rho2[j+1],rho2[j+2],rho2[j+3]);
  }

  writeparam(param, fp_out, fp_param);

  fp_log = fopen(logfile, "a");  
  fprintf(fp_log,"   ================================================\n");
  fclose(fp_log);
  
  return 0;
}

