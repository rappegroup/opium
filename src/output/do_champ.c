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

/****************************************************************************
 * generate *.champ output (Champ compatible) *
 *****************************************************************************
 * INPUT:  param_t structure, *.ps, *.loc
 * OUTPUT: *.champ
 ****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#ifndef M_PI 
#define M_PI 3.14159265358979323846
#endif

#include "parameter.h"        /* defines structure: 'param_t' */
#include "cdim.h"        /* fortran code parameters */
#include "do_champ.h"          /* the module's own header */
#include "nlm.h"
#include "common_blocks.h"

void interp2_(double *,double *);
void writeparam(param_t *param, FILE *fp, FILE *fp_param);
void readPS(param_t *param);
int do_champ(param_t *param, FILE *fp_param, char *logfile){

  int i,j,k,ic,icount,ncore,ii;
  char filename[180];
  int ilk[10];
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

  static double rnl[N0][NPDM];

  xao[0]='s';
  xao[1]='p';
  xao[2]='d';
  xao[3]='f';
  xao[4]='g';

  for (i=0;i<10;i++)
    ilk[i]=0;

  ncore=param->norb-param->nval;

  fp_log = fopen(logfile, "a");
  fprintf(fp_log,"<<<do_champ -- champ>>>\n");  
  fclose(fp_log);

  zeff = param->z;
  for (i=0; i<param->norb - param->nval; i++)
    zeff -= param->wnl[i];

  sprintf(filename, "%s.loc", param->name);
  fp_file = fopen(filename, "rb");
  fread(rvloc, sizeof(double), param->ngrid, fp_file);
  fclose(fp_file);
  interp2_(rvloc,rvloc2);

  for (i=0; i<param->nll; i++) {
    sprintf(filename, "%s.psi.ps.l=%d", param->name,i);
    fp_file = fopen(filename, "rb");
    fread(rnl[i], sizeof(double), param->ngrid, fp_file);
    fclose(fp_file);
  }

  for (i=0; i<param->nll; i++) {
    sprintf(filename, "%s.pot.ps.l=%d", param->name,i);
    fp_file = fopen(filename, "rb");
    fread(dumm, sizeof(double), param->ngrid, fp_file);
    interp2_(dumm,dumm2);
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
  
  if (param->nboxes > 0) {
      lmax++;
      lloc=lmax;
      
      fp_log = fopen(logfile, "a");
      fprintf(fp_log," Making l+1 the local potential %d\n",lloc);
      fclose(fp_log);
      
      for (i=0; i<ngg; i++){
	  rvv2[icount][i] = rvloc2[i];
      }
  } 
  
  /* End new section */  
  
  if (param->rpcc > 0.){
      sprintf(filename, "%s.rho_pcore", param->name);
      fp_file = fopen(filename, "rb");
      fread(rscore, sizeof(double), param->ngrid, fp_file);
      interp2_(rscore,rscore2);
      fclose(fp_file);
  }
  
  sprintf(filename, "%s.rho_nl", param->name);
  fp_file = fopen(filename, "rb");
  fread(rho, sizeof(double), param->ngrid, fp_file);
  fclose(fp_file);
  interp2_(rho,rho2);

  /* section 1 : header */

  sprintf(filename, "%s.champ", param->name);
  fp_out = fopen(filename, "w");

  /* 1st line of champ format: The atom symbol, icorr, irel, nicore */
 
  switch (param->ixc) {
    
  case 0 : sprintf(xstr,"%s","ca"); break;
  case 1 : sprintf(xstr,"%s","pw"); break;
  case 2 : sprintf(xstr,"%s","pb"); break;

  default :
    fprintf(stderr,"Can not determine XC string for Champ format, using ca \n");
    sprintf(xstr,"%s","ca");
  }


  fprintf(fp_out,"# OPIUM generated potential for %2s\n", param->symbol); 

  fprintf(fp_out,"# %3d%7.2f%5d\n",lmax+1, zeff, 10);

  /* 3rd line of champ format: text : grid_type=3 ngrid a b  */
  fprintf(fp_out,"# %3d%5d%20.12E%20.12E%20.12E\n",3,ngg-1,rgrid_.a,rgrid_.b,exp(rgrid_.b)); 
  
  fprintf(fp_out," %3d%3d%5d\n",lmax+1,3,ngg-1); 

  /* section 2: radial grid */
  fprintf(fp_out,"Rydberg   units   \n"); 

  fprintf(fp_out,"r   V  xais yaxis   \n"); 
 

  for (j=1; j<ngg; j++) {
       fprintf(fp_out, "%#20.12E\t", rgrid_.r[j]);
       for (i=0; i<lmax+1; i++){
	    fprintf(fp_out, "%#20.12E\t",rvv2[i][j]/rgrid_.r[j]);
       }
       fprintf(fp_out, "\n");
  }
  fprintf(fp_out, "\n");

  writeparam(param, fp_out, fp_param);

/* added wissam */
  sprintf(filename, "%s.psi.champ", param->name);
  fp_file = fopen(filename, "w");
  fprintf(fp_file, "%5d%5d%5d", param->nll , 3 ,param->ngrid);
  fprintf(fp_file, "    OPIUM generated %s potential\n", param->symbol);

  for (i=0;i<param->ngrid;i++) {
    fprintf(fp_file, "%#20.12E\t",  rgrid_.r[j]);
    for (kk=0; kk<param->nll;kk++){
           fprintf(fp_file, "%#20.12E\t", rnl[kk][i]/rgrid_.r[j]);
       }
       fprintf(fp_file, "\n");
  }
  fprintf(fp_file, "\n");
  writeparam(param, fp_file, fp_param);

/* wissam end */

  fp_log = fopen(logfile, "a");  
  fprintf(fp_log,"   ================================================\n");
  fclose(fp_log);
  
  return 0;
}

