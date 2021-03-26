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
void nrelsproj(param_t *param, char *);
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
  nrelsproj(param,logfile);

  fp_log = fopen(logfile, "a");
  fprintf(fp_log,"<<<do_champ -- champ>>>\n");  
  fclose(fp_log);

  zeff = param->z;
  for (i=0; i<param->norb - param->nval; i++)
    zeff -= param->wnl[i];
  
  
  /*  sprintf(filename, "%s.loc", param->name);
  if (fp_file = fopen(filename, "rb")) {
    fread(rvloc, sizeof(double), param->ngrid, fp_file);
    fclose(fp_file);
  } else {
    fp_log = fopen(logfile, "a");
    fprintf(fp_log,"Looks like you never ran nl yet you have augmentation functions :( --EXIT!\n");
    printf("Looks like you never ran nl yet you have augmentation functions :( --EXIT!\n");
    fclose(fp_log);
    exit(1);
    }*/


  for (i=0; i<param->nll; i++) {
    sprintf(filename, "%s.psi.ps.l=%d", param->name,i);
    fp_file = fopen(filename, "rb");
    fread(rnl[i], sizeof(double), param->ngrid, fp_file);
    fclose(fp_file);
  }

  for (i=0; i<param->nll; i++) {
    sprintf(filename, "%s.pot.ps.l=%d", param->name,i);
    fp_file = fopen(filename, "rb");
    fread(rvv[i], sizeof(double), param->ngrid, fp_file);
    fclose(fp_file);
  }

  lmax=param->nll-1;
  lloc=nlm_label(param->nlm[param->localind+ncore]).l;


  /* End new section */  
  
  if (param->rpcc > 0.){
      sprintf(filename, "%s.rho_pcore", param->name);
      fp_file = fopen(filename, "rb");
      fread(rscore, sizeof(double), param->ngrid, fp_file);
      fclose(fp_file);
  }


  /*  
  sprintf(filename, "%s.rho_nl", param->name);
  fp_file = fopen(filename, "rb");
  fread(rho, sizeof(double), param->ngrid, fp_file);
  fclose(fp_file);
  */


  /* section 1 : header */

  sprintf(filename, "%s.champ", param->name);
  fp_out = fopen(filename, "w");
 
  fprintf(fp_out,"# OPIUM generated potential for %2s (Hartree units)\n", param->symbol); 
  fprintf(fp_out,"%3d%7.2f%5d : npot, z_val, r_asymp   \n",lmax+1, zeff, 0.0);

  fprintf(fp_out,"%3d%5d%20.12E%20.12E : igrid_ps, nr_ps, r0_ps, h_ps \n",2,param->ngrid,grid_.r[0],grid_.h); 
  fprintf(fp_out, "%#20.12E\t", grid_.r[0]);
  for (i=0; i<lmax+1; i++){
     fprintf(fp_out, "%#20.12E\t",0.5*rvv[i][0]/grid_.r[0]);
  }
  fprintf(fp_out, "\n");


  for (j=1; j<param->ngrid; j++) {
    fprintf(fp_out, "%#20.12E\t", grid_.r[j]);
    for (i=0; i<lmax+1; i++){
      fprintf(fp_out, "%#20.12E\t",0.5*rvv[i][j]/grid_.r[j]);
    }
    fprintf(fp_out, "\n");
  }
  fprintf(fp_out, "\n");


  writeparam(param, fp_out, fp_param);

/* added wissam */
  sprintf(filename, "%s.psi.champ", param->name);
  fp_file = fopen(filename, "w");
  fprintf(fp_file, "%5d%5d%5d", param->nll , 2 ,param->ngrid);
  fprintf(fp_file, "    OPIUM generated %s potential\n", param->symbol);

  for (i=0;i<param->ngrid;i++) {
    fprintf(fp_file, "%#20.12E\t",  grid_.r[j]);
    for (kk=0; kk<param->nll;kk++){
           fprintf(fp_file, "%#20.12E\t", rnl[kk][i]/grid_.r[j]);
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

