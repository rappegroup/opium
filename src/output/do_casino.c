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
 * generate *.casino output (Casino semi-local potential)                *
 *****************************************************************************
 * INPUT:  param_t structure, *.ps, *.loc
 * OUTPUT: *.casino
 ****************************************************************************/

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#include "parameter.h"        /* defines structure: 'param_t' */
#include "cdim.h"        /* fortran code parameters */
#include "do_casino.h"          /* the module's own header */
#include "nlm.h"
#include "common_blocks.h"

void interp2_(int *, double *,double *, double *);
void writeparam(param_t *param, FILE *fp, FILE *fp_param);
void nrelsproj(param_t *param, char *);
int do_casino(param_t *param, FILE *fp_param, char *logfile){

  int i,j,k,ic,icount,ncore;
  char filename[180];
  FILE *fp;
  FILE *fp_log;
  double zeff;	                
  char xctype[20];              
  double r_1;
  double zz=0;
  int ngg=nrgrid_.nr-10;
  
  static double rvcore[N0][NPDM];
  static double r[NPDM];
  static double dumm[NPDM];
  static double dumm2[NPDM];
  static double rvv[N0][NPDM];

  fp_log = fopen(logfile, "a");
  fprintf(fp_log,"<<<do_casino>>>\n");  
  fclose(fp_log);

  /* section 0 : check if this potential format is supported */
  /* casino doesn't support DNL potentials */

  nrelsproj(param,logfile);
  /* new routine to set the arrays for semicore states */

  if (param->nboxes > 0) {
    fp_log = fopen(logfile, "a");
    fprintf(fp_log," casino format does not support the use of augmentation operators\n");
    fclose(fp_log);
    return 1;
  }

  if (param->rpcc > 1e-12) {
    fp_log = fopen(logfile, "a");
    fprintf(fp_log," casino format does not support the use of a partical core correction\n");
    fclose(fp_log);
    return 1;
  }

  for (i=0; i<param->nll; i++) {
    sprintf(filename, "%s.pot.ps.l=%d", param->name,i);

    fp = fopen(filename, "rb");
    fread(dumm, sizeof(double), param->ngrid, fp);    
    fseek(fp,sizeof(double) ,param->ngrid);
    interp2_(&nrgrid_.nr,dumm,dumm2,rgrid_.r);
    for (j=0; j<ngg; j++){
      rvv[i][j]=dumm2[j];
    }
    fclose(fp);
  }

  /* compute zeff */
  zeff = param->z;
  for (i=0; i<param->norb - param->nval; i++)
    zeff -= param->wnl[i];
  
  /* section 1 : header */

  sprintf(filename, "%s.casino", param->name);
  fp = fopen(filename, "w");

  /* Casino header */

  switch (param->ixc) {
  case -1 : 
    switch ((!strcmp(param->reltype, "nrl")) ? 0:1){
    case 0 : sprintf(xctype,"%s","HF"); break;
    case 1 : sprintf(xctype,"%s","DF"); break;
    default :
      printf(" Problem with Casino header \n");
      exit(1);
    };break;
  
  case 0 :     
    switch ((!strcmp(param->reltype, "nrl")) ? 0:1){
    case 0 : sprintf(xctype,"%s","PZLDA-NREL"); break;
    case 1 : sprintf(xctype,"%s","PZLDA-REL"); break;
    default :
      printf(" Problem with Casino header \n");
      exit(1);
    };break;

  case 1 :     
    switch ((!strcmp(param->reltype, "nrl")) ? 0:1){
    case 0 : sprintf(xctype,"%s","PWLDA-NREL"); break;
    case 1 : sprintf(xctype,"%s","PWLDA-REL"); break;
    default :
      printf(" Problem with Casino header \n");
      exit(1);
    };break;

  case 2 :     
    switch ((!strcmp(param->reltype, "nrl")) ? 0:1){
    case 0 : sprintf(xctype,"%s","PBEGGA-NREL"); break;
    case 1 : sprintf(xctype,"%s","PBEGGA-REL"); break;
    default :
      printf(" Problem with Casino header \n");
      exit(1);
    };break;

  case 3 :     
    switch ((!strcmp(param->reltype, "nrl")) ? 0:1){
    case 0 : sprintf(xctype,"%s","PW91GGA-NREL"); break;
    case 1 : sprintf(xctype,"%s","PW91GGA-REL"); break;
    default :
      printf(" Problem with Casino header \n");
      exit(1);
    };break;

  case 4 :     
    switch ((!strcmp(param->reltype, "nrl")) ? 0:1){
    case 0 : sprintf(xctype,"%s","WCGGA-NREL"); break;
    case 1 : sprintf(xctype,"%s","WCGGA-REL"); break;
    default :
      printf(" Problem with Casino header \n");
      exit(1);
    };break;

  case 5 :
    switch ((!strcmp(param->reltype, "nrl")) ? 0:1){
    case 0 : sprintf(xctype,"%s","PBEsolGGA-NREL"); break;
    case 1 : sprintf(xctype,"%s","PBEsolGGA-REL"); break;
    default :
      printf(" Problem with Casino header \n");
      exit(1);
    };break;

  default :
    printf(" Problem with Casino header \n");
    exit(1);
  }

  

  fprintf(fp,"%s Opium generated real space pseudopotential for %s \n",xctype,param->name);  
  fprintf(fp,"Atomic number and pseudo-charge \n");  
  fprintf(fp,"%-1.0f %-2.2f \n",param->z,zeff);
  fprintf(fp,"Energy units (rydberg/hartree/ev): \n");  
  fprintf(fp,"rydberg \n");  
  fprintf(fp,"Angular momentum of local component (0=s,1=p,2=d..)\n");  
  fprintf(fp,"%d \n",param->localind);    
  fprintf(fp,"NLRULE override (1) VMC/DMC (2) config gen (0 ==> input/default value) \n");
  fprintf(fp,"0 0 \n");  
  fprintf(fp,"Number of grid points \n");
  fprintf(fp,"%d \n",ngg);  
  fprintf(fp,"R(i) in atomic units\n");

  for (i=0; i<ngg ; i++) {    
    fprintf(fp,"    % 2.15E \n",rgrid_.r[i]);
  }
  
  for (j=0;j<param->nll ; j++) {
    fprintf(fp,"r*potential (L=%d) in Ry \n",j);
    for (i=0; i<ngg ; i++) {    
      fprintf(fp,"    % 2.15E \n",rvv[j][i]);
    }
  }
  
  writeparam(param, fp, fp_param);  

  fp_log = fopen(logfile, "a");
  fprintf(fp_log, "   ================================================\n");
  
  fclose(fp_log);
  return 0;
}

