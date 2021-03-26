/*
 * Copyright (c) 1998-2011 The OPIUM Group
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
 * generate *.teter output (for ABINIT)
 ******************************************************************************/
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <stdlib.h>

#ifndef M_PI 
#define M_PI 3.14159265358979323846
#endif

#include "parameter.h"        /* defines structure: 'param_t' */
#include "cdim.h"        /* fortran code parameters */
#include "do_teter.h"          /* the module's own header */
#include "common_blocks.h"
#include "nlm.h"

#define MAX(a, b)   (((a) > (b)) ? (a):(b))
#define streq(a,b) (!strcasecmp(a,b))

void writeparam(param_t *param, FILE *fp, FILE *fp_param);
void nrelsproj(param_t *param, char *);
int do_teter(param_t *param, FILE *fp_param, char *logfile){

  int i,j,k,ic,ncore;
  char filename[180];
  FILE *fp,*fp2;
  FILE *fp_log;
  int kk=2;
  time_t t;
  double zeff;
  double ztt;              
  double xmin;             
  double rmax;
  char datestring[7];
  int nn,ll,lixc;
  int config=-1;
  char lc[4];
  double occ;
  int nnl,ng;
  double dij[N0];
  int nind[N0];
  char sgn='+';
  int lloc=0;
  int lmax=0;

  static double rscore[NPDM];
  static double nlcore[NPDM];
  static double rvcore[N0][NPDM];
  static double rnl[N0][NPDM];
  static double rr[NPDM],rvloc[NPDM];
  static double rab[NPDM];
  static double beta[NPDM];
  static double ddd[NPDM];

  fp_log = fopen(logfile, "a");
  fprintf(fp_log,"<<<do_teter>>>\n");  
  fclose(fp_log);

  lc[0]='S';
  lc[1]='P';
  lc[2]='D';
  lc[3]='F';

  if ((!strcmp(param->reltype, "nrl")) || (!strcmp(param->reltype, "srl"))){
    /*    nrelorbnl(param,config,logfile);
	  nrelsproj(param,logfile);*/
    printf("!!ERROR!!: Currently the Teter format can only be used for fully relativistic psps \n");  
    fp_log = fopen(logfile, "a");
    fprintf(fp_log,"!!ERROR!!: Currently the Teter format can only be used for fully relativistic psps \n");
    fclose(fp_log);
    return 1;
  } else {
    relorbnl(param,config,logfile);
    /*relsproj(param,logfile);*/
  }

  ncore=aorb_.norb-aorb_.nval;

  for (i=0;i<param->ngrid;i++)
    rscore[i]=0.0;

  ztt = pow(param->z,2./3.);
  xmin = log(ztt * param->a);
  rmax = exp(param->b*(param->ngrid-1)+xmin)/param->z;
  lloc=param->localind;
  for (i=0;i<param->ngrid;i++){
    rr[i]=exp(xmin+i*param->b)/param->z;
    rab[i]=param->b*rr[i];
  }
  
  for (j=0;j<param->nll;j++)
    for (i=0;i<param->ngrid;i++)
      if (aval_.rcall[j] < rr[i]) nind[j]=i;

  /* compute zeff */
  zeff = param->z;
  for (i=0; i<param->norb - param->nval; i++)
    zeff -= param->wnl[i];

  for (i=0; i<param->nll; i++) {
    if (!strcmp(param->reltype, "frl")) {
      sgn=(adat_.so[i] > 0) ? '+' : '-' ; 
      sprintf(filename, "%s.pot.ps.l=%d%c", param->name,aorb_.lo[i],sgn);
    }else{
      sprintf(filename, "%s.pot.ps.l=%d", param->name,aorb_.lo[i]);
    }
    fp = fopen(filename, "rb");
    fread(rvcore[i], sizeof(double), param->ngrid, fp);
    fseek(fp,sizeof(double) ,param->ngrid);
    fclose(fp);
  }

  for (i=0; i<param->nll; i++) {
    if (!strcmp(param->reltype, "frl")) {
      sgn=(adat_.so[i] > 0) ? '+' : '-' ; 
      sprintf(filename, "%s.psi.ps.l=%d%c", param->name,aorb_.lo[i],sgn);

    }else{
      sprintf(filename, "%s.psi.ps.l=%d", param->name,aorb_.lo[i]);
    }      
    fp = fopen(filename, "rb");
    fread(rnl[i], sizeof(double), param->ngrid, fp);
    fclose(fp);
  }

  if (param->rpcc > 0.){
    sprintf(filename, "%s.rho_pcore", param->name);
    fp = fopen(filename, "rb");
    fread(rscore, sizeof(double), param->ngrid, fp);
    fclose(fp);
  }
  /* read in the local potential  */

  lmax=aorb_.lo[param->nll-1];
  lloc=param->localind;

  if (param->nboxes > 0) {
    /*    fprintf(fp_log," fhi format does not support the use of augmentation operators\n");*/
    
    sprintf(filename, "%s.loc", param->name);
    if (fp = fopen(filename, "rb")) {
      fread(rvloc, sizeof(double), param->ngrid, fp);
      fclose(fp);
    } else {
      fp_log = fopen(logfile, "a");
      fprintf(fp_log,"Looks like you never ran nl yet you have augmentation functions :( --EXIT!\n");
      printf("Looks like you never ran nl yet you have augmentation functions :( --EXIT!\n");
      fclose(fp_log);
      exit(1);
    }

    fp_log = fopen(logfile, "a");
    fprintf(fp_log," Making l+1 the local potential %d\n",kk);
    fclose(fp_log);

    for (i=0; i<param->ngrid; i++){
      rvcore[param->nll][i]=rvloc[i];
      rvcore[param->nll+1][i]=rvloc[i];
    }
    lloc = aorb_.lo[param->nll-1]+1;
    lmax=lloc;
    aorb_.lo[param->nll]=lmax;
    aorb_.lo[param->nll+1]=lmax;
    param->nll++;
    param->nll++;
  }

  if (param->ixc == 0) lixc=2;
  if (param->ixc < 0) lixc=99;
  if (param->ixc == 1) lixc=7;
  if (param->ixc == 2) lixc=11;
  if (param->ixc == 3) lixc=99;
  if (param->ixc == 4) lixc=23;
  if (param->ixc == 5) lixc=24;

  /* start header */
  sprintf(filename, "%s.teter", param->name);
  fp = fopen(filename, "w");
  
  fprintf(fp,"OPIUM generated psp \n");
  time(&t);
  strftime(datestring, 7, "%y%m%d", localtime(&t));
  fprintf(fp, "%10.5f%10.5f  %s\t\t\tzatom,zion,pspdat\n",param->z,zeff,datestring);
  fprintf(fp, "%5d%5d%5d%5d%10d%10.5f\tpspcod,pspxc,lmax,lloc,mmax,r2well\n",5,
	  lixc,lmax,lloc,param->ngrid,0.);
  fprintf(fp," %22.15e %10.5f 2 r1,al,pspso\n",rr[0],param->b);
  for (j=0;j<param->nll;j++){  
    ll=aorb_.lo[j];
    if (ll == 0) {
      fprintf(fp," %d 0 0 1 0.00 l,e99.0,e99.9,nproj,rcpsp \n",ll);    
    }else{
      fprintf(fp," %d 0 0 2 0.00 l,e99.0,e99.9,nproj,rcpsp \n",ll);    
    }
    fprintf(fp,"     .00000000    .0000000000    .0000000000    .00000000   rms,ekb1,ekb2,epsatm\n");
  }
  fprintf(fp," 0   0   0   rchrg,fchrg,qchrg\n");

  for (j=0;j<param->nll;j++){    
    ll=aorb_.lo[j];
    fprintf(fp," %d =l for pseudopotential\n",ll);
    for (i=0;i<param->ngrid;i++){
      fprintf(fp,"%21.13e ",rvcore[j][i]/(2.0*rr[i]));
      if( (i+1)%3 == 0 ) fprintf(fp,"\n");
    }
    if(i%3 !=0 ) fprintf(fp,"\n");
  }    

  for (j=0;j<param->nll;j++){    
    ll=aorb_.lo[j];
    fprintf(fp," %d =l for wavefunctions\n",ll);
    for (i=0;i<param->ngrid;i++){
      fprintf(fp,"%21.13e ",rnl[j][i]);
      if( (i+1)%3 == 0 ) fprintf(fp,"\n");
    }
    if(i%3 !=0 ) fprintf(fp,"\n");
  }    
  
  fprintf(fp,"OPIUM generated pseudopotential\n");
  writeparam(param, fp, fp_param);  




  return 0;
}



