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
 * generate *.upf output (QE-PWSCF Norm-Conserving compatible)                *
 ******************************************************************************/
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#ifndef M_PI 
#define M_PI 3.14159265358979323846
#endif

#include "parameter.h"        /* defines structure: 'param_t' */
#include "cdim.h"        /* fortran code parameters */
#include "do_ncpp.h"          /* the module's own header */
#include "common_blocks.h"
#include "nlm.h"

#define MAX(a, b)   (((a) > (b)) ? (a):(b))
#define streq(a,b) (!strcasecmp(a,b))

void writeparam(param_t *param, FILE *fp, FILE *fp_param);
void nrelsproj(param_t *param, char *);
double simpson(double *, double *, int); 
int do_qeupf(param_t *param, FILE *fp_param, char *logfile){

  int i,j,k,ic,ncore;
  char filename[180];
  FILE *fp,*fp2;
  FILE *fp_log;
  int kk=2;
  double zeff;
  double ztt;              
  double xmin;             
  double rmax;
  int nn,ll;
  int config=-1;
  char lc[4];
  double occ;
  int nnl,ng;
  double dij[N0];
  int nind[N0];
  char sgn='+';

  static double rscore[NPDM];
  static double nlcore[NPDM];
  static double rvcore[N0][NPDM];
  static double rnl[N0][NPDM];
  static double rr[NPDM];
  static double rab[NPDM];
  static double beta[NPDM];
  static double ddd[NPDM];

  fp_log = fopen(logfile, "a");
  fprintf(fp_log,"<<<do_upf>>>\n");  
  fclose(fp_log);

  lc[0]='S';
  lc[1]='P';
  lc[2]='D';
  lc[3]='F';

  if ((!strcmp(param->reltype, "nrl")) || (!strcmp(param->reltype, "srl"))) {
    nrelorbnl(param,config,logfile);
    nrelsproj(param,logfile);

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


  /* read in the local potential  */
  sprintf(filename, "%s.loc", param->name);
  fp = fopen(filename, "rb");
  fread(nlcore, sizeof(double), param->ngrid, fp);
  fclose(fp);


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
  /* read in pcc */
  if (param->rpcc > 0.){
    sprintf(filename, "%s.rho_pcore", param->name);
    fp = fopen(filename, "rb");
    fread(rscore, sizeof(double), param->ngrid, fp);
    fclose(fp);
  }
  

  /* First put in the param file: */

  sprintf(filename, "%s.upf", param->name);
  fp = fopen(filename, "w");

  fprintf(fp,"  <PP_INFO>\n");
  writeparam(param, fp, fp_param);  

  /* !! writeparam closes but doesnt open?*/
  fp = fopen(filename, "a");
  fprintf(fp,"  </PP_INFO>\n");

  /* start header */
  
  fprintf(fp,"\n\n\n<PP_HEADER>\n");
  fprintf(fp,"   0         Version Number\n");
  fprintf(fp,"   %s        Element\n",param->symbol);  
  fprintf(fp,"   NC        Norm - Conserving pseudopotential\n");
  if (param->rpcc > 1e-6){
    fprintf(fp,"    T      Nonlinear Core Correction\n");
  }else{
    fprintf(fp,"    F      Nonlinear Core Correction\n");
  }
  if (param->ixc == 0) {
    fprintf(fp,"SLA  PZ   NOGX NOGC    PZ   Exchange-Correlation functional\n");
  } else if (param->ixc == 2) {
    fprintf(fp,"SLA  PW   PBE  PBE     PBE  Exchange-Correlation functional\n");
  } else {
    printf("!!WARNING!!: Your choice of XC functional is not currently supported for this type of output (UPF), will print LDA Perdew-Zunger\n");
    fprintf(fp,"SLA  PZ   NOGX NOGC    PZ   Exchange-Correlation functional\n");
  }
  fprintf(fp," %lg          Z valence\n",zeff);
  fprintf(fp," %lg          Total energy\n",0.0);
  fprintf(fp," %f   %f     Suggested cutoff for wfc and rho\n",0.0,0.0);
  fprintf(fp," %d           Max angular momentum component\n",aorb_.lo[param->nll-1]);  
  fprintf(fp," %d           Number of points in mesh\n",param->ngrid);  

  if (param->nboxes > 0) {
    fprintf(fp," %d  %d     Number of Wavefuncitons, Number of Projectors\n",param->nll,param->nll);
  }else{
    fprintf(fp," %d  %d     Number of Wavefuncitons, Number of Projectors\n",param->nll,param->nll-1);
  }
  fprintf(fp," Wavefunctions         nl  l   occ\n");

  for (j=0;j<param->nll;j++) {
    nn=aorb_.no[j];
    ll=aorb_.lo[j];
    fprintf(fp,"                       %c  %d  %f\n",lc[ll],ll,adat_.wnl[j]);    
  }
  fprintf(fp,"</PP_HEADER>\n");
  
  fprintf(fp,"<PP_MESH>\n");
  fprintf(fp,"  <PP_R>\n");
  for (i=0;i<param->ngrid;i++){
    fprintf(fp," %19.16le ",rr[i]);
    if( (i+1)%4 == 0 ) fprintf(fp,"\n");
  }
  if(i%4 !=0 ) fprintf(fp,"\n");
  fprintf(fp,"  </PP_R>\n");
  fprintf(fp,"  <PP_RAB>\n");
 
  for (i=0;i<param->ngrid;i++){
    fprintf(fp," %19.16le ",rab[i]);
    if( (i+1)%4 == 0 ) fprintf(fp,"\n");
  }
  if(i%4 !=0 ) fprintf(fp,"\n");

  fprintf(fp,"  </PP_RAB>\n");
  fprintf(fp,"</PP_MESH>\n");

  if (param->rpcc > 1e-6){
    
    fprintf(fp,"\n\n<PP_NLCC>\n");
    for (i=0;i<param->ngrid;i++){
      fprintf(fp,"%19.16le ",rscore[i]/rr[i]/rr[i]/4.0/M_PI);
      if( (i+1)%4 == 0 ) fprintf(fp,"\n");
    }
    if(i%4 !=0 ) fprintf(fp,"\n");
    fprintf(fp,"</PP_NLCC>\n");
  }

  fprintf(fp,"\n\n<PP_LOCAL>\n");
  for (i=0;i<param->ngrid;i++){
    fprintf(fp,"%19.16le ",nlcore[i]/rr[i]);
    if( (i+1)%4 == 0 ) fprintf(fp,"\n");
  }
  if(i%4 !=0 ) fprintf(fp,"\n");
  fprintf(fp,"</PP_LOCAL>\n");

  sprintf(filename, "%s.beta", param->name);
  fp2 = fopen(filename, "w");
  
  fprintf(fp,"\n\n<PP_NONLOCAL>\n");  
  nnl=0;
  for (j=0;j<param->nll;j++){  
    nn=aorb_.no[j];
    ll=aorb_.lo[j];

    for (i=0;i<param->ngrid;i++){  
      fprintf(fp2," %19.16le  %19.16le %d \n",rr[i],rnl[j][i],j);	
    }
    fprintf(fp2,"&\n");

    if ((param->localind != j) || (param->nboxes > 0)) {
      nnl+=1; 
      fprintf(fp,"  <PP_BETA>\n");
      
      for (i=0;i<param->ngrid;i++){  
	beta[i]=rnl[j][i]*(rvcore[j][i]-nlcore[i])/rr[i];
	ddd[i]=beta[i]*rnl[j][i];    
      }
      
      fprintf(fp,"  %d  %d      Beta    L\n",nnl,ll);
      fprintf(fp,"  %d  \n",nind[j]);

      for (i=0;i<param->ngrid;i++){
	fprintf(fp," %19.16le ",beta[i]);
	if( (i+1)%4 == 0 ) fprintf(fp,"\n");
      }
      if(i%4 !=0 ) fprintf(fp,"\n");
      
      fprintf(fp,"  </PP_BETA>\n");
      ng=param->ngrid;
      if (ng%2 !=0) ng-=1;
      dij[nnl-1]=1.0/simpson(ddd,rab,ng);

    }else{
    }
  }
  fprintf(fp,"  <PP_DIJ>\n");
  fprintf(fp," %d   Number of nonzero Dij\n",nnl);
  
  for (i=0;i<nnl;i++){
    fprintf(fp," %d %d %19.16le \n",i+1,i+1,dij[i]);
  }
  fprintf(fp,"  </PP_DIJ>\n");
  fprintf(fp,"</PP_NONLOCAL>\n");  
  fclose(fp2);



  fprintf(fp,"<PP_PSWFC>\n");  
  for (j=0;j<param->nll;j++){  
    nn=aorb_.no[j];
    ll=aorb_.lo[j];
    occ=adat_.wnl[j];
    occ = MAX(occ,0);
    fprintf(fp,"%c    %d  %5.2f          Wavefunction\n",lc[ll],ll,occ);
    for (i=0;i<param->ngrid;i++){
      fprintf(fp," %19.16le ",rnl[j][i]);
      rscore[i]+=occ*rnl[j][i]*rnl[j][i];
      if( (i+1)%4 == 0 ) fprintf(fp,"\n");
    }
    if(i%4 !=0 ) fprintf(fp,"\n");
  }
  fprintf(fp,"</PP_PSWFC>\n");  

  fprintf(fp,"\n\n<PP_RHOATOM>\n");
  for (i=0;i<param->ngrid;i++){
    fprintf(fp," %19.16le ",rscore[i]);
    if( (i+1)%4 == 0 ) fprintf(fp,"\n");
  }
  if(i%4 !=0 ) fprintf(fp,"\n");
  fprintf(fp,"</PP_RHOATOM>\n");

  if (!strcmp(param->reltype, "frl")) {
    fprintf(fp,"\n\n<PP_ADDINFO>\n");
    for (j=0;j<param->nll;j++){  
      nn=aorb_.no[j];
      ll=aorb_.lo[j];
      occ=adat_.wnl[j];
      occ = MAX(occ,0);
      fprintf(fp," %d%c  %d  %d  %5.2f %5.2f \n",nn,lc[ll],nn,ll,ll+adat_.so[j],occ);
    }    
    for (j=0;j<param->nll;j++){  
      if ((param->localind != j) || (param->nboxes > 0)) {      
        nn=aorb_.no[j];
	ll=aorb_.lo[j];
	occ=adat_.wnl[j];
	occ = MAX(occ,0);
	fprintf(fp,"  %d  %5.2f \n",ll,ll+adat_.so[j]);
      }
    }
    fprintf(fp,"%12.8f,  %12.8f,  %12.8f,  %12.8f \n",xmin,rmax,param->z,param->b);
    fprintf(fp,"</PP_ADDINFO>\n");
  }
  fclose(fp);
  return 0;
}

double simpson(double *fn, double *rab, int ng){

  double val=0.0;
  double r12=1.0/12.0;
  double p1,p2,p3;
  int i;

  p3=r12*fn[0]*rab[0];
  for (i=1;i<ng-1;i+=2) {
    p1=p3;
    p2=r12*fn[i]*rab[i];
    p3=r12*fn[i+1]*rab[i+1];
    val+=(p3+p1)*4.0+16.0*p2;
  }
  return val;
}


