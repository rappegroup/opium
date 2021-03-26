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

#ifndef M_PI 
#define M_PI 3.14159265358979323846
#endif

#include "parameter.h"        /* defines structure: 'param_t' */
#include "cdim.h"        /* fortran code parameters */
#include "do_ncpp.h"          /* the module's own header */
#include "nlm.h"

#define MAX(a, b)   (((a) > (b)) ? (a):(b))
#define streq(a,b) (!strcasecmp(a,b))

void writeparam(param_t *param, FILE *fp, FILE *fp_param);
void nrelsproj(param_t *param, char *);
double simpson(double *, double *, int); 
int do_qeupf(param_t *param, FILE *fp_param, char *logfile){

  int i,j,k,ic,ncore;
  char filename[180];
  FILE *fp;
  FILE *fp_log;
  int kk=2;
  double zeff;
  double ztt;              
  double xmin;             
  double rmax;
  int nn,ll;
  char lc[4];
  char upfversion[7];
  double occ;
  int nnl,ng;
  double dij[N0];
  int nind[N0];

  static double rscore[NPDM];
  static double nlcore[NPDM];
  static double rvcore[N0][NPDM];
  static double rnl[N0][NPDM];
  static double rr[NPDM];
  static double rab[NPDM];
  static double beta[NPDM];
  static double ddd[NPDM];

  sprintf(upfversion,"%s","\'2.0.0\'");
  fp_log = fopen(logfile, "a");
  fprintf(fp_log,"<<<do_upf>>>\n");  
  fclose(fp_log);

  nrelsproj(param,logfile);

  lc[0]='S';
  lc[1]='P';
  lc[2]='D';
  lc[3]='F';

  ncore=param->norb-param->nval;

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
      if (param->rc[j] < rr[i]) nind[j]=i;

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
    sprintf(filename, "%s.pot.ps.l=%d", param->name,i);
    fp = fopen(filename, "rb");
    fread(rvcore[i], sizeof(double), param->ngrid, fp);
    fseek(fp,sizeof(double) ,param->ngrid);
    fclose(fp);
  }

  for (i=0; i<param->nll; i++) {
    sprintf(filename, "%s.psi.ps.l=%d", param->name,i);
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



  /* First put in the param file: */

  sprintf(filename, "%s.upf", param->name);
  fp = fopen(filename, "w");

  fprintf(fp,"<UPF version=%s>\n",upfversion);
  fprintf(fp,"  <PP_INFO>\n");

  writeparam(param, fp, fp_param);  

  /* !! writeparam closes but doesnt open?*/
  fp = fopen(filename, "a");
  fprintf(fp,"  </PP_INFO>\n");

  /* start header */
  fprintf(fp,"  <!--                               -->\n");
  fprintf(fp,"  <!-- END OF HUMAN READABLE SECTION -->\n");
  fprintf(fp,"  <!--                               -->\n");


  fprintf(fp,"  <PP_HEADER generated='Generated using \"OPIUM\"'\n");
  fprintf(fp,"             author=\"anonymous\" \n");
  fprintf(fp,"             comment=\"\"\n");
  fprintf(fp,"             element=\"%s\"\n",param->symbol);
  fprintf(fp,"             pseudo_type=\"NC\"\n");

  if (streq(param->reltype,"nrl")){
    fprintf(fp,"             relativistic=\"no\"\n");
  }
  if (streq(param->reltype,"srl")){
    fprintf(fp,"             relativistic=\"scalar\"\n");
  }
  if (streq(param->reltype,"frl")){
    fprintf(fp,"             relativistic=\"full\"\n");
  }

  fprintf(fp,"             is_ultrasoft=\"F\"\n");
  fprintf(fp,"             is_paw=\"F\"\n");
  fprintf(fp,"             is_coulomb=\"F\"\n");
  fprintf(fp,"             has_so=\"F\"\n");
  fprintf(fp,"             has_wfc=\"T\"\n");
  fprintf(fp,"             has_gipaw=\"F\"\n");

  if (param->rpcc > 1e-6){
      fprintf(fp,"             core_correction=\"T\"\n");
  }else{
      fprintf(fp,"             core_correction=\"F\"\n");
  }

  fprintf(fp,"             functional=\"PBE\"\n");
  fprintf(fp,"             z_valence=\"%18.16le\"\n",zeff);
  fprintf(fp,"             total_psenergy=\"%18.16le\"\n",0.0);
  fprintf(fp,"             wfc_cutoff=\"%18.16le\" \n",0.0);
  fprintf(fp,"             rho_cutoff=\"%18.16le\" \n",0.0);
  fprintf(fp,"             l_max=\"%d\"\n",param->nll-1);
  fprintf(fp,"             l_max_rho=\"%d\"\n",2*(param->nll-1));
  fprintf(fp,"             l_local=\"%d\"\n",nlm_label(param->nlm[param->localind+ncore]).l);
  fprintf(fp,"             mesh_size=\"%d\"\n",param->ngrid);
  fprintf(fp,"             number_of_wfc=\"%d\"\n",param->nll);

  if (param->nboxes > 0) {
    fprintf(fp,"             number_of_proj=\"%d\"/>\n",param->nll);
  }else{
    fprintf(fp,"             number_of_proj=\"%d\"/>\n",param->nll);
  }
  
  fprintf(fp,"  <PP_MESH dx=\"%18.16le\" mesh=\"%d\" xmin=\"%18.16le\" rmax=\"%18.16le\" \nzmesh=\"%18.16le\">\n",param->b,param->ngrid,xmin,rmax,param->z);

  fprintf(fp,"    <PP_R type=\"real\" size=\"%d\" columns=\"4\">\n",param->ngrid);
  for (i=0;i<param->ngrid;i++){
    fprintf(fp," %19.16le ",rr[i]);
    if( (i+1)%4 == 0 ) fprintf(fp,"\n");
  }
  if(i%4 !=0 ) fprintf(fp,"\n");
  fprintf(fp,"    </PP_R>\n");

  fprintf(fp,"    <PP_RAB type=\"real\" size=\"%d\" columns=\"4\">\n",param->ngrid);
  for (i=0;i<param->ngrid;i++){
    fprintf(fp," %19.16le ",rab[i]);
    if( (i+1)%4 == 0 ) fprintf(fp,"\n");
  }
  if(i%4 !=0 ) fprintf(fp,"\n");

  fprintf(fp,"    </PP_RAB>\n");
  fprintf(fp,"  </PP_MESH>\n");

  fprintf(fp,"  <PP_LOCAL type=\"real\" size=\"%d\" columns=\"4\">\n",param->ngrid);
  for (i=0;i<param->ngrid;i++){
    fprintf(fp,"%19.16le ",nlcore[i]/rr[i]);
    if( (i+1)%4 == 0 ) fprintf(fp,"\n");
  }
  if(i%4 !=0 ) fprintf(fp,"\n");
  fprintf(fp,"  </PP_LOCAL>\n");


  fprintf(fp,"  <PP_NONLOCAL>\n");  
  nnl=0;
  for (j=0;j<param->nll;j++){  
    nn=nlm_label(param->nlm[ncore+j]).n;
    ll=nlm_label(param->nlm[ncore+j]).l;
    printf(" nn ll local %d %d %d \n",nn,ll,param->localind);
    if ((param->localind != j) || (param->nboxes > 0)) {
      nnl+=1; 
      fprintf(fp,"    <PP_BETA.%d type=\"real\" size=\"%d\" columns=\"4\" index=\"%d\" label=\"%d%c\" angular_momentum=\"%d\" cutoff_radius_index=\"%d\"\ncutoff_radius=\"%f\" norm_conserving_radius=\"%f\">\n",nnl,param->ngrid,j+1,nn,lc[ll],ll,nind[j],param->rc[j],param->rc[j]);
      
      for (i=0;i<param->ngrid;i++){  
	beta[i]=rnl[j][i]*(rvcore[j][i]-nlcore[i])/rr[i];
	ddd[i]=beta[i]*rnl[j][i];    
      }
      
      for (i=0;i<param->ngrid;i++){
	fprintf(fp," %19.16le ",beta[i]);
	if( (i+1)%4 == 0 ) fprintf(fp,"\n");
      }
      if(i%4 !=0 ) fprintf(fp,"\n");
      
      fprintf(fp,"    </PP_BETA.%d>\n",nnl);
      ng=param->ngrid;
      if (ng%2 !=0) ng-=1;
      dij[j]=simpson(ddd,rab,ng);

    }else{
    }
  }
  fprintf(fp,"    <PP_DIJ type=\"real\" size=\"%d\" columns=\"4\">\n",nnl);
  
  for (i=0;i<nnl;i++){
    fprintf(fp," %19.16le ",dij[i]);
    if( (i+1)%4 == 0 ) fprintf(fp,"\n");
  }
  if(i%4 !=0 ) fprintf(fp,"\n");
  fprintf(fp,"    </PP_DIJ>\n");
  fprintf(fp,"  </PP_NONLOCAL>\n");  

  fprintf(fp,"  <PP_PSPWFC>\n");  
  for (j=0;j<param->nll;j++){  
    nn=nlm_label(param->nlm[ncore+j]).n;
    ll=nlm_label(param->nlm[ncore+j]).l;
    occ=param->wnl[ncore+j];
    occ = MAX(occ,0);
    fprintf(fp,"    <PP_CHI.%d type=\"real\" size=\"%d\" columns=\"4\" index=\"%d\" label=\"%d%c\" l=\"%d\" occupation=\"%f\" n=\"%d\"\npseudo_energy=\" \" cutoff_radius=\"%f\" ultrasoft_cutoff_radius=\"%f\">\n",j+1,param->ngrid,j+1,nn,lc[ll],ll,occ,nn,param->rc[j],param->rc[j]);
  for (i=0;i<param->ngrid;i++){
    fprintf(fp," %19.16le ",rnl[j][i]);
    rscore[i]+=occ*rnl[j][i]*rnl[j][i];
    if( (i+1)%4 == 0 ) fprintf(fp,"\n");
  }
  if(i%4 !=0 ) fprintf(fp,"\n");
  fprintf(fp,"    </PP_CHI.%d>\n",j+1);
  }
  fprintf(fp,"  </PP_PSPWFC>\n");  

  fprintf(fp,"  <PP_RHOATOM type=\"real\" size=\"%d\" columns=\"4\">\n",param->ngrid);
  for (i=0;i<param->ngrid;i++){
    fprintf(fp," %19.16le ",rscore[i]);
    if( (i+1)%4 == 0 ) fprintf(fp,"\n");
  }
  if(i%4 !=0 ) fprintf(fp,"\n");
  fprintf(fp,"  </PP_RHOATOM>\n");

  fprintf(fp,"</UPF>\n");

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


