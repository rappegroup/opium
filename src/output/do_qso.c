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
 * generate *.qso.xml output (for QBOX)
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
#include "do_qso.h"          /* the module's own header */
#include "common_blocks.h"
#include "nlm.h"

#define MAX(a, b)   (((a) > (b)) ? (a):(b))
#define streq(a,b) (!strcasecmp(a,b))

void writeparam(param_t *param, FILE *fp, FILE *fp_param);
void nrelsproj(param_t *param, char *);
void interp2_(int *,double *,double *, double *);

int do_qso(param_t *param, FILE *fp_param, char *logfile){

  int i,j,k,ic,ncore,nmesh;
  char filename[180];
  FILE *fp,*fp2,*fp3;
  FILE *fp_log;
  int kk=2;
  time_t t;
  double zeff;
  double ztt;              
  double xmin;             
  double rmax;
  double mesh;
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
  double rmmax;

  static double rv[NPDM],rw[NPDM];
  static double dumm[NPDM];
  static double qsmesh[NPDM];

  fp_log = fopen(logfile, "a");
  fprintf(fp_log,"<<<do_qso>>>\n");  
  fclose(fp_log);

  if ((!strcmp(param->reltype, "nrl")) || (!strcmp(param->reltype, "srl"))){
    nrelorbnl(param,config,logfile);
    nrelsproj(param,logfile);
  } else {
    printf("!!ERROR!!: The qso can not do fully relativistic psps\n");
    fp_log = fopen(logfile, "a");
    fprintf(fp_log,"!!ERROR!!: The qso can not do fully relativistic psps\n");
    fclose(fp_log);
    return 1;
  }

  if (param->rpcc > 0.){
    printf(" qso can not do pcc currently\n");
    fp_log = fopen(logfile, "a");
    fprintf(fp_log," qso can not do pcc currently\n");
    fclose(fp_log);
    return 1;
  }

  ncore=aorb_.norb-aorb_.nval;
  ztt = pow(param->z,2./3.);
  xmin = log(ztt * param->a);
  rmax = exp(param->b*(param->ngrid-1)+xmin)/param->z;
  lloc=param->localind;

  zeff = param->z;
  for (i=0; i<param->norb - param->nval; i++)
    zeff -= param->wnl[i];
  
  lmax=aorb_.lo[param->nll-1];
  /*  mesh=0.002;
      nmesh=10000;*/
  nmesh=param->ngridl;
  mesh=param->lspc;
  rmmax=mesh*nmesh;
  
  for (j=0; j<nmesh; j++) {
    qsmesh[j]=(double)j *mesh;
  }

  sprintf(filename, "%s.qso.xml", param->name);
  fp = fopen(filename, "w");
  fprintf(fp,"<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");
  fprintf(fp,"<fpmd:species xmlns:fpmd=\"http://www.quantum-simulation.org/ns/fpmd/fpmd-1.0\"\n"); 
  fprintf(fp,"xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\"\n");
  fprintf(fp,"xsi:schemaLocation=\"http://www.quantum-simulation.org/ns/fpmd/fpmd-1.0\n");
  fprintf(fp,"species.xsd\">\n");
  fprintf(fp,"<description>\n");


  fprintf(fp,"OPIUM generated pseudopotential \n");
  time(&t);
  strftime(datestring, 7, "%y%m%d", localtime(&t));
  writeparam(param, fp, fp_param);  

  sprintf(filename, "%s.qso.xml", param->name);
  fp = fopen(filename, "a");
  fprintf(fp,"</description>\n");
  
  fprintf(fp,"<symbol>%s</symbol>\n",param->symbol);  
  fprintf(fp,"<atomic_number>%3.0f</atomic_number>\n",param->z);  
  fprintf(fp,"<mass>%4.0f</mass>\n",param->mass);  
  fprintf(fp,"<norm_conserving_pseudopotential>\n");  
  fprintf(fp,"<valence_charge>%lg</valence_charge>\n",zeff);    
  fprintf(fp,"<lmax>%d</lmax>\n",lmax);    
  fprintf(fp,"<llocal>%d</llocal>\n",lloc);    
  fprintf(fp,"<nquad>0</nquad>\n");    
  fprintf(fp,"<rquad>0</rquad>\n");    
  fprintf(fp,"<mesh_spacing>%lg</mesh_spacing>\n",mesh);    

  for (i=0; i<param->nll; i++) {

    sprintf(filename, "%s.pot.ps.l=%d", param->name,aorb_.lo[i]);
    fp2 = fopen(filename, "rb");
    sprintf(filename, "%s.psi.ps.l=%d", param->name,aorb_.lo[i]);
    fp3 = fopen(filename, "rb");
    
    fprintf(fp,"<projector l=\"%d\" size=\"%d\">\n",i,nmesh);  
    fprintf(fp,"<radial_potential>\n");

    fread(rv, sizeof(double), param->ngrid, fp2);
    fseek(fp2,sizeof(double) ,param->ngrid);

    for (j=0; j<param->ngrid; j++) {    
      /*printf(" rv: %d %20.10f  %20.10f\n", j,rv[j],grid_.r[j]);*/
      rv[j]=rv[j]/(2*grid_.r[j]);
    }
    interp2_(&nmesh,rv,dumm,qsmesh);

    for (j=0; j<nmesh; j++) {
      fprintf(fp,"%20.10f  \n",dumm[j]);      
    }

    fprintf(fp,"</radial_potential>\n");
    fprintf(fp,"<radial_function>\n");

    fread(rw, sizeof(double), param->ngrid, fp3);

    for (j=0; j<param->ngrid; j++) {    
      rw[j]=rw[j]/(grid_.r[j]);
    }

    interp2_(&nmesh,rw,dumm,qsmesh);

    for (j=0; j<nmesh; j++) {
      fprintf(fp,"%20.10f  \n",dumm[j]);      
    }
    
    fprintf(fp,"</radial_function>\n");
    fprintf(fp,"</projector>\n");
  } 
 fprintf(fp,"</norm_conserving_pseudopotential>\n");  
  fprintf(fp,"</fpmd:species>\n");  

  fclose(fp3);
  fclose(fp2);
  fclose(fp);
  return 0;
}





