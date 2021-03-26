/*
 * Copyright (c) 1998-2004 The OPIUM Group
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
/* 
 * $Id: do_qplot.c,v 1.7 2004/10/02 18:34:49 ewalter Exp $
 */


/* Bessel transform of psi and V_i */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>


#include "parameter.h"        /* defines structure: 'param_t' */
#include "fortparam.h"
#include "do_qplot.h"           /* the module's own header */
#include "common_blocks.h"
#include "nlm.h"

void btrans_(double *, char *);
void readPS(param_t *param);

int do_qplot(param_t *param, char *logfile){

  int i;
  int scount=0;
  int pcount=0;
  int dcount=0;
  int fcount=0;
  int lcolor=0;
  int lsty=0;
  int ncore;
  int ic;

  char filename[160];
  FILE *parm;
  FILE *fp;
  char *comm;
  char lc = 0;
  double zeff;

  comm= (char *) malloc(120*sizeof(char));

  readPS(param);

  ncore=param->norb-param->nval;

  /*ic=0;
  sprintf(filename, "%s.pot_ps", param->name);
  fp = fopen(filename, "rb");
  for (i=0; i<param->nval; i++) {
    if (i == param->ipot[ic]) {
      fread(totpot_.rvcore[ic], sizeof(double), param->ngrid, fp);
      ic++;
    }else{
      fseek(fp,sizeof(double)*param->ngrid,1);
    }
  }
  fclose(fp);*/

  /*
  sprintf(filename, "%s.eig_ae", param->name);
  fp = fopen(filename, "rb");
  for (i=0; i<param->nval; i++){
    fread(&atomic_.en[i], sizeof(double),1, fp);
    fread(&atomic_.wnl[i], sizeof(double),1, fp);
    fread(&atomic_.nlm[i], sizeof(int),1, fp);
    fread(&nmax_.nmax[i], sizeof(int),1, fp);
    fread(&nmax_.maxim, sizeof(int),1, fp);
    fread(&atomic_.xion, sizeof(double),1, fp);
  }
  fclose(fp); */

  sprintf(filename, "%s.loc", param->name);
  fp = fopen(filename, "rb");
  fread(nlcore_.rvloc, sizeof(double), param->ngrid, fp);
  fclose(fp);

  sprintf(filename, "%s.psi_nl", param->name);
  fp = fopen(filename, "rb");
  for (i=0; i<param->nval; i++) {
    fread(atomic_.rnl[i], sizeof(double), param->ngrid, fp);
  }
  fclose(fp);

  zeff=atomic_.xion;
  for (i=0; i<param->nval; i++) {
    zeff +=param->wnl[i+ncore];
  }
  nll_.nll=param->nll;

  btrans_(&zeff,param->name);

  if ((parm = fopen("vq.par","w")) != NULL) {
    
    fprintf(parm,"# q-potential param file for xmgrace\n");
    fprintf(parm,"g0 on\n");
    fprintf(parm,"with g0\n");
    fprintf(parm,"world xmin 0\n");
    fprintf(parm,"world xmax 100\n");
    fprintf(parm,"title \"Bessel transform of V_ionic: %s\"\n",param->symbol);
    fprintf(parm,"title font 0\n");
    fprintf(parm,"title size 1.500000\n");
    fprintf(parm,"title color 1\n");
    fprintf(parm,"subtitle \"atom name: %s\"\n",param->name);
    fprintf(parm,"subtitle font 0\n");
    fprintf(parm,"subtitle size 1.000000\n");
    fprintf(parm,"subtitle color 1\n");
    fprintf(parm,"xaxis on\n");
    fprintf(parm,"xaxis tick major 10\n");
    fprintf(parm,"xaxis tick minor 5\n");
    fprintf(parm,"xaxis label \"q\"\n");
    fprintf(parm,"yaxis on\n");
    fprintf(parm,"yaxis label \"V_ion(q)\"\n");
    fprintf(parm,"legend on\n");
    fprintf(parm,"legend loctype view\n");
    fprintf(parm,"legend 0.85, 0.8\n");
    
    for (i=0; i<param->nll;i++){
      
      if (nlm_label(param->nlm[param->ipot[i]+ncore]).l == 0) {
	scount++;
	lcolor=1;
	lsty=scount;
	lc='s';
      }else if (nlm_label(param->nlm[param->ipot[i]+ncore]).l == 1) {
	pcount++;
	lcolor=2;
	lsty=pcount;
	lc='p';
      }else if (nlm_label(param->nlm[param->ipot[i]+ncore]).l == 2) {
	dcount++;
	lcolor=3;
	lsty=dcount;
	lc='d';
      }else if (nlm_label(param->nlm[param->ipot[i]+ncore]).l == 3) {
	fcount++;
	lcolor=4;
	lsty=fcount;
	lc='f';
      }
      fprintf(parm," s%d hidden false \n",i);
      fprintf(parm," s%d type xy \n",i);
      fprintf(parm," s%d symbol 0 \n",i);
      fprintf(parm," s%d line type 1 \n",i);
      fprintf(parm," s%d line linestyle %d \n",i,lsty);
      fprintf(parm," s%d line linewidth 2.0 \n",i);
      fprintf(parm," s%d line color %d \n",i,lcolor);
      fprintf(parm," s%d legend \"V_%d%c\" \n",i,
	      nlm_label(param->nlm[param->ipot[i]+ncore]).n,lc);
    }
    
    fprintf(parm," s%d hidden false \n",i);
    fprintf(parm," s%d type xy \n",i);
    fprintf(parm," s%d symbol 0 \n",i);
    fprintf(parm," s%d line type 1 \n",i);
    fprintf(parm," s%d line linestyle %d \n",i,3);
    fprintf(parm," s%d line linewidth 3.0 \n",i);
    fprintf(parm," s%d line color %d \n",i,14);
    fprintf(parm," s%d legend \"V_loc\"\n" ,i);
    

    fprintf(parm," s%d hidden false \n",i+1);
    fprintf(parm," s%d type xy \n",i+1);
    fprintf(parm," s%d symbol 0 \n",i+1);
    fprintf(parm," s%d line type 1 \n",i+1);
    fprintf(parm," s%d line linestyle %d \n",i+1,3);
    fprintf(parm," s%d line linewidth 2.5 \n",i+1);
    fprintf(parm," s%d line color %d \n",i+1,13);
    fprintf(parm," s%d legend \"Rho_pcc\"\n" ,i+1);
    
    
    fclose(parm);
    
  } 
  
  sprintf(comm, "xmgrace %s.plt_vq -autoscale y -p vq.par  -saveall %s_viq.agr & ", param->name,param->name);
  system(comm);

  scount=0;
  pcount=0;
  dcount=0;
  fcount=0;
  lcolor=0;
  lsty=0;
 
  if ((parm = fopen("psiq.par","w")) != NULL) {
    
    fprintf(parm,"# q-wfn param file for xmgrace\n");
    fprintf(parm,"g0 on\n");
    fprintf(parm,"with g0\n");
    fprintf(parm,"world xmin 0\n");
    fprintf(parm,"world xmax 10\n");
    fprintf(parm,"title \"Bessel transform of Psi: %s\"\n",param->symbol);
    fprintf(parm,"title font 0\n");
    fprintf(parm,"title size 1.500000\n");
    fprintf(parm,"title color 1\n");
    fprintf(parm,"subtitle \"atom name: %s\"\n",param->name);
    fprintf(parm,"subtitle font 0\n");
    fprintf(parm,"subtitle size 1.000000\n");
    fprintf(parm,"subtitle color 1\n");
    fprintf(parm,"xaxis on\n");
    fprintf(parm,"xaxis tick major 1\n");
    fprintf(parm,"xaxis tick minor 0.5\n");
    fprintf(parm,"xaxis label \"q\"\n");
    fprintf(parm,"yaxis on\n");
    fprintf(parm,"yaxis label \"Psi(q)\"\n");
    fprintf(parm,"legend on\n");
    fprintf(parm,"legend loctype view\n");
    fprintf(parm,"legend 0.85, 0.8\n");
    
    for (i=0; i<param->nval;i++){
      
      if (nlm_label(param->nlm[i+ncore]).l == 0) {
	scount++;
	lcolor=1;
	lsty=scount;
	lc='s';
      }else if (nlm_label(param->nlm[i+ncore]).l == 1) {
	pcount++;
	lcolor=2;
	lsty=pcount;
	lc='p';
      }else if (nlm_label(param->nlm[i+ncore]).l == 2) {
	dcount++;
	lcolor=3;
	lsty=dcount;
	lc='d';
      }else if (nlm_label(param->nlm[i+ncore]).l == 3) {
	fcount++;
	lcolor=4;
	lsty=fcount;
	lc='f';
      }

      fprintf(parm," s%d hidden false \n",i);
      fprintf(parm," s%d type xy \n",i);
      fprintf(parm," s%d symbol 0 \n",i);
      fprintf(parm," s%d line type 1 \n",i);
      fprintf(parm," s%d line linestyle %d \n",i,lsty);
      fprintf(parm," s%d line linewidth 2.0 \n",i);
      fprintf(parm," s%d line color %d \n",i,lcolor);
      fprintf(parm," s%d legend \"Psi_%d%c\" \n",i,
	      nlm_label(param->nlm[i+ncore]).n,lc);
    }
    
    fclose(parm);
    
  } 
    
  sprintf(comm, "xmgrace %s.plt_wq -autoscale y -p psiq.par  -saveall %s_psiq.agr & ",
		  param->name,param->name);
  system(comm);
  free(comm);

  return 0;
  }
  
