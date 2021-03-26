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

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>


#include "parameter.h"
#include "cdim.h"
#include "common_blocks.h"
#include "do_logplot.h"
#include "do_tc.h"
#include "nlm.h"

#define streq(a,b) (!strcasecmp(a,b))

int do_logplot(param_t *param, char *logfile){

  int i,k;
  int scount=0;
  int pcount=0;
  int dcount=0;
  int fcount=0;
  int scount2=0;
  int pcount2=0;
  int dcount2=0;
  int fcount2=0;
  int lcolor=0;
  int lsty=0;
  int ncore;
  int donl=1;
  static char filename[160];

  double eae[10],enl[10],lae[10],lnl[10];
  int npot[10];
  int ntpot;
  int doifc=0;

  FILE *parm,*fp;
  char *comm;
  char lc = 0;
  

  #define comm_size 240
  comm= (char *) malloc(comm_size*sizeof(char));

  ncore=param->norb-param->nval;
  if (param->ilogder == -67) {
    fprintf(stderr," Cannot plot log. der. without setting Loginfo keyblock \n");
    return 1;
  } else {
    /*    param->ilogder = 0; */
  }

  if (param->ilogder > param->nconfigs) {
    fprintf(stderr," You want the log. deriv. of configuration %d, but there are only %d configurations! \n  ABORT \n",
	    param->ilogder,param->nconfigs);
    return 1;
  } else {
    do_tc(param,logfile,param->ilogder,doifc,donl);
  }
  sprintf(filename, "%s.logdeAE", param->name);
  fp = fopen(filename, "rb");
  for (i=0; i<param->nval; i++) {
    fread(&eae[i], sizeof(double), 1, fp);
    fread(&lae[i], sizeof(double), 1, fp);
  }
  fclose(fp);
  sprintf(filename, "%s.logdeNL", param->name);
  fp = fopen(filename, "rb");
  for (i=0; i<param->nval; i++) {
    fread(&enl[i], sizeof(double), 1, fp);
    fread(&lnl[i], sizeof(double), 1, fp);
  }
  fclose(fp);

  if ((parm = fopen("logd.par","w")) != NULL) {
    
    fprintf(parm,"# Log Deriv par file for xmgrace\n");
    fprintf(parm,"g0 on\n");
    fprintf(parm,"with g0\n");
    fprintf(parm,"title \"Log_derivs for %s at r\\slog\\N=%6.2f\"\n",param->symbol,param->rphas);
    fprintf(parm,"title font 0\n");
    fprintf(parm,"title size 1.500000\n");
    fprintf(parm,"title color 1\n");
    if (param->psmeth=='o') {
      fprintf(parm,"subtitle \"Optimized Pseudopotential Method\"\n");
    }else if (param->psmeth=='k') {
      fprintf(parm,"subtitle \"Kerker Pseudopotential Method\"\n");
    }else if (param->psmeth=='t') {
      fprintf(parm,"subtitle \"Troullier-Martins Pseudopotential Method\"\n");
    }
    fprintf(parm,"subtitle font 0\n");
    fprintf(parm,"subtitle size 1.000000\n");
    fprintf(parm,"subtitle color 1\n");
    fprintf(parm,"xaxis on\n");
    fprintf(parm,"xaxis label \"E (Ry)\"\n");
    fprintf(parm,"yaxis on\n");
    fprintf(parm,"yaxis label \" r * d[ln\\xy\\f{}]/dr \"\n");
    fprintf(parm,"world ymin -10\n");
    fprintf(parm,"world ymax 10\n");
    fprintf(parm,"yaxis tick major 10 \n");
    fprintf(parm,"yaxis tick minor 2 \n");
    fprintf(parm,"legend on\n");
    fprintf(parm,"legend loctype view\n");
    fprintf(parm,"legend 0.85, 0.8\n");
    
    for (i=0; i<param->nll;i++){
      
      if (nlm_label(param->nlm[i+ncore]).l == 0) {
	scount++;
	lcolor=1;
	lsty=1;
	lc='s';

      }else if (nlm_label(param->nlm[i+ncore]).l == 1) {
	pcount++;
	lcolor=2;
	lsty=1;
	lc='p';

      }else if (nlm_label(param->nlm[i+ncore]).l == 2) {
	dcount++;
	lcolor=3;
	lsty=1;
	lc='d';

      }else if (nlm_label(param->nlm[i+ncore]).l == 3) {
	fcount++;
	lcolor=4;
	lsty=1;
	lc='f';

      }

      fprintf(parm," s%d hidden false \n",i);
      fprintf(parm," s%d type xy \n",i);
      fprintf(parm," s%d symbol 0 \n",i);
      fprintf(parm," s%d line type 1 \n",i);
      fprintf(parm," s%d line linestyle %d \n",i,lsty);
      fprintf(parm," s%d line linewidth 2.0 \n",i);
      fprintf(parm," s%d line color %d \n",i,lcolor);
      fprintf(parm," s%d legend \"%c AE\" \n",i,lc);
    }
    
    k=i;
    
    for (i=0; i<param->nll;i++){
      
      if (nlm_label(param->nlm[i+ncore]).l == 0) {
	scount++;
	scount2++;
	lcolor=1;
	lsty=3;
	lc='s';

      }else if (nlm_label(param->nlm[i+ncore]).l == 1) {
	pcount++;
	pcount2++;
	lcolor=2;
	lsty=3;
	lc='p';

      }else if (nlm_label(param->nlm[i+ncore]).l == 2) {
	dcount++;
	dcount2++;
	lcolor=3;
	lsty=3;
	lc='d';

      }else if (nlm_label(param->nlm[i+ncore]).l == 3) {
	fcount++;
	fcount2++;
	lcolor=4;
	lsty=3;
	lc='f';
	
      }
      fprintf(parm," s%d hidden false \n",i+k);
      fprintf(parm," s%d type xy \n",i+k);
      fprintf(parm," s%d symbol 0 \n",i+k);
      fprintf(parm," s%d line type 1 \n",i+k);
      fprintf(parm," s%d line linestyle %d \n",i+k,lsty);
      fprintf(parm," s%d line linewidth 4.0 \n",i+k);
      fprintf(parm," s%d line color %d \n",i+k,lcolor);
      fprintf(parm," s%d legend \"%c NL\" \n",i+k,lc);
    }

    for (i=0; i<param->nval;i++){
      
      fprintf(parm,"with ellipse \n");
      fprintf(parm,"ellipse on \n");
      fprintf(parm,"ellipse loctype world \n");
      fprintf(parm,"ellipse g0 \n");
      fprintf(parm,"ellipse %lg , %lg , %lg , %lg \n",eae[i]-(param->emax-param->emin)*0.03*0.5,lae[i]-20.0*0.035*0.5,
	      eae[i]+(param->emax-param->emin)*0.03*0.5,lae[i]+20.0*0.035*0.5);

      fprintf(parm,"ellipse linestyle 1 \n");
      fprintf(parm,"ellipse linewidth 2.0 \n");
      fprintf(parm,"ellipse color %d \n",nlm_label(param->nlm[i+ncore]).l+1);
      fprintf(parm,"ellipse fill color %d\n",0);
      fprintf(parm,"ellipse fill pattern 4\n");
      fprintf(parm,"ellipse def \n");

    }
    for (i=0; i<param->nval;i++){
      
      fprintf(parm,"with ellipse \n");
      fprintf(parm,"ellipse on \n");
      fprintf(parm,"ellipse loctype world \n");
      fprintf(parm,"ellipse g0 \n");
      fprintf(parm,"ellipse %lg , %lg , %lg , %lg \n",enl[i]-(param->emax-param->emin)*0.03*0.5,lnl[i]-20.0*0.035*0.5,
	      enl[i]+(param->emax-param->emin)*0.03*0.5,lnl[i]+20.0*0.035*0.5);
      fprintf(parm,"ellipse linestyle 3 \n");
      fprintf(parm,"ellipse linewidth 2.0 \n");
      fprintf(parm,"ellipse color %d \n",nlm_label(param->nlm[i+ncore]).l+1);
      fprintf(parm,"ellipse fill color %d\n",nlm_label(param->nlm[i+ncore]).l+1);
      fprintf(parm,"ellipse fill pattern 0\n");
      fprintf(parm,"ellipse def \n");

    }
  }
  fclose(parm);
  
  snprintf(comm,comm_size,"xmgrace $XMGRACE_OPTS %s.logd_plt -p logd.par -autoscale xy -saveall %s.logd_agr & ",
	   param->name,param->name);
  system(comm);
  free(comm);
  
  return 0;
}

