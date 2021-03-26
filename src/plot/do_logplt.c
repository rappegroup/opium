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
 * $Id: do_logplt.c,v 1.13 2004/10/10 22:10:06 ewalter Exp $
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>


#include "parameter.h"
#include "fortparam.h"
#include "common_blocks.h"
#include "do_logplt.h"
#include "do_tc.h"
#include "nlm.h"

#define streq(a,b) (!strcasecmp(a,b))

int do_logplt(param_t *param, char *logfile){

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
  static char filename[160];

  double eae[10],enl[10];
  int npot[10];
  int ntpot;

  FILE *parm,*fp;
  char *comm;
  char lc = 0;

  comm= (char *) malloc(120*sizeof(char));

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
    do_tc(param,logfile,param->ilogder);
  }
  sprintf(filename, "%s.logdeAE", param->name);
  fp = fopen(filename, "rb");
  for (i=0; i<param->nll; i++)
    fread(&eae[i], sizeof(double), 1, fp);
  fclose(fp);
  sprintf(filename, "%s.logdeNL", param->name);
  fp = fopen(filename, "rb");
  for (i=0; i<param->nll; i++)
    fread(&enl[i], sizeof(double), 1, fp);
  fclose(fp);

  if ((parm = fopen("logd.par","w")) != NULL) {
    
    fprintf(parm,"# Log Deriv par file for xmgrace\n");
    fprintf(parm,"g0 on\n");
    fprintf(parm,"with g0\n");
    fprintf(parm,"title \"Log_derivs: %s\"\n",param->symbol);
    fprintf(parm,"title font 0\n");
    fprintf(parm,"title size 1.500000\n");
    fprintf(parm,"title color 1\n");
    fprintf(parm,"subtitle \"atom name: %s, radius=%6.2f\"\n",param->name,param->rphas);
    fprintf(parm,"subtitle font 0\n");
    fprintf(parm,"subtitle size 1.000000\n");
    fprintf(parm,"subtitle color 1\n");
    fprintf(parm,"xaxis on\n");
    fprintf(parm,"xaxis label \"E (Ry)\"\n");
    fprintf(parm,"yaxis on\n");
    fprintf(parm,"yaxis label \" r * d[log(psi)]/dr \"\n");
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
	lcolor=1+(scount-1)*4;
	lsty=1;
	lc='s';
	fprintf(parm,"with line \n");
	fprintf(parm,"line on \n");
	fprintf(parm,"line loctype world \n");
	fprintf(parm,"line g0 \n");
	fprintf(parm,"line %lg , %lg , %lg , %lg \n",eae[i],10.0,eae[i],-10.0);
	fprintf(parm,"line linewidth 2.0 \n");
	fprintf(parm,"line color %d\n",lcolor);
	fprintf(parm,"line def \n");

      }else if (nlm_label(param->nlm[i+ncore]).l == 1) {
	pcount++;
	lcolor=2+(pcount-1)*4;
	lsty=1;
	lc='p';
	fprintf(parm,"with line \n");
	fprintf(parm,"line on \n");
	fprintf(parm,"line loctype world \n");
	fprintf(parm,"line g0 \n");
	fprintf(parm,"line %lg , %lg , %lg , %lg \n",eae[i],10.0,eae[i],-10.0);
	fprintf(parm,"line linewidth 2.0 \n");
	fprintf(parm,"line color %d\n",lcolor);
	fprintf(parm,"line def \n");

      }else if (nlm_label(param->nlm[i+ncore]).l == 2) {
	dcount++;
	lcolor=3+(dcount-1)*4;
	lsty=1;
	lc='d';
	fprintf(parm,"with line \n");
	fprintf(parm,"line on \n");
	fprintf(parm,"line loctype world \n");
	fprintf(parm,"line g0 \n");
	fprintf(parm,"line %lg , %lg , %lg , %lg \n",eae[i],10.0,eae[i],-10.0);
	fprintf(parm,"line linewidth 2.0 \n");
	fprintf(parm,"line color %d\n",lcolor);
	fprintf(parm,"line def \n");

      }else if (nlm_label(param->nlm[i+ncore]).l == 3) {
	fcount++;
	lcolor=4+(fcount-1)*4;
	lsty=1;
	lc='f';
	fprintf(parm,"with line \n");
	fprintf(parm,"line on \n");
	fprintf(parm,"line loctype world \n");
	fprintf(parm,"line g0 \n");
	fprintf(parm,"line %lg , %lg , %lg , %lg \n",eae[i],10.0,eae[i],-10.0);
	fprintf(parm,"line linewidth 2.0 \n");
	fprintf(parm,"line color %d\n",lcolor);
	fprintf(parm,"line def \n");

      }

      fprintf(parm," s%d hidden false \n",i);
      fprintf(parm," s%d type xy \n",i);
      fprintf(parm," s%d symbol 0 \n",i);
      fprintf(parm," s%d line type 1 \n",i);
      fprintf(parm," s%d line linestyle %d \n",i,lsty);
      fprintf(parm," s%d line linewidth 2.0 \n",i);
      fprintf(parm," s%d line color %d \n",i,lcolor);
      fprintf(parm," s%d legend \"%d%c AE\" \n",i,
	      nlm_label(param->nlm[i+ncore]).n,lc);
    }
    
    k=i;
    
    for (i=0; i<param->nll;i++){
      
      if (nlm_label(param->nlm[i+ncore]).l == 0) {
	scount++;
	scount2++;
	lcolor=1+(scount2-1)*4;
	lsty=3;
	lc='s';
	fprintf(parm,"with line \n");
	fprintf(parm,"line on \n");
	fprintf(parm,"line loctype world \n");
	fprintf(parm,"line g0 \n");
	fprintf(parm,"line %lg , %lg , %lg , %lg \n",enl[i],10.0,enl[i],-10.0);
	fprintf(parm,"line linewidth 2.0 \n");
	fprintf(parm,"line linestyle %d\n",2);
	fprintf(parm,"line color %d\n",lcolor);
	fprintf(parm,"line def \n");
      }else if (nlm_label(param->nlm[i+ncore]).l == 1) {
	pcount++;
	pcount2++;
	lcolor=2+(pcount2-1)*4;
	lsty=3;
	lc='p';
	fprintf(parm,"with line \n");
	fprintf(parm,"line on \n");
	fprintf(parm,"line loctype world \n");
	fprintf(parm,"line g0 \n");
	fprintf(parm,"line %lg , %lg , %lg , %lg \n",enl[i],10.0,enl[i],-10.0);
	fprintf(parm,"line linewidth 2.0 \n");
	fprintf(parm,"line linestyle %d\n",2);
	fprintf(parm,"line color %d\n",lcolor);
	fprintf(parm,"line def \n");
      }else if (nlm_label(param->nlm[i+ncore]).l == 2) {
	dcount++;
	dcount2++;
	lcolor=3+(dcount2-1)*4;
	lsty=3;
	lc='d';
	fprintf(parm,"with line \n");
	fprintf(parm,"line on \n");
	fprintf(parm,"line loctype world \n");
	fprintf(parm,"line g0 \n");
	fprintf(parm,"line %lg , %lg , %lg , %lg \n",enl[i],10.0,enl[i],-10.0);
	fprintf(parm,"line linewidth 2.0 \n");
	fprintf(parm,"line linestyle %d\n",2);
	fprintf(parm,"line color %d\n",lcolor);
	fprintf(parm,"line def \n");
      }else if (nlm_label(param->nlm[i+ncore]).l == 3) {
	fcount++;
	fcount2++;
	lcolor=4+(fcount2-1)*4;
	lsty=3;
	lc='f';
	fprintf(parm,"with line \n");
	fprintf(parm,"line on \n");
	fprintf(parm,"line loctype world \n");
	fprintf(parm,"line g0 \n");
	fprintf(parm,"line %lg , %lg , %lg , %lg \n",enl[i],10.0,enl[i],-10.0);
	fprintf(parm,"line linewidth 2.0 \n");
	fprintf(parm,"line linestyle %d\n",2);
	fprintf(parm,"line color %d\n",lcolor);
	fprintf(parm,"line def \n");
      }
      fprintf(parm," s%d hidden false \n",i+k);
      fprintf(parm," s%d type xy \n",i+k);
      fprintf(parm," s%d symbol 0 \n",i+k);
      fprintf(parm," s%d line type 1 \n",i+k);
      fprintf(parm," s%d line linestyle %d \n",i+k,lsty);
      fprintf(parm," s%d line linewidth 4.0 \n",i+k);
      fprintf(parm," s%d line color %d \n",i+k,lcolor);
      fprintf(parm," s%d legend \"%d%c NL\" \n",i+k,
	      nlm_label(param->nlm[i+ncore]).n,lc);
    }
  }
  fclose(parm);
  
  sprintf(comm, "xmgrace %s.plt_logd -p logd.par -autoscale xy -saveall %s.logd_agr & ",
		  param->name,param->name);
  system(comm);
  free(comm);
  
  return 0;
}

