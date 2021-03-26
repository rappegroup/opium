/*
 * Copyright (c) 1998-2005 The OPIUM Group
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
 * $Id: do_wplot.c,v 1.8 2004/10/02 18:34:49 ewalter Exp $
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "parameter.h"        /* defines structure: 'param_t' */
#include "cdim.h"
#include "common_blocks.h"
#include "do_wplot.h"           /* the module's own header */
#include "nlm.h"

#define streq(a,b) (!strcasecmp(a,b))

int do_wplot(param_t *param, char *logfile, char *pltyp){

  int i,k=0;
  int scount=0;
  int pcount=0;
  int dcount=0;
  int fcount=0;
  int lcolor=0;
  int lsty=0;
  int ncore;

  FILE *parm;
  char *comm;
  char lc=0,sc=0;
  
  #define comm_size 240
  comm= (char *) malloc(comm_size*sizeof(char));

  ncore=param->norb-param->nval;
  
  if ((parm = fopen("w.par","w")) != NULL) {
      
    fprintf(parm,"# wavefunction param file for xmgrace\n");
    
    if (streq(pltyp,"n")) {
      fprintf(parm,"title \"Valence wavefunctions for %s\"\n",param->symbol);
      if (param->psmeth=='o') {
	fprintf(parm,"subtitle \"Optimized Pseudopotential Method\"\n");
      }else if (param->psmeth=='k') {
	fprintf(parm,"subtitle \"Kerker Pseudopotential Method\"\n");
      }else if (param->psmeth=='t') {
	fprintf(parm,"subtitle \"Troullier-Martins Pseudopotential Method\"\n");
      }
    } else {
      fprintf(parm,"title \"All electron wavefunctions for %s\"\n",param->symbol);
    }
    fprintf(parm,"g0 on\n");
    fprintf(parm,"with g0\n");
    fprintf(parm,"world xmin 0\n");
    fprintf(parm,"world xmax 5\n");
    fprintf(parm,"title font 0\n");
    fprintf(parm,"title size 1.500000\n");
    fprintf(parm,"title color 1\n");
    fprintf(parm,"subtitle font 0\n");
    fprintf(parm,"subtitle size 1.000000\n");
    fprintf(parm,"subtitle color 1\n");
    fprintf(parm,"xaxis  on\n");
    fprintf(parm,"xaxis  tick major 1\n");
    fprintf(parm,"xaxis  tick minor 0.5\n");
    fprintf(parm,"xaxis  label \"r (a.u.)\"\n");
    fprintf(parm,"yaxis  on\n");
    fprintf(parm,"yaxis  label \"\\+r \\xy\\f{}(r)\"\n");
    fprintf(parm,"legend on\n");
    fprintf(parm,"legend loctype view\n");
    fprintf(parm,"legend 0.85, 0.8\n");
    
    if (streq(param->reltype, "nrl")){  
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
	fprintf(parm," s%d line linewidth 1.5 \n",i);
	fprintf(parm," s%d line color %d \n",i,lcolor);
	fprintf(parm," s%d legend \"\\xy\\f{}\\s%d%c\\N\\S\\-AE\\N\" \n",i,
		nlm_label(param->nlm[i+ncore]).n,lc);
      }
      k=i;
    } else {
      
      for (i=ncore; i<param->norb;i++){
	
	if (nlm_label(param->nlm[i]).l == 0) {
	  scount++;
	  lcolor=1;
	  lsty=scount;
	  lc='s';
	}else if (nlm_label(param->nlm[i]).l == 1) {
	  pcount++;
	  lcolor=2;
	  lsty=pcount;
	  lc='p';
	}else if (nlm_label(param->nlm[i]).l == 2) {
	  dcount++;
	  lcolor=3;
	  lsty=dcount;
	  lc='d';
	}else if (nlm_label(param->nlm[i]).l == 3) {
	  fcount++;
	  lcolor=4;
	  lsty=fcount;
	  lc='f';
	}
	k=i-ncore;
	fprintf(parm," s%d hidden false \n",k);
	fprintf(parm," s%d type xy \n",k);
	fprintf(parm," s%d symbol 0 \n",k);
	fprintf(parm," s%d line type 1 \n",k);
	fprintf(parm," s%d line linestyle %d \n",k,lsty+3);
	fprintf(parm," s%d line linewidth 5.0 \n",k);
	fprintf(parm," s%d line color %d \n",k,lcolor);
	fprintf(parm," s%d legend \"\\xy\\f{}\\s%d%c\\S\\-avg\\N\" \n",k,
		nlm_label(param->nlm[i]).n,lc);
	
      }
      k++;
      for (i=0; i<atmwave_.numorb;i++){
	
	if (atmwave_.llo[i] == 0) {
	  scount++;
	  lcolor=1;
	  lsty=scount;
	  lc='s';
	  sc=' ';
	} else if (atmwave_.llo[i] == 1) {
	  pcount++;
	  lcolor=2;
	  lsty=pcount;
	  lc='p';
	  sc='-';
	  if (atmwave_.sso[i] < 0.0) sc='+';
	} else if (atmwave_.llo[i] == 2) {	
	  dcount++;
	  lcolor=3;
	  lsty=dcount;
	  lc='d';
	  sc='-';
	  if (atmwave_.sso[i] < 0.0) sc='+';
	}else if (atmwave_.llo[i] == 3) {	
	  fcount++;
	  lcolor=4;
	  lsty=fcount;
	  lc='f';
	  sc='-';
	  if (atmwave_.sso[i] < 0.0) sc='+';
	}
	fprintf(parm," s%d hidden false \n",i+k);
	fprintf(parm," s%d type xy \n",i+k);
	fprintf(parm," s%d symbol 0 \n",i+k);
	fprintf(parm," s%d line type 1 \n",i+k);
	fprintf(parm," s%d line linestyle %d \n",i+k,lsty-1);
	fprintf(parm," s%d line linewidth 1.0 \n",i+k);
	fprintf(parm," s%d line color %d \n",i+k,lcolor);
	fprintf(parm," s%d legend \"\\xy\\f{}\\s%d%c%c\\N\\S\\-AE\\N\" \n",i+k,atmwave_.nno[i],lc,sc);
      }

      k=i+k;
    }
    
    if (streq(pltyp,"n")) {

      
      for (i=0; i<param->nval;i++) {
	
	if (nlm_label(param->nlm[i+ncore]).l == 0) {
	  scount++;
	  lcolor=1;
	  lsty=scount+1;
	  lc='s';
	  /*	  fprintf(parm,"with line \n");
	  fprintf(parm,"line on \n");
	  fprintf(parm,"line loctype world \n");
	  fprintf(parm,"line g0 \n");
	  fprintf(parm,"line %lg , %lg , %lg , %lg \n",param->rc[i],1.0,param->rc[i],0.0);
	  fprintf(parm,"line linewidth 1.0 \n");
	  fprintf(parm,"line color %d\n",lcolor);
	  fprintf(parm,"line def \n");*/

	}else if (nlm_label(param->nlm[i+ncore]).l == 1) {
	  pcount++;
	  lcolor=2;
	  lsty=pcount+1;
	  lc='p';
	  /*fprintf(parm,"with line \n");
	  fprintf(parm,"line on \n");
	  fprintf(parm,"line loctype world \n");
	  fprintf(parm,"line g0 \n");
	  fprintf(parm,"line %lg , %lg , %lg , %lg \n",param->rc[i],1.0,param->rc[i],0.0);
	  fprintf(parm,"line linewidth 2.0 \n");
	  fprintf(parm,"line color %d\n",lcolor);
	  fprintf(parm,"line def \n");*/

	}else if (nlm_label(param->nlm[i+ncore]).l == 2) {
	  dcount++;
	  lcolor=3;
	  lsty=dcount+1;
	  lc='d';
	  /*fprintf(parm,"with line \n");
	  fprintf(parm,"line on \n");
	  fprintf(parm,"line loctype world \n");
	  fprintf(parm,"line g0 \n");
	  fprintf(parm,"line %lg , %lg , %lg , %lg \n",param->rc[i],1.0,param->rc[i],0.0);
	  fprintf(parm,"line linewidth 3.0 \n");
	  fprintf(parm,"line color %d\n",lcolor);
	  fprintf(parm,"line def \n");*/

	}else if (nlm_label(param->nlm[i+ncore]).l == 3) {
	  fcount++;
	  lcolor=4;
	  lsty=fcount+1;
	  lc='f';
	  /*fprintf(parm,"with line \n");
	  fprintf(parm,"line on \n");
	  fprintf(parm,"line loctype world \n");
	  fprintf(parm,"line g0 \n");
	  fprintf(parm,"line %lg , %lg , %lg , %lg \n",param->rc[i],1.0,param->rc[i],0.0);
	  fprintf(parm,"line linewidth 3.5 \n");
	  fprintf(parm,"line color %d\n",lcolor);
	  fprintf(parm,"line def \n");*/

	}
	
	fprintf(parm," s%d hidden false \n",k+i);
	fprintf(parm," s%d type xydx \n",k+i);
	fprintf(parm," s%d symbol 0 \n",k+i);
	fprintf(parm," s%d line type 1 \n",k+i);
	fprintf(parm," s%d line linestyle %d\n",k+i,lsty);
	fprintf(parm," s%d line linewidth 2.5 \n",k+i);
	fprintf(parm," s%d line color %d \n",k+i,lcolor);
	fprintf(parm," s%d legend \"\\xy\\f{}\\s%d%c\\N\\S\\-NL\\N\" \n",i+k,nlm_label(param->nlm[i+ncore]).n,lc);
	fprintf(parm," s%d errorbar on\n",k+i);
	fprintf(parm," s%d errorbar place both\n",k+i);
	fprintf(parm," s%d errorbar color %d\n",k+i,lcolor);
	fprintf(parm," s%d errorbar pattern 1\n",k+i);
	fprintf(parm," s%d errorbar size 10.000000\n",k+i);
	fprintf(parm," s%d errorbar linewidth 2.0\n",k+i);
	fprintf(parm," s%d errorbar linestyle 1\n",k+i);
      }
    }
    
    fclose(parm);
      
  } 

  if (streq(pltyp,"n")) {    
    snprintf(comm, comm_size,"xmgrace -timestamp $XMGRACE_OPTS %s.ae_plt -settype xydx %s.nl_plt -p w.par -autoscale y -saveall %s_nl.agr & ", 
	    param->name,param->name,param->name);
  } else {

    snprintf(comm, comm_size,"xmgrace -timestamp $XMGRACE_OPTS %s.ae_plt -p w.par -autoscale y -saveall %s_ae.agr & ", 
	    param->name,param->name);
  }

  system(comm);
  free(comm);

  return 0;
}

