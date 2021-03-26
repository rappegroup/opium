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
 * $Id: do_vplot.c,v 1.9 2004/10/02 18:34:49 ewalter Exp $
 */

/* screened pseudo potential plotting with xmgrace */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>


#include "parameter.h"
#include "cdim.h"
#include "common_blocks.h"
#include "do_vplot.h"
#include "do_ps.h"
#include "nlm.h"

#define streq(a,b) (!strcasecmp(a,b))

int do_vplot(param_t *param, char *logfile, char *pltyp){

  int i;
  int scount=0;
  int pcount=0;
  int dcount=0;
  int fcount=0;
  int lcolor=0;
  int lsty=0;
  int ncore;

  FILE *parm;
  char *comm;
  char lc=0;

  comm= (char *) malloc(120*sizeof(char));

  ncore=param->norb-param->nval;

  readAE(param);

  /* check for semicore */

  if (streq(pltyp,"s")) {

    parm = fopen("vs.par","w");

    if (!parm) {
      printf("Could not open grace param file 'vs.par': %s\n",
	  strerror(errno));
      return 1;
    }

    fprintf(parm,"# screened potential param file for xmgrace\n");
    fprintf(parm,"g0 on\n");
    fprintf(parm,"with g0\n");
    fprintf(parm,"world xmin 0\n");
    fprintf(parm,"world xmax 5\n");
    fprintf(parm,"title \"Screened pseudopotential for %s\"\n",param->symbol);
    fprintf(parm,"world ymin -50\n");
    fprintf(parm,"world ymax 50\n");
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
    fprintf(parm,"xaxis tick major 1\n");
    fprintf(parm,"xaxis tick minor 0.5\n");
    fprintf(parm,"xaxis label \"r (a.u.)\"\n");
    fprintf(parm,"yaxis on\n");
    fprintf(parm,"yaxis tick major 10 \n");
    fprintf(parm,"yaxis tick minor 2 \n");
    fprintf(parm,"yaxis label \"V\\sscr\\N (Ry)\"\n");
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
      fprintf(parm," s%d type xydx \n",i);
      fprintf(parm," s%d symbol 0 \n",i);
      fprintf(parm," s%d line type 1 \n",i);
      fprintf(parm," s%d line linestyle %d \n",i,lsty);
      fprintf(parm," s%d line linewidth 2.0 \n",i);
      fprintf(parm," s%d line color %d \n",i,lcolor);
      fprintf(parm," s%d legend \"V\\s%d%c\\N\" \n",i,
	  nlm_label(param->nlm[param->ipot[i]+ncore]).n,lc);
      fprintf(parm," s%d errorbar on\n",i);
      fprintf(parm," s%d errorbar place both\n",i);
      fprintf(parm," s%d errorbar color %d\n",i,lcolor);
      fprintf(parm," s%d errorbar pattern 1\n",i);
      fprintf(parm," s%d errorbar size 10.000000\n",i);
      fprintf(parm," s%d errorbar linewidth 2.0\n",i);
      fprintf(parm," s%d errorbar linestyle 1\n",i);
    }
    i=param->nll;
    fprintf(parm," s%d hidden false \n",i);
    fprintf(parm," s%d type xydx \n",i);
    fprintf(parm," s%d symbol 0 \n",i);
    fprintf(parm," s%d line type 1 \n",i);
    fprintf(parm," s%d line linestyle %d \n",i,3);
    fprintf(parm," s%d line linewidth 3.0 \n",i);
    fprintf(parm," s%d line color %d \n",i,14);
    fprintf(parm," s%d legend \"V\\sloc\\N\"\n" ,i);

    fclose(parm);

    sprintf(comm, "xmgrace -timestamp $XMG_OPTS -settype xydx %s.vs_plt -p vs.par -autoscale y -saveall %s_vs.agr & ",
	param->name,param->name);
    system(comm);

  } else {

    parm = fopen("vi.par","w");

    if (!parm) {
      printf("Could not open grace param file 'vi.par': %s\n",
	  strerror(errno));
      return 1;
    }

    fprintf(parm,"# IPS par file for xmgrace\n");
    fprintf(parm,"g0 on\n");
    fprintf(parm,"with g0\n");
    fprintf(parm,"world xmin 0\n");
    fprintf(parm,"world xmax 5\n");
    fprintf(parm,"world ymin -50\n");
    fprintf(parm,"world ymax 50\n");
    fprintf(parm,"title \"Ionic pseudopotential for %s\"\n",param->symbol);
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
    fprintf(parm,"xaxis  on\n");
    fprintf(parm,"xaxis  tick major 1\n");
    fprintf(parm,"xaxis  tick minor 0.5\n");
    fprintf(parm,"xaxis  label \"r (a.u.)\"\n");
    fprintf(parm,"yaxis  on\n");
    fprintf(parm,"yaxis  tick major 10 \n");
    fprintf(parm,"yaxis  tick minor 2 \n");
    fprintf(parm,"yaxis  label \"V\\sion\\N (Ry)\"\n");
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
      fprintf(parm," s%d type xydx \n",i);
      fprintf(parm," s%d symbol 0 \n",i);
      fprintf(parm," s%d line type 1 \n",i);
      fprintf(parm," s%d line linestyle %d \n",i,lsty);
      fprintf(parm," s%d line linewidth 2.0 \n",i);
      fprintf(parm," s%d line color %d \n",i,lcolor);
      fprintf(parm," s%d legend \"V\\s%d%c\\N\" \n",i,
	  nlm_label(param->nlm[param->ipot[i]+ncore]).n,lc);
      fprintf(parm," s%d errorbar on\n",i);
      fprintf(parm," s%d errorbar place both\n",i);
      fprintf(parm," s%d errorbar color %d\n",i,lcolor);
      fprintf(parm," s%d errorbar pattern 1\n",i);
      fprintf(parm," s%d errorbar size 10.000000\n",i);
      fprintf(parm," s%d errorbar linewidth 2.0\n",i);
      fprintf(parm," s%d errorbar linestyle 1\n",i);
    }
    /* add the AE ionic potential */
    i=param->nll;
    fprintf(parm," s%d hidden false \n",i);
    fprintf(parm," s%d type xydx \n",i);
    fprintf(parm," s%d symbol 0 \n",i);
    fprintf(parm," s%d line type 1 \n",i);
    fprintf(parm," s%d line linestyle %d \n",i,3);
    fprintf(parm," s%d line linewidth 2.0 \n",i);
    fprintf(parm," s%d line color %d \n",i,15);
    fprintf(parm," s%d legend \"2Z\\seff\\N /r\"\n" ,i);

    fprintf(parm," s%d hidden false \n",i+1);
    fprintf(parm," s%d type xydx \n",i+1);
    fprintf(parm," s%d symbol 0 \n",i+1);
    fprintf(parm," s%d line type 1 \n",i+1);
    fprintf(parm," s%d line linestyle %d \n",i+1,3);
    fprintf(parm," s%d line linewidth 3.0 \n",i+1);
    fprintf(parm," s%d line color %d \n",i+1,14);
    fprintf(parm," s%d legend \"V\\sloc\\N\"\n" ,i+1);

    fclose(parm);

    sprintf(comm, "xmgrace -timestamp $XMGRACE_OPTS -settype xydx %s.vi_plt -p vi.par -autoscale y -saveall %s_vi.agr & ",
	param->name,param->name);

    system(comm);
  }

  free(comm);
  return 0;
}

/* vim: cindent sw=2 showmatch
 */
