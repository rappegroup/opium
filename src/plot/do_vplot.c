/*
 * $Id: do_vplot.c,v 1.6 2004/06/16 21:25:54 mbarnes Exp $
 */

/* screened pseudo potential plotting with xmgrace */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>


#include "parameter.h"
#include "fortparam.h"
#include "common_blocks.h"
#include "do_vplot.h"
#include "do_ps.h"
#include "nlm.h"

#define streq(a,b) (*a==*b && !strcmp(a+1,b+1))

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

    parm = fopen("sps.par","w");

    if (!parm) {
      printf("Could not open grace param file 'sps.par': %s\n",
	  strerror(errno));
      return 1;
    }

    fprintf(parm,"# screened potential param file for xmgrace\n");
    fprintf(parm,"g0 on\n");
    fprintf(parm,"with g0\n");
    fprintf(parm,"world xmin 0\n");
    fprintf(parm,"world xmax 5\n");
    fprintf(parm,"title \"Screened pseudopotential: %s\"\n",param->symbol);
    fprintf(parm,"world ymin -50\n");
    fprintf(parm,"world ymax 50\n");
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
    fprintf(parm,"xaxis label \"R (Bohr)\"\n");
    fprintf(parm,"yaxis on\n");
    fprintf(parm,"yaxis tick major 10 \n");
    fprintf(parm,"yaxis tick minor 2 \n");
    fprintf(parm,"yaxis label \"V_scr (Ry)\"\n");
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
      fprintf(parm," s%d legend \"%d%c\" \n",i,
	  nlm_label(param->nlm[param->ipot[i]+ncore]).n,lc);
    }

    fprintf(parm," s%d hidden false \n",i+1);
    fprintf(parm," s%d type xy \n",i+1);
    fprintf(parm," s%d symbol 0 \n",i+1);
    fprintf(parm," s%d line type 1 \n",i+1);
    fprintf(parm," s%d line linestyle %d \n",i+1,3);
    fprintf(parm," s%d line linewidth 3.0 \n",i+1);
    fprintf(parm," s%d line color %d \n",i+1,14);
    fprintf(parm," s%d legend \"V_loc\"\n" ,i+1);

    fclose(parm);

    sprintf(comm, "xmgrace %s.plt_sps -p sps.par -autoscale y -saveall %s_sps.agr & ",
	param->name,param->name);
    system(comm);

  } else {

    parm = fopen("ips.par","w");

    if (!parm) {
      printf("Could not open grace param file 'sps.par': %s\n",
	  strerror(errno));
      return 1;
    }

    fprintf(parm,"# IPS par file for xmgrace\n");
    fprintf(parm,"g0 on\n");
    fprintf(parm,"with g0\n");
    fprintf(parm,"    world xmin 0\n");
    fprintf(parm,"    world xmax 5\n");
    fprintf(parm,"    world ymin -50\n");
    fprintf(parm,"    world ymax 50\n");
    fprintf(parm,"    title \"Ionic pseudopotential: %s\"\n",param->symbol);
    fprintf(parm,"    title font 0\n");
    fprintf(parm,"    title size 1.500000\n");
    fprintf(parm,"    title color 1\n");
    fprintf(parm,"    subtitle \"atom name: %s\"\n",param->name);
    fprintf(parm,"    subtitle font 0\n");
    fprintf(parm,"    subtitle size 1.000000\n");
    fprintf(parm,"    subtitle color 1\n");
    fprintf(parm,"    xaxis  on\n");
    fprintf(parm,"    xaxis  tick major 1\n");
    fprintf(parm,"    xaxis  tick minor 0.5\n");
    fprintf(parm,"    xaxis  label \"R (Bohr)\"\n");
    fprintf(parm,"    yaxis  on\n");
    fprintf(parm,"    yaxis  tick major 10 \n");
    fprintf(parm,"    yaxis  tick minor 2 \n");
    fprintf(parm,"    yaxis  label \"V_ion (Ryd)\"\n");
    fprintf(parm,"    legend on\n");
    fprintf(parm,"    legend loctype view\n");
    fprintf(parm,"    legend 0.85, 0.8\n");

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
      fprintf(parm," s%d legend \"%d%c\" \n",i,
	  nlm_label(param->nlm[param->ipot[i]+ncore]).n,lc);
    }
    /* add the AE ionic potential */

    fprintf(parm," s%d hidden false \n",i+1);
    fprintf(parm," s%d type xy \n",i+1);
    fprintf(parm," s%d symbol 0 \n",i+1);
    fprintf(parm," s%d line type 1 \n",i+1);
    fprintf(parm," s%d line linestyle %d \n",i+1,3);
    fprintf(parm," s%d line linewidth 2.0 \n",i+1);
    fprintf(parm," s%d line color %d \n",i+1,15);
    fprintf(parm," s%d legend \"2Z_eff/r\"\n" ,i+1);

    fprintf(parm," s%d hidden false \n",i+2);
    fprintf(parm," s%d type xy \n",i+2);
    fprintf(parm," s%d symbol 0 \n",i+2);
    fprintf(parm," s%d line type 1 \n",i+2);
    fprintf(parm," s%d line linestyle %d \n",i+2,3);
    fprintf(parm," s%d line linewidth 3.0 \n",i+2);
    fprintf(parm," s%d line color %d \n",i+2,14);
    fprintf(parm," s%d legend \"V_loc\"\n" ,i+2);

    fclose(parm);

    sprintf(comm, "xmgrace %s.plt_ips -p ips.par -autoscale y -saveall %s_ips.agr & ",
	param->name,param->name);

    system(comm);
  }

  free(comm);
  return 0;
}

/* vim: cindent sw=2 showmatch
 */
