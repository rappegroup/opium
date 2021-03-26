/*
 * $Id: do_logplt.c,v 1.6 2004/06/16 21:25:54 mbarnes Exp $
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

#define streq(a,b) (*a==*b && !strcmp(a+1,b+1))

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

  int npot[10];
  int ntpot;

  FILE *parm;
  char *comm;
  char lc = 0;

  comm= (char *) malloc(120*sizeof(char));

  /*  if ((streq(param->reltype, "frl")||(streq(param->reltype, "srl")))){
      fprintf(fp_log,"No log. deriv. plotting for srl or frl modes yet :( ");
      return 1;
      }  
      if (param->ilogder != 1) {
      fprintf(fp_log,"No log. derivs. computed, check the [LogInfo] keyblock ");
      return 1;
      }*/


  for (i=0; i<10; i++){
    npot[i]=0;
  }

  ntpot=0;

  for (i=0; i<param->nval; i++){
    npot[nlm_label(atomic_.nlm[i]).l]++;
    if (npot[nlm_label(atomic_.nlm[i]).l] == 1) {
      ntpot++;
    }
  }
  param->nll=ntpot;

  do_tc(param,logfile,param->ilogder);
  
  ncore=param->norb-param->nval;
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
    fprintf(parm,"world ymin -50\n");
    fprintf(parm,"world ymax 50\n");
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
      }else if (nlm_label(param->nlm[i+ncore]).l == 1) {
	pcount++;
	lcolor=2+(pcount-1)*4;
	lsty=1;
	lc='p';
      }else if (nlm_label(param->nlm[i+ncore]).l == 2) {
	dcount++;
	lcolor=3+(dcount-1)*4;
	lsty=1;
	lc='d';
      }else if (nlm_label(param->nlm[i+ncore]).l == 3) {
	fcount++;
	lcolor=4+(fcount-1)*4;
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
      }else if (nlm_label(param->nlm[i+ncore]).l == 1) {
	pcount++;
	pcount2++;
	lcolor=2+(pcount2-1)*4;
	lsty=3;
	lc='p';
      }else if (nlm_label(param->nlm[i+ncore]).l == 2) {
	dcount++;
	dcount2++;
	lcolor=3+(dcount2-1)*4;
	lsty=3;
	lc='d';
      }else if (nlm_label(param->nlm[i+ncore]).l == 3) {
	fcount++;
	fcount2++;
	lcolor=4+(fcount2-1)*4;
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

