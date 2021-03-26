/*
 * $Id: do_spsplt.c,v 1.3 2004/06/16 21:25:54 mbarnes Exp $
 */


/* screened pseudo potential plotting with xmgrace */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>


#include "parameter.h"        /* defines structure: 'param_t' */
#include "fortparam.h"
#include "common_blocks.h"
#include "do_spsplt.h"           /* the module's own header */
#include "nlm.h"

int do_spsplt(param_t *param, char *logfile){

  int i,j,k;
  int scount=0;
  int pcount=0;
  int dcount=0;
  int fcount=0;
  int lcolor=0;
  int lsty=0;
  int nc;
  int ncore;

  FILE *parm;
  char *comm;
  char lc,sc;

  comm= (char *) malloc(120*sizeof(char));

  ncore=param->norb-param->nval;
  if ((parm = fopen("sps.par","w")) != NULL) {
    
    fprintf(parm,"# SPS par file for xmgrace\n");
    fprintf(parm,"g0 on\n");
    fprintf(parm,"with g0\n");
    fprintf(parm,"    world xmin 0\n");
    fprintf(parm,"    world xmax 5\n");
    fprintf(parm,"    world ymin -50\n");
    fprintf(parm,"    world ymax 0\n");
    fprintf(parm,"    title \"Screened pseudopotential: %s\"\n",param->symbol);
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
    fprintf(parm,"    yaxis  tick major 5 \n");
    fprintf(parm,"    yaxis  tick minor 1 \n");
    fprintf(parm,"    yaxis  label \"r * Psi\"\n");
    fprintf(parm,"    legend on\n");
    fprintf(parm,"    legend loctype view\n");
    fprintf(parm,"    legend 0.85, 0.8\n");

    for (i=0; i<param->nll;i++){
      
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
      k=i;
      fprintf(parm," s%d hidden false \n",k);
      fprintf(parm," s%d type xy \n",k);
      fprintf(parm," s%d symbol 0 \n",k);
      fprintf(parm," s%d line type 1 \n",k);
      fprintf(parm," s%d line linestyle %d \n",k,lsty);
      fprintf(parm," s%d line linewidth 2.0 \n",k);
      fprintf(parm," s%d line color %d \n",k,lcolor);
      fprintf(parm," s%d legend \"%d%c\" \n",k,
	      nlm_label(param->nlm[i+ncore]).n,lc);
    }

    fprintf(parm," s%d hidden false \n",k+1);
    fprintf(parm," s%d type xy \n",k+1);
    fprintf(parm," s%d symbol 0 \n",k+1);
    fprintf(parm," s%d line type 1 \n",k+1);
    fprintf(parm," s%d line linestyle %d \n",k+1,3);
    fprintf(parm," s%d line linewidth 3.0 \n",k+1);
    fprintf(parm," s%d line color %d \n",k+1,14);
    fprintf(parm," s%d legend \"V_loc\"\n" ,k+1);
    

    fclose(parm);

  } 
    
  sprintf(comm, "xmgrace %s.plt_sps -p sps.par -saveall %s_sps.agr & ", param->name,param->name);
  system(comm);
  free(comm);
  
  return 0;
}

