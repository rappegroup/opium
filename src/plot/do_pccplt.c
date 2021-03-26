/*
 * $Id: do_pccplt.c,v 1.5 2004/06/16 21:25:54 mbarnes Exp $
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>


#include "parameter.h"        /* defines structure: 'param_t' */
#include "fortparam.h"
#include "common_blocks.h"
#include "do_pccplt.h"
#include "nlm.h"

int do_pccplt(param_t *param, char *logfile){

  int ncore;

  FILE *parm;
  char *comm;

  comm= (char *) malloc(120*sizeof(char));

  ncore=param->norb-param->nval;

  if ((parm = fopen("pcc.par","w")) != NULL) {
    
    fprintf(parm,"# PCC par file for xmgrace\n");
    fprintf(parm,"g0 on\n");
    fprintf(parm,"with g0\n");
    fprintf(parm,"world xmin 0\n");
    fprintf(parm,"world xmax 5\n");
    fprintf(parm,"title \"core density: %s\"\n",param->symbol);
    fprintf(parm,"title font 0\n");
    fprintf(parm,"title size 1.500000\n");
    fprintf(parm,"title color 1\n");
    fprintf(parm,"subtitle \"atom name: %s\"\n",param->name);
    fprintf(parm,"subtitle font 0\n");
    fprintf(parm,"subtitle size 1.000000\n");
    fprintf(parm,"subtitle color 1\n");
    fprintf(parm,"xaxis  on\n");
    fprintf(parm,"xaxis  tick major 1\n");
    fprintf(parm,"xaxis  tick minor 0.5\n");
    fprintf(parm,"xaxis  label \"r (Bohr)\"\n");
    fprintf(parm,"yaxis  on\n");
    fprintf(parm,"yaxis  label \"4*pi*r^2 rho\"\n");
    fprintf(parm,"legend on\n");
    fprintf(parm,"legend loctype view\n");
    fprintf(parm,"legend 0.85, 0.8\n");

    /* full core */
    
    fprintf(parm," s0 hidden false \n");
    fprintf(parm," s0 type xy \n");
    fprintf(parm," s0 symbol 0 \n");
    fprintf(parm," s0 line type 1 \n");
    fprintf(parm," s0 line linestyle 1 \n");
    fprintf(parm," s0 line linewidth 2.0 \n");
    fprintf(parm," s0 line color 1 \n");
    fprintf(parm," s0 legend \" core density \" \n");


    if (param->rpcc > 1e-12) {
      /* partial core */
      
      fprintf(parm," s1 hidden false \n");
      fprintf(parm," s1 type xy \n");
      fprintf(parm," s1 symbol 0 \n");
      fprintf(parm," s1 line type 1 \n");
      fprintf(parm," s1 line linestyle 1 \n");
      fprintf(parm," s1 line linewidth 2.0 \n");
      fprintf(parm," s1 line color 2 \n");
      fprintf(parm," s1 legend \" partial core \" \n");

      /*      fprintf(parm," s2 hidden false \n");
	      fprintf(parm," s2 type xy \n");
	      fprintf(parm," s2 symbol 0 \n");
	      fprintf(parm," s2 line type 1 \n");
	      fprintf(parm," s2 line linestyle 1 \n");
	      fprintf(parm," s2 line linewidth 1.0 \n");
	      fprintf(parm," s2 line color 3 \n");
	      fprintf(parm," s2 legend \" partial core 1st deriv.\" \n");
	      
	      fprintf(parm," s3 hidden false \n");
	      fprintf(parm," s3 type xy \n");
	      fprintf(parm," s3 symbol 0 \n");
	      fprintf(parm," s3 line type 1 \n");
	      fprintf(parm," s3 line linestyle 1 \n");
	      fprintf(parm," s3 line linewidth 1.0 \n");
	      fprintf(parm," s3 line color 4 \n");
	      fprintf(parm," s3 legend \" partial core 2nd deriv.\" \n"); */

      /* valence density */
      
      fprintf(parm," s2 hidden false \n");
      fprintf(parm," s2 type xy \n");
      fprintf(parm," s2 symbol 0 \n");
      fprintf(parm," s2 line type 1 \n");
      fprintf(parm," s2 line linestyle 1 \n");
      fprintf(parm," s2 line linewidth 2.0 \n");
      fprintf(parm," s2 line color 3 \n");
      fprintf(parm," s2 legend \" valence density \" \n");
      
    } else {
      /* valence density */
      
      fprintf(parm," s1 hidden false \n");
      fprintf(parm," s1 type xy \n");
      fprintf(parm," s1 symbol 0 \n");
      fprintf(parm," s1 line type 1 \n");
      fprintf(parm," s1 line linestyle 1 \n");
      fprintf(parm," s1 line linewidth 2.0 \n");
      fprintf(parm," s1 line color 2 \n");
      fprintf(parm," s1 legend \" valence density \" \n");
    }
    fclose(parm);
  } 

  sprintf(comm, "xmgrace %s.plt_pcc -p pcc.par -saveall %s.pcc_agr & ", param->name,param->name);
  system(comm);
  free(comm);
  
  return 0;
}

