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
#include <string.h>
#include <math.h>
#include <stdlib.h>

#include "parameter.h"
#include "cdim.h"
#include "nlm.h"
#include "common_blocks.h"   

void nrelorbae(param_t *param, int config, char *logfile) {
  
  int i,j,ic,jc,ii,jj;
  FILE *fp_log;
  int npot[10];
  
  for (i=0; i<10; i++){
    npot[i]=0;
  }

  aorb_.norb = param->norb;
  aorb_.nval = param->nval;
  aorb_.ncore = param->norb-param->nval;

  adat_.xion = param->z;

  for (i=0; i<aorb_.ncore; i++){
    aorb_.nlm[i] = param->nlm[i];
    adat_.wnl[i] = param->wnl[i];
    adat_.xion -= param->wnl[i];
    adat_.en[i] = param->en[i];
    aorb_.no[i] = nlm_label(param->nlm[i]).n;
    aorb_.lo[i] = nlm_label(param->nlm[i]).l;
  }

  for (i=aorb_.ncore; i<aorb_.norb; i++){

    if (config<0) {  
      aorb_.nlm[i] = param->nlm[i];
      adat_.wnl[i] = param->wnl[i];
      aorb_.no[i] = nlm_label(param->nlm[i]).n;
      aorb_.lo[i] = nlm_label(param->nlm[i]).l;
      if (config==-1) adat_.en[i] = param->en[i];
    }else{
      aorb_.nlm[i] = param->nlm_conf[config][i-aorb_.ncore];
      adat_.wnl[i] = param->wnl_conf[config][i-aorb_.ncore];
      aorb_.no[i] = nlm_label(param->nlm_conf[config][i-aorb_.ncore]).n;
      aorb_.lo[i] = nlm_label(param->nlm_conf[config][i-aorb_.ncore]).l;
      adat_.en[i] = param->en_conf[config][i-aorb_.ncore];
    }
  }

  for (i=aorb_.ncore; i<aorb_.norb; i++){

    ic=i-aorb_.ncore;
    aval_.rcall[ic] = param->rc[ic];   
    aval_.ibd[ic] = 1;
    if (adat_.wnl[i]<0) aval_.ibd[ic] = 0;
    if (aval_.ibd[ic]==0)  adat_.wnl[i]=0.0;

    adat_.xion -= adat_.wnl[i];
  }

    
  /*     fp_log = fopen(logfile, "a");
   fprintf(fp_log,"  Core orbitals - \n");
   for (i=0; i<aorb_.ncore; i++){
    fprintf(fp_log," |%d%d0>  %6.3f  %6.3f \n",aorb_.no[i],aorb_.lo[i],adat_.wnl[i],adat_.en[i]);
  }
  fprintf(fp_log,"  Valence orbitals - \n");
  for (i=aorb_.ncore; i<aorb_.norb; i++){
    ic=i-aorb_.ncore;
    fprintf(fp_log," |%d%d0>  %6.3f %6.3f %d  \n",aorb_.no[i],aorb_.lo[i],adat_.wnl[i],adat_.en[i],aval_.ibd[ic]);
    }
    fclose(fp_log);*/
}


