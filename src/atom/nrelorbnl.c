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

void nrelorbnl(param_t *param, int config, char *logfile) {

  int i,ii,j,k,l;
  int ncore;
  FILE *fp_log;

  ncore=param->norb - param->nval;
  aorb_.norb = param->nval;
  aorb_.nval = param->nval;
  aorb_.ncore = 0;
  adat_.xion = param->z;

  for (i=0; i<ncore; i++)
    adat_.xion -= param->wnl[i];

  /* map the AE |nlm> labels onto the NL |nlm> labels */
  for (i=0; i<param->nval; i++){
    if (config<0) {
      aorb_.nlm[i] = param->nlm[ncore + i];
      adat_.wnl[i] = param->wnl[ncore + i];
      /*      adat_.en[i] = param->en[ncore + i];*/
    }else{
      aorb_.nlm[i] = param->nlm_conf[config][i];
      adat_.wnl[i] = param->wnl_conf[config][i];
      /*      adat_.en[i] = param->en_conf[config][i];*/
    } 
    aval_.ibd[i] = 1;
    if (adat_.wnl[i]<0) aval_.ibd[i] = 0;
    if (aval_.ibd[i]==0)  adat_.wnl[i]=0.0;

    aval_.rcall[i]=param->rc[i];

    adat_.xion -= adat_.wnl[i];
    k = 0;
    for (j=0; j<i; j++)
	if (nlm_label(aorb_.nlm[j]).l == nlm_label(aorb_.nlm[i]).l) k++;
    aorb_.lo[i] = nlm_label(aorb_.nlm[i]).l;    
    aorb_.no[i] = k+1+aorb_.lo[i];
    aorb_.nlm[i] = (aorb_.no[i])*100 + aorb_.lo[i]*10
	+ nlm_label(aorb_.nlm[i]).m;
  }
  /*  for (i=0; i<param->nval; i++){
    printf(" %d %d %d %f %d\n",aorb_.no[i],aorb_.lo[i],aorb_.nlm[i],adat_.wnl[i],aval_.ibd[i]);
  }
  
  for (i=0; i<aorb_.nval; i++){
    printf("nrelorbnl  i,en,lo,so %d %f %d %f\n",i,adat_.en[i],aorb_.lo[i],adat_.so[i]);
    }*/


}
