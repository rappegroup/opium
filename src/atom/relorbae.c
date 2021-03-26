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

void relorbae(param_t *param, int config, char *logfile) {
  
  int i,j,ic,jc,ncore,ii,jj;
  double sc;
  FILE *fp_log;
  int npot[10];
  int ni,nj,nto,mto;

  /* important:  all orbital counting should use the non-rel norb,ncore,nval since the 
     param file uses these */

  /* add etrial struct for rel sols */

  for (i=0; i<10; i++){
    npot[i]=0;
  }
  param->isemi=0;

  ncore = param->norb-param->nval;
  adat_.xion = param->z;
  
  aorb_.norb=0;
  aorb_.ncore=0;
  aorb_.nval=0;

  /* populate core orbitals */  
  sc=-0.5;
  for (i=0; i<ncore; i++){
    for (j=0; j<2; j++){
      aorb_.no[aorb_.norb] = nlm_label(param->nlm[i]).n;
      aorb_.lo[aorb_.norb] = nlm_label(param->nlm[i]).l;
      adat_.so[aorb_.norb] = sc;
      aorb_.nlm[aorb_.norb]=param->nlm[i];
      adat_.en[aorb_.norb]=param->en[i];
      adat_.wnl[aorb_.norb] = 2.0*(nlm_label(param->nlm[i]).l+sc)+1.0;
      adat_.xion -= adat_.wnl[aorb_.norb];
      aval_.ibd[aorb_.norb]=param->ibound[i];

      if (fabs(adat_.wnl[aorb_.norb])>0.1) aorb_.norb++;
      sc=-1.0*sc;
    }
    
  }
  sc=-0.5;
  ic=0;
  aorb_.ncore = aorb_.norb;

  for (i=ncore; i<param->norb; i++){

    for (j=0; j<2; j++){
      if (config <0) {
	aorb_.no[aorb_.norb] = nlm_label(param->nlm[i]).n;
	aorb_.lo[aorb_.norb] = nlm_label(param->nlm[i]).l;
	adat_.so[aorb_.norb] = sc;
	adat_.en[aorb_.norb] = param->en[i];
	  
	if (param->wnl[i] < 0) {
	  adat_.wnl[aorb_.norb]=0.0;
	  aval_.ibd[ic]=0;
	}else{
	  adat_.wnl[aorb_.norb] = param->wnl[i]*(2*(aorb_.lo[aorb_.norb]+adat_.so[aorb_.norb])+1)/(4*aorb_.lo[aorb_.norb]+2);
	  aval_.ibd[ic]=1;
	}

      }else{
	aorb_.no[aorb_.norb] = nlm_label(param->nlm_conf[config][i-ncore]).n;
	aorb_.lo[aorb_.norb] = nlm_label(param->nlm_conf[config][i-ncore]).l;
	adat_.so[aorb_.norb] = sc;
	adat_.en[aorb_.norb] = param->en_conf[config][i-ncore];

	if (param->wnl_conf[config][i-ncore] < 0) {
	  adat_.wnl[aorb_.norb]=0.0;
	  aval_.ibd[ic]=0;
	}else{
	  adat_.wnl[aorb_.norb] = param->wnl_conf[config][i-ncore]*(2*(aorb_.lo[aorb_.norb]
				      +adat_.so[aorb_.norb])+1)/(4*aorb_.lo[aorb_.norb]+2);
	  aval_.ibd[ic]=1;
	}

      }
      /*      optparam_.qcl[aorb_.norb-aorb_.ncore]=param->qc[i-ncore];
	      optparam_.nbl[aorb_.norb-aorb_.ncore]=param->nb[i-ncore];*/
      aval_.rcall[aorb_.norb-aorb_.ncore]=param->rc[i-ncore];
      /*      printf(" aorb_.norb-aorb_.ncore,i-ncore,rc %d %d %lg \n",aorb_.norb-aorb_.ncore,i-ncore,aval_.rcall[aorb_.norb-aorb_.ncore]);*/
      aorb_.nlm[aorb_.norb]=param->nlm[i];
      adat_.xion -= adat_.wnl[aorb_.norb];
      if (aorb_.lo[aorb_.norb]+adat_.so[aorb_.norb]>0.0) {
	aorb_.norb++;
	ic++;
      }
      sc=-1.0*sc;


    }
  }

  aorb_.nval =  aorb_.norb - aorb_.ncore;
  
  for (i=aorb_.ncore; i<aorb_.norb; i++)
    for (j=i+1; j<aorb_.norb; j++) 
      if ((aorb_.lo[j]==aorb_.lo[i])&&(aorb_.no[j]!=aorb_.no[i])) param->isemi=1;
  

  /*    fp_log = fopen(logfile, "a");
  fprintf(fp_log,"  Core orbitals - nlm s j k occ etrial\n");
  for (i=0; i<aorb_.ncore; i++){
    fprintf(fp_log," |%d%d0> %+1.0f/2 %+1.0f/2 %+1.0f %6.3f %6.3f  \n",aorb_.no[i],aorb_.lo[i],2*adat_.so[i],2*(aorb_.lo[i]+adat_.so[i]),-2*adat_.so[i]*(aorb_.lo[i]+adat_.so[i]+0.5),adat_.wnl[i],adat_.en[i]);
  }
  
  fprintf(fp_log,"  Valence orbitals - nlm s j k occ etrial\n");
  for (i=aorb_.ncore; i<aorb_.norb; i++){
    fprintf(fp_log," |%d%d0> %+1.0f/2 %+1.0f/2 %+1.0f %6.3f %6.3f \n",aorb_.no[i],aorb_.lo[i],2*adat_.so[i],2*(aorb_.lo[i]+adat_.so[i]),-2*adat_.so[i]*(aorb_.lo[i]+adat_.so[i]+0.5),adat_.wnl[i],adat_.en[i]);
  }
  fclose(fp_log);
  */
  
  for (i=0; i< aorb_.norb; i++) {
    for (j=0; j< aorb_.norb; j++) {
      if (i==j) continue;
      if (aorb_.no[i] != aorb_.no[j]) continue;
      if (aorb_.lo[i] != aorb_.lo[j]) continue;
      if (fabs(adat_.so[i] - adat_.so[j])>0.001) continue;
      printf("%d %d two orbitals the same! -- abort! \n",i,j);
      exit(1);
    }
  }
}


