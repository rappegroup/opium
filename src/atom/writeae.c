/*
 * Copyright (c) 1998-2008 The OPIUM Group
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

/* standard libraries */
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "parameter.h"        /* defines structure: 'param_t' */
#include "nlm.h"              /* nlm_label call */
#include "cdim.h"             /* fortran code parameters */
#include "common_blocks.h"    /* fortran common blocks */
#include "energy.h"           /* this is for saving the energies */

void writeAE(param_t *param) {

  int i,j,k,ncore,mcore;
  FILE *fp;
  char filename[160];
  
  /* New plot section */
  
  sprintf(filename, "%s.pcc_plt", param->name);
  fp = fopen(filename, "w");

  
  if (param->rpcc>1e-12){
    for (j=0;j<param->ngrid;j++)
      fprintf(fp,"%lg %lg \n",grid_.r[j],rscore_.rscoretot[j]);
    fprintf(fp,"@ \n");
  }else{  
    for (j=0;j<param->ngrid;j++)
      fprintf(fp,"%lg %lg \n",grid_.r[j],rscore_.rscore[j]);
    fprintf(fp,"@ \n");
  }

  if (param->rpcc>1e-12){
    for (j=0;j<param->ngrid;j++)
      fprintf(fp,"%lg %lg \n",grid_.r[j],rscore_.rscore[j]);
    fprintf(fp,"@ \n");
  }
  fclose(fp);


  /* AE PLOT INFO */
  sprintf(filename, "%s.ae_plt", param->name);
  fp = fopen(filename, "w");

  ncore = param->norb-param->nval;  


  if (!strcmp(param->reltype, "nrl")){

    for (i=0; i<aorb_.nval;i++) {

      if (aval_.ibd[i]==1) {
	for (j=0;j<param->ngrid;j++)
	   fprintf(fp,"%lg %lg \n",grid_.r[j],wfn_.rnl[i+ncore][j]);
      } else {
	j=0;
	while(grid_.r[j] < param->rc[i]*2) {
	  fprintf(fp,"%lg %lg \n",grid_.r[j],wfn_.rnl[i+ncore][j]);
	  j++;

	}
      }
      fprintf(fp,"@ \n");
    }
  } else {
    k=0;
    for (i=0; i<param->nvalrel; i++) {
      if (aval_.ibd[k]==1) {
	for (j=0;j<param->ngrid;j++)
	  fprintf(fp,"%lg %lg \n",grid_.r[j],wfnrel_.rnla[i+param->ncorerel][j]);
	fprintf(fp,"@ \n");
      } else {
	j=0;
	while(grid_.r[j] < param->rc[k]*2) {
	  fprintf(fp,"%lg %lg \n",grid_.r[j],wfnrel_.rnla[i+param->ncorerel][j]);
	  j++;
	}
	fprintf(fp,"@ \n");
      }
      if (aorb_.lo[i+param->ncorerel]==0) {
	k++;
      } else {
	if ((aorb_.lo[i+param->ncorerel] == aorb_.lo[i-1+param->ncorerel])) {
	  k++;
	}
      }
    }
    
    if (param->ixc >= 0) {

      for (i=0; i<param->nval; i++) {
	if (aval_.ibd[i]==1) {
	  for (j=0;j<param->ngrid;j++)
	    if (param->rc[i] < grid_.r[j]) fprintf(fp,"%lg %lg \n",grid_.r[j],wfn_.rnl[i][j]);
	  fprintf(fp,"@ \n");
	} else {
	  j=0;
	  while(grid_.r[j] < param->rc[i]*2) {
	    if (param->rc[i] < grid_.r[j]) fprintf(fp,"%lg %lg \n",grid_.r[j],wfn_.rnl[i][j]);
	    j++;
	  }
	  fprintf(fp,"@ \n");
	}	  
      }
      
    }
  }
  
  fclose(fp);



  /* now dump the all-electron wave function into a binary file */
    
  sprintf(filename, "%s.psi_ae", param->name);
  fp = fopen(filename, "wb");
  ncore = aorb_.norb - aorb_.nval;
  if ((!strcmp(param->reltype, "nrl"))||(param->ixc < 0)){
    for (i=0; i<aorb_.nval; i++) 
      fwrite(wfn_.rnl[ncore+i], sizeof(double), param->ngrid, fp);
  }else{
    for (i=0; i<aorb_.nval; i++)
      fwrite(wfn_.rnl[i], sizeof(double), param->ngrid, fp);
  }
  fclose(fp);



  sprintf(filename, "%s.psi_ae_core", param->name);
  fp = fopen(filename, "wb");
  if ((!strcmp(param->reltype, "nrl")||(param->ixc < 0))){
    ncore = aorb_.norb - aorb_.nval;
    for (i=0; i<ncore; i++) 
      fwrite(wfn_.rnl[i], sizeof(double), param->ngrid, fp);
  }else{
    ncore = 0;
    for (i=0; i<ncore; i++)
      fwrite(wfn_.rnl[i], sizeof(double), param->ngrid, fp);
    /* fix  for srl */
  }
  fclose(fp);
  
  /* now dump the all-electron potential into a binary file */
  
  sprintf(filename, "%s.pot_ae", param->name);
  fp = fopen(filename, "wb");
  if ((!strcmp(param->reltype, "nrl"))||(param->ixc < 0)){
    ncore = aorb_.norb-aorb_.nval;
    mcore=ncore;
  }else{ 
    ncore = aorb_.norb-aorb_.nval;
    mcore=0;
  }    
  for (i=0; i<aorb_.nval; i++) {
    fwrite(totpot_.rvcore[i+mcore], sizeof(double), param->ngrid, fp);
    fwrite(totpot_.rvps[i+mcore], sizeof(double), param->ngrid, fp);
  }
  fclose(fp);

  sprintf(filename, "%s.eng", param->name);
  fp = fopen(filename, "wb");
  fwrite(&aval_.etot, sizeof(double),1, fp);
  fclose(fp);

  sprintf(filename, "%s.eig_ae", param->name);
  fp = fopen(filename, "wb");

  for (i=0; i<aorb_.nval; i++) {
    if ((!strcmp(param->reltype, "nrl"))||(param->ixc < 0)){
      ncore = aorb_.norb-aorb_.nval;
      mcore=ncore;
    }else{ 
      ncore = aorb_.norb-aorb_.nval;
      mcore=0;
    }    
    fwrite(&adat_.en[mcore+i], sizeof(double),1, fp);
    fwrite(&adat_.wnl[ncore+i], sizeof(double),1, fp);
    fwrite(&aorb_.nlm[ncore+i], sizeof(int),1, fp);
    fwrite(&aorb_.no[ncore+i], sizeof(int),1, fp);
    fwrite(&aorb_.lo[ncore+i], sizeof(int),1, fp);
    fwrite(&adat_.so[ncore+i], sizeof(double),1, fp);
    fwrite(&aorb_.nmax[mcore+i], sizeof(int),1, fp);
    fwrite(&aorb_.maxim, sizeof(int),1, fp);
    fwrite(&adat_.xion, sizeof(double),1, fp);
    fwrite(&aval_.ibd[i], sizeof(int),1, fp);
  }
  fclose(fp);
    
  ncore = param->norb-param->nval;

  sprintf(filename, "%s.eig_ae_core", param->name);
  fp = fopen(filename, "wb");

  ncore = aorb_.norb-aorb_.nval;
  for (i=0; i<ncore; i++) {
    fwrite(&adat_.en[i], sizeof(double),1, fp);
    fwrite(&adat_.wnl[i], sizeof(double),1, fp);
    fwrite(&aorb_.nlm[i], sizeof(int),1, fp);
    fwrite(&aorb_.nmax[i], sizeof(int),1, fp);
    fwrite(&aorb_.maxim, sizeof(int),1, fp);
    fwrite(&adat_.xion, sizeof(double),1, fp);
    fwrite(&aval_.ibd[i], sizeof(int),1, fp);
  }
  fclose(fp);
  
  /* now dump the PCC into a binary file if NLCC requested */

  if (param->rpcc>1e-12){
    sprintf(filename, "%s.rho_pcore", param->name);
    fp = fopen(filename, "wb");
    fwrite(rscore_.rscore, sizeof(double), param->ngrid, fp);
    fwrite(rscore_.rdd, sizeof(double), param->ngrid, fp);
    fwrite(rscore_.rddd, sizeof(double), param->ngrid, fp);
    fclose(fp);

    sprintf(filename, "%s.rho_fcore", param->name);
    fp = fopen(filename, "wb");
    fwrite(rscore_.rscoretot, sizeof(double), param->ngrid, fp);
    fclose(fp);
  }else{
    sprintf(filename, "%s.rho_fcore", param->name);
    fp = fopen(filename, "wb");
    fwrite(rscore_.rscore, sizeof(double), param->ngrid, fp);
    fclose(fp);
  }

}

