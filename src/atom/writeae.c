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
  
  /* First write out all AE things that we want to plot */

  /* first the partial core correction */  
  /* if there is a core correction, rscore=partial core , rscoretot=real core */
  /* if not, rscore=real core */

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


  /* NRL */
  if (!strcmp(param->reltype, "nrl")){

    for (i=0; i<aorb_.nval;i++) {

      if (aval_.ibd[i]==1) {
	for (j=0;j<param->ngrid;j++)
	   fprintf(fp,"%20.10lg %20.10lg \n",grid_.r[j],wfn_.rnl[i+ncore][j]);
      } else {
	j=0;
	while(grid_.r[j] < param->rc[i]*2) {
	  fprintf(fp,"%20.10lg %20.10lg \n",grid_.r[j],wfn_.rnl[i+ncore][j]);
	  j++;

	}
      }
      fprintf(fp,"@ \n");
    }
  } else {

 /* SRL or FRL */
    k=0;
    for (i=0; i<param->nvalrel; i++) {
      if (aval_.ibd[k]==1) {
	for (j=0;j<param->ngrid;j++)
	  fprintf(fp,"%20.10lg %20.10lg \n",grid_.r[j],wfnrel_.rnla[i+param->ncorerel][j]);
	fprintf(fp,"@ \n");
      } else {
	j=0;
	while(grid_.r[j] < param->rc[k]*2) {
	  fprintf(fp,"%20.10lg %20.10lg \n",grid_.r[j],wfnrel_.rnla[i+param->ncorerel][j]);
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
    
    /* If it is DFT and SRL we should have an avergaged wfn by now */
    if ((param->ixc >= 0)&&(!strcmp(param->reltype, "srl"))) {

      for (i=0; i<param->nval; i++) {
	if (aval_.ibd[i]==1) {
	  for (j=0;j<param->ngrid;j++)
	    if (param->rc[i] < grid_.r[j]) fprintf(fp,"%20.10lg %20.10lg \n",grid_.r[j],wfn_.rnl[i][j]);
	  fprintf(fp,"@ \n");
	} else {
	  j=0;
	  while(grid_.r[j] < param->rc[i]*2) {
	    if (param->rc[i] < grid_.r[j]) fprintf(fp,"%20.10lg %20.10lg \n",grid_.r[j],wfn_.rnl[i][j]);
	    j++;
	  }
	  fprintf(fp,"@ \n");
	}	  
      }
      
    }
  }
  
  fclose(fp);
  /* done with plots */


  /* now dump the all-electron wave function into a binary file */
    
  sprintf(filename, "%s.psi_ae", param->name);
  fp = fopen(filename, "wb");

  /* we are still carrying around the core unless we are dft+srl */
  if ((!strcmp(param->reltype, "nrl"))||(param->ixc < 0)||(!strcmp(param->reltype, "frl"))) {
    ncore = aorb_.norb - aorb_.nval;
  }else{
    ncore = 0;
  }

  for (i=0; i<aorb_.nval; i++) {
    fwrite(wfn_.rnl[ncore+i], sizeof(double), param->ngrid, fp);
  }
  fclose(fp);  

  sprintf(filename, "%s.psi_ae_ascii", param->name);
  fp = fopen(filename, "w");
  for (i=0; i<aorb_.nval; i++) {
    for (j=0;j<param->ngrid;j++)
      fprintf(fp,"%20.10lg, %20.10lg \n ",grid_.r[j],wfn_.rnl[ncore+i][j]);
    fprintf(fp,"\n");
  }
  fclose(fp);

  /* now dump the all-electron potential into a binary file */
  sprintf(filename, "%s.pot_ae", param->name);
  fp = fopen(filename, "wb");

  for (i=0; i<aorb_.nval; i++) {
    fwrite(totpot_.rvcore[i+ncore], sizeof(double), param->ngrid, fp);
    fwrite(totpot_.rvps[i+ncore], sizeof(double), param->ngrid, fp);
  }
  fclose(fp);

  /* dump energy for later */
  sprintf(filename, "%s.eng", param->name);
  fp = fopen(filename, "wb");
  fwrite(&aval_.etot, sizeof(double),1, fp);
  fclose(fp);

  /* this is sort of a generic catch-all structure */
  sprintf(filename, "%s.eig_ae", param->name);
  fp = fopen(filename, "wb");

  for (i=0; i<aorb_.nval; i++) {

  /* we are still carrying around the core unless we are dft+srl */
    ncore = aorb_.norb-aorb_.nval;
    if ((!strcmp(param->reltype, "nrl"))||(param->ixc < 0)||(!strcmp(param->reltype, "frl"))) {
      mcore=ncore;
    }else{ 
      mcore=0;
    }    
    /* ncore : these variables are NOT changed in 'average.f' */ 
    /* mcore : these variables are changed 'average.f' (mcore=0 since only the valence part remains) */
    /* ibd : is a avalence only variable */

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



  /*  Fixed core stuff
  sprintf(filename, "%s.psi_ae_core", param->name);
  fp = fopen(filename, "wb");

  if ((!strcmp(param->reltype, "nrl"))||(param->ixc < 0)||(!strcmp(param->reltype, "frl"))) {
    ncore = aorb_.norb - aorb_.nval;
    for (i=0; i<ncore; i++) 
      fwrite(wfn_.rnl[i], sizeof(double), param->ngrid, fp);
  }else{
    ncore = 0;
    for (i=0; i<ncore; i++)
      fwrite(wfn_.rnl[i], sizeof(double), param->ngrid, fp);
  }
  fclose(fp);
  
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
  */


}

