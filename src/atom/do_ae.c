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
#include "nlm.h"
#include "cdim.h"
#include "do_ae.h"
#include "common_blocks.h"   
#include "energy.h"

#define streq(a,b) (!strcasecmp(a,b))

int do_ae(param_t *param, char *logfile){
  
  int i,j,k;
  int irel=0;
  int irelxc=0;
  int ncore;
  int mcore;
  FILE *fp_log;
  char *rp = report;
  int ifrl;
  int config=-1;
  int ipsp=0;
  int ifc=0;
  int iexit=0;
  int iprint=1;
  double exccut_temp;
  double temp_eigen[10];
  double temp_norm[10];
  int ikstor[3];
  double rkstor[3];
  char filename[80];

  /* clear report buffer */
  *rp = 0;

  /* set log file name */
  sprintf(filenames_.file_log, "%s", logfile);

  fp_log=fopen(logfile,"a");

  fprintf(fp_log,"\n ======================================================================== \n");
  fprintf(fp_log," Begin AE calculation\n");
  fprintf(fp_log,  " ======================================================================== \n");

  if (param->exccut < 0) {
    exccut_temp=0.0;
  } else {
    exccut_temp=param->exccut;
  }

  if (!strcmp(param->reltype, "nrl")){
    
    /* do the NRL type AE solve */

    fprintf(fp_log," Performing non-relativistic AE calculation...  \n" );
    fclose(fp_log);

    /* set up the nlm's ; config=-1 is the reference state */
    nrelorbae(param,config,logfile);

    /* starting guess for the AE potential */
    startae(param, aorb_.norb);

    /* find the SCF solution */
    irel=0;
    irelxc=0;    
    if (param->ixc < 0) {
      hfsolve_(&param->z,&param->ixc,&exccut_temp,&ipsp,&ifc,&iexit,&irel,&iprint);
    }else {
      dftsolve_(&param->z,&param->ixc,&exccut_temp,&ipsp,&ifc,&iexit,&irel,&irelxc,&iprint);
    }
    if (iexit) {
      printf("Terminal error in: scpot <-- do_ae\n EXITING OPIUM \n");
      exit(1);
    }

  } else {

    fprintf(fp_log," Performing relativistic AE calculation...  \n" );
    fclose(fp_log);

    /* turn the l-orbitals into j-orbitals */
    relorbae(param,config,logfile);

    /* starting guess for the AE potential */
    startae(param, aorb_.norb);

    irel=1;
    if (streq(param->relxc,"rxc")){ 
      irelxc=1;
    }else{
      irelxc=0;
    }

    if (param->ixc < 0) {
      dfsolve_(&param->z,&param->ixc,&exccut_temp,&ipsp,&ifc,&iexit,&irel,&iprint);
    }else {
      dftsolve_(&param->z,&param->ixc,&exccut_temp,&ipsp,&ifc,&iexit,&irel,&irelxc,&iprint);
    }

    if (iexit) {
      printf("Terminal error in: atm <-- do_ae\n EXITING OPIUM \n");
      exit(1);
    }
    
    config=-2;

    param->nvalrel=aorb_.nval;
    param->ncorerel=aorb_.ncore;

    if (!strcmp(param->reltype, "srl")) {
      if (param->ixc < 0) {
	/*printf(" nval= %d \n ", aorb_.nval);*/
      }else{
	nrelorbae(param,config,logfile);
	average_(&param->ixc,&param->exccut, &iexit); 
      }
    }
    if (iexit) {
      printf("Terminal error in: average <-- do_ae\n EXITING OPIUM \n");
      exit(1);
    }
    
  }

  enae[0] = aval_.etot;

  if (param->rpcc>1e-12) {
    if (param->ixc < 0) {
      fp_log = fopen(logfile, "a");
      fprintf(fp_log, "  Hartree-Fock exchange is incompatible with a partial core correction  -- STOP\n");
      printf("  Hartree-Fock exchange is incompatible with a partial core correction -- STOP\n");
      fclose(fp_log);
      exit(1);
    }

    fp_log = fopen(logfile, "a");
    fprintf(fp_log, "   ================================================\n");
    fprintf(fp_log,"<<<pseudizing core>>>\n");
    fclose(fp_log);
    getpcc_(&iexit);

    if (iexit) {
      printf("Terminal error in: getpcc <-- do_ae\n EXITING OPIUM \n");
      exit(1);
    }

    fp_log = fopen(logfile, "a");
    fprintf(fp_log, "   ================================================\n");
    fclose(fp_log);
  }
  /* record data in report */
  rp = write_reportae(param,rp,config,temp_eigen,temp_norm,aorb_.norb);
  
  writeAE(param);
  
  fp_log = fopen(logfile, "a");
  fprintf(fp_log,"\n ======================================================================== \n");
  fprintf(fp_log," End AE calculation\n");
  fprintf(fp_log," ======================================================================== \n");

  fclose(fp_log);

  return 0;
}

/* report section */

void do_ae_report(FILE *fp){
  fprintf(fp, "%s", report);
}

char * write_reportae(param_t *param, char *rp,int config, double temp_eigen[], double temp_norm[], int nol) {
  
  int i,j,l,ncore;
  double e_rel,n_rel;
  char sgn='+';

  ncore=aorb_.norb - aorb_.nval;

  rp += sprintf(rp, 
		" AE:Orbital    Filling       Eigenvalues[Ry]          Norm\n"
		"    ----------------------------------------------------------\n");

  if (!strcmp(param->reltype, "nrl")||!strcmp(param->reltype, "frl")){
    for (i=0; i<nol; i++)
      if (i>ncore-1) {
	if (!strcmp(param->reltype, "frl")) {
	  sgn=(adat_.so[i] > 0) ? '+' : '-' ; 
	  rp += sprintf(rp, "\t%3d%c\t%6.3f\t%19.10f\t%14.10f\n",
			aorb_.nlm[i],sgn, adat_.wnl[i], adat_.en[i],aval_.rnorm[i]);
	}else{
	  rp += sprintf(rp, "\t%3d\t%6.3f\t%19.10f\t%14.10f\n",
			aorb_.nlm[i], adat_.wnl[i], adat_.en[i],aval_.rnorm[i]);
	}	  
	temp_eigen[i-ncore]=adat_.en[i];
	temp_norm[i-ncore]=aval_.rnorm[i];
      }else{
	if (!strcmp(param->reltype, "frl")) {
	  sgn=(adat_.so[i] > 0) ? '+' : '-' ; 
	  rp += sprintf(rp, "\t%3d%c\t%6.3f\t%19.10f\n",
			aorb_.nlm[i],sgn,adat_.wnl[i], adat_.en[i]);
	}else{
	  rp += sprintf(rp, "\t%3d\t%6.3f\t%19.10f\n",
			aorb_.nlm[i], adat_.wnl[i], adat_.en[i]);
	}

      }

  } else if (!strcmp(param->reltype, "srl")){
    j=aorb_.ncore;

    if ((config<0)&&(param->ixc>=0)) {

      for (i=0; i<param->nval; i++)
	rp += sprintf(rp, "\t%3d\t%6.3f\t%19.10f\t%14.10f\n",
		      aorb_.nlm[i+ncore], adat_.wnl[i+ncore], adat_.en[i],
		      aval_.rnorm[i]);
    }else{
      for (i=0; i<param->nval; i++){

	l=aorb_.lo[j];

	if (l != 0){
	  e_rel = (l*adat_.en[j] + (l+1)*adat_.en[j+1])/ (double)(2.*l+1.);
	  n_rel = (l*aval_.rnorm[j] + (l+1)*aval_.rnorm[j+1]) / (double)(2.*l+1.);
	  ++j;
	}else{
	  e_rel = adat_.en[j];
	  n_rel = aval_.rnorm[j];
	}
	temp_eigen[i]=e_rel;
	temp_norm[i]=n_rel;
	
	ncore=param->norb-param->nval;
	if (config<0) {
	  rp+=sprintf(rp, "\t%3d\t%6.3f\t%19.10f\t%14.10f\n",
		      param->nlm[i+ncore], param->wnl[i+ncore],e_rel, n_rel);
	}else {
	  rp+=sprintf(rp, "\t%3d\t%6.3f\t%19.10f\t%14.10f\n",
		      param->nlm_conf[config][i], param->wnl_conf[config][i],e_rel, n_rel);
	}
	j++;
      }
    }
  }
  
  rp += sprintf(rp, "\n      E_tot = %19.10f Ry\n", aval_.etot);

  return rp;
  
}

