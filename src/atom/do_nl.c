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
#include <stdlib.h>

#include "parameter.h"        /* defines structure: 'param_t' */
#include "nlm.h"              /* nlm_label call */
#include "cdim.h"        /* fortran code parameters */
#include "do_nl.h"            /* the module's own header */
#include "common_blocks.h"    /* fortran common blocks */
#include "energy.h"           /* this is for saving the energies */

/* report feature */
//comment by JY hfsolve_ defined here
static char report[8000];
void nrelorbnl(param_t *param, int, char *); 
char * write_reportnl(param_t *param, char *rp,int,double temp_eigen[], double temp_norm[], int, int);
void dftsolve_(double  *, int * ,double *, int *, int *, int *, int *, int *, int *);
void readPS(param_t *param);
void writeNL(param_t *param);
void writeSL(param_t *param);
void hfsolve_(double  *, int * ,double *, int * , int *, int *, int *, int *);
void nrelsproj(param_t *param, char *);
void relsproj(param_t *param, char *);

int do_nl(param_t *param, char *logfile, int insl ){

  int i;  
  FILE *fp_log;
  FILE *fp;
  double zeff;
  int config=-1;
  int iprint=0;
  int ipsp=1;
  int ifc=0;
  int irel=0;
  int irelxc=0;
  int iexit=0;
  char *rp=report;
  char filename[160];
  double temp_eigen[10];
  double temp_norm[10];
  double exccut_temp;

  /* set the log file */
  sprintf(filenames_.file_log, "%s", logfile);

  fp_log = fopen(logfile, "a");

  fprintf(fp_log,"\n ======================================================================== \n");
  if (insl == 0) fprintf(fp_log," Begin SL calculation\n");
  if (insl == 1) fprintf(fp_log," Begin NL calculation\n");
  fprintf(fp_log," ======================================================================== \n");
  fclose(fp_log);

  if (param->exccut < 0) {
    exccut_temp=0.0;
  } else {
    exccut_temp=param->exccut;
  }
  
  aorb_.norb = param->nval;
  aorb_.nval = param->nval;
  aorb_.ncore = 0; 
  ipsp = 1;   
  nlpot2_.inl = insl;  

  if (!strcmp(param->reltype, "frl")) { 
    relsproj(param,logfile);
    relorbnl(param,config,logfile);
  }else{
    nrelsproj(param,logfile);
    nrelorbnl(param,config,logfile);
  }

  readPS(param);

  zeff=adat_.xion;
  for (i=0; i<param->nll; i++) {
    zeff +=adat_.wnl[i];


  }
  irel=0;
  irelxc=0;
  iprint=1;
  if (param->ixc < 0 || param->ixc == 7) {
    nlpot2_.inl = 0;  
    if (insl != 0 ){
      printf("!!NOTE!!: All Hartree-Fock and hybrid functional tests are done with the semi-local form of the potential \n");
      printf("!!NOTE!!: Use the 'sl' command instead of the 'nl' command to remove this warning...\n");

      fp_log = fopen(logfile, "a");

      fprintf(fp_log,"!!NOTE!!: All Hartree-Fock and hybrid functional tests are done with the semi-local form of the potential \n");
      fprintf(fp_log,"!!NOTE!!: Use the 'sl' command instead of the 'nl' command to remove this warning...\n");

      fclose(fp_log);

    }
    insl=0;
    //printf("zeff=%lf ixc=%i exccut_temp=%lf ipsp=%i ifc=%i iexit=%i irel=%i iprint=%i",zeff,param->ixc,exccut_temp,ipsp,ifc,iexit,irel,iprint);
    hfsolve_(&zeff,&param->ixc,&exccut_temp,&ipsp,&ifc,&iexit,&irel,&iprint);
  }else {
    //printf("do nonlocal \n");
    //printf("zeff=%lf ixc=%i exccut_temp=%lf ipsp=%i ifc=%i iexit=%i irel=%i irelxc=%i iprint=%i",zeff,param->ixc,exccut_temp,ipsp,ifc,iexit,irel,irelxc,iprint);
    dftsolve_(&zeff,&param->ixc,&exccut_temp,&ipsp,&ifc,&iexit,&irel,&irelxc,&iprint);
  }

  sprintf(filename, "%s.psi_last", param->name);
  fp = fopen(filename, "wb");
  for (i=0; i<param->nval; i++){
    fwrite(&aorb_.nlm[i], sizeof(int), 1, fp);
    fwrite(&adat_.wnl[i], sizeof(double), 1, fp);
    fwrite(&adat_.en[i], sizeof(double), 1, fp);
    fwrite(wfn_.rnl[i], sizeof(double), param->ngrid, fp);
  }
  fclose(fp);
  
  if (iexit) {
    if (insl == 0) printf("Terminal error in: scpot <-- do_sl\n EXITING OPIUM \n");
    if (insl == 1) printf("Terminal error in: scpot <-- do_nl\n EXITING OPIUM \n");
    exit(1);
  }

  if (insl != 0) {
    writeNL(param);
  } else {
    writeSL(param);
    if (param->nboxes != 0) {
      fp_log = fopen(logfile, "a");
      fprintf(fp_log,"You MUST run nl before writing out your psp file since you are \n");
      fprintf(fp_log,"using augmentation functions in your [KBdesign] keyblock.  \n\n");
      printf("You MUST run nl before writing out your psp file since you are \n");
      printf("using augmentation functions in your [KBdesign] keyblock.  \n\n");
    }
  }
  rp = write_reportnl(param,rp,config,temp_eigen,temp_norm,insl,aorb_.nval);

  ennl[0] = aval_.etot;

  fp_log = fopen(logfile, "a");

  fprintf(fp_log,"\n ======================================================================== \n");
  if (insl == 0) fprintf(fp_log," End SL calculation\n");
  if (insl == 1) fprintf(fp_log," End NL calculation\n");
  fprintf(fp_log," ======================================================================== \n");
  fclose(fp_log);
  return 0;
}

/* report section */

void do_nl_report(FILE *fp){
  fprintf(fp, "%s", report);
}

char * write_reportnl(param_t *param, char *rp, int config, double temp_eigen[], double temp_norm[], int insl, int nol) {
  
  int i;
  double temp_eigen_tot, temp_norm_tot;
  int igh=0;
  char sgn=' ';

  if (param->ixc<0) 
    rp+=sprintf(rp, "\n\n NOTICE!! :Ghost testing not done for HF psps yet, sorry :( \n");

  if (config<0) {

    if (insl==0) {
      rp+=sprintf(rp, 
		  "\n SL:Orbital    Filling       Eigenvalues[Ry]          Norm       Ghost\n"
		  "    ------------------------------------------------------------------\n");
    }else{
      rp+=sprintf(rp, 
		  "\n NL:Orbital    Filling       Eigenvalues[Ry]          Norm       Ghost\n"
		  "    ------------------------------------------------------------------\n");
    }      
    

    for (i=0; i<nol; i++){
      if (!strcmp(param->reltype, "frl")) sgn=(adat_.so[i] > 0) ? '+' : '-' ; 
      if (param->ixc<0) {
	rp+=sprintf(rp, "\t%3d%c\t%6.3f\t%19.10f\t%14.10f\t%6s\n",
		    aorb_.nlm[i],sgn,adat_.wnl[i], adat_.en[i],
		    aval_.rnorm[i],"???");

      }else if (aval_.ibd[i]==1) {
	rp+=sprintf(rp, "\t%3d%c\t%6.3f\t%19.10f\t%14.10f\t%6s\n",
		    aorb_.nlm[i],sgn,adat_.wnl[i], adat_.en[i],
		    aval_.rnorm[i],(local_.nlghost[i]>0)?"yes":((local_.nlghost[i]<0)?"?":"no"));
      }else{
	rp+=sprintf(rp, "\t%3d%c\t%6.3f\t%19.10f\t\t\t%6s\n",
		    aorb_.nlm[i],sgn,adat_.wnl[i], adat_.en[i],
		    (local_.nlghost[i]>0)?"yes":((local_.nlghost[i]<0)?"?":"no"));
	
      }
      igh+=local_.nlghost[i];
      
    }
    if (igh > 0) {
      rp+=sprintf(rp, "\n");
      rp+=sprintf(rp, "  !!ERROR!! Ghosts are present in pseudopotential \n");
      rp+=sprintf(rp, "  !!ERROR!! See log file for more information \n");
    }else{
      if (param->ixc >= 0) rp+=sprintf(rp, "\n\t ========== No ghosts in potential!!========== \n");
    }

  }else{
    if (insl==0) {
      rp+=sprintf(rp, 
		  "\n SL:Orbital    Filling       Eigenvalues[Ry]          Norm       Ghost\n"
		  "    ------------------------------------------------------------------\n");
    }else{
      rp+=sprintf(rp, 
		  "\n NL:Orbital    Filling       Eigenvalues[Ry]          Norm       Ghost\n"
		  "    ------------------------------------------------------------------\n");
    }
    sgn=' ';
    for (i=0; i<nol; i++){
      if (!strcmp(param->reltype, "frl")) sgn=(adat_.so[i] > 0) ? '+' : '-' ; 
      if (param->ixc<0) {
	rp+=sprintf(rp, "\t%3d%c\t%6.3f\t%19.10f\t%14.10f\t%6s\n",
		    aorb_.nlm[i], sgn,adat_.wnl[i], adat_.en[i],
		    aval_.rnorm[i],"???");
	
	
      }else if (aval_.ibd[i]==1) {
	rp+=sprintf(rp, "\t%3d%c\t%6.3f\t%19.10f\t%14.10f\t%6s\n",
		    aorb_.nlm[i], sgn,adat_.wnl[i], adat_.en[i],
		    aval_.rnorm[i],(local_.nlghost[i]>0)?"yes":((local_.nlghost[i]<0)?"?":"no"));
      }else{
	rp+=sprintf(rp, "\t%3d%c\t%6.3f\t%19.10f\t\t\t%6s\n",
		    aorb_.nlm[i],sgn, adat_.wnl[i], adat_.en[i],
		    (local_.nlghost[i]>0)?"yes":((local_.nlghost[i]<0)?"?":"no"));
	
      }
    }
  }
  rp+=sprintf(rp, "\n      E_tot = %19.10f Ry\n",
	      aval_.etot);
  
  if (config>=0) {
    if (insl==0) {
      rp+=sprintf(rp, 
		  "\n AE-SL:Orbital Filling       Eigenvalues[mRy]         Norm[1e-3] \n"
		  " AE-SL- --------------------------------------------------------------\n");
    }else{
      rp+=sprintf(rp, 
		  "\n AE-NL:Orbital Filling       Eigenvalues[mRy]         Norm[1e-3] \n"
		  " AE-NL- --------------------------------------------------------------\n");
    }      

    temp_eigen_tot = temp_norm_tot = 0.;
    
    sgn=' ';
    for (i=0; i<nol; i++) {
      if (!strcmp(param->reltype, "frl")) sgn=(adat_.so[i] > 0) ? '+' : '-' ; 
      if (aval_.ibd[i]==1) {
	if (insl==0) {
	  rp+=sprintf(rp, " AE-SL- %3d%c\t%6.3f\t%19.10f\t%14.10f\t\n",
		      aorb_.nlm[i],sgn, adat_.wnl[i], 1000.0*(temp_eigen[i]-adat_.en[i]),
		      1000.0*(temp_norm[i]-aval_.rnorm[i]));
	}else{
	  rp+=sprintf(rp, " AE-NL- %3d%c\t%6.3f\t%19.10f\t%14.10f\t\n",
		      aorb_.nlm[i], sgn,adat_.wnl[i], 1000.0*(temp_eigen[i]-adat_.en[i]),
		      1000.0*(temp_norm[i]-aval_.rnorm[i]));
	}
	temp_eigen_tot += fabs(temp_eigen[i]-adat_.en[i]);
	temp_norm_tot += fabs(temp_norm[i]-aval_.rnorm[i]);
      }
    }
    
    if (insl==0) {
      rp+=sprintf(rp, " AE-SL-  total error =\t%19.10f\t%14.10f\n\n",
		  1000*temp_eigen_tot, 1000*temp_norm_tot);
    }else{
      rp+=sprintf(rp, " AE-NL-  total error =\t%19.10f\t%14.10f\n\n",
		  1000*temp_eigen_tot, 1000*temp_norm_tot);
    }      
    rp+=sprintf(rp," =====================================================================\n");
  }    
  
  return rp;
}

void readPS(param_t *param) {
  int i,j,junk,ncore;
  FILE *fp;
  char filename[160];
  double junk2;
  char sgn='+';

  ncore=param->norb-param->nval;

  for (j=0; j<aorb_.nval; j++) {
    if (!strcmp(param->reltype, "frl")) {
      sgn=(adat_.so[j] > 0) ? '+' : '-' ; 
      sprintf(filename, "%s.psi.ps.l=%d%c", param->name,nlm_label(aorb_.nlm[j]).l,sgn);
    }else{
      sprintf(filename, "%s.psi.ps.l=%d", param->name,nlm_label(aorb_.nlm[j]).l);
    }
    fp = fopen(filename, "rb");
    fread(wfn_.rnl[j], sizeof(double), param->ngrid, fp);
    fclose(fp);
  }

  for (j=0; j<aorb_.nval; j++) {
    if (!strcmp(param->reltype, "frl")) {
      sgn=(adat_.so[j] > 0) ? '+' : '-' ; 
      sprintf(filename, "%s.pot.ps.l=%d%c", param->name,nlm_label(aorb_.nlm[j]).l,sgn);
    }else{
      sprintf(filename, "%s.pot.ps.l=%d", param->name,nlm_label(aorb_.nlm[j]).l);
    }
    
    fp = fopen(filename, "rb");
    fread(totpot_.rvcore[j], sizeof(double), param->ngrid, fp);
    fread(totpot_.rvps[j], sizeof(double), param->ngrid, fp);
    fclose(fp);
  }
  
  for (i=0; i<param->nval; i++) 
    for (j=0; j<param->ngrid; j++) 
      totpot_.rvcoul[j]=totpot_.rvps[i][j]-totpot_.rvcore[i][j] ;

  sprintf(filename, "%s.eig_ae", param->name);
  fp = fopen(filename, "rb");
  for (i=0; i<aorb_.nval; i++){
    fread(&adat_.en[i], sizeof(double),1, fp);
    fread(&junk2, sizeof(double),1, fp);
    fread(&junk, sizeof(int),1, fp);
    fread(&junk, sizeof(int),1, fp);
    fread(&junk, sizeof(int),1, fp);
    fread(&junk2, sizeof(double),1, fp);
    fread(&junk, sizeof(int),1, fp);
    fread(&junk, sizeof(int),1, fp);
    fread(&junk2, sizeof(double),1, fp);
    fread(&junk, sizeof(int),1, fp);
  }
  fclose(fp); 

  if (param->rpcc>1e-12){
    sprintf(filename, "%s.rho_pcore", param->name);
    fp = fopen(filename, "rb");
    fread(rscore_.rscore, sizeof(double), param->ngrid, fp);
    fread(rscore_.rdd, sizeof(double), param->ngrid, fp);
    fread(rscore_.rddd, sizeof(double), param->ngrid, fp);
    fclose(fp);

    sprintf(filename, "%s.rho_fcore", param->name);
    fp = fopen(filename, "wb");
    fread(rscore_.rscoretot, sizeof(double), param->ngrid, fp);
    fclose(fp);
  }else{
    sprintf(filename, "%s.rho_fcore", param->name);
    fp = fopen(filename, "wb");
    fread(rscore_.rscore, sizeof(double), param->ngrid, fp);
    fclose(fp);
  }
  
  /*  for (i=0; i<aorb_.nval; i++){
    printf("READPS  i,en,lo,so,nmax %d %f %d %f %d \n",i,adat_.en[i],aorb_.lo[i],adat_.so[i],aorb_.nmax[i]);
    }*/


    }

void writeNL(param_t *param) {

  int i,j;
  int iset=0;
  FILE *fp;
  char filename[160];

  /* now dump the local potential into a binary file */
  sprintf(filename, "%s.loc", param->name);
  fp = fopen(filename, "wb");
  fwrite(nlcore_.rvloc, sizeof(double), param->ngrid, fp);
  fclose(fp);
  
  sprintf(filename, "%s.psi_nl", param->name);
  fp = fopen(filename, "wb");
  for (i=0; i<param->nval; i++)
    fwrite(wfn_.rnl[i], sizeof(double), param->ngrid, fp);
  fclose(fp);

  sprintf(filename, "%s.rho_nl", param->name);
  fp = fopen(filename, "wb");
  fwrite(rscore_.rsval, sizeof(double), param->ngrid, fp);
  fclose(fp);

  sprintf(filename, "%s.pcc_plt", param->name);
  fp = fopen(filename, "a");
  for (j=0;j<param->ngrid;j++)
    fprintf(fp,"%20.10lg %20.10lg \n",grid_.r[j],rscore_.rsval[j]);
  fprintf(fp,"@ \n");
  fclose(fp);

  sprintf(filename, "%s.nl_plt", param->name);
  fp = fopen(filename, "w");
  for (i=0; i<aorb_.nval;i++) {
    if (aval_.ibd[i]==1) {
      for (j=0;j<param->ngrid;j++)
	if ((grid_.r[j]<param->rc[i])||(iset)) {
	  fprintf(fp,"%20.10lg %20.10lg 0.0\n",grid_.r[j],wfn_.rnl[i][j]);
	}else{
	  fprintf(fp,"%20.10lg %20.10lg 1e-8\n",grid_.r[j],wfn_.rnl[i][j]);
	  iset=1;
	}
    } else {
      j=0;
      while(grid_.r[j] < param->rc[i]) {
	fprintf(fp,"%20.10lg %20.10lg 0.0\n",grid_.r[j],wfn_.rnl[i][j]);
	j++;
      }
    }
    fprintf(fp,"@ \n");
    iset=0;
  }
  fclose(fp);

  if (param->ixc >= 0 && param->ixc != 7) {
    sprintf(filename, "%s.vi_plt", param->name);
    fp = fopen(filename, "a");
    for (j=0;j<param->ngrid;j++){
      fprintf(fp,"%20.10lg %20.10lg 0.0\n",grid_.r[j],nlcore_.rvloc[j]/grid_.r[j]);
    }
    fprintf(fp,"@ \n");
    fclose(fp);
    
    sprintf(filename, "%s.vs_plt", param->name);
    fp = fopen(filename, "a");
    for (j=0;j<param->ngrid;j++){
      fprintf(fp,"%20.10lg %20.10lg 0.0\n",grid_.r[j],(nlcore_.rvloc[j]+totpot_.rvcoul[j])/grid_.r[j]);
    }
    fprintf(fp,"@ \n");
    fclose(fp);
  }
}
void writeSL(param_t *param) {

  int i,j;
  int iset=0;
  FILE *fp;
  char filename[160];

  if (param->nboxes == 0) {
    sprintf(filename, "%s.loc", param->name);
    fp = fopen(filename, "wb");
    fwrite(nlcore_.rvloc, sizeof(double), param->ngrid, fp);
    fclose(fp);
  } 

  sprintf(filename, "%s.psi_sl", param->name);
  fp = fopen(filename, "wb");
  for (i=0; i<param->nval; i++)
    fwrite(wfn_.rnl[i], sizeof(double), param->ngrid, fp);
  fclose(fp);

  sprintf(filename, "%s.rho_sl", param->name);
  fp = fopen(filename, "wb");
  fwrite(rscore_.rsval, sizeof(double), param->ngrid, fp);
  fclose(fp);

  sprintf(filename, "%s.pcc_plt", param->name);
  fp = fopen(filename, "a");
  for (j=0;j<param->ngrid;j++)
    fprintf(fp,"%20.10lg %20.10lg \n",grid_.r[j],rscore_.rsval[j]);
  fprintf(fp,"@ \n");
  fclose(fp);

  sprintf(filename, "%s.nl_plt", param->name);
  fp = fopen(filename, "w");
  for (i=0; i<param->nval;i++) {
    if (aval_.ibd[i]==1) {
      for (j=0;j<param->ngrid;j++)
	if ((grid_.r[j]<param->rc[i])||(iset)) {
	  fprintf(fp,"%20.10lg %20.10lg 0.0\n",grid_.r[j],wfn_.rnl[i][j]);
	}else{
	  fprintf(fp,"%20.10lg %20.10lg 1e-8\n",grid_.r[j],wfn_.rnl[i][j]);
	  iset=1;
	}
    } else {
      j=0;
      while(grid_.r[j] < param->rc[i]) {
	fprintf(fp,"%20.10lg %20.10lg 0.0\n",grid_.r[j],wfn_.rnl[i][j]);
	j++;
      }
    }
    fprintf(fp,"@ \n");
    iset=0;
  }
  fclose(fp);

  if (param->ixc >= 0 && param->ixc != 7) {
    sprintf(filename, "%s.vi_plt", param->name);
    fp = fopen(filename, "a");
    for (j=0;j<param->ngrid;j++){
      fprintf(fp,"%20.10lg %20.10lg 0.0\n",grid_.r[j],nlcore_.rvloc[j]/grid_.r[j]);
    }
    fprintf(fp,"@ \n");
    fclose(fp);
    
    sprintf(filename, "%s.vs_plt", param->name);
    fp = fopen(filename, "a");
    for (j=0;j<param->ngrid;j++){
      fprintf(fp,"%20.10lg %20.10lg 0.0\n",grid_.r[j],(nlcore_.rvloc[j]+totpot_.rvcoul[j])/grid_.r[j]);
    }
    fprintf(fp,"@ \n");
    fclose(fp);
  }
}
