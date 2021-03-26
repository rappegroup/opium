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
#include <stdlib.h>
#include <string.h>

#include "parameter.h"
#include "cdim.h"        /* fortran code parameters */
#include "do_ke.h"
#include "common_blocks.h"
#include "nlm.h"

void kcomp_(int *, char *, double[10][N0] , int[10][N0]);
static char report[8000];
void nrelorbnl(param_t *param, int, char *logfile); 
void readAE(param_t *param);
void readPS(param_t *param);
char * write_reportke(param_t *param , char *rp, double[10][N0] , int[10][N0], int);
int do_ke(param_t *param, char *logfile){

  int i;
  int irel;
  int config=-1;
  int ikstor[10][N0];
  double rkstor[10][N0];
  FILE *fp_log,*fp;
  char *rp=report;
  char filename[80];

  fp_log = fopen(logfile, "a");

  fprintf(fp_log,"\n ======================================================================== \n");
  fprintf(fp_log," Begin KE convergence testing\n");
  fprintf(fp_log," ======================================================================== \n");
  fclose(fp_log);  

  /* set the log file */
  sprintf(filenames_.file_log, "%s", logfile);

  /*  aorb_.norb=param->nll;
      aorb_.nval=param->nll;*/

  if (!strcmp(param->reltype, "frl")) { 
    relsproj(param,logfile);
    relorbnl(param,config,logfile);
    irel=1;
  }else{
    nrelsproj(param,logfile);
    nrelorbnl(param,config,logfile);
    irel=0;
  }
  readPS(param);

  sprintf(filename, "%s.kedat", param->name);
  kcomp_(&irel,filename,rkstor,ikstor);

  rp=write_reportke(param,rp,rkstor,ikstor,aorb_.nval);
  fp_log = fopen(logfile, "a");

  fprintf(fp_log,"\n ======================================================================== \n");
  fprintf(fp_log," End KE convergence testing\n");
  fprintf(fp_log," ======================================================================== \n");

  fclose(fp_log);
  return 0;
}


/* report section */

void do_ke_report(FILE *fp){
  fprintf(fp, "%s", report);
}

char * write_reportke(param_t *param , char *rp, double rkstor[10][N0], int ikstor[10][N0], int nol) {
  
  int i,j;
  char sgn=' ';
  
  j=0;
  rp+=sprintf(rp, "\t ===   Ecut necessary for ~1    eV convergence error / electron === \n");
  rp+=sprintf(rp, "\t --------------------------------------------------------------- \n");
  rp+=sprintf(rp, "\t\t \t Ecut[Ry] \t error [meV/e] \n");
  for (i=0; i<nol; i++){
    if (!strcmp(param->reltype, "frl")) sgn=(adat_.so[i] > 0) ? '+' : '-' ; 
    if (aval_.ibd[i]==1) {
      if (ikstor[j][i] < 0) {
	rp+=sprintf(rp, "\t%3d%c \t\t   -----  >400 Ry   -----\n",aorb_.nlm[i],sgn);
      }else{
	rp+=sprintf(rp, "\t%3d%c \t\t   %d \t %16.3f \n",aorb_.nlm[i],sgn,ikstor[j][i],rkstor[j][i]);
      }
    }
  }     
  j=1;
  rp+=sprintf(rp, "\t ===   Ecut necessary for ~100 meV convergence error / electron === \n");
  rp+=sprintf(rp, "\t --------------------------------------------------------------- \n");
  rp+=sprintf(rp, "\t\t \t Ecut[Ry] \t error [meV/e] \n");
  for (i=0; i<nol; i++){
    if (!strcmp(param->reltype, "frl")) sgn=(adat_.so[i] > 0) ? '+' : '-' ; 
    if (aval_.ibd[i]==1) {
      if (ikstor[j][i] < 0) {
	rp+=sprintf(rp, "\t%3d%c \t\t   -----  >400 Ry   -----\n",aorb_.nlm[i],sgn);
      }else{
	rp+=sprintf(rp, "\t%3d%c \t\t   %d \t %16.3f \n",aorb_.nlm[i],sgn,ikstor[j][i],rkstor[j][i]);
      }
    }
  }     
  rp+=sprintf(rp, "\n");
  j=2;
  rp+=sprintf(rp, "\t ===   Ecut necessary for  ~10 meV convergence error / electron === \n");
  rp+=sprintf(rp, "\t --------------------------------------------------------------- \n");
  rp+=sprintf(rp, "\t\t \t Ecut[Ry] \t error [meV/e] \n");
  for (i=0; i<nol; i++){
    if (!strcmp(param->reltype, "frl")) sgn=(adat_.so[i] > 0) ? '+' : '-' ; 
    if (aval_.ibd[i]==1) {
      if (ikstor[j][i] < 0) {
	rp+=sprintf(rp, "\t%3d%c \t\t   -----  >400 Ry   -----\n",aorb_.nlm[i],sgn);
      }else{
	rp+=sprintf(rp, "\t%3d%c \t\t   %d \t %16.3f \n",aorb_.nlm[i],sgn,ikstor[j][i],rkstor[j][i]);
      }
    }
  }     
  rp+=sprintf(rp, "\n");
  j=3;
  rp+=sprintf(rp, "\t ===   Ecut necessary for   ~1 meV convergence error / electron === \n");
  rp+=sprintf(rp, "\t --------------------------------------------------------------- \n");
  rp+=sprintf(rp, "\t\t \t Ecut[Ry] \t error [meV/e] \n");
  for (i=0; i<nol; i++){
    if (!strcmp(param->reltype, "frl")) sgn=(adat_.so[i] > 0) ? '+' : '-' ; 
    if (aval_.ibd[i]==1) {
      if (ikstor[j][i] < 0) {
	rp+=sprintf(rp, "\t%3d%c \t\t   -----  >400 Ry   -----\n",aorb_.nlm[i],sgn);
      }else{
	rp+=sprintf(rp, "\t%3d%c \t\t   %d \t %16.3f \n",aorb_.nlm[i],sgn,ikstor[j][i],rkstor[j][i]);
      }
    }
  }     

  return rp;
}

