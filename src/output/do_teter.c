/*
 * Copyright (c) 1998-2010 The OPIUM Group
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
#include <math.h>
#include <time.h>

#include "parameter.h"
#include "do_teter.h"
#include "uniPPlib.h"
#include "cdim.h"        
#include "nlm.h"

void writeparam(param_t *param, FILE *fp, FILE *fp_param);
void nrelsproj(param_t *param, char *);

int do_teter(param_t *param, FILE *fp_param, char *logfile){

  int i, l, k;
  uniPP unipp;
  char filename[180];
  FILE *fp;
  FILE *fp_log;
  time_t t;
  int icount,lixc;
  char datestring[7];
  double rpccmax;   
  int ill[N0];
  int ncore,c;
  int kk=0;
  double z_ion,m_mesh,a_mesh;
  int l_max,l_loc,nlcc;
    
  static double rscore[NPDM],rdd[NPDM],rddd[NPDM];
  static double rvcore[N0][NPDM],rvloc[NPDM];
  static double rnl[N0][NPDM];

  fp_log = fopen(logfile, "a");
  fprintf(fp_log,"<<<do_teter>>>\n");
  fclose(fp_log);
}
  /*  ncore=param->norb - param->nval;
  

  if ((!strcmp(param->reltype, "nrl")) || (!strcmp(param->reltype, "srl")) && (param->ixc >= 0)) {
    nrelorbnl(param,config,logfile);
    nrelsproj(param,logfile);

  } else {
    relorbnl(param,config,logfile);
    relsproj(param,logfile);

  }

  for (i=0; i<param->nll; i++) {
    if (!strcmp(param->reltype, "frl")) {
      sgn=(adat_.so[i] > 0) ? '+' : '-' ; 
      printf(" i lo %d %d \n",i,aorb_.lo[i]);
      sprintf(filename, "%s.pot.ps.l=%d%c", param->name,aorb_.lo[i],sgn);
    }else{
      sprintf(filename, "%s.pot.ps.l=%d", param->name,aorb_.lo[i]);
    }
    fp = fopen(filename, "rb");
    fread(rvcore[i], sizeof(double), param->ngrid, fp);
    fseek(fp,sizeof(double) ,param->ngrid);
    fclose(fp);
  }

  for (i=0; i<param->nll; i++) {
    if (!strcmp(param->reltype, "frl")) {
      sgn=(adat_.so[i] > 0) ? '+' : '-' ; 
      sprintf(filename, "%s.psi.ps.l=%d%c", param->name,aorb_.lo[i],sgn);
      printf(" i lo %d %d \n",i,aorb_.lo[i]);

    }else{
      sprintf(filename, "%s.psi.ps.l=%d", param->name,aorb_.lo[i]);
    }      
    fp = fopen(filename, "rb");
    fread(rnl[i], sizeof(double), param->ngrid, fp);
    fclose(fp);
  }


  ncore=param->norb - param->nval;

  z_ion=param->z;
  for (i=0; i<param->norb - param->nval; i++)
    z_ion -= param->wnl[i];
  l_max=param->nll;
  if (param->nboxes > 0)
    l_loc=param->nll;
  else
    l_loc = nlm_label(param->nlm[param->localind+ncore]).l;
  if (param->rpcc > 0.)
    nlcc=1;
  else
    nlcc=0;
  m_mesh=param->ngrid;
  a_mesh=exp(param->b);



  if (param->nboxes > 0) {
    fprintf(fp_log," teter format does not support the use of augmentation operators for now\n");
    printf(" teter format does not support the use of augmentation operators for now\n");
    exit(1);
  }
*/  /*    sprintf(filename, "%s.loc", param->name);
    if (fp = fopen(filename, "rb")) {
      fread(rvloc, sizeof(double), param->ngrid, fp);
      fclose(fp);
    } else {
      fp_log = fopen(logfile, "a");
      fprintf(fp_log,"Looks like you never ran nl yet you have augmentation functions :( --EXIT!\n");
      printf("Looks like you never ran nl yet you have augmentation functions :( --EXIT!\n");
      fclose(fp_log);
      exit(1);
      }

    fp_log = fopen(logfile, "a");
    fprintf(fp_log," Making l+1 the local potential %d\n",kk);
    fclose(fp_log);

    }  */

/*if (param->ixc == 0) lixc=2;
  if (param->ixc < 0) lixc=99;
  if (param->ixc == 1) lixc=7;
  if (param->ixc == 2) lixc=11;
  if (param->ixc == 3) lixc=99;
  if (param->ixc == 4) lixc=23;
  if (param->ixc == 5) lixc=24;

  
  sprintf(filename, "%s.teter", param->name);
  fp = fopen(filename, "w");
  fprintf(fp, "OPIUM generated %s potential\n", param->symbol);
  time(&t);
  strftime(datestring, 7, "%y%m%d", localtime(&t));
  fprintf(fp, "%10.5f%10.5f  %s\t\t\tzatom,zion,pspdat\n",
    param->z, z_ion, datestring);
  fprintf(fp, "%5d%5d%5d%5d%10d%10.5f\tpspcod,pspxc,lmax,lloc,mmax,r2well\n",
    5, lixc, l_max-1, l_loc,m_mesh, 0.);

  for (j=0;j<param->nll;j++) {
    nn=aorb_.no[j];
    ll=aorb_.lo[j];
    rc=aval_.rcall[j];
    


  fprintf(fp, "%10.5f%10.5f%10.5f\t\t\trchrg,fchrg,qchrg\n",
	  nlcc?rpccmax:0., nlcc?1.:0., 0.);

  writeparam(param, fp, fp_param);

  fp_log = fopen(logfile, "a");  
  fprintf(fp_log,"   ================================================\n");
  fclose(fp_log);
  
  return 0;
}
*/
