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

/* screened or ionic pseudo potential plotting with xmgrace */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>


#include "parameter.h"
#include "cdim.h"
#include "common_blocks.h"
#include "do_vplot.h"
#include "do_ps.h"
#include "nlm.h"

#define streq(a,b) (!strcasecmp(a,b))
#define MAX(a, b)   (((a) > (b)) ? (a):(b))
void nrelsproj(param_t *param, char *);
int do_vplot(param_t *param, char *logfile, char *pltyp){

  int i,j,k=0,kk=0;
  int scount=0;
  int pcount=0;
  int dcount=0;
  int fcount=0;
  int lcolor=0;
  int lsty=0;
  int lwid=2.0;
  int ncore;

  FILE *parm;
  char *comm;
  char *xc;
  char *met;
  char lc=0,sc=0;

  #define comm_size 240
  #define met_size 45
  #define xc_size 45
  comm= (char *) malloc(comm_size*sizeof(char));
  met= (char *) malloc(met_size*sizeof(char));
  xc= (char *) malloc(xc_size*sizeof(char));

  if ((!strcmp(param->reltype, "nrl")) || (!strcmp(param->reltype, "srl")) && (param->ixc >= 0)) {
    nrelsproj(param,logfile);
  } else {
    relsproj(param,logfile);
  }
  ncore=param->norb-param->nval;

  readAE(param);

  if (streq(pltyp,"s")) {

  /* PLOT THE SCREENED POTENTIAL */

    parm = fopen("vs.par","w");

    if (!parm) {
      printf("Could not open grace param file 'vs.par': %s\n",
	  strerror(errno));
      return 1;
    }

    fprintf(parm,"# screened potential param file for xmgrace\n");
    fprintf(parm,"g0 on\n");
    fprintf(parm,"with g0\n");
    fprintf(parm,"world xmin 0\n");
    fprintf(parm,"world xmax 5\n");
    fprintf(parm,"title \"Screened pseudopotential for %s\"\n",param->symbol);
    fprintf(parm,"world ymin -50\n");
    fprintf(parm,"world ymax 50\n");
    fprintf(parm,"title font 0\n");
    fprintf(parm,"title size 1.500000\n");
    fprintf(parm,"title color 1\n");


    if (param->psmeth=='o') {
      snprintf(met,met_size,"Optimized Pseudopotential Method");
    }else if (param->psmeth=='k') {
      snprintf(met,met_size,"Kerker Pseudopotential Method");
    }else if (param->psmeth=='t') {
      snprintf(met,met_size,"Troullier-Martins Pseudopotential Method");
    }
    
    if (param->ixc == 0) {
      snprintf(xc,xc_size,"XC=Perdew-Zunger LDA");
    }else if (param->ixc == 1) {
      snprintf(xc,xc_size,"XC=Perdew-Wang LDA");
    }else if (param->ixc == 2) {
      snprintf(xc,xc_size,"XC=Perdew-Burke-Ernzerhof GGA");
    }else if (param->ixc == 3) {
      snprintf(xc,xc_size,"XC=Perdew-Wang GGA");
    }else if (param->ixc == 4) {    
      snprintf(xc,xc_size,"XC=Wu-Cohen GGA");
    }else if (param->ixc == 5) {    
      snprintf(xc,xc_size,"XC=PBESol GGA");
    }else if (param->ixc == 6) {    
      snprintf(xc,xc_size,"XC=VWN5 LDA");
    }else if (param->ixc == -1) {
      snprintf(xc,xc_size,"Hartree-Fock Exchange");
    }

    fprintf(parm,"subtitle \"%s   %s   %s\"\n",xc,met,param->reltype);

    fprintf(parm,"subtitle font 0\n");
    fprintf(parm,"subtitle size 1.000000\n");
    fprintf(parm,"subtitle color 1\n");
    fprintf(parm,"xaxis on\n");
    fprintf(parm,"xaxis tick major 1\n");
    fprintf(parm,"xaxis tick minor 0.5\n");
    fprintf(parm,"xaxis label \"r (a.u.)\"\n");
    fprintf(parm,"yaxis on\n");
    fprintf(parm,"yaxis tick major 10 \n");
    fprintf(parm,"yaxis tick minor 2 \n");
    fprintf(parm,"yaxis label \"V\\sscr\\N (Ry)\"\n");
    fprintf(parm,"legend on\n");
    fprintf(parm,"legend loctype view\n");
    fprintf(parm,"legend 0.85, 0.8\n");

    if (!streq(param->reltype, "frl")){  
      for (i=0; i<param->nval;i++){
	
	if ((param->lpot[i]) == 0) {
	  scount++;
	  lcolor=1;
	  lsty=scount;
	  lc='s';
	}else if ((param->lpot[i]) == 1) {
	  pcount++;
	  lcolor=2;
	  lsty=pcount;
	  lc='p';
	}else if ((param->lpot[i]) == 2) {
	  dcount++;
	  lcolor=3;
	  lsty=dcount;
	  lc='d';
	}else if ((param->lpot[i]) == 3) {
	  fcount++;
	  lcolor=4;
	  lsty=fcount;
	  lc='f';
	}
	fprintf(parm," s%d hidden false \n",i);
	fprintf(parm," s%d type xydx \n",i);
	fprintf(parm," s%d symbol 0 \n",i);
	fprintf(parm," s%d line type 1 \n",i);
	fprintf(parm," s%d line linestyle %d \n",i,lsty);
	fprintf(parm," s%d line linewidth 2.0 \n",i);
	fprintf(parm," s%d line color %d \n",i,lcolor);
	fprintf(parm," s%d legend \"V\\s%c\\N\" \n",i,lc);
	fprintf(parm," s%d errorbar on\n",i);
	fprintf(parm," s%d errorbar place both\n",i);
	fprintf(parm," s%d errorbar color %d\n",i,lcolor);
	fprintf(parm," s%d errorbar pattern 1\n",i);
	fprintf(parm," s%d errorbar size 10.000000\n",i);
	fprintf(parm," s%d errorbar linewidth 2.0\n",i);
	fprintf(parm," s%d errorbar linestyle 1\n",i);
      }
    }else{
      for (i=0; i<param->nval;i++) {
	if (nlm_label(param->nlm[i+ncore]).l == 0) {
	  scount++;
	  lcolor=1;
	  lsty=scount;
	  lc='s';
	  sc=' ';
	  fprintf(parm," s%d hidden false \n",k);
	  fprintf(parm," s%d type xydx \n",k);
	  fprintf(parm," s%d symbol 0 \n",k);
	  fprintf(parm," s%d line type 1 \n",k);
	  fprintf(parm," s%d line linestyle %d \n",k,lsty);
	  fprintf(parm," s%d line linewidth 2.0 \n",k);
	  fprintf(parm," s%d line color %d \n",k,lcolor);
	  fprintf(parm," s%d legend \"V\\s%c\\N\" \n",k,lc);
	  fprintf(parm," s%d errorbar on\n",k);
	  fprintf(parm," s%d errorbar place both\n",k);
	  fprintf(parm," s%d errorbar color %d\n",k,lcolor);
	  fprintf(parm," s%d errorbar pattern 1\n",k);
	  fprintf(parm," s%d errorbar size 10.000000\n",k);
	  fprintf(parm," s%d errorbar linewidth 2.0\n",k);
	  fprintf(parm," s%d errorbar linestyle 1\n",k);
	  k++;
	}else{
	  for (j=0;j<2;j++) {
	    if (j == 0) sc='-';
	    if (j == 1) sc='+';
	    if (nlm_label(param->nlm[i+ncore]).l == 1) {	  
	      pcount++;
	      lcolor=2;
	      lwid=2.0;
	      lsty=1;
	      if (j==1) {
		lwid=4.0;
		lsty=3;
	      }
	      lc='p';
	      
	    }else if (nlm_label(param->nlm[i+ncore]).l == 2) {	  
	      dcount++;
	      lcolor=3;
	      lwid=2.0;
	      lsty=1;
	      if (j==1) {
		lwid=4.0;
		lsty=3;
	      }
	      lc='d';
	      
	    }else if (nlm_label(param->nlm[i+ncore]).l == 3) {	  
	      fcount++;
	      lcolor=4;
	      lwid=2.0;
	      lsty=1;
	      if (j==1) {
		lwid=4.0;
		lsty=3;
	      }
	      lc='f';
	      
	    }
	    fprintf(parm," s%d hidden false \n",k);
	    fprintf(parm," s%d type xydx \n",k);
	    fprintf(parm," s%d symbol 0 \n",k);
	    fprintf(parm," s%d line type 1 \n",k);
	    fprintf(parm," s%d line linestyle %d \n",k,lsty);
	    fprintf(parm," s%d line linewidth %d \n",k,lwid);
	    fprintf(parm," s%d line color %d \n",k,lcolor);
	    fprintf(parm," s%d legend \"V\\s%c%c\\N\" \n",k,lc,sc);
	    fprintf(parm," s%d errorbar on\n",k);
	    fprintf(parm," s%d errorbar place both\n",k);
	    fprintf(parm," s%d errorbar color %d\n",k,lcolor);
	    fprintf(parm," s%d errorbar pattern 1\n",k);
	    fprintf(parm," s%d errorbar size 10.000000\n",k);
	    fprintf(parm," s%d errorbar linewidth 2.0\n",k);
	    fprintf(parm," s%d errorbar linestyle 1\n",k);
	    k++;
	  }
	}
      }
    }
    
    if (param->ixc >= 0) {
      kk=param->nll;
      fprintf(parm," s%d hidden false \n",kk);
      fprintf(parm," s%d type xydx \n",kk);
      fprintf(parm," s%d symbol 0 \n",kk);
      fprintf(parm," s%d line type 1 \n",kk);
      fprintf(parm," s%d line linestyle %d \n",kk,3);
      fprintf(parm," s%d line linewidth 3.0 \n",kk);
      fprintf(parm," s%d line color %d \n",kk,9);
      fprintf(parm," s%d legend \"V\\sloc\\N\"\n" ,kk);
    }
    
    fclose(parm);

    snprintf(comm,comm_size,"xmgrace -timestamp $XMG_OPTS -settype xydx %s.vs_plt -p vs.par -autoscale y -saveall %s_vs.agr & ",
	param->name,param->name);
    system(comm);

  } else {

  /* PLOT THE DESCREENED (i.e. IONIC) POTENTIAL */

    parm = fopen("vi.par","w");

    if (!parm) {
      printf("Could not open grace param file 'vi.par': %s\n",
	  strerror(errno));
      return 1;
    }

    fprintf(parm,"# IPS par file for xmgrace\n");
    fprintf(parm,"g0 on\n");
    fprintf(parm,"with g0\n");
    fprintf(parm,"world xmin 0\n");
    fprintf(parm,"world xmax 5\n");
    fprintf(parm,"world ymin -50\n");
    fprintf(parm,"world ymax 50\n");
    fprintf(parm,"title \"Ionic pseudopotential for %s\"\n",param->symbol);
    fprintf(parm,"title font 0\n");
    fprintf(parm,"title size 1.500000\n");
    fprintf(parm,"title color 1\n");


    if (param->psmeth=='o') {
      snprintf(met,met_size,"Optimized Pseudopotential Method");
    }else if (param->psmeth=='k') {
      snprintf(met,met_size,"Kerker Pseudopotential Method");
    }else if (param->psmeth=='t') {
      snprintf(met,met_size,"Troullier-Martins Pseudopotential Method");
    }
    
    if (param->ixc == 0) {
      snprintf(xc,xc_size,"XC=Perdew-Zunger LDA");
    }else if (param->ixc == 1) {
      snprintf(xc,xc_size,"XC=Perdew-Wang LDA");
    }else if (param->ixc == 2) {
      snprintf(xc,xc_size,"XC=Perdew-Burke-Ernzerhof GGA");
    }else if (param->ixc == 3) {
      snprintf(xc,xc_size,"XC=Perdew-Wang GGA");
    }else if (param->ixc == 4) {    
      snprintf(xc,xc_size,"XC=Wu-Cohen GGA");
    }else if (param->ixc == 5) {    
      snprintf(xc,xc_size,"XC=PBESol GGA");
    }else if (param->ixc == 6) {    
      snprintf(xc,xc_size,"XC=VWN5 LDA");
    }else if (param->ixc == -1) {
      snprintf(xc,xc_size,"Hartree-Fock Exchange");
    }

    fprintf(parm,"subtitle \"%s   %s   %s\"\n",xc,met,param->reltype);

    fprintf(parm,"subtitle font 0\n");
    fprintf(parm,"subtitle size 1.000000\n");
    fprintf(parm,"subtitle color 1\n");
    fprintf(parm,"xaxis  on\n");
    fprintf(parm,"xaxis  tick major 1\n");
    fprintf(parm,"xaxis  tick minor 0.5\n");
    fprintf(parm,"xaxis  label \"r (a.u.)\"\n");
    fprintf(parm,"yaxis  on\n");
    fprintf(parm,"yaxis  tick major 10 \n");
    fprintf(parm,"yaxis  tick minor 2 \n");
    fprintf(parm,"yaxis  label \"V\\sion\\N (Ry)\"\n");
    fprintf(parm,"legend on\n");
    fprintf(parm,"legend loctype view\n");
    fprintf(parm,"legend 0.85, 0.8\n");

    if (!streq(param->reltype, "frl")){  
      for (i=0; i<param->nval;i++){
	
	if (param->lpot[i] == 0) {
	  scount++;
	  lcolor=1;
	  lsty=scount;
	  lc='s';
	}else if (param->lpot[i] == 1) {
	  pcount++;
	  lcolor=2;
	  lsty=pcount;
	  lc='p';
	}else if (param->lpot[i] == 2) {
	  dcount++;
	  lcolor=3;
	  lsty=dcount;
	  lc='d';
	}else if (param->lpot[i] == 3) {
	  fcount++;
	  lcolor=4;
	  lsty=fcount;
	  lc='f';
	}
	fprintf(parm," s%d hidden false \n",i);
	fprintf(parm," s%d type xydx \n",i);
	fprintf(parm," s%d symbol 0 \n",i);
	fprintf(parm," s%d line type 1 \n",i);
	fprintf(parm," s%d line linestyle %d \n",i,lsty);
	fprintf(parm," s%d line linewidth 2.0 \n",i);
	fprintf(parm," s%d line color %d \n",i,lcolor);
	fprintf(parm," s%d legend \"V\\s%c\\N\" \n",i,lc);
	fprintf(parm," s%d errorbar on\n",i);
	fprintf(parm," s%d errorbar place both\n",i);
	fprintf(parm," s%d errorbar color %d\n",i,lcolor);
	fprintf(parm," s%d errorbar pattern 1\n",i);
	fprintf(parm," s%d errorbar size 10.000000\n",i);
	fprintf(parm," s%d errorbar linewidth 2.0\n",i);
	fprintf(parm," s%d errorbar linestyle 1\n",i);
      }
    }else{
      for (i=0; i<param->nval;i++) {

	if (nlm_label(param->nlm[i+ncore]).l == 0) {
	  scount++;
	  lcolor=1;
	  lsty=scount;
	  lc='s';
	  sc=' ';
	  fprintf(parm," s%d hidden false \n",k);
	  fprintf(parm," s%d type xydx \n",k);
	  fprintf(parm," s%d symbol 0 \n",k);
	  fprintf(parm," s%d line type 1 \n",k);
	  fprintf(parm," s%d line linestyle %d \n",k,lsty);
	  fprintf(parm," s%d line linewidth 2.0 \n",k);
	  fprintf(parm," s%d line color %d \n",k,lcolor);
	  fprintf(parm," s%d legend \"V\\s%c\\N\" \n",k,lc);
	  fprintf(parm," s%d errorbar on\n",k);
	  fprintf(parm," s%d errorbar place both\n",k);
	  fprintf(parm," s%d errorbar color %d\n",k,lcolor);
	  fprintf(parm," s%d errorbar pattern 1\n",k);
	  fprintf(parm," s%d errorbar size 10.000000\n",k);
	  fprintf(parm," s%d errorbar linewidth 2.0\n",k);
	  fprintf(parm," s%d errorbar linestyle 1\n",k);
	  k++;
	}else{
	  for (j=0;j<2;j++) {
	    if (j == 0) sc='-';
	    if (j == 1) sc='+';
	    if (nlm_label(param->nlm[i+ncore]).l == 1) {	  
	      pcount++;
	      lcolor=2;
	      lwid=2.0;
	      lsty=1;
	      if (j==1) {
		lwid=4.0;
		lsty=3;
	      }
	      lc='p';
	      
	    }else if (nlm_label(param->nlm[i+ncore]).l == 2) {	  
	      dcount++;
	      lcolor=3;
	      lwid=2.0;
	      lsty=1;
	      if (j==1) {
		lwid=4.0;
		lsty=3;
	      }
	      lc='d';
	      
	    }else if (nlm_label(param->nlm[i+ncore]).l == 3) {	  
	      fcount++;
	      lcolor=4;
	      lwid=2.0;
	      lsty=1;
	      if (j==1) {
		lwid=4.0;
		lsty=3;
	      }
	      lc='f';
	      
	    }
	    
	    fprintf(parm," s%d hidden false \n",k);
	    fprintf(parm," s%d type xydx \n",k);
	    fprintf(parm," s%d symbol 0 \n",k);
	    fprintf(parm," s%d line type 1 \n",k);
	    fprintf(parm," s%d line linestyle %d \n",k,lsty);
	    fprintf(parm," s%d line linewidth %d \n",k,lwid);
	    fprintf(parm," s%d line color %d \n",k,lcolor);
	    fprintf(parm," s%d legend \"V\\s%c%c\\N\" \n",k,lc,sc);
	    fprintf(parm," s%d errorbar on\n",k);
	    fprintf(parm," s%d errorbar place both\n",k);
	    fprintf(parm," s%d errorbar color %d\n",k,lcolor);
	    fprintf(parm," s%d errorbar pattern 1\n",k);
	    fprintf(parm," s%d errorbar size 10.000000\n",k);
	    fprintf(parm," s%d errorbar linewidth 2.0\n",k);
	    fprintf(parm," s%d errorbar linestyle 1\n",k);
	    k++;

	  }
	}
      }
    }
    /* add the AE ionic potential */

    kk=param->nll;
    fprintf(parm," s%d hidden false \n",kk);
    fprintf(parm," s%d type xydx \n",kk);
    fprintf(parm," s%d symbol 0 \n",kk);
    fprintf(parm," s%d line type 1 \n",kk);
    fprintf(parm," s%d line linestyle %d \n",kk,3);
    fprintf(parm," s%d line linewidth 2.0 \n",kk);
    fprintf(parm," s%d line color %d \n",kk,15);
    fprintf(parm," s%d legend \"2Z\\seff\\N /r\"\n" ,kk);

    if (param->ixc >= 0) {
      fprintf(parm," s%d hidden false \n",kk+1);
      fprintf(parm," s%d type xydx \n",kk+1);
      fprintf(parm," s%d symbol 0 \n",kk+1);
      fprintf(parm," s%d line type 1 \n",kk+1);
      fprintf(parm," s%d line linestyle %d \n",kk+1,3);
      fprintf(parm," s%d line linewidth 3.0 \n",kk+1);
      fprintf(parm," s%d line color %d \n",kk+1,9);
      fprintf(parm," s%d legend \"V\\sloc\\N\"\n" ,kk+1);
    }
    fclose(parm);

    snprintf(comm,comm_size, "xmgrace -timestamp $XMGRACE_OPTS -settype xydx %s.vi_plt -p vi.par -autoscale y -saveall %s_vi.agr & ",
	param->name,param->name);

    system(comm);
  }

  free(comm);
  return 0;
}

/* vim: cindent sw=2 showmatch
 */
