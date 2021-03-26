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
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "parameter.h"        /* defines structure: 'param_t' */
#include "cdim.h"
#include "common_blocks.h"
#include "do_wplot.h"           /* the module's own header */
#include "nlm.h"

#define streq(a,b) (!strcasecmp(a,b))

int do_wplot(param_t *param, char *logfile, char *pltyp){

  int i,j,k=0;
  int scount=0;
  int pcount=0;
  int dcount=0;
  int fcount=0;
  int lcolor=0;
  int lsty=0;
  int ncore;
  int lwid=2;

  FILE *parm;
  char *comm;
  char lc=0,sc=0;
  char *xc;
  char *met;
  
  #define comm_size 240
  #define met_size 45
  #define xc_size 45
  comm= (char *) malloc(comm_size*sizeof(char));
  met= (char *) malloc(met_size*sizeof(char));
  xc= (char *) malloc(xc_size*sizeof(char));

  ncore=param->norb-param->nval;
  
  if ((parm = fopen("w.par","w")) != NULL) {
      
    fprintf(parm,"# wavefunction param file for xmgrace\n");
    
    if (streq(pltyp,"n")) {
      fprintf(parm,"title \"Valence wavefunctions for %s\"\n",param->symbol);
      if (param->psmeth=='o') {
	fprintf(parm,"subtitle \"Optimized Pseudopotential Method\"\n");
      }else if (param->psmeth=='k') {
	fprintf(parm,"subtitle \"Kerker Pseudopotential Method\"\n");
      }else if (param->psmeth=='t') {
	fprintf(parm,"subtitle \"Troullier-Martins Pseudopotential Method\"\n");
      }
    } else {
      fprintf(parm,"title \"All electron wavefunctions for %s\"\n",param->symbol);
    }
    if (streq(pltyp,"n")) {
      if (param->psmeth=='o') {
	snprintf(met,met_size,"Optimized Pseudopotential Method");
      }else if (param->psmeth=='k') {
	snprintf(met,met_size,"Kerker Pseudopotential Method");
      }else if (param->psmeth=='t') {
	snprintf(met,met_size,"Troullier-Martins Pseudopotential Method");
      }
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

    fprintf(parm,"g0 on\n");
    fprintf(parm,"with g0\n");
    fprintf(parm,"world xmin 0\n");
    fprintf(parm,"world xmax 5\n");
    fprintf(parm,"title font 0\n");
    fprintf(parm,"title size 1.500000\n");
    fprintf(parm,"title color 1\n");
    fprintf(parm,"subtitle font 0\n");
    fprintf(parm,"subtitle size 1.000000\n");
    fprintf(parm,"subtitle color 1\n");
    fprintf(parm,"xaxis  on\n");
    fprintf(parm,"xaxis  tick major 1\n");
    fprintf(parm,"xaxis  tick minor 0.5\n");
    fprintf(parm,"xaxis  label \"r (a.u.)\"\n");
    fprintf(parm,"yaxis  on\n");
    fprintf(parm,"yaxis  label \"\\+r \\xy\\f{}(r)\"\n");
    fprintf(parm,"legend on\n");
    fprintf(parm,"legend loctype view\n");
    fprintf(parm,"legend 0.85, 0.8\n");

    /* Do NRL AE ; same steps if HF or DFT */    
    if (streq(param->reltype, "nrl")){  
      for (i=0; i<param->nval;i++){
	if (nlm_label(param->nlm[i+ncore]).l == 0) {
	  scount++;
	  lcolor=1;
	  lsty=scount;
	  lc='s';
	}else if (nlm_label(param->nlm[i+ncore]).l == 1) {
	  pcount++;
	  lcolor=2;
	  lsty=pcount;
	  lc='p';
	}else if (nlm_label(param->nlm[i+ncore]).l == 2) {
	  dcount++;
	  lcolor=3;
	  lsty=dcount;
	  lc='d';
	}else if (nlm_label(param->nlm[i+ncore]).l == 3) {
	  fcount++;
	  lcolor=4;
	  lsty=fcount;
	  lc='f';
	}
	fprintf(parm," s%d hidden false \n",i);
	fprintf(parm," s%d type xy \n",i);
	fprintf(parm," s%d symbol 0 \n",i);
	fprintf(parm," s%d line type 1 \n",i);
	fprintf(parm," s%d line linestyle %d \n",i,lsty);
	fprintf(parm," s%d line linewidth 2.0 \n",i);
	fprintf(parm," s%d line color %d \n",i,lcolor);
	fprintf(parm," s%d legend \"\\xy\\f{}\\s%d%c\\N\\S\\-AE\\N\" \n",i,
		nlm_label(param->nlm[i+ncore]).n,lc);
      }
      k=i;
    } else {

      /* Now plot Rel AE wfns */
      for (i=0;i<param->nval;i++){
	if (nlm_label(param->nlm[i+ncore]).l == 0) {
	  scount++;
	  lcolor=1;
	  lsty=scount;
	  lc='s';
	  sc=' ';
	  fprintf(parm," s%d hidden false \n",k);
	  fprintf(parm," s%d type xy \n",k);
	  fprintf(parm," s%d symbol 0 \n",k);
	  fprintf(parm," s%d line type 1 \n",k);
	  fprintf(parm," s%d line linestyle %d \n",k,lsty);
	  fprintf(parm," s%d line linewidth 2.0 \n",k);
	  fprintf(parm," s%d line color %d \n",k,lcolor);
	  fprintf(parm," s%d legend \"\\xy\\f{}\\s%d%c%c\\N\\S\\-AE\\N\" \n",
		  k,nlm_label(param->nlm[i+ncore]).n,lc,sc);
	  k++;
	  
	} else {
	  for (j=0;j<2;j++) {
	    if (j == 0) sc='-';
	    if (j == 1) sc='+';


	    if (nlm_label(param->nlm[i+ncore]).l == 1) {	  
	      pcount++;
	      lcolor=2;
	      lsty=1;
	      lwid=2.0;
	      if (j==1) {
		lwid=4.0;
		lsty=3;
	      }

	      lc='p';
	      
	    }else if (nlm_label(param->nlm[i+ncore]).l == 2) {	  
	      dcount++;
	      lcolor=3;
	      lsty=1;
	      lwid=2.0;
	      if (j==1) {
		lwid=4.0;
		lsty=3;
	      }
	      lc='d';

	    }else if (nlm_label(param->nlm[i+ncore]).l == 3) {	  
	      fcount++;
	      lcolor=4;
	      lsty=1;
	      lwid=2.0;
	      if (j==1) {
		lwid=4.0;
		lsty=3;
	      }
	      lc='f';

	    }

	    fprintf(parm," s%d hidden false \n",k);
	    fprintf(parm," s%d type xy \n",k);
	    fprintf(parm," s%d symbol 0 \n",k);
	    fprintf(parm," s%d line type 1 \n",k);
	    fprintf(parm," s%d line linestyle %d \n",k,lsty);
	    fprintf(parm," s%d line linewidth %d\n",k,lwid);
	    fprintf(parm," s%d line color %d \n",k,lcolor);
	    fprintf(parm," s%d legend \"\\xy\\f{}\\s%d%c%c\\N\\S\\-AE\\N\" \n",
		    k,nlm_label(param->nlm[i+ncore]).n,lc,sc);
	    k++;
	  }
	}

      }

      if ((param->ixc >= 0 )&&(streq(param->reltype, "srl"))){	  
	/* if DFT, there is an average wfn to plot */
	for (i=0; i<param->nval;i++){
	  if (nlm_label(param->nlm[i+ncore]).l == 0) {
	    scount++;
	    lcolor=1;
	    lsty=scount;
	    lc='s';
	  }else if (nlm_label(param->nlm[i+ncore]).l == 1) {
	    pcount++;
	    lcolor=2;
	    lsty=pcount;
	    lc='p';
	  }else if (nlm_label(param->nlm[i+ncore]).l == 2) {
	    dcount++;
	    lcolor=3;
	    lsty=dcount;
	    lc='d';
	  }else if (nlm_label(param->nlm[i+ncore]).l == 3) {
	    fcount++;
	    lcolor=4;
	    lsty=fcount;
	    lc='f';
	  }
	  fprintf(parm," s%d hidden false \n",k);
	  fprintf(parm," s%d type xy \n",k);
	  fprintf(parm," s%d symbol 0 \n",k);
	  fprintf(parm," s%d line type 1 \n",k);
	  fprintf(parm," s%d line linestyle %d \n",k,lsty+3);
	  fprintf(parm," s%d line linewidth 5.0 \n",k);
	  fprintf(parm," s%d line color %d \n",k,lcolor);
	  fprintf(parm," s%d legend \"\\xy\\f{}\\s%d%c\\S\\-avg\\N\" \n",k,
		  nlm_label(param->nlm[i+ncore]).n,lc);
	  k++;
      	}
      }
    }

    /* now do the pseudowfns if asked for */    
    if (streq(pltyp,"n")) {
      if (!streq(param->reltype, "frl")){  
	
	for (i=0; i<param->nval;i++) {
	  
	  if (nlm_label(param->nlm[i+ncore]).l == 0) {
	    scount++;
	    lcolor=1;
	    lsty=4;
	    lc='s';
	  }else if (nlm_label(param->nlm[i+ncore]).l == 1) {
	    pcount++;
	    lcolor=2;
	    lsty=4;
	    lc='p';
	  }else if (nlm_label(param->nlm[i+ncore]).l == 2) {
	    dcount++;
	    lcolor=3;
	    lsty=4;
	    lc='d';
	  }else if (nlm_label(param->nlm[i+ncore]).l == 3) {
	    fcount++;
	    lcolor=4;
	    lsty=4;
	    lc='f';
	  }
	  
	  fprintf(parm," s%d hidden false \n",k+i);
	  fprintf(parm," s%d type xydx \n",k+i);
	  fprintf(parm," s%d symbol 0 \n",k+i);
	  fprintf(parm," s%d line type 1 \n",k+i);
	  fprintf(parm," s%d line linestyle %d\n",k+i,lsty);
	  fprintf(parm," s%d line linewidth 2.0 \n",k+i);
	  fprintf(parm," s%d line color %d \n",k+i,lcolor);
	  fprintf(parm," s%d legend \"\\xy\\f{}\\s%d%c\\N\\S\\-NL\\N\" \n",i+k,nlm_label(param->nlm[i+ncore]).n,lc);
	  fprintf(parm," s%d errorbar on\n",k+i);
	  fprintf(parm," s%d errorbar place both\n",k+i);
	  fprintf(parm," s%d errorbar color %d\n",k+i,lcolor);
	  fprintf(parm," s%d errorbar pattern 1\n",k+i);
	  fprintf(parm," s%d errorbar size 10.000000\n",k+i);
	  fprintf(parm," s%d errorbar linewidth 2.0\n",k+i);
	  fprintf(parm," s%d errorbar linestyle 1\n",k+i);
	}
      }else{
	for (i=0; i<param->nval;i++) {
	  if (nlm_label(param->nlm[i+ncore]).l == 0) {
	    scount++;
	    lcolor=1;
	    lsty=7;
	    lc='s';
	    sc=' ';
	    fprintf(parm," s%d hidden false \n",k);
	    fprintf(parm," s%d type xy \n",k);
	    fprintf(parm," s%d symbol 0 \n",k);
	    fprintf(parm," s%d line type 1 \n",k);
	    fprintf(parm," s%d line linestyle %d \n",k,lsty);
	    fprintf(parm," s%d line linewidth 3.0 \n",k);
	    fprintf(parm," s%d line color %d \n",k,lcolor);
	    fprintf(parm," s%d legend \"\\xy\\f{}\\s%d%c%c\\N\\S\\-NL\\N\" \n",
		    k,nlm_label(param->nlm[i+ncore]).n,lc,sc);
	    k++;
	    
	  } else {
	    
	    for (j=0;j<2;j++) {
	      if (j == 0) sc='-';
	      if (j == 1) sc='+';
	      
	      if (nlm_label(param->nlm[i+ncore]).l == 1) {	  

		pcount++;
		lcolor=2;
		lsty=7;
		lwid=3.0;
		if (j==1) {
		  lwid=5.0;
		  lsty=3;
		}
		lc='p';
	      
	      }else if (nlm_label(param->nlm[i+ncore]).l == 2) {	  
		
		dcount++;
		lcolor=3;
		lsty=7;
		lwid=3.0;
		if (j==1) {
		  lwid=5.0;
		  lsty=3;
		}
		lc='d';
		
	      }else if (nlm_label(param->nlm[i+ncore]).l == 3) {	  

		fcount++;
		lcolor=4;
		lsty=7;
		lwid=3.0;
		if (j==1) {
		  lwid=5.0;
		  lsty=3;
		}
		lc='f';

	      }
	      
	      fprintf(parm," s%d hidden false \n",k);
	      fprintf(parm," s%d type xy \n",k);
	      fprintf(parm," s%d symbol 0 \n",k);
	      fprintf(parm," s%d line type 1 \n",k);
	      fprintf(parm," s%d line linestyle %d \n",k,lsty);
	      fprintf(parm," s%d line linewidth 3.0 \n",k);
	      fprintf(parm," s%d line color %d \n",k,lcolor);
	      fprintf(parm," s%d legend \"\\xy\\f{}\\s%d%c%c\\N\\S\\-NL\\N\" \n",
		      k,nlm_label(param->nlm[i+ncore]).n,lc,sc);
	      k++;
	    }
	  }
	  
	}
      }
      
      fclose(parm);
      
    }
  }

  if (streq(pltyp,"n")) {    
    snprintf(comm, comm_size,"xmgrace -timestamp $XMGRACE_OPTS %s.ae_plt -settype xydx %s.nl_plt -p w.par -autoscale y -saveall %s_nl.agr & ", 
	    param->name,param->name,param->name);
  } else {

    snprintf(comm, comm_size,"xmgrace -timestamp $XMGRACE_OPTS %s.ae_plt -p w.par -autoscale y -saveall %s_ae.agr & ", 
	    param->name,param->name);
  }

  system(comm);
  free(comm);

  return 0;
}

