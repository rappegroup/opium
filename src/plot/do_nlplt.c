/*
 * $Id: do_nlplt.c,v 1.3 2004/06/16 21:25:54 mbarnes Exp $
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>


#include "parameter.h"        /* defines structure: 'param_t' */
#include "fortparam.h"
#include "common_blocks.h"
#include "do_nlplt.h"           /* the module's own header */
#include "nlm.h"

#define streq(a,b) (*a==*b && !strcmp(a+1,b+1))

int do_nlplt(param_t *param, char *logfile){

  int i,j,k;
  int scount=0;
  int pcount=0;
  int dcount=0;
  int fcount=0;
  int lcolor=0;
  int lsty=0;
  int nc;
  int ncore;

  FILE *parm;
  char *comm;
  char lc,sc;

  comm= (char *) malloc(120*sizeof(char));

  ncore=param->norb-param->nval;
  if ((parm = fopen("nl.par","w")) != NULL) {
    
    fprintf(parm,"# NL par file for xmgrace\n");
    fprintf(parm,"g0 on\n");
    fprintf(parm,"with g0\n");
    fprintf(parm,"    world xmin 0\n");
    fprintf(parm,"    world xmax 5\n");
    fprintf(parm,"    world ymin -3\n");
    fprintf(parm,"    world ymax 3\n");
    fprintf(parm,"    title \"Valence wavefunctions: %s\"\n",param->symbol);
    fprintf(parm,"    title font 0\n");
    fprintf(parm,"    title size 1.500000\n");
    fprintf(parm,"    title color 1\n");
    fprintf(parm,"    subtitle \"atom name: %s\"\n",param->name);
    fprintf(parm,"    subtitle font 0\n");
    fprintf(parm,"    subtitle size 1.000000\n");
    fprintf(parm,"    subtitle color 1\n");
    fprintf(parm,"    xaxis  on\n");
    fprintf(parm,"    xaxis  tick major 1\n");
    fprintf(parm,"    xaxis  tick minor 0.5\n");
    fprintf(parm,"    xaxis  label \"R (Bohr)\"\n");
    fprintf(parm,"    yaxis  on\n");
    fprintf(parm,"    yaxis  tick major 1 \n");
    fprintf(parm,"    yaxis  tick minor 0.5 \n");
    fprintf(parm,"    yaxis  label \"r * Psi\"\n");
    fprintf(parm,"    legend on\n");
    fprintf(parm,"    legend loctype view\n");
    fprintf(parm,"    legend 0.85, 0.8\n");
    
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
	fprintf(parm," s%d legend \"%d%c AE\" \n",i,
		nlm_label(param->nlm[i+ncore]).n,lc);
      }
    } else {
      
      for (i=ncore; i<param->norb;i++){
	
	if (nlm_label(param->nlm[i]).l == 0) {
	  scount++;
	  lcolor=1;
	  lsty=scount;
	  lc='s';
	}else if (nlm_label(param->nlm[i]).l == 1) {
	  pcount++;
	  lcolor=2;
	  lsty=pcount;
	  lc='p';
	}else if (nlm_label(param->nlm[i]).l == 2) {
	  dcount++;
	  lcolor=3;
	  lsty=dcount;
	  lc='d';
	}else if (nlm_label(param->nlm[i]).l == 3) {
	  fcount++;
	  lcolor=4;
	  lsty=fcount;
	  lc='f';
	}
	k=i-ncore;
	fprintf(parm," s%d hidden false \n",k);
	fprintf(parm," s%d type xy \n",k);
	fprintf(parm," s%d symbol 0 \n",k);
	fprintf(parm," s%d line type 1 \n",k);
	fprintf(parm," s%d line linestyle %d \n",k,lsty+5);
	fprintf(parm," s%d line linewidth 3.0 \n",k);
	fprintf(parm," s%d line color %d \n",k,lcolor);
	fprintf(parm," s%d legend \"%d%c Rel_ave.\" \n",k,
		nlm_label(param->nlm[i]).n,lc);

      }
      k++;
      for (i=0; i<atmwave_.numorb;i++){
	
	if (atmwave_.llo[i] == 0) {
	  scount++;
	  lcolor=1;
	  lsty=scount;
	  lc='s';
	  sc=' ';
	} else if (atmwave_.llo[i] == 1) {
	  pcount++;
	  lcolor=2;
	  lsty=pcount;
	  lc='p';
	  sc='-';
	  if (atmwave_.sso[i] < 0.0) sc='+';
	} else if (atmwave_.llo[i] == 2) {	
	  dcount++;
	  lcolor=3;
	  lsty=dcount;
	  lc='d';
	  sc='-';
	  if (atmwave_.sso[i] < 0.0) sc='+';
	}else if (atmwave_.llo[i] == 3) {	
	  fcount++;
	  lcolor=4;
	  lsty=fcount;
	  lc='f';
	  sc='-';
	  if (atmwave_.sso[i] < 0.0) sc='+';
	}
	fprintf(parm," s%d hidden false \n",i+k);
	fprintf(parm," s%d type xy \n",i+k);
	fprintf(parm," s%d symbol 0 \n",i+k);
	fprintf(parm," s%d line type 1 \n",i+k);
	fprintf(parm," s%d line linestyle %d \n",i+k,lsty);
	fprintf(parm," s%d line linewidth 2.0 \n",i+k);
	fprintf(parm," s%d line color %d \n",i+k,lcolor);
	fprintf(parm," s%d legend \"%d%c%c AE\" \n",i+k,
		atmwave_.nno[i],lc,sc);
      }
      
      /*
	nc=i;
	
	for (i=0; i<atmwave_.numorb;i++){
	
	if (atmwave_.llo[i] == 0) {
	scount++;
	lcolor=1;
	lsty=scount;
	lc='s';
	sc=' ';
	} else if (atmwave_.llo[i] == 1) {
	pcount++;
	lcolor=2;
	lsty=pcount;
	lc='p';
	sc='-';
	if (atmwave_.sso[i] < 0.0) sc='+';
	} else if (atmwave_.llo[i] == 2) {	
	dcount++;
	lcolor=3;
	lsty=dcount;
	lc='d';
	sc='-';
	if (atmwave_.sso[i] < 0.0) sc='+';
	}else if (atmwave_.llo[i] == 3) {	
	fcount++;
	lcolor=4;
	lsty=fcount;
	lc='f';
	sc='-';
	if (atmwave_.sso[i] < 0.0) sc='+';
	}
	fprintf(parm," s%d hidden false \n",i+nc);
	fprintf(parm," s%d type xy \n",i+nc);
	fprintf(parm," s%d symbol 0 \n",i+nc);
	fprintf(parm," s%d line type 1 \n",i+nc);
	fprintf(parm," s%d line linestyle %d \n",i+nc,lsty);
	fprintf(parm," s%d line linewidth 2.0 \n",i+nc);
	fprintf(parm," s%d line color %d \n",i+nc,lcolor);
	fprintf(parm," s%d legend \"%d%c%c\" \n",i+nc,
	atmwave_.nno[i],lc,sc);
	}
      */ 
    }
    
    k=i;
    
    for (i=0; i<param->nval;i++) {
      
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

      fprintf(parm," s%d hidden false \n",k+i);
      fprintf(parm," s%d type xy \n",k+i);
      fprintf(parm," s%d symbol 0 \n",k+i);
      fprintf(parm," s%d line type 1 \n",k+i);
      fprintf(parm," s%d line linestyle %d \n",k+i,lsty);
      fprintf(parm," s%d line linewidth 2.0 \n",k+i);
      fprintf(parm," s%d line color %d \n",k+i,lcolor);
      fprintf(parm," s%d legend \"%d%c NL\" \n",k+i,
	      nlm_label(param->nlm[i+ncore]).n,lc);
    }


    fclose(parm);
    
  } 
    
  sprintf(comm, "xmgrace %s.plt_ae %s.plt_nl -p nl.par -saveall %s.nl_agr & ", 
	  param->name,param->name,param->name);
  system(comm);
  free(comm);
  
  return 0;
}

