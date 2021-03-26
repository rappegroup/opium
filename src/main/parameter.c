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
/****************************************************************************
 * Read the parameter settings from 'fp' using 'FlexiLib' and set default    *
 * values.                                                                   *
 ****************************************************************************/

/* standard libraries */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

/* FlexiLib */
#include "flexi.h"
#include "cdim.h"
#include "nlm.h"
#include "common_blocks.h"

/* module's own header */
#include "parameter.h"

#define streq(a,b) (!strcasecmp(a,b))

double symtoz(char *sym , char *longname, double *mass);
int read_param(param_t *param, FILE *fp, FILE *fp_log){

  int i, j, k,kk;
  int n,ii,jj,mm;        
  int lc[7];
  char **ens,**etemp;                /* temp string array for eigenvalue guesses */
  char ***ensc,***etemp2; /* temp string array for eigenvalue guesses of test configs */
  double *b_start, *b_end;  /* temp arrays for box start and end values */
  char **b_unit;            /* temp string array to read the box unit */
  double r_l, r_r;          /* temp radius variables [a.u.] */
  char la[]={"spdf"};
  int ncore; 
  char loc='s';
  char pccmeth='l';
  char avgtype='b';
  int loctemp=0;
  int npot[10];
  int ntemp,nbtemp,nt1,nt2;
  double wtemp,rtemp,qtemp;
  int *ntemp2;
  double *wtemp2;

  /* Prepare for 1st flexi_gather_keys() */  
  
  /* [Atom] */
  param->symbol = (char *) malloc(160*sizeof(char));
  param->longname = (char *) malloc(160*sizeof(char));
  flexi_request_key("Atom",1,"%s %d", param->symbol, &param->norb);

  /* [Configs] */
  flexi_request_key("Configs",0,"%d", &param->nconfigs);
  /* default */
  param->nconfigs = 0;

  /* [EVTol] */
  flexi_request_key("Tol",0,"%lg %lg", &consts_.etol,&consts_.vtol);
  /* default */
  consts_.etol = 1.e-8;   
  consts_.vtol = 1.e-6;   
  consts_.maxit = 100;      

  /* [Pseudo] */
  flexi_request_key("Pseudo",1,"%d", &param->nval);

  /* [Pcc] */
  flexi_request_key("Pcc",0,"%lg", &param->rpcc);
  flexi_request_key("Pcc",0,"%c",&pccmeth);

  /* default */
  param->rpcc=0.0;
  
  /* [LogInfo] */
  flexi_request_key("Loginfo",0,"%d %lg %lg %lg", &param->ilogder,
		    &param->rphas, &param->emin, &param->emax);
  /* default */
  param->ilogder = -67;
  param->rphas = 2.00;
  param->emin = -5.;
  param->emax = 0.;
  
  /* [KBdesign] */
  flexi_request_key("KBdesign",0,"%c %d", &loc, &param->nboxes);
  /* default */
  param->local = 0;
  param->nboxes = 0;
  
  /* [XC] */
  param->xcparam = (char *) malloc(160*sizeof(char));
  flexi_request_key("XC",0,"%s",param->xcparam);
  flexi_request_key("XC",0,"%lg",&param->exccut);
  /* default */
  param->exccut=-600;
  strcpy(param->xcparam, "lda");
  
  /* [Conmax] */
  flexi_request_key("Conmax",0,"%d %lg %d", &param->switch1, &param->encsm,
		    &param->limsm);
  flexi_request_key("Conmax",0,"%d %d", &param->switch2, &param->nstep);
  flexi_request_key("Conmax",0,"%d %lg", &param->switch3, &param->tolcon);
  /* default */
  param->switch1 = 1;
  param->switch2 = 1;
  param->switch3 = 1;
  param->encsm = 1e-5;
  param->limsm = 5;
  param->nstep = 1;
  param->tolcon=1e-5;

  /* 1st flexi_gather_keys() */  
  if (flexi_gather_keys(fp)) return 1;

  /* get derived values from 1st set of keys */

  /* element name and symbol from Z */

  rpcc_.rpcc = param->rpcc;
  ipccmeth_.ipccmeth=-1;
  /*  printf("  pccmeth= %c \n",pccmeth); */
  if ((pccmeth == 'l') || (pccmeth == 'L')) {
    ipccmeth_.ipccmeth=0; 
  }else if ((pccmeth=='f') || (pccmeth=='F')) {
    ipccmeth_.ipccmeth=1;
  }else {
    printf("Can not determine partial core correction method ; must be either lfc or fuchs \n");
    exit(1);
  }
 
  param->mass=-1.0;
  param->z = symtoz(param->symbol,param->longname,&param->mass);

  param->ixc = -100;
  if (streq(param->xcparam, "pbesol")) param->ixc = 5;
  if (streq(param->xcparam, "wcgga")) param->ixc = 4;
  if (streq(param->xcparam, "pw91gga")) param->ixc = 3;
  if (streq(param->xcparam, "gga")) param->ixc = 2;
  if (streq(param->xcparam, "pbegga")) param->ixc = 2;
  if (streq(param->xcparam, "pbe")) param->ixc = 2;
  if (streq(param->xcparam, "pwlda")) param->ixc = 1;
  if (streq(param->xcparam, "pzlda")) param->ixc = 0;
  if (streq(param->xcparam, "lda")) param->ixc = 0;
  /*  if (streq(param->xcparam, "vwn5lda")) param->ixc = 6;
      if (streq(param->xcparam, "vwn5")) param->ixc = 6;*/
  if (streq(param->xcparam, "hf")) param->ixc = -1;

  if (param->ixc == -100) {
    printf(" \nXC functional not understood \n ");
    printf("Please set the [XC] keyblock in your param file\n ");
    printf("like this: \n\n ");
    printf("  [XC] \n ");
    printf("  pzlda \n\n\n ");
    printf("acceptable values are:  pzlda, pwlda, pw91gga, pbegga, wcgga, pbesol, and hf \n"); 

    exit(1);
  }    

  logarith_.rphas = param->rphas;
  logarith_.elogmin = param->emin;
  logarith_.elogmax = param->emax;

  /* conmax parameters passed to subroutine scr */
  scrconmax_.isw1 = param->switch1;
  scrconmax_.isw2 = param->switch2;
  scrconmax_.isw3 = param->switch3;
  scrconmax_.enc = param->encsm;
  scrconmax_.lim = param->limsm;
  scrconmax_.nstp = param->nstep;
  scrconmax_.tolc = param->tolcon;

  /* Prepare for 2nd flexi_gather_keys() */
  flexi_clear_keys();
  rewind(fp);

  /*EEE*/

  /* [Atom] */
  flexi_request_key("Atom",1,"%s %d", param->symbol, &param->norb);
  param->ibound = (int *)malloc(param->norb*sizeof(int));
  param->etrial = (double *)malloc(param->norb*sizeof(double));
  param->nlm = (int *)malloc(param->norb*sizeof(int));
  param->wnl = (double *)malloc(param->norb*sizeof(double));
  param->en = (double *)malloc(param->norb*sizeof(double));
  ens = (char **)malloc(param->norb*sizeof(char*));
  etemp = (char **)malloc(param->norb*sizeof(char*));
  for (i=0; i<param->norb; i++){
    ens[i] = (char *)malloc(20*sizeof(char));
    etemp[i] = (char *)malloc(20*sizeof(char));
    flexi_request_key("Atom",1,"%d %lg %s", &param->nlm[i],
		      &param->wnl[i], ens[i]);
  }

  if (flexi_gather_keys(fp)) return 1;
  
  /* Prepare for 3rd flexi_gather_keys() */
  flexi_clear_keys();
  rewind(fp);
  
  
  /* [Configs] */
  flexi_request_key("Configs",0,"%d", &param->nconfigs);
  param->nlm_conf = (int **)malloc(param->nconfigs*sizeof(int *));
  param->wnl_conf = (double **)malloc(param->nconfigs*sizeof(double *));
  param->en_conf = (double **)malloc(param->nconfigs*sizeof(double *)); 
  ensc = (char ***)malloc(param->nconfigs*sizeof(char **));
  etemp2 = (char ***)malloc(param->nconfigs*sizeof(char **));
  ntemp2 = (int *)malloc(param->nconfigs*sizeof(int));
  wtemp2 = (double *)malloc(param->nconfigs*sizeof(double));
  for (i=0; i<param->nconfigs; i++){
    param->nlm_conf[i] = (int *)malloc(param->nval*sizeof(int));
    param->wnl_conf[i] = (double *)malloc(param->nval*sizeof(double));
    param->en_conf[i] = (double *)malloc(param->nval*sizeof(double));
    ensc[i] = (char **)malloc(param->nval*sizeof(char *));
    etemp2[i] = (char **)malloc(param->nval*sizeof(char *));
    for (j=0; j<param->nval; j++){
      ensc[i][j] = (char *)malloc(20*sizeof(char));
      etemp2[i][j] = (char *)malloc(20*sizeof(char));
      k=flexi_request_key("Configs",1,"%d %lg %s", &param->nlm_conf[i][j],
			&param->wnl_conf[i][j], ensc[i][j]);
    }
  }

  if (flexi_gather_keys(fp)) return 1;

  /* Prepare for 4th flexi_gather_keys() */
  flexi_clear_keys();
  rewind(fp);

  /* [Pseudo] */
  flexi_request_key("Pseudo",1,"%d ", &param->nval);

  param->rc = (double *)malloc(param->nval*sizeof(double));
  param->qc = (double *)malloc(param->nval*sizeof(double));
  param->nb = (int *)malloc(param->nval*sizeof(int));

  for (i=0; i<param->nval; i++) {
    flexi_request_key("Pseudo",1,"%lg",&param->rc[i]);
  }
  
  /* meth is a character for method type (k or o) */
  flexi_request_key("Pseudo",1,"%c",&param->psmeth);

  for (i=0; i<param->nval; i++) {
    param->qc[i]=0.0;
    param->nb[i]=0;
  }

  if (flexi_gather_keys(fp)) return 1;

  /* Prepare for 5th flexi_gather_keys() */
  flexi_clear_keys();
  rewind(fp);

  /* [HFsmooth] */  
  param->qpopt=0;
  param->qptol=1e-6;
  param->rlocalr = (double *)malloc(param->norb*sizeof(double));  
  flexi_request_key("HFsmooth",0,"%d %lg", &param->qpopt,&param->qptol);
  for (i=0; i<param->nval; i++) {
    param->rlocalr[i]=0.0;
    flexi_request_key("HFsmooth",0,"%lg", &param->rlocalr[i]);
  }
  /* [KBdesign] */
  flexi_request_key("KBdesign",0,"%c %d", &loc, &param->nboxes);
  param->box_start = (int *)malloc(param->nboxes*sizeof(int));
  param->box_end = (int *)malloc(param->nboxes*sizeof(int));
  param->box_height = (double *)malloc(param->nboxes*sizeof(double));
  b_unit = (char **)malloc(param->nboxes*sizeof(char *));
  b_start = (double *)malloc(param->nboxes*sizeof(double));
  b_end = (double *)malloc(param->nboxes*sizeof(double));
  for (i=0; i<param->nboxes; i++){
    b_unit[i] = (char *)malloc(3*sizeof(char));
    flexi_request_key("KBdesign",0,"%s %lg %lg %lg", 
		      b_unit[i], &b_start[i], &b_end[i], &param->box_height[i]);
  }

  if (flexi_gather_keys(fp)) return 1;

  /* Prepare for 6th flexi_gather_keys() */
  flexi_clear_keys();
  rewind(fp);

  /*[Optinfo]*/
  /* Do we need to get qc and nb? (are we optimized)? */
  if (param->psmeth=='o') {

    for (i=0; i<param->nval; i++) {
      flexi_request_key("Optinfo",1,"%lg %d",&param->qc[i],&param->nb[i]);
    }
    flexi_request_key("Optinfo",0,"%c",&param->optmeth);
    param->optmeth='n';
  }

  if (flexi_gather_keys(fp)) return 1;


  /* Prepare for 7th flexi_gather_keys() */
  flexi_clear_keys();
  rewind(fp);

  /* [Grid] */
  flexi_request_key("Grid",0,"%d %lg %lg", &param->ngrid, 
		    &param->a, &param->b);
  /* default */
  param->a = .0001;
  param->b = .013;
  param->ngrid = 1201;

  /* [LinearGrid] */
  flexi_request_key("QSOMesh",0,"%d %lg", &param->ngridl, 
		    &param->lspc);
  
  param->ngridl=1000;
  param->lspc=0.02;

  /* [RelGrid] */
  flexi_request_key("Relgrid",0,"%d %lg %lg", &param->ngrid2, 
		    &param->a2, &param->b2);
  /* default */
  param->b2=1.0/70.0;
  param->a2=exp(-4.0*2.3025851);
  param->ngrid2 = 2201;

  if (flexi_gather_keys(fp)) return 1;

  /* Prepare for 8th flexi_gather_keys() */
  flexi_clear_keys();
  rewind(fp);

  /* [Relativity] */
  param->reltype = (char *) malloc(80*sizeof(char));
  param->relxc = (char *) malloc(80*sizeof(char));
  flexi_request_key("Relativity",0,"%s %s",param->reltype,param->relxc);
  /* default */
  strcpy(param->reltype, "nrl");

  if (flexi_gather_keys(fp)) return 1;

  /* Prepare for 9th flexi_gather_keys() */
  flexi_clear_keys();
  rewind(fp);

  /* [Relativity] */
  flexi_request_key("Average",0,"%c",&avgtype);
  /* default */
  avgtype='b';
  if (streq(param->reltype,"srl")) {
    if ((param->ixc >= 2)&&(param->ixc < 6)) {   /* if GGA */
      avgtype='m';
    }
  }
  if (flexi_gather_keys(fp)) return 1;

  /*  printf("avgtype  %c \n", avgtype);*/
  if (streq(param->reltype,"srl")) {
    if ((avgtype == 'b') || (avgtype == 'B')) {
      iavgtype_.iavgtype=0; 
      if ((param->ixc >= 2)&&(param->ixc < 6)) {   /* if GGA */
	printf("NOTE: averaging \"byrc\" is not correct when using a GGA XC functional ; use caution! \n");
	fprintf(fp_log,"NOTE: averaging \"byrc\" is not correct when using a GGA XC functional ; use caution! \n");
      }
    }else if ((avgtype=='m') || (avgtype=='M')) {
      iavgtype_.iavgtype=1; 
    }else {
      printf("WARNING: Can not determine rel. average type ; must be either \"byrc\" or \"minrc\"  \n");
      fprintf(fp_log,"WARNING: Can not determine rel. average type ; must be either \"byrc\" or \"minrc\" \n");
      exit(1);
    }
  }

  /* Reorder ATOM, PSEUDO, and CONFIGS blocks to be in s p d f s p d f order */

  ncore = param->norb - param->nval;  
  param->lpot = (int *)malloc(80*sizeof(int));
  param->npot = (int *)malloc(80*sizeof(int));


  /* get the npot value, = 0 for lowest n for each l, =1 for next etc. */
  
  for (i=0 ; i<param->nval ; i++) {
    ii=i+ncore;
    param->lpot[i]=nlm_label(param->nlm[ii]).l;
    lc[i]=0;
  }
  for (i=0;i<param->nval;i++) { 
    for (j=0;j<4;j++) {
      if (param->lpot[i]==j) {
	param->npot[i]=lc[j];
	lc[j]++;
      }
    }
  }
  
  for (i=0 ; i<param->nval ; i++) {
    ii=i+ncore;
    mm = ii;
    for (j=i+1 ; j<param->nval ; j++) {
      
      jj=j+ncore;
      
      if (param->npot[j] < param->npot[mm-ncore]) {
	mm=jj;
      }else{
	if (param->npot[j] == param->npot[mm-ncore]) {
	  if ((nlm_label(param->nlm[jj]).l < nlm_label(param->nlm[mm]).l)) {
	    mm = jj;
	  }
	}
      }
    }

    ntemp=param->nlm[ii];
    wtemp=param->wnl[ii];
    
    nt1=param->npot[i];
    nt2=param->lpot[i];
    
    strcpy(etemp[i],ens[ii]);
    for (k=0 ; k<param->nconfigs ; k++) {
      ntemp2[k]=param->nlm_conf[k][i];
      wtemp2[k]=param->wnl_conf[k][i];
      strcpy(etemp2[k][i],ensc[k][i]);
    }

    rtemp=param->rc[i];    
    if (param->psmeth=='o') {
      qtemp=param->qc[i];
      nbtemp=param->nb[i];
    }
    
    param->nlm[ii]=param->nlm[mm];
    param->wnl[ii]=param->wnl[mm];
    
    param->npot[i]=param->npot[mm-ncore];
    param->lpot[i]=param->lpot[mm-ncore];
    
    strcpy(ens[ii],ens[mm]);
    for (k=0 ; k<param->nconfigs ; k++) {
      param->nlm_conf[k][i]=param->nlm_conf[k][mm-ncore];
      param->wnl_conf[k][i]=param->wnl_conf[k][mm-ncore];
      strcpy(ensc[k][i],ensc[k][mm-ncore]);
    }
    param->rc[i]=param->rc[mm-ncore];
    if (param->psmeth=='o') {
      param->qc[i]=param->qc[mm-ncore];
      param->nb[i]=param->nb[mm-ncore];
    }
    
    param->nlm[mm]=ntemp;
    param->wnl[mm]=wtemp;
    param->npot[mm-ncore]=nt1;
    param->lpot[mm-ncore]=nt2;
    
    strcpy(ens[mm],etemp[i]);
    for (k=0 ; k<param->nconfigs ; k++) {
      param->nlm_conf[k][mm-ncore]=ntemp2[k];
      param->wnl_conf[k][mm-ncore]=wtemp2[k];
      strcpy(ensc[k][mm-ncore],etemp2[k][i]);   
    }
    
    param->rc[mm-ncore]=rtemp;
    if (param->psmeth=='o') {
      param->qc[mm-ncore]=qtemp;
      param->nb[mm-ncore]=nbtemp;
    }
  }


  /*for (i=0 ; i<param->nval ; i++) {
    ii=i+ncore;
    printf("i nlm %d %d %lg %s %lg %lg %d \n", ii,param->nlm[ii],param->wnl[ii],ens[ii],param->rc[i],param->qc[i],param->nb[i]);
  }
  
    for (k=0 ; k<param->nconfigs ; k++) 
    for (i=0 ; i<param->nval ; i++) {  
      printf("k i ensc %d %d %s \n ",k,i,ensc[k][i]);
    }
    exit(1);*/
  

  nlpot2_.inl = 0;          /* set AE mode flag */
  for (i=0; i<param->nval; i++) {
    aval_.rcall[i] = param->rc[i];
    /*    optparam_.qcl[i] = param->qc[i];
	  optparam_.nbl[i] = param->nb[i];*/
  }
  
  /* ADD CHECK for semicore violations by user */

  /* Cleanup */
  flexi_clear_keys();
  rewind(fp);

  /*  param->ipot = (int *)malloc(80*sizeof(int));*/

  /* DONE WITH INPUT VALUES */

  /* Post processing section */

  /* set up semicore info */

  /*  
  param->ipot = (int *)malloc(80*sizeof(int));

  for (i=0; i<10; i++){
    npot[i]=0;
  }
  for (i=0; i<20; i++){
    param->ipot[i]=0;
  }
  param->nll=0;

  ncore = param->norb - param->nval;
  for (i=ncore; i<param->norb; i++){
    npot[nlm_label(param->nlm[i]).l]++;
    if (npot[nlm_label(param->nlm[i]).l] == 1) {
      param->ipot[param->nll]=i-ncore;
      aval_.rcall[param->nll] = param->rc[i-ncore];
      optparam_.qcl[param->nll] = param->qc[i-ncore];
      optparam_.nbl[param->nll] = param->nb[i-ncore];
      param->nll++;
    }
  }

  aorb_.nval = param->nval;
  aorb_.ncore = param->norb - param->nval;
  psdat_.nll = param->nll;
  */
  /* Section to make KB local idiot proof */

  param->local=-67;
  switch (loc) {
    
  case '0' : 
  case 's' : 
  case 'S' : loctemp=0 ; break; 

  case '1' : 
  case 'p' : 
  case 'P' : loctemp=1 ; break; 

  case '2' : 
  case 'd' : 
  case 'D' : loctemp=2 ; break; 

  case '3' : 
  case 'f' : 
  case 'F' : loctemp=3 ; break; 
    
  default : 
    fprintf(stderr," Can't determine local potential, check [KBdesign] key-block \n"); 
    exit(1); 
  }
  ncore=param->norb-param->nval;
  for (i=ncore;i<param->norb;i++) {
    if (loctemp == nlm_label(param->nlm[i]).l) {
      
      /*EJW  param->local is the l value of the local pot
	param->localind is the valence index (index-ncore) of the local pot */
      
      param->local=loctemp;
      param->localind=i-ncore;
      break;
    } 
  }

  if (param->local == -67) {
    fprintf(stderr," Can't determine local potential, check [KBdesign] key-block \n"); 
    exit(1); 
  }

  /* This is where initial guesses are chosen */
  /* etrial is used as the default eigenvalue if an unbound state is found */

  /* process reference config eigenvalue guesses and free arrays */
  for (i=0; i<param->norb; i++){
    
    if (i >= ncore) param->ibound[i-ncore]=1;
    k = 0;
    j = 0;
    
    while ( k<strlen(ens[i]) && ens[i][k]=='-' ) ++k;
    if (k==strlen(ens[i])){
      /* user wants opium to guess a good eigenvalue */
      n = param->nlm[i]/100;
      param->etrial[i] = -1. * (param->z * param->z) / (double)(n*n*n*n);
      param->en[i]=param->etrial[i];
      if ((i >=ncore)&&(param->wnl[i]< -1e-12)) {
	printf(" You must supply an eigenvalue for an unbound state! \n Problem with state: %d \n",param->nlm[i]);
	fprintf(fp_log," You must supply an eigenvalue for an unbound state! \n Problem with state: %d \n",param->nlm[i]);
	exit(1);
      }
    }else {
      /* user provides an eigenvalue guess */
      param->etrial[i]= atof(ens[i]);
      param->en[i]=param->etrial[i];
    }

    if (param->wnl[i] < -1e-12) {
      if (i>ncore) {
	/*	param->ibound[i-ncore]=0;
	  param->wnl[i]=0.0;*/
      } else {
	printf(" !ERROR! A core state can not be unbound! \n Problem with state: %d \n",param->nlm[i]);
	fprintf(fp_log, " !ERROR! A core state can not be unbound! \n Problem with state: %d \n",param->nlm[i]);
	exit(1);
	}
    }
    
    free(ens[i]);
  }
  free(ens);
  
  /* process test config eigenvalue guesses and free arrays */
  for (i=0; i<param->nconfigs; i++){
    for (j=0; j<param->nval; j++){
      k = 0;

      while ( k<strlen(ensc[i][j]) && ensc[i][j][k]=='-' ) ++k;
      if (k==strlen(ensc[i][j])){
        /* user wants opium to guess a good eigenvalue */
        n = param->nlm_conf[i][j]/100;
        param->en_conf[i][j] = -1. * (param->z * param->z) / (double)(n*n*n*n);

      }else{
        /* user provides an eigenvalue guess */
        param->en_conf[i][j] = atof(ensc[i][j]);
      }

      free(ensc[i][j]);
    }
    free(ensc[i]);
  }
  free(ensc);

  /* Grid set up */
  
  grid_.r1=param->a;
  grid_.h=param->b;
  grid_.z=param->z;
  
  grid_.np=0;

  /* if (streq(param->reltype,"nrl")) { */
    grid_.r[0] = param->a * pow(param->z, -1./3.);
    /*  } else {
    grid_.r[0] = param->a / param->z;
    }*/
  for (i=1; i<NPDM; i++) {
    grid_.r[i]=grid_.r[0] * exp(param->b * (i));
    if (grid_.r[i]>120.0) break;
    grid_.np++;
  }
  if (grid_.np > param->ngrid) {
    printf("check grid; # points needed %d > max np %d \n",grid_.np,param->ngrid);
    return 1;
  }
  if (grid_.r[grid_.np]<30.0) {
    printf("check grid; grid extends to less than 30 Bohr \n");
    return 1;
  }

  param->ngrid=grid_.np;
  

  /* non-rel grid ok, how about rel-grid */
  
  rgrid_.b=param->b2;
  rgrid_.a=param->a2 / (param->b2*param->z);
  
  nrgrid_.nr=0;
  for (i=0; i<NPDM; i++){
    rgrid_.r[i]=rgrid_.a * (exp(rgrid_.b*(i))-1.0);
    rgrid_.rab[i] =(rgrid_.r[i]+rgrid_.a)*rgrid_.b;
    if (rgrid_.r[i]>80.0) break;
    nrgrid_.nr++;
  }

  /* process the box limits in the KBdesign key-block and free temp arrays */
  for (i=0; i<param->nboxes; i++){
    if (streq(b_unit[i], "au") || streq(b_unit[i], "a.u.") 
	|| streq(b_unit[i], "AU") || streq(b_unit[i], "A.U.")){
      
      /* need to convert box limits to grid units */
      r_l = param->a * pow(param->z, -1./3.);
      r_r = r_l * exp(param->b * (param->ngrid-1));
      
      /* outside grid test */
      if(b_start[i]<r_l)
        param->box_start[i] = 1;            /* Fortran array index conv. */
      else if(b_start[i]>r_r)
        param->box_start[i] = param->ngrid; /* Fortran array index conv. */
      if(b_end[i]<r_l)
        param->box_end[i] = 1;              /* Fortran array index conv. */
      else if(b_end[i]>r_r)
        param->box_end[i] = param->ngrid;   /* Fortran array index conv. */

      /* inside grid test */
      for (j=1; j<param->ngrid; j++){
        r_r = r_l * exp(param->b);
        /* check start limit */
        if ((r_l<b_start[i]) && (b_start[i]<r_r)){
          if (fabs(b_start[i]-r_l) < fabs(b_start[i]-r_r))
            param->box_start[i] = j;    /* notice Fortran array index conv. */
          else
            param->box_start[i] = j+1;  /* notice Fortran array index conv. */ 
        }
        /* check end limit */
        if ((r_l<b_end[i]) && (b_end[i]<r_r)){
          if (fabs(b_end[i]-r_l) < fabs(b_end[i]-r_r))
            param->box_end[i] = j;      /* notice Fortran array index conv. */
          else
            param->box_end[i] = j+1;    /* notice Fortran array index conv. */
        }
        r_l = r_r;
      }
    }else{
      /* limits are assumed to be in grid units */
      param->box_start[i] = (int)b_start[i];
      param->box_end[i] = (int)b_end[i];
    }
    free(b_unit[i]);
  }
  free(b_unit);
  free(b_start);
  free(b_end);  


  /* enforce GGA damping from 0.0 -> 0.001 au */
  if ((param->ixc >= 2)&&(param->ixc < 6)) {   /* if GGA */
    if (fabs(param->exccut + 600) < 1e-6 ) {   /* if exccut not set */
      param->exccut=1e-3;     
    }
  }else{
    param->exccut=0.0;
  }

  /*  printf(" exccut= %lg \n",param->exccut);*/

  /* set up DNL stuff */
   
  local_.iloc = param->localind+1;
  box_.numbox = param->nboxes;

  /*  for (i=0 ; i<4; i++) {*/
  for (i=0; i<N0; i++) { /* BST 09/08/10: N0 (defined in cdim.h) replaces hard coded value of 4 here. */
    if (i < param->nboxes){
      box_.iboxstart[i] = param->box_start[i];
      box_.iboxend[i]   = param->box_end[i];
      box_.boxheight[i] = param->box_height[i];
    }else{
      box_.iboxstart[i] = 1;
      box_.iboxend[i]   = 2;
      box_.boxheight[i] = 0.;
    }
  }

  /* malloc the strings for psp output info, the execstring, compile info, and execution info: */
  param->execstring = (char *) malloc(160*sizeof(char));
  param->version = (char *) malloc(160*sizeof(char));
  param->chost = (char *) malloc(160*sizeof(char));
  param->cdate = (char *) malloc(160*sizeof(char));
  param->csys = (char *) malloc(160*sizeof(char));
  param->ehost = (char *) malloc(160*sizeof(char));
  param->edate = (char *) malloc(160*sizeof(char));


  
  if (streq(param->reltype,"nrl")) {
    if (streq(param->relxc,"rxc")) { 
      strcpy(param->relxc, "nxc");
      fprintf(fp_log, " WARNING!! Can't use rel corrections for XC in nrl calculation ; rel corrections removed \n");
      printf(" WARNING!! Can't use rel corrections for XC in nrl calculation ; rel corrections removed \n");
    }else if (!streq(param->relxc,"nxc")) {
      strcpy(param->relxc, "nxc");
    }
  }else {
    if (  (!streq(param->relxc,"nxc"))&&(!streq(param->relxc,"rxc"))) {
      strcpy(param->relxc, "rxc");
    }
  }

/* Finally! Write the log file header */
  
  fprintf(fp_log, " File prefix             : %s \n",param->name);
  fprintf(fp_log, " Element                 : %s(%s) \n",param->longname,param->symbol);
  fprintf(fp_log, " Z                       : %lg \n",param->z);
  fprintf(fp_log, " \n");
  fprintf(fp_log, " Number of all-electron orbitals   : %d \n",param->norb);
  fprintf(fp_log, " Number of pseudo       orbitals   : %d \n",param->nval);
  fprintf(fp_log, " \n");
  
  if (streq(param->reltype,"nrl")) {
    fprintf(fp_log, " Pseudopotential is non-relativistic \n");
  }else if (streq(param->reltype,"srl")){
    fprintf(fp_log, " Pseudopotential is scalar-relativistic \n");         
    if ( (streq(param->relxc,"rxc")) && (!streq(param->reltype,"nrl")) ) {
      fprintf(fp_log, " Relativistic corrections applied to XC potential \n");         
    }else{
      fprintf(fp_log, " Relativistic corrections NOT applied to XC potential \n");         
    }
    if (iavgtype_.iavgtype == 0) {
      fprintf(fp_log, " Using individual cut-off radii for srl averaging \n");         
    }else if (iavgtype_.iavgtype == 1) {
      fprintf(fp_log, " Using the minimum cut-off radius for srl averaging \n");         
    }
  }else if (streq(param->reltype,"frl")){
    fprintf(fp_log, " Fully relativistic pseudopotential \n");
    if ( (streq(param->relxc,"rxc")) && (!streq(param->reltype,"nrl")) ) {
      fprintf(fp_log, " Relativistic corrections applied to XC potential \n");         
    }else{
      fprintf(fp_log, " Relativistic corrections NOT applied to XC potential \n");         
    }
  }else{
    fprintf(fp_log, " !!ERROR!! Can not determine relativity \n");
    printf(" !!ERROR!! Can not determine relativity \n");
    exit(1);
  }
    
  if (streq(param->xcparam,"lda")|| streq(param->xcparam,"pzlda")) {
    fprintf(fp_log, " Exchage-correlation functional is : Perdew-Zunger LDA \n");
  }else if (streq(param->xcparam,"pwlda")){
    fprintf(fp_log, " Exchage-correlation functional is : Perdew-Wang LDA \n");
  }else if (streq(param->xcparam,"vwn5")|| streq(param->xcparam,"vwn5lda")) {
    fprintf(fp_log, " Exchage-correlation functional is : Vosko, Wilk, Nusair (V) LDA  \n");
  }else if (streq(param->xcparam,"hf")){
    fprintf(fp_log, " Using Hartree-Fock Exchange (Beta!!)                \n");
    if (param->qpopt == 0) {
      fprintf(fp_log, " HF pseudopotential NOT being smoothed! \n");
    }else if (param->qpopt != 0) {
      if (param->qpopt == 1) 
	fprintf(fp_log, " HF smoothing method: Trail & Needs for Troullier Martins  \n");
      if (param->qpopt == 2) 
	fprintf(fp_log, " HF smoothing method: AWR#2 for Optimized/RRKJ \n");
      if (param->qpopt == 3) 
	fprintf(fp_log, " HF smoothing method: AWR#3 for Optimized/RRKJ \n");
      if (param->qpopt == 4) 
	fprintf(fp_log, " HF smoothing method: AWR#4 for Optimized/RRKJ \n");
      fprintf(fp_log, " HF smoothing eigenvalue tol: %10.3e \n",param->qptol);
    }else{
      fprintf(fp_log, " !!WARNING!! Unknown HF smoothing method! \n");
    }
    
  }else if (streq(param->xcparam,"gga")|| streq(param->xcparam,"pbegga")) {
    fprintf(fp_log, " Exchage-correlation functional is : Perdew, Burke, Ernzerhof GGA \n");
  }else if (streq(param->xcparam,"pw91gga")){
    fprintf(fp_log, " Exchage-correlation functional is : Perdew Wang 91 GGA \n");
  }else if (streq(param->xcparam,"wcgga")){
    fprintf(fp_log, " Exchage-correlation functional is : Wu-Cohen GGA \n");
  }else if (streq(param->xcparam,"pbesol")){
    fprintf(fp_log, " Exchage-correlation functional is : PBEsol \n");
  }
  if (param->exccut > 0.0) {
    fprintf(fp_log, " GGA smoothed from 0 to %f bohr using switch to LDA \n",fabs(param->exccut));
    /*fprintf(fp_log, " GGA smoothed in all steps of calculation \n");*/
  } else if (param->exccut < 0.0) {
    fprintf(fp_log, " GGA smoothed from 0 to %f bohr using quadratic extrapolation\n",fabs(param->exccut));
    fprintf(fp_log, " GGA smoothed in descreening step only \n");
  }

  ncore=param->norb - param->nval;

  fprintf(fp_log, " The %c potential is used for the KB construction \n",
	  la[nlm_label(param->nlm[param->localind+ncore]).l]);

  if (param->psmeth=='o') {
    fprintf(fp_log, " Optimized (RRKJ) pseudopotential method will be used \n");
    fprintf(fp_log, "\n");
    fprintf(fp_log, " nl    cutoff radii    q-max      # bessel functions\n");
    for (i=0; i<param->nval; i++){
      fprintf(fp_log, " %d%c       %5.3f        %5.3f              %d \n", 
	      nlm_label(param->nlm[i+ncore]).n,la[nlm_label(param->nlm[i+ncore]).l],
	      param->rc[i],param->qc[i],param->nb[i]);
    }
  }else if (param->psmeth=='k') {
    fprintf(fp_log, " Kerker pseudopotential method will be used \n");
    fprintf(fp_log, "\n");
    fprintf(fp_log, " ang. mom.   cutoff radii \n");
    n=param->norb - param->nval;
    for (i=0; i<param->nval; i++){
      fprintf(fp_log, "  %d            %5.3f \n", nlm_label(param->nlm[i+ncore]).l,param->rc[i]);
    }
  }else if (param->psmeth=='t') {
    fprintf(fp_log, " Troullier-Martins pseudopotential method will be used \n");
    fprintf(fp_log, "\n");
    fprintf(fp_log, " ang. mom.   cutoff radii \n");
    n=param->norb - param->nval;
    for (i=0; i<param->nval; i++){
      fprintf(fp_log, "   %d            %5.3f \n", nlm_label(param->nlm[i+ncore]).l,param->rc[i]);
    }
  }else{
    fprintf(fp_log, " Cannot determine pseudopotential method \n");
    exit(1);
  }

  if (param->rpcc > 1e-6) {
    fprintf(fp_log, " Partial core radius: %lg \n",param->rpcc);
    if (ipccmeth_.ipccmeth == 0) {
      fprintf(fp_log, " Using the Louie-Froyen-Cohen type pcc \n");
    }else{
      fprintf(fp_log, " Using the Fuchs-Scheffler type pcc \n");
    }
  }

  fprintf(fp_log, "\n");

  fprintf(fp_log, " Reference Configuration: core  <-/-> valence \n");

  for (i=0; i<ncore; i++){
    fprintf(fp_log," %d%c%lg ", nlm_label(param->nlm[i]).n,la[nlm_label(param->nlm[i]).l],
	    param->wnl[i]);
  }
  fprintf(fp_log," <-/->");
  for (i=ncore; i<param->norb; i++){
    fprintf(fp_log," %d%c%lg ",nlm_label(param->nlm[i]).n,la[nlm_label(param->nlm[i]).l],
	    param->wnl[i]);
  }
  fprintf(fp_log, "\n");
  fprintf(fp_log, "\n");

  fprintf(fp_log," Grid definitions: \n");

  /*  if (streq(param->reltype,"nrl")){*/
  fprintf(fp_log," a_grid=%-11.2e b_grid=%-11.2e # points=%d \n r(1)=%-11.2e   r(np)=%-11.2e \n",
	  param->a,param->b,grid_.np,grid_.r[0],grid_.r[grid_.np-1]);
  /* }
  
 if (streq(param->reltype,"srl")){
    fprintf(fp_log," Non-relativitic grid:  a_grid=%-11.2e b_grid=%-11.2e # points=%d \n  r(1)=%-11.2e r(np)=%-11.2e \n",
	    param->a,param->b,grid_.np,grid_.r[0],grid_.r[grid_.np-1]);
    
    fprintf(fp_log," Relativistic grid   :  a_grid=%-11.2e b_grid=%-11.2e # points=%d \n  r(1)=%-11.2e r(np)=%-11.2e \n",
      param->a2,param->b2,nrgrid_.nr,rgrid_.r[0],rgrid_.r[nrgrid_.nr-1]);
  }*/
  
  fprintf(fp_log, "\n");
  
  fprintf(fp_log, " dEmax tolerance:  %-11.2e   dVmax tolerance:  %-11.2e  \n",consts_.etol,consts_.vtol);
  
  fprintf(fp_log, "\n");
  
  if (param->nboxes > 0) {
    for (i=0; i<param->nboxes; i++){
      fprintf(fp_log, " Box %d  %lg --> %lg  height: %lg \n",i+1,grid_.r[param->box_start[i]],grid_.r[param->box_end[i]],param->box_height[i]);
    }
  }
  
  return 0;
  
}  

double symtoz(char *sym, char *longname, double *mass) {
  
  double val=0.0;

  if (streq(sym,"H")){
    strcpy(longname,"Hydrogen   ");
    val =  1;
    *mass = 1.008;

  }else if (streq(sym,"He")){
    strcpy(longname,"Helium     ");
    val =  2;
    *mass = 4.002602;

  }else if (streq(sym,"Li")){
    strcpy(longname,"Lithium    ");
    val =  3;
    *mass = 6.94;

  }else if (streq(sym,"Be")){
    strcpy(longname,"Beryllium  ");
    val =  4;
    *mass = 9.012182;

  }else if (streq(sym,"B")){
    strcpy(longname,"Boron      ");
    val =  5;
    *mass =     10.81;

  }else if (streq(sym,"C")){
    strcpy(longname,"Carbon     ");
    val =  6;
    *mass =     12.011;

  }else if (streq(sym,"N")){
    strcpy(longname,"Nitrogen   ");
    val =  7;
    *mass =     14.007;

  }else if (streq(sym,"O")){
    strcpy(longname,"Oxygen     ");
    val =  8;
    *mass =     15.999;

  }else if (streq(sym,"F")){
    strcpy(longname,"Fluorine   ");
    val =  9;
    *mass =     18.9984032;

  }else if (streq(sym,"Ne")){
    strcpy(longname,"Neon       ");
    val = 10;
    *mass =     20.1797;

  }else if (streq(sym,"Na")){
    strcpy(longname,"Sodium     ");
    val = 11;
    *mass =     22.98976928;

  }else if (streq(sym,"Mg")){
    strcpy(longname,"Magnesium  ");
    val = 12;
    *mass =     24.3050;

  }else if (streq(sym,"Al")){
    strcpy(longname,"Aluminum   ");
    val = 13;
    *mass =     26.9815386;

  }else if (streq(sym,"Si")){
    strcpy(longname,"Silicon    ");
    val = 14;
    *mass =     28.085;

  }else if (streq(sym,"P")){
    strcpy(longname,"Phosphorus ");
    val = 15;
    *mass =     30.973762;

  }else if (streq(sym,"S")){
    strcpy(longname,"Sulfur     ");
    val = 16;
    *mass =     32.06;

  }else if (streq(sym,"Cl")){
    strcpy(longname,"Chlorine   ");
    val = 17;
    *mass =     35.45;

  }else if (streq(sym,"Ar")){
    strcpy(longname,"Argon      ");
    val = 18;
    *mass =     39.948;

  }else if (streq(sym,"K")){
    strcpy(longname,"Potassium  ");
    val = 19;
    *mass =     39.0983;

  }else if (streq(sym,"Ca")){
    strcpy(longname,"Calcium    ");
    val = 20;
    *mass =     40.078;

  }else if (streq(sym,"Sc")){
    strcpy(longname,"Scandium   ");
    val = 21;
    *mass =     44.955912;

  }else if (streq(sym,"Ti")){
    strcpy(longname,"Titanium   ");
    val = 22;
    *mass =     47.867;

  }else if (streq(sym,"V")){
    strcpy(longname,"Vanadium   ");
    val = 23;
    *mass =     50.9415;

  }else if (streq(sym,"Cr")){
    strcpy(longname,"Chromium   ");
    val = 24;
    *mass =     51.9961;

  }else if (streq(sym,"Mn")){
    strcpy(longname,"Manganese  ");
    val = 25;
    *mass =     54.938045;

  }else if (streq(sym,"Fe")){
    strcpy(longname,"Iron       ");
    val = 26;
    *mass =     55.845;

  }else if (streq(sym,"Co")){
    strcpy(longname,"Cobalt     ");
    val = 27;
    *mass =     58.933195;

  }else if (streq(sym,"Ni")){
    strcpy(longname,"Nickel     ");
    val = 28;
    *mass =     58.6934;

  }else if (streq(sym,"Cu")){
    strcpy(longname,"Copper     ");
    val = 29;
    *mass =     63.546;

  }else if (streq(sym,"Zn")){
    strcpy(longname,"Zinc       ");
    val = 30;
    *mass =     65.38;

  }else if (streq(sym,"Ga")){
    strcpy(longname,"Gallium    ");
    val = 31;
    *mass =     69.723;

  }else if (streq(sym,"Ge")){
    strcpy(longname,"Germanium  ");
    val = 32;
    *mass =     72.63;

  }else if (streq(sym,"As")){
    strcpy(longname,"Arsenic    ");
    val = 33;
    *mass =     74.92160;

  }else if (streq(sym,"Se")){
    strcpy(longname,"Selenium   ");
    val = 34;
    *mass =     78.96;

  }else if (streq(sym,"Br")){
    strcpy(longname,"Bromine    ");
    val = 35;
    *mass =     79.904;

  }else if (streq(sym,"Kr")){
    strcpy(longname,"Krypton    ");
    val = 36;
    *mass =     83.798;

  }else if (streq(sym,"Rb")){
    strcpy(longname,"Rubidium   ");
    val = 37;
    *mass =     85.4678;

  }else if (streq(sym,"Sr")){
    strcpy(longname,"Strontium  ");
    val = 38;
    *mass =     87.62;

  }else if (streq(sym,"Y")){
    strcpy(longname,"Yttrium    ");
    val = 39;
    *mass =     88.90585;

  }else if (streq(sym,"Zr")){
    strcpy(longname,"Zirconium  ");
    val = 40;
    *mass =     91.224;

  }else if (streq(sym,"Nb")){
    strcpy(longname,"Niobium    ");
    val = 41;
    *mass =     92.90638;

  }else if (streq(sym,"Mo")){
    strcpy(longname,"Molybdenum ");
    val = 42;
    *mass =     95.96;

  }else if (streq(sym,"Tc")){
    strcpy(longname,"Technetium ");
    val = 43;
    *mass =98;

  }else if (streq(sym,"Ru")){
    strcpy(longname,"Ruthenium  ");
    val = 44;
    *mass =     101.07;

  }else if (streq(sym,"Rh")){
    strcpy(longname,"Rhodium    ");
    val = 45;
    *mass =     102.90550;

  }else if (streq(sym,"Pd")){
    strcpy(longname,"Palladium  ");
    val = 46;
    *mass =     106.42;

  }else if (streq(sym,"Ag")){
    strcpy(longname,"Silver     ");
    val = 47;
    *mass =     107.8682;

  }else if (streq(sym,"Cd")){
    strcpy(longname,"Cadmium    ");
    val = 48;
    *mass =     112.411;

  }else if (streq(sym,"In")){
    strcpy(longname,"Indium     ");
    val = 49;
    *mass =     114.818;

  }else if (streq(sym,"Sn")){
    strcpy(longname,"Tin        ");
    val = 50;
    *mass =     118.710;

  }else if (streq(sym,"Sb")){
    strcpy(longname,"Antimony   ");
    val = 51;
    *mass =     121.760;

  }else if (streq(sym,"Te")){
    strcpy(longname,"Tellurium  ");
    val = 52;
    *mass =     127.60;

  }else if (streq(sym,"I")){
    strcpy(longname,"Iodine     ");
    val = 53;
    *mass =     126.90447;

  }else if (streq(sym,"Xe")){
    strcpy(longname,"Xenon      ");
    val = 54;
    *mass =     131.293;

  }else if (streq(sym,"Cs")){
    strcpy(longname,"Cesium     ");
    val = 55;
    *mass =     132.9054519;

  }else if (streq(sym,"Ba")){
    strcpy(longname,"Barium     ");
    val = 56;
    *mass =     137.327;

  }else if (streq(sym,"La")){
    strcpy(longname,"Lanthanum  ");
    val = 57;
    *mass =     138.90547;

  }else if (streq(sym,"Ce")){
    strcpy(longname,"Cerium     ");
    val = 58;
    *mass =     140.116;

  }else if (streq(sym,"Pr")){
    strcpy(longname,"Praseodymium");
    val = 59;
    *mass =     140.90765;

  }else if (streq(sym,"Nd")){
    strcpy(longname,"Neodymium  ");
    val = 60;
    *mass =     144.242;

  }else if (streq(sym,"Pm")){
    strcpy(longname,"Promethium ");
    val = 61;
    *mass =    145;

  }else if (streq(sym,"Sm")){
    strcpy(longname,"Samarium   ");
    val = 62;
    *mass =     150.36;

  }else if (streq(sym,"Eu")){
    strcpy(longname,"Europium   ");
    val = 63;
    *mass =     151.964;

  }else if (streq(sym,"Gd")){
    strcpy(longname,"Gadolinium ");
    val = 64;
    *mass =     157.25;

  }else if (streq(sym,"Tb")){
    strcpy(longname,"Terbium    ");
    val = 65;
    *mass =     158.92535;

  }else if (streq(sym,"Dy")){
    strcpy(longname,"Dysprosium ");
    val = 66;
    *mass =     162.500;

  }else if (streq(sym,"Ho")){
    strcpy(longname,"Holmium    ");
    val = 67;
    *mass =     164.93032;

  }else if (streq(sym,"Er")){
    strcpy(longname,"Erbium     ");
    val = 68;
    *mass =     167.259;

  }else if (streq(sym,"Tm")){
    strcpy(longname,"Thulium    ");
    val = 69;
    *mass =     168.93421;

  }else if (streq(sym,"Yb")){
    strcpy(longname,"Ytterbium  ");
    val = 70;
    *mass =     173.054;

  }else if (streq(sym,"Lu")){
    strcpy(longname,"Lutetium   ");
    val = 71;
    *mass =     174.9668;

  }else if (streq(sym,"Hf")){
    strcpy(longname,"Hafnium    ");
    val = 72;
    *mass =     178.49;

  }else if (streq(sym,"Ta")){
    strcpy(longname,"Tantalum   ");
    val = 73;
    *mass =     180.94788;

  }else if (streq(sym,"W")){
    strcpy(longname,"Tungsten   ");
    val = 74;
    *mass =     183.84;

  }else if (streq(sym,"Re")){
    strcpy(longname,"Rhenium    ");
    val = 75;
    *mass =     186.207;

  }else if (streq(sym,"Os")){
    strcpy(longname,"Osmium     ");
    val = 76;
    *mass =     190.23;

  }else if (streq(sym,"Ir")){
    strcpy(longname,"Iridium    ");
    val = 77;
    *mass =     192.217;

  }else if (streq(sym,"Pt")){
    strcpy(longname,"Platinum   ");
    val = 78;
    *mass =     195.084;

  }else if (streq(sym,"Au")){
    strcpy(longname,"Gold       ");
    val = 79;
    *mass =     196.966569;

  }else if (streq(sym,"Hg")){
    strcpy(longname,"Mercury    ");
    val = 80;
    *mass =     200.59;

  }else if (streq(sym,"Tl")){
    strcpy(longname,"Thallium   ");
    val = 81;
    *mass =     204.38;

  }else if (streq(sym,"Pb")){
    strcpy(longname,"Lead       ");
    val = 82;
    *mass =     207.2;

  }else if (streq(sym,"Bi")){
    strcpy(longname,"Bismuth    ");
    val = 83;
    *mass =     208.98040;

  }else if (streq(sym,"Po")){
    strcpy(longname,"Polonium   ");
    val = 84;
    *mass =     209;

  }else if (streq(sym,"At")){
    strcpy(longname,"Astatine   ");
    val = 85;
    *mass =     210;

  }else if (streq(sym,"Rn")){
    strcpy(longname,"Radon      ");
    val = 86;
    *mass =     222;

  }else if (streq(sym,"Fr")){
    strcpy(longname,"Francium   ");
    val = 87;
    *mass =     223;

  }else if (streq(sym,"Ra")){
    strcpy(longname,"Radium     ");
    val = 88;
    *mass =     226;

  }else if (streq(sym,"Ac")){
    strcpy(longname,"Actinium   ");
    val = 89;
    *mass =     227;

  }else if (streq(sym,"Th")){
    strcpy(longname,"Thorium    ");
    val = 90;
    *mass =     232.03806;

  }else if (streq(sym,"Pa")){
    strcpy(longname,"Protactinium");
    val = 91;
    *mass =     231.03588;

  }else if (streq(sym,"U")){
    strcpy(longname,"Uranium    ");
    val = 92;
    *mass =     238.02891;

  }else if (streq(sym,"Np")){
    strcpy(longname,"Neptunium  ");
    val = 93;
    *mass =     237;

  }else if (streq(sym,"Pu")){
    strcpy(longname,"Plutonium  ");
    val = 94;
    *mass =     244;

  }else if (streq(sym,"Am")){
    strcpy(longname,"Americium  ");
    val = 95;
    *mass =     243;

  }else if (streq(sym,"Cm")){
    strcpy(longname,"Curium     ");
    val = 96;
    *mass =     247;

  }else if (streq(sym,"Bk")){
    strcpy(longname,"Berkelium  ");
    val = 97;
    *mass =     247;

  }else if (streq(sym,"Cf")){
    strcpy(longname,"Californium");
    val = 98;
    *mass =     251;

  }else if (streq(sym,"Es")){
    strcpy(longname,"Einsteinium");
    val = 99;
    *mass =     252;

  }else if (streq(sym,"Fm")){
    strcpy(longname,"Fermium    ");
    val =100;
    *mass =     257;

  }else if (streq(sym,"Md")){
    strcpy(longname,"Mendelevium");
    val =101;
    *mass =     258;

  }else if (streq(sym,"No")){
    strcpy(longname,"Nobelium   ");
    val =102;
    *mass =     259;

  }else if (streq(sym,"Lr")){
    strcpy(longname,"Lawrencium ");
    val =103;
    *mass =     262;

  }else if (streq(sym,"Rf")){
    strcpy(longname,"Rutherfordium");
    val =104;
    *mass =     265;

  }else if (streq(sym,"Db")){
    strcpy(longname,"Dubnium    ");
    val =105;
    *mass =     268;

  }else if (streq(sym,"Sg")){
    strcpy(longname,"Seaborgium ");
    val =106;
    *mass =     271;

  }else if (streq(sym,"Bh")){
    strcpy(longname,"Bohrium    ");
    val =107;
    *mass =     270;

  }else if (streq(sym,"Hs")){
    strcpy(longname,"Hassium    ");
    val =108;
    *mass =     277;

  }else if (streq(sym,"Mt")){
    strcpy(longname,"Meitnerium ");
    val =109;
    *mass = 276;
  }
  return val;
}
  
