/* 
 * $Id: parameter.c,v 1.8 2004/06/16 21:25:54 mbarnes Exp $
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
#include "fortparam.h"
#include "nlm.h"
#include "common_blocks.h"

/* module's own header */
#include "parameter.h"

#define streq(a,b) (*a==*b && !strcmp(a+1,b+1))

double symtoz(char *sym , char *longname);
int read_param(param_t *param, FILE *fp, FILE *fp_log){

  int i, j, k;
  int n;        
  char **ens;                /* temp string array for eigenvalue guesses */
  char ***ensc; /* temp string array for eigenvalue guesses of test configs */
  double *b_start, *b_end;  /* temp arrays for box start and end values */
  char **b_unit;            /* temp string array to read the box unit */
  double r_l, r_r;          /* temp radius variables [a.u.] */
  char la[]={"spdf"};
  int ncore; 
  char loc='s';
  int loctemp=0;
  int npot[10];

  /* Prepare for 1st flexi_gather_keys() */  

  /* [Atom] */
  param->symbol = (char *) malloc(160*sizeof(char));
  param->longname = (char *) malloc(160*sizeof(char));
  flexi_request_key("Atom",1,"%s %d", param->symbol, &param->norb);

  /* [Relativity] */
  param->reltype = (char *) malloc(80*sizeof(char));
  flexi_request_key("Relativity",0,"%s",param->reltype);
  /* default */
  strcpy(param->reltype, "nrl");

  /* [Grid] */
  flexi_request_key("Grid",0,"%d %lg %lg", &param->ngrid, 
		    &param->a, &param->b);
  /* default */
  param->a = .0001;
  param->b = .013;
  param->ngrid = 1201;
  
  /* [RelGrid] */
  flexi_request_key("Relgrid",0,"%d %lg %lg", &param->ngrid2, 
		    &param->a2, &param->b2);
  /* default */
  param->b2=1.0/70.0;
  param->a2=exp(-4.0*2.3025851);
  param->ngrid2 = 2201;

  /* [Configs] */
  flexi_request_key("Configs",0,"%d", &param->nconfigs);
  /* default */
  param->nconfigs = 0;

  /* [EVTol] */
  flexi_request_key("Tol",0,"%lg %lg", &consts_.etol,&consts_.vtol);
  /* default */
  consts_.etol = 1.e-8;   
  consts_.vtol = 1.e-6;   
  consts_.maxit = 800;      

  /* [Pseudo] */
  flexi_request_key("Pseudo",1,"%d", &param->nval);

  /* [Pcc] */
  flexi_request_key("Pcc",0,"%lg", &param->rpcc);
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
  flexi_request_key("XC",0,"%lg",&param->rxccut);
  /* default */
  strcpy(param->xcparam, "lda");
  param->rxccut=0.01;
  
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

  param->z = symtoz(param->symbol,param->longname);
  rxccut_.rxccut=param->rxccut;

  param->ixc = 0;
  if (streq(param->xcparam, "gga")) param->ixc = 2;
  if (streq(param->xcparam, "pwlda")) param->ixc = 1;
  if (streq(param->xcparam, "pzlda")) param->ixc = 0;
  if (streq(param->xcparam, "lda")) param->ixc = 0;

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

  /* [Atom] */
  flexi_request_key("Atom",1,"%s %d", param->symbol, &param->norb);
  param->nlm = (int *)malloc(param->norb*sizeof(int));
  param->wnl = (double *)malloc(param->norb*sizeof(double));
  param->en = (double *)malloc(param->norb*sizeof(double));
  ens = (char **)malloc(param->norb*sizeof(char*));
  for (i=0; i<param->norb; i++){
    ens[i] = (char *)malloc(20*sizeof(char));
    flexi_request_key("Atom",1,"%d %lg %s", &param->nlm[i],
		      &param->wnl[i], ens[i]);
  }
  
  /* [Configs] */
  flexi_request_key("Configs",0,"%d", &param->nconfigs);
  param->nlm_conf = (int **)malloc(param->nconfigs*sizeof(int *));
  param->wnl_conf = (double **)malloc(param->nconfigs*sizeof(double *));
  param->en_conf = (double **)malloc(param->nconfigs*sizeof(double *)); 
  param->ensave_conf = (double **)malloc(param->nconfigs*sizeof(double *)); 
  ensc = (char ***)malloc(param->nconfigs*sizeof(char **));
  for (i=0; i<param->nconfigs; i++){
    param->nlm_conf[i] = (int *)malloc(param->nval*sizeof(int));
    param->wnl_conf[i] = (double *)malloc(param->nval*sizeof(double));
    param->en_conf[i] = (double *)malloc(param->nval*sizeof(double));
    param->ensave_conf[i] = (double *)malloc(param->nval*sizeof(double)); 
    ensc[i] = (char **)malloc(param->nval*sizeof(char *));
    for (j=0; j<param->nval; j++){
      ensc[i][j] = (char *)malloc(20*sizeof(char));
      flexi_request_key("Configs",0,"%d %lg %s", &param->nlm_conf[i][j],
			&param->wnl_conf[i][j], ensc[i][j]);
    }
  }

  /* [Pseudo] */
  flexi_request_key("Pseudo",1,"%d ", &param->nval);

  param->rc = (double *)malloc(param->nval*sizeof(double));
  param->qc = (double *)malloc(param->nval*sizeof(double));
  param->nb = (int *)malloc(param->nval*sizeof(int));

  for (i=0; i<param->nval; i++) {
    flexi_request_key("Pseudo",1,"%lg",&param->rc[i]);
  }

  /* mth is a character for method type (k or o) */
  flexi_request_key("Pseudo",1,"%c",&param->psmeth);

  for (i=0; i<param->nval; i++) {
    param->qc[i]=0.0;
    param->nb[i]=0;
  }
  
  /* [KBdesign] */
  flexi_request_key("KBdesign",0,"%s %d", &param->local, &param->nboxes);
  param->box_start = (int *)malloc(param->nboxes*sizeof(int));
  param->box_end = (int *)malloc(param->nboxes*sizeof(int));
  param->box_height = (double *)malloc(param->nboxes*sizeof(double));
  b_unit = (char **)malloc(param->nboxes*sizeof(char *));
  b_start = (double *)malloc(param->nboxes*sizeof(double));
  b_end = (double *)malloc(param->nboxes*sizeof(double));
  for (i=0; i<param->nboxes; i++){
    b_unit[i] = (char *)malloc(20*sizeof(char));
    flexi_request_key("KBdesign",0,"%s %lg %lg %lg", 
		      b_unit[i], &b_start[i], &b_end[i], &param->box_height[i]);
  }

  /* 2nd flexi_gather_keys() */
  if (flexi_gather_keys(fp)) return 1;

  /* Prepare for 3rd flexi_gather_keys() */
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

  /* 3rd flexi_gather_keys() */
  if (flexi_gather_keys(fp)) return 1;

  nlpot2_.inl = 0;          /* set AE mode flag */
  for (i=0; i<param->nval; i++) {
    rcall_.rcall[i] = param->rc[i];
    lparam_.qcl[i] = param->qc[i];
    lparam_.nbl[i] = param->nb[i];
  }

  /* Cleanup */
  flexi_clear_keys();
  rewind(fp);

  /* DONE WITH INPUT VALUES */

  /* Post processing section */

  /* set up semicore info */

  param->ipot = (int *)malloc(20*sizeof(int));

  for (i=0; i<4; i++){
    npot[i]=0;
  }
  for (i=0; i<20; i++){
    param->ipot[i]=0;
  }
  param->nll=0;
  
  ncore=param->norb-param->nval;
  for (i=ncore; i<param->norb; i++){
    npot[nlm_label(param->nlm[i]).l]++;
    if (npot[nlm_label(param->nlm[i]).l] == 1) {
      param->ipot[param->nll]=i-ncore;
      rcall_.rcall[param->nll] = param->rc[i-ncore];
      lparam_.qcl[param->nll] = param->qc[i-ncore];
      lparam_.nbl[param->nll] = param->nb[i-ncore];
      param->nll++;
    }
  }

  npm_.nvales = param->nval;
  npm_.ncores = param->norb - param->nval;
  ncore = param->norb - param->nval;
  nll_.nll = param->nll;
  rpcc_.rpcc = param->rpcc;

  /* Section to make KB local idiot proof */

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

  for (i=ncore;i<param->norb;i++) {
    if (loctemp == nlm_label(param->nlm[i]).l) param->local=loctemp;
  }

  /* This is where initial guesses are chosen */
  /* ensave is used as the default eigenvalue if an unbound state is found */

  param->ibound = (int *)malloc(param->norb*sizeof(int));
  param->ensave = (double *)malloc(param->norb*sizeof(double));
  /* process reference config eigenvalue guesses and free arrays */
  for (i=0; i<param->norb; i++){
    param->ibound[i]=1;
    k = 0;
    j = 0;
    
    while ( k<strlen(ens[i]) && ens[i][k]=='-' ) ++k;
    
    if (k==strlen(ens[i])){
      /* user wants opium to guess a good eigenvalue */
      n = param->nlm[i]/100;
      param->en[i] = -1. * (param->z * param->z) / (double)(n*n*n*n);
      param->ensave[i]=0.0;
      
    }else {
      /* user provides an eigenvalue guess */
      param->en[i] = atof(ens[i]);
      param->ensave[i]=param->en[i];
      if (param->en[i] > 0.0) param->ibound[i]=0;
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
	param->ensave_conf[i][j]=0.0;
      }else{
        /* user provides an eigenvalue guess */
        param->en_conf[i][j] = atof(ensc[i][j]);
	param->ensave_conf[i][j]=param->en_conf[i][j];
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
  grid_.r[0] = param->a * pow(param->z, -1./3.);
  for (i=1; i<NPDM; i++) {
    grid_.r[i]=grid_.r[0] * exp(param->b * (i));
    if (grid_.r[i]>80.0) break;
    grid_.np++;
  }
  
  if (grid_.np > param->ngrid) {
    printf("check Non-rel grid; # points needed > max np \n");
    return 1;
  }
  if (grid_.r[grid_.np]<30.0) {
    printf("check Non-rel grid; grid extends to less than 30 Bohr \n");
    return 1;
  }

  param->ngrid=grid_.np;
  
  /* non-rel grid ok, how about rel-grid */
  
  rgrid_.b=param->b2;
  rgrid_.a=param->a2 / (param->b2*param->z);
  
  nrgrid_.nr=0;
  for (i=0; i<NRMAX; i++){
    rgrid_.r[i]=rgrid_.a * (exp(rgrid_.b*(i))-1.0);
    rgrid_.rab[i] =(rgrid_.r[i]+rgrid_.a)*rgrid_.b;
    if (rgrid_.r[i]>80.0) break;
    nrgrid_.nr++;
  }

  if (nrgrid_.nr > param->ngrid2) {
    printf("check Rel grid; # points needed > max np \n");
    return 1;
  }
  if (rgrid_.r[nrgrid_.nr]<30.0) {
    printf("check Rel grid; grid extends to less than 30 Bohr %d \n",nrgrid_.nr);
    return 1;
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


  /* set up DNL stuff */

  local_.iloc = param->local + 1;
  box_.numbox = param->nboxes;

  for (i=0; i<4; i++){
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
  }else if (streq(param->reltype,"frl")){
    fprintf(fp_log, " Fully relativistic pseudopotentials are not implemented \n");
    exit(1);
  }
  if (streq(param->xcparam,"lda")|| streq(param->xcparam,"pzlda")) {
    fprintf(fp_log, " Exchage-correlation functional is : Perdew-Zunger LDA \n");
  }else if (streq(param->xcparam,"pwlda")){
    fprintf(fp_log, " Exchage-correlation functional is : Perdew-Wang LDA \n");
  }else if (streq(param->xcparam,"gga")){
    fprintf(fp_log, " Exchage-correlation functional is : Perdew, Burke, Ernzerhof GGA \n");
  }

  ncore=param->norb - param->nval;

  fprintf(fp_log, " The %c potential is used for the KB construction \n",
	  la[nlm_label(param->nlm[param->local+ncore]).l]);

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
    fprintf(fp_log, " LFC partial core radius: %lg \n",param->rpcc);
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

  if (streq(param->reltype,"nrl")){
    fprintf(fp_log," a_grid=%-11.2e b_grid=%-11.2e # points=%d \n r(1)=%-11.2e   r(np)=%-11.2e \n",
	    param->a,param->b,grid_.np,grid_.r[0],grid_.r[grid_.np]);
  }
  
  if (streq(param->reltype,"srl")){
    fprintf(fp_log," Non-relativitic grid:  a_grid=%-11.2e b_grid=%-11.2e # points=%d \n  r(1)=%-11.2e r(np)=%-11.2e \n",
	    param->a,param->b,grid_.np,grid_.r[0],grid_.r[grid_.np]);
    
    fprintf(fp_log," Relativistic grid   :  a_grid=%-11.2e b_grid=%-11.2e # points=%d \n  r(1)=%-11.2e r(np)=%-11.2e \n",
	    param->a2,param->b2,nrgrid_.nr,rgrid_.r[0],rgrid_.r[nrgrid_.nr-1]);
  }

  fprintf(fp_log, "\n");

  fprintf(fp_log, " dEmax tolerance:  %-11.2e   dVmax tolerance:  %-11.2e  \n",consts_.etol,consts_.vtol);

  fprintf(fp_log, "\n");

  return 0;

}  

double symtoz(char *sym, char *longname) {
  
  double val=0.0;

  if (streq(sym,"H")){
    strcpy(longname,"Hydrogen   ");
    val =  1;

  }else if (streq(sym,"He")){
    strcpy(longname,"Helium     ");
    val =  2;
    
  }else if (streq(sym,"Li")){
    strcpy(longname,"Lithium    ");
    val =  3;
    
  }else if (streq(sym,"Be")){
    strcpy(longname,"Beryllium  ");
    val =  4;
    
  }else if (streq(sym,"B")){
    strcpy(longname,"Boron      ");
    val =  5;
    
  }else if (streq(sym,"C")){
    strcpy(longname,"Carbon     ");
    val =  6;
    
  }else if (streq(sym,"N")){
    strcpy(longname,"Nitrogen   ");
    val =  7;
    
  }else if (streq(sym,"O")){
    strcpy(longname,"Oxygen     ");
    val =  8;
    
  }else if (streq(sym,"F")){
    strcpy(longname,"Fluorine   ");
    val =  9;
    
  }else if (streq(sym,"Ne")){
    strcpy(longname,"Neon       ");
    val = 10;
    
  }else if (streq(sym,"Na")){
    strcpy(longname,"Sodium     ");
    val = 11;
    
  }else if (streq(sym,"Mg")){
    strcpy(longname,"Magnesium  ");
    val = 12;
    
  }else if (streq(sym,"Al")){
    strcpy(longname,"Aluminum   ");
    val = 13;
    
  }else if (streq(sym,"Si")){
    strcpy(longname,"Silicon    ");
    val = 14;
    
  }else if (streq(sym,"P")){
    strcpy(longname,"Phosphorus ");
    val = 15;
    
  }else if (streq(sym,"S")){
    strcpy(longname,"Sulfur     ");
    val = 16;
    
  }else if (streq(sym,"Cl")){
    strcpy(longname,"Chlorine   ");
    val = 17;
    
  }else if (streq(sym,"Ar")){
    strcpy(longname,"Argon      ");
    val = 18;
    
  }else if (streq(sym,"K")){
    strcpy(longname,"Potassium  ");
    val = 19;
    
  }else if (streq(sym,"Ca")){
    strcpy(longname,"Calcium    ");
    val = 20;
    
  }else if (streq(sym,"Sc")){
    strcpy(longname,"Scandium   ");
    val = 21;
    
  }else if (streq(sym,"Ti")){
    strcpy(longname,"Titanium   ");
    val = 22;
    
  }else if (streq(sym,"V")){
    strcpy(longname,"Vanadium   ");
    val = 23;
    
  }else if (streq(sym,"Cr")){
    strcpy(longname,"Chromium   ");
    val = 24;
    
  }else if (streq(sym,"Mn")){
    strcpy(longname,"Manganese  ");
    val = 25;
    
  }else if (streq(sym,"Fe")){
    strcpy(longname,"Iron       ");
    val = 26;
    
  }else if (streq(sym,"Co")){
    strcpy(longname,"Cobalt     ");
    val = 27;
    
  }else if (streq(sym,"Ni")){
    strcpy(longname,"Nickel     ");
    val = 28;
    
  }else if (streq(sym,"Cu")){
    strcpy(longname,"Copper     ");
    val = 29;
    
  }else if (streq(sym,"Zn")){
    strcpy(longname,"Zinc       ");
    val = 30;
    
  }else if (streq(sym,"Ga")){
    strcpy(longname,"Gallium    ");
    val = 31;
    
  }else if (streq(sym,"Ge")){
    strcpy(longname,"Germanium  ");
    val = 32;
    
  }else if (streq(sym,"As")){
    strcpy(longname,"Arsenic    ");
    val = 33;
    
  }else if (streq(sym,"Se")){
    strcpy(longname,"Selenium   ");
    val = 34;
    
  }else if (streq(sym,"Br")){
    strcpy(longname,"Bromine    ");
    val = 35;
    
  }else if (streq(sym,"Kr")){
    strcpy(longname,"Krypton    ");
    val = 36;
    
  }else if (streq(sym,"Rb")){
    strcpy(longname,"Rubidium   ");
    val = 37;
    
  }else if (streq(sym,"Sr")){
    strcpy(longname,"Strontium  ");
    val = 38;
    
  }else if (streq(sym,"Y")){
    strcpy(longname,"Yttrium    ");
    val = 39;
    
  }else if (streq(sym,"Zr")){
    strcpy(longname,"Zirconium  ");
    val = 40;
    
  }else if (streq(sym,"Nb")){
    strcpy(longname,"Niobium    ");
    val = 41;
    
  }else if (streq(sym,"Mo")){
    strcpy(longname,"Molybdenum ");
    val = 42;
    
  }else if (streq(sym,"Tc")){
    strcpy(longname,"Technetium ");
    val = 43;
    
  }else if (streq(sym,"Ru")){
    strcpy(longname,"Ruthenium  ");
    val = 44;
    
  }else if (streq(sym,"Rh")){
    strcpy(longname,"Rhodium    ");
    val = 45;
    
  }else if (streq(sym,"Pd")){
    strcpy(longname,"Palladium  ");
    val = 46;
    
  }else if (streq(sym,"Ag")){
    strcpy(longname,"Silver     ");
    val = 47;
    
  }else if (streq(sym,"Cd")){
    strcpy(longname,"Cadmium    ");
    val = 48;
    
  }else if (streq(sym,"In")){
    strcpy(longname,"Indium     ");
    val = 49;
    
  }else if (streq(sym,"Sn")){
    strcpy(longname,"Tin        ");
    val = 50;
    
  }else if (streq(sym,"Sb")){
    strcpy(longname,"Antimony   ");
    val = 51;
    
  }else if (streq(sym,"Te")){
    strcpy(longname,"Tellurium  ");
    val = 52;
    
  }else if (streq(sym,"I")){
    strcpy(longname,"Iodine     ");
    val = 53;
    
  }else if (streq(sym,"Xe")){
    strcpy(longname,"Xenon      ");
    val = 54;
    
  }else if (streq(sym,"Ce")){
    strcpy(longname,"Cesium     ");
    val = 55;
    
  }else if (streq(sym,"Ba")){
    strcpy(longname,"Barium     ");
    val = 56;
    
  }else if (streq(sym,"La")){
    strcpy(longname,"Lanthanum  ");
    val = 57;
    
  }else if (streq(sym,"Ce")){
    strcpy(longname,"Cerium     ");
    val = 58;
    
  }else if (streq(sym,"Pr")){
    strcpy(longname,"Praseodymium");
    val = 59;
    
  }else if (streq(sym,"Nd")){
    strcpy(longname,"Neodymium  ");
    val = 60;
    
  }else if (streq(sym,"Pm")){
    strcpy(longname,"Promethium ");
    val = 61;
    
  }else if (streq(sym,"Sm")){
    strcpy(longname,"Samarium   ");
    val = 62;
    
  }else if (streq(sym,"Eu")){
    strcpy(longname,"Europium   ");
    val = 63;
    
  }else if (streq(sym,"Gd")){
    strcpy(longname,"Gadolinium ");
    val = 64;
    
  }else if (streq(sym,"Tb")){
    strcpy(longname,"Terbium    ");
    val = 65;
    
  }else if (streq(sym,"Dy")){
    strcpy(longname,"Dysprosium ");
    val = 66;
    
  }else if (streq(sym,"Ho")){
    strcpy(longname,"Holmium    ");
    val = 67;
    
  }else if (streq(sym,"Er")){
    strcpy(longname,"Erbium     ");
    val = 68;
    
  }else if (streq(sym,"Tm")){
    strcpy(longname,"Thulium    ");
    val = 69;
    
  }else if (streq(sym,"Yb")){
    strcpy(longname,"Ytterbium  ");
    val = 70;
    
  }else if (streq(sym,"Lu")){
    strcpy(longname,"Lutetium   ");
    val = 71;
    
  }else if (streq(sym,"Hf")){
    strcpy(longname,"Hafnium    ");
    val = 72;
    
  }else if (streq(sym,"Ta")){
    strcpy(longname,"Tantalum   ");
    val = 73;
    
  }else if (streq(sym,"W")){
    strcpy(longname,"Tungsten   ");
    val = 74;
    
  }else if (streq(sym,"Re")){
    strcpy(longname,"Rhenium    ");
    val = 75;
    
  }else if (streq(sym,"Os")){
    strcpy(longname,"Osmium     ");
    val = 76;
    
  }else if (streq(sym,"Ir")){
    strcpy(longname,"Iridium    ");
    val = 77;
    
  }else if (streq(sym,"Pt")){
    strcpy(longname,"Platinum   ");
    val = 78;
    
  }else if (streq(sym,"Au")){
    strcpy(longname,"Gold       ");
    val = 79;
    
  }else if (streq(sym,"Hg")){
    strcpy(longname,"Mercury    ");
    val = 80;
    
  }else if (streq(sym,"Ti")){
    strcpy(longname,"Thallium   ");
    val = 81;
    
  }else if (streq(sym,"Pb")){
    strcpy(longname,"Lead       ");
    val = 82;
    
  }else if (streq(sym,"Bi")){
    strcpy(longname,"Bismuth    ");
    val = 83;
    
  }else if (streq(sym,"Po")){
    strcpy(longname,"Polonium   ");
    val = 84;
    
  }else if (streq(sym,"At")){
    strcpy(longname,"Astatine   ");
    val = 85;
    
  }else if (streq(sym,"Rn")){
    strcpy(longname,"Radon      ");
    val = 86;
    
  }else if (streq(sym,"Fr")){
    strcpy(longname,"Francium   ");
    val = 87;
    
  }else if (streq(sym,"Ra")){
    strcpy(longname,"Radium     ");
    val = 88;
    
  }else if (streq(sym,"Ac")){
    strcpy(longname,"Actinium   ");
    val = 89;
    
  }else if (streq(sym,"Th")){
    strcpy(longname,"Thorium    ");
    val = 90;
    
  }else if (streq(sym,"Pa")){
    strcpy(longname,"Protactinium");
    val = 91;
    
  }else if (streq(sym,"U")){
    strcpy(longname,"Uranium    ");
    val = 92;
    
  }else if (streq(sym,"Np")){
    strcpy(longname,"Neptunium  ");
    val = 93;
    
  }else if (streq(sym,"Pu")){
    strcpy(longname,"Plutonium  ");
    val = 94;
    
  }else if (streq(sym,"Am")){
    strcpy(longname,"Americium  ");
    val = 95;
    
  }else if (streq(sym,"Cm")){
    strcpy(longname,"Curium     ");
    val = 96;
    
  }else if (streq(sym,"Bk")){
    strcpy(longname,"Berkelium  ");
    val = 97;
    
  }else if (streq(sym,"Cf")){
    strcpy(longname,"Californium");
    val = 98;
    
  }else if (streq(sym,"Es")){
    strcpy(longname,"Einsteinium");
    val = 99;
    
  }else if (streq(sym,"Fm")){
    strcpy(longname,"Fermium    ");
    val =100;
    
  }else if (streq(sym,"Md")){
    strcpy(longname,"Mendelevium");
    val =101;
    
  }else if (streq(sym,"No")){
    strcpy(longname,"Nobelium   ");
    val =102;
    
  }else if (streq(sym,"Lr")){
    strcpy(longname,"Lawrencium ");
    val =103;
    
  }else if (streq(sym,"Rf")){
    strcpy(longname,"Rutherfordium");
    val =104;
    
  }else if (streq(sym,"Db")){
    strcpy(longname,"Dubnium    ");
    val =105;
    
  }else if (streq(sym,"Sg")){
    strcpy(longname,"Seaborgium ");
    val =106;
    
  }else if (streq(sym,"Bh")){
    strcpy(longname,"Bohrium    ");
    val =107;
    
  }else if (streq(sym,"Hs")){
    strcpy(longname,"Hassium    ");
    val =108;
    
  }else if (streq(sym,"Mt")){
    strcpy(longname,"Meitnerium ");
    val =109;
  }
  return val;
}
  
