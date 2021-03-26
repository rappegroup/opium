/*
 * $Id: do_ncpp.c,v 1.4 2004/07/23 14:52:33 ewalter Exp $
 */

/****************************************************************************
 * generate *.ncpp output (PWSCF Norm-Conserving compatible)                *
 *****************************************************************************
 * INPUT:  param_t structure, *.ps, *.loc
 * OUTPUT: *.ncpp
 ****************************************************************************/

#include <stdio.h>
#include <math.h>
#include <string.h>

#ifndef M_PI 
#define M_PI 3.14159265358979323846
#endif

#include "parameter.h"        /* defines structure: 'param_t' */
#include "fortparam.h"        /* fortran code parameters */
#include "do_ncpp.h"          /* the module's own header */

int do_ncpp(param_t *param, char *logfile){

  int i,j,k,ic;
  char filename[180];
  FILE *fp;
  FILE *fp_log;
  double zeff;	                /* effective Z */
  int lmax;                     /* max l value */   
  char lcore;                   /* T/F core correction? */ 
  char xctype[5];               /* pz or pbe for now */
  double ztt;                   /* store z^2/3 */
  double xmin;                  /* ncpp grid */
  double r_1;
  int lchi;
  double locc ;                 /* l and occ for wavefunctions */
  
  static double rscore[NPDM];
  static double nlcore[NPDM];
  static double rvcore[NVALE0+1][NPDM];
  static double rnl[N0][NPDM];
  static double r[NPDM];

  fp_log = fopen(logfile, "a");
  fprintf(fp_log,"<<<do_ncpp>>>\n");  

  /* section 0 : check if this potential format is supported */
  /* ncpp doesn't support DNL potentials */

  if (param->nboxes > 0) {
    fprintf(fp_log," ncpp format does not support the use of augmentation operators\n");
    return 1;
  }

  /* read in the local potential from a binary file created by do_nl() */
  sprintf(filename, "%s.loc", param->name);
  fp = fopen(filename, "rb");
  fread(nlcore, sizeof(double), param->ngrid, fp);
  fclose(fp);
  
  sprintf(filename, "%s.psi_nl", param->name);
  fp = fopen(filename, "rb");
  for (i=0; i<param->nval; i++)
    fread(rnl[i], sizeof(double), param->ngrid, fp);
  fclose(fp);

  sprintf(filename, "%s.pot_ps", param->name);
  fp = fopen(filename, "rb");
  for (i=0; i<param->nval; i++) {
    fread(rvcore[i], sizeof(double), param->ngrid, fp);
  }	
  fclose(fp);

  sprintf(filename, "%s.rho_core", param->name);
  fp = fopen(filename, "rb");
  fread(rscore, sizeof(double), param->ngrid, fp);
  fclose(fp);

  /* section 1 : header */

  sprintf(filename, "%s.ncpp", param->name);
  fp = fopen(filename, "w");

  /* 1st line of ncpp format: the XC type and a comment */  

  if ((!strcmp(param->xcparam, "lda") || (!strcmp(param->xcparam, "lda"))))
    sprintf(xctype,"%s"," \'pz\'");
  else  
    sprintf(xctype,"%s","\'pbe\'");
  fprintf(fp,"%s %s --Opium generated potential--\n",xctype,param->name);  
  
  /* 2nd line of ncpp format: symbol, zeff,lmax,0 (for analytic),*/
  /* 0 (for analytic), .T./.F. core-correction?,l-local,.F. (bhstype)*/

  /* compute zeff */
  zeff = param->z;
  for (i=0; i<param->norb - param->nval; i++)
    zeff -= param->wnl[i];

  /* lmax is total num ang mom - 1 */
  lmax = param->nll - 1;

  /* check if core-correction */
  if (param->rpcc > 1e-12)
    lcore = 'T';
  else
    lcore = 'F';

  fprintf(fp,"\'%s\',%f,%d,0,0,.%c.,%d,.F. \n",param->symbol,zeff,lmax,lcore,param->local);  

  /* 3rd line of ncpp format: Z used in mesh, xmin , dx ,np , */
  /* num wavefunctions */

  /* The grid specification is different: */
  /* our grid  r_i = (a * z^-1/3) * exp((i-1)*b) */ 
  /* (a and b used to be called r1 and h) */
  /* ncpp grid r_n = (exp( (xmin + (n-1)*b)) / z */
  /*               = (exp(xmin)/z * exp((n-1)*b) */
  /* therefore the b's are the same and, to get the same grid*/
  /* a * z^-1/3 = exp(xmin)/z  or */
  /* ln(z^2/3 * a) = xmin  --EJW */
  
  ztt = pow(param->z,2./3.);
  xmin = log(ztt * param->a);

  fprintf(fp,"%f %f %f %d %d\n",param->z,xmin,param->b,param->ngrid,param->nval);  

  /* section 2: V_nl */

   r_1 = param->a * pow(param->z,-1./3.);
   for (k=0; k<param->ngrid ; k++) {
     r[k] = r_1 * exp(param->b * k);
   }
   for (j=0; j<lmax+1; j++) {
     fprintf(fp," Pseudo l=%d\n",j);
     for (k=0; k<param->ngrid ; k++) {
       fprintf(fp, "%1.15e  ", rvcore[j][k]/r[k]);
       if (!((k+1)%4)) fprintf(fp, "\n");
     }
     if (k%4) fprintf(fp, "\n");
   }

  /* section 3: partial core */
  if (param->rpcc > 1e-12) {
    for (k=0; k<param->ngrid; k++) {
      fprintf(fp, "%1.15e  ", rscore[k]/(r[k]*r[k]*4*M_PI));
      if (!((k+1)%4)) fprintf(fp, "\n");
    }
     if (k%4) fprintf(fp, "\n");
  }

  /* section 4: Wavefunctions */
  for (j=0; j<param->nval; j++) {
    i=param->norb - param->nval + j;
    lchi = (param->nlm[i]%100)/10;
    locc = param->wnl[i]; 
    fprintf(fp," Wavefunction %d\n",j+1);
    fprintf(fp," %d  %f\n",lchi,locc);
    for (k=0; k<param->ngrid; k++) {
      fprintf(fp, "%1.15e  ", rnl[j][k]);
      if (!((k+1)%4)) fprintf(fp, "\n");
    }
    if (k%4) fprintf(fp, "\n");
  }

  fclose(fp);
  fclose(fp_log);
  return 0;
}

