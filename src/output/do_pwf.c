/*
 * Copyright (c) 1998-2005 The OPIUM Group
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
/* 
 * $Id: do_pwf.c,v 1.6 2004/10/02 18:34:49 ewalter Exp $
 */

/*                                                                          */
/*generate *.pwf output                                                     */
/****************************************************************************/
/*INPUT:  param_t structure, *.ps, *.loc                                    */
/*OUTPUT: *.pwf                                                             */
/****************************************************************************/

/* standard libraries */
#include <stdio.h>

#include "parameter.h"        /* defines structure: 'param_t' */
#include "cdim.h"        /* fortran code parameters */
#include "do_pwf.h"           /* the module's own header */
#include "common_blocks.h"    /* fortran common blocks */

/* from atom/do_nl.c */
void nrelorbnl(param_t *param, int);
/* creative passing of the 1st argument to follow */
void writepwf_(void *, char *);

int do_pwf(param_t *param, FILE *fp_param, char *logfile){

  int i;              /* loop counter */
  char filename[180]; /* filename */
  int c,config;              /* dummy character */
  FILE *fp;           /* file pointer */
  FILE *fp_log; 

  fp_log = fopen(logfile, "a");
  fprintf(fp_log,"<<<do_pwf>>>\n");
  
  /* set the log file */
  sprintf(filenames_.file_log, "%s", logfile);
  
  /* set the pwf common block */
  /*  sprintf(filenames_.file_ps, "%s.ps", param->name);*/
  pwf_.zeff = param->z;
  for (i=0; i<param->norb-param->nval; i++)
    pwf_.zeff -= param->wnl[i];
  pwf_.rpcc = param->rpcc;  
  pwf_.igrid = param->ngrid;
  pwf_.spacing = param->b;
  pwf_.first = param->a;
  pwf_.nval = param->nll;
  if (param->nboxes > 0)
    pwf_.inl = 1;
  else
    pwf_.inl = 0;

  if (param->nll == 1)
    param->ist1 = 3;
  else
    param->ist1 = 2;

  pwf_.ist1 = param->ist1;
  pwf_.ilocal = param->local;
  pwf_.ilocalind = param->localind;

  param->ncut = param->ngrid;
  pwf_.ncut = grid_.np;

  np_.nvales=param->nval;
  np_.ncores=param->norb-param->nval;
  npp_.npots=param->nll;

  config=-1;
  nrelorbnl(param,config);

  /* read in the local potential from a binary file created by do_nl() */
  sprintf(filename, "%s.loc", param->name);
  fp = fopen(filename, "rb");
  fread(nlcore_.rvloc, sizeof(double), param->ngrid, fp);
  fclose(fp);

  sprintf(filename, "%s.psi_nl", param->name);
  fp = fopen(filename, "rb");
  for (i=0; i<param->nval; i++)
    fread(atomic_.rnl[i], sizeof(double), param->ngrid, fp);
  fclose(fp);

  sprintf(filename, "%s.pot_ps", param->name);
  fp = fopen(filename, "rb");
  for (i=0; i<param->nval; i++) {
    fread(totpot_.rvcore[i], sizeof(double), param->ngrid, fp);
  }	
  fclose(fp);

  if (param->rpcc > 0.){
    sprintf(filename, "%s.rho_pcore", param->name);
    fp = fopen(filename, "rb");
    fread(rscore_.rscore, sizeof(double), param->ngrid, fp);
    fclose(fp);
  }

  fprintf(fp_log,"<<calling: writepwf>>\n");
  fclose(fp_log);
  sprintf(filename, "%s.pwf", param->name);
  fp = fopen(filename, "w");
  writepwf_(fp, param->symbol);
  
  /* append the parameter file to the end of the pwf file */
  fprintf(fp, "\n");
  fprintf(fp, "############################################################\n");
  fprintf(fp, "#    Opium Parameter File                                  #\n");
  fprintf(fp, "############################################################\n");
  fprintf(fp, "\n");
  rewind(fp_param);
  while((c=fgetc(fp_param)) != EOF) fputc(c, fp);
  fprintf(fp, "############################################################\n");
  fclose(fp);
  
  fp_log = fopen(logfile, "a");
  fprintf(fp_log, "   ================================================\n");
  fclose(fp_log);
  return 0;
}
