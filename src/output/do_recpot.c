/*
 * Copyright (c) 1998-2004 The OPIUM Group
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
 * $Id: do_recpot.c,v 1.5 2004/10/10 22:10:06 ewalter Exp $
 */

/*                                                                          */
/*generate *.recpot output for CASTEP                                       */
/****************************************************************************/
/*INPUT:  param_t structure, *.ps, *.loc                                    */
/*OUTPUT: *.pwf                                                             */
/****************************************************************************/

/* standard libraries */
#include <stdio.h>

#include "parameter.h"        /* defines structure: 'param_t' */
#include "fortparam.h"        /* fortran code parameters */
#include "do_recpot.h"           /* the module's own header */
#include "common_blocks.h"    /* fortran common blocks */

/* from atom/do_nl.c */
void nrelorbnl(param_t *param, int);
/* creative passing of the 1st argument to follow */
void writerecpot_(void *,void *,  char *);

int do_recpot(param_t *param, FILE *fp_param, char *logfile){

  int i;              /* loop counter */
  char filename[180]; /* filename */
  int c,config;              /* dummy character */
  FILE *fp;           /* file pointer */
  FILE *fp_log; 

  fp_log = fopen(logfile, "a");
  fprintf(fp_log,"<<<do_recpot>>>\n");
  if (param->rpcc > 1e-12) {
    fprintf(fp_log,"The recpot format does not support a partial core correction for now.  Sorry!\n");
    fprintf(fp_log, "   ================================================\n");
    fclose(fp_log);
    return 1;
  }
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

  if (param->nval == 1)
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

  sprintf(filename, "%s.rho_core", param->name);
  fp = fopen(filename, "rb");
  fread(rscore_.rscore, sizeof(double), param->ngrid, fp);
  fclose(fp);
      
  fprintf(fp_log,"<<calling: writerecpot>>\n");
  fclose(fp_log);
  sprintf(filename, "%s.recpot", param->name);
  fp = fopen(filename, "w");

  writerecpot_(fp, fp_param, &param->psmeth);
  fclose(fp);

  fp_log = fopen(logfile, "a");
  fprintf(fp_log, "   ================================================\n");
  fclose(fp_log);
  return 0;
}

void writeparam_(FILE *fp, FILE *fp_param, double *ecutev) {
  int c;
  int cut;

  fprintf(fp, "START COMMENT \n");

  cut = 0.8 * *ecutev;
  fprintf(fp, "%-4d  COARSE \n",cut);

  cut = 0.9 * *ecutev;
  fprintf(fp, "%-4d  MEDIUM\n",cut);

  cut = *ecutev;
  fprintf(fp, "%-4d  FINE\n",cut);

  fprintf(fp, "Pseudopotential generated by OPIUM \n");

  fprintf(fp, "############################################################\n");

  /* append the parameter file to the header of the recpot file */
  fprintf(fp, "\n");
  fprintf(fp, "############################################################\n");
  fprintf(fp, "#    Opium Parameter File                                  #\n");
  fprintf(fp, "############################################################\n");
  fprintf(fp, "\n");
  rewind(fp_param);
  while((c=fgetc(fp_param)) != EOF) fputc(c, fp);


  fprintf(fp, "END COMMENT \n");
  fprintf(fp, " 3 5 \n");

}