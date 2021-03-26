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
 * $Id: opium.c,v 1.14 2004/10/12 20:22:11 ewalter Exp $
 */

/****************************************************************************
* This is OPIUM's frontend                                                  *
*****************************************************************************
*                                                                           *
* MAIN FILE                                                                 *
*                                                                           *
****************************************************************************/

#define VERSION "1.0.3"

/* standard libraries */
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <stdlib.h>
#include <errno.h>
/*const unsigned short int *ctype_b; */

#include "parameter.h"        /* defines structure: 'param_t'   */
#include "nlm.h"              /* nlm_label call */

/* other program modules */
#include "do_ae.h"
#include "do_ps.h"
#include "do_nl.h"
#include "do_tc.h"
#include "do_pwf.h"
#include "do_recpot.h"
#include "do_upf.h"
#include "do_fhi.h"
#include "do_plot.h"
#include "do_ncpp.h"
#include "do_logplt.h"

#define streq(a,b) (!strcasecmp(a,b))

static void read_in_param(param_t *param, char *name, char *logfile, char *parafile);
static int read_stringline(FILE *fp, char *a);
static void do_command(param_t *param, char *paramfile, char *logfile, 
  char *command);
static void do_help();
static void do_chelp();
static void do_phelp();
static void do_khelp();

static int verbosity;

/****************************************************************************
* MAIN                                                                      *
****************************************************************************/

int main(int argc, char *argv[]){

  /**************************************************************************
  * define local variables                                                  *
  **************************************************************************/
  
  time_t t0 = time(&t0);       /* current starting time */
  param_t param;               /* parameters from parameter file */  
  int i;    
  int off;
  int quit;
  int doplot;
  char interactline[160]; /* line string for interactive mode */
  char lastline[160];     /* last line string for interactive mode */
  char command[40];
  FILE *fp_log;           /* log file */
  char *paramfile;

  
  /**************************************************************************
  * check command line options given                                        *
  **************************************************************************/
  
  if ((argc == 2 && (
	   streq(argv[1], "commands") || streq(argv[1], "command") ||
	   streq(argv[1], "-command") || streq(argv[1], "--command") ||
	   streq(argv[1], "comm") || streq(argv[1], "com") || streq(argv[1], "-c") ||
	   streq(argv[1], "--c"))))
  {
    do_chelp();
    exit(0);
  }

  if ((argc == 2 && (
	streq(argv[1], "plotting") || streq(argv[1], "--plot") || streq(argv[1], "-plot") ||
	streq(argv[1], "plot") || streq(argv[1], "plt") || streq(argv[1], "-p") ||
	streq(argv[1], "--p")))){
    do_phelp();
    exit(0);
  }

  if ((argc == 2 && (
        streq(argv[1], "keyblock") || streq(argv[1], "-key") || streq(argv[1], "keys") || streq(argv[1], "--key") ||
	streq(argv[1], "--keys") || streq(argv[1], "-k") || streq(argv[1], "--k"))))
  {
    do_khelp();
    exit(0);
  }

  if ((argc < 3) || (argc > 1 && (streq(argv[1], "help") || streq(argv[1], "-help") ||
	  streq(argv[1], "--help") || streq(argv[1], "-h") || streq(argv[1], "--h"))))
  {
    do_help();
    exit(0);
  }

  /**************************************************************************
  * Open the log file                                                       *
  **************************************************************************/
  
  fp_log = fopen(argv[2], "w");

  /**************************************************************************
  * Print greetings                                                         *
  **************************************************************************/
  
  fprintf(fp_log,
      "\n ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
  fprintf(fp_log,
      "           OPIUM  Version: %10s                             \n", VERSION);
  fprintf(fp_log,
      "           =====================================            \n");
  fprintf(fp_log,
      " See http://opium.sourceforge.net for help and information                 \n");  
  fprintf(fp_log,
      " Copyright 2004 : The OPIUM project                         \n");

  fprintf(fp_log," time of execution: %s\n", ctime(&t0));  
  fprintf(fp_log,
      " ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
  fclose(fp_log);
  
  /**************************************************************************
  * Subroutine call structure                                               *
  **************************************************************************/

  paramfile = (char *) malloc(50*sizeof(char));  
  if (argc>3){  /* command line options */

    /* first read in the parameter file */

    read_in_param(&param, argv[1], argv[2], paramfile);

    /* then process all command line arguments */

    for (i=3; i<argc; i++){
      if (streq(argv[i], "plot")) {
	if (i<argc-1) {
	  do_plot(&param, argv[2], argv[i+1]);
	} else {
	  do_phelp();
	}
	i++;
      } else { 
	do_command(&param, paramfile, argv[2], argv[i]);
      }
    }
  }else{        /* enter interactive mode */
    printf("------------------------------------------------\n");
    printf("------- Welcome to OPIUM version %8s ------\n", VERSION);
    printf("------- opening interactive opium session ------\n");
    printf("------------------------------------------------\n");
    fp_log = fopen(argv[2], "a");
    fprintf(fp_log,"------- opening interactive opium session ------\n");
    fclose(fp_log);
    /* "endless" loop of interactive input */
    quit=0;
    for(;;){
      /* prompt and input */
      printf("opium[param=%s|log=%s|verb=%d]>> ", argv[1], argv[2], verbosity);
      read_stringline(stdin, interactline);
      /* read in the parameter file (again) */
      read_in_param(&param, argv[1], argv[2], paramfile);

      /* check whether '!' was given */
      if (strchr(interactline, '!'))
        strcpy(interactline, lastline);
      /* dissect the input */
      off = 0;
      doplot=0;

      while(off<strlen(interactline) && 
            sscanf(interactline+off, "%s", command)==1){
        off += strlen(command)+1;
	
        if (streq(command, "quit") || streq(command, "exit")){
          quit=1;
          break;
        }
	
        if (streq(command, "plt") || streq(command, "plot")){
	  sscanf(interactline+off,"%s",command);
	  off+=strlen(command)+1;
	  if (do_plot(&param, argv[2], command)==-1) 
	    do_phelp();
	}else{
	  do_command(&param, paramfile, argv[2], command);
	}
      }
      if (quit) break;
      if (strlen(interactline))
        strcpy(lastline, interactline);  /* store the processed line */
    }
    printf("------------------------------------------------\n");
    printf("------- closing interactive opium session ------\n");
    printf("------------------------------------------------\n");
  }
  
  fp_log = fopen(argv[2], "a");
  fprintf(fp_log,"------------- closing opium session ------------\n");
  fclose(fp_log);
  
  return 0;
}  


static void read_in_param(param_t *param, char *name, char *logfile, char *paramfile){
  
  FILE *fp_param; 
  FILE *fp_log = fopen(logfile, "a");
  char *tokp;

  tokp=strtok(name, ".");

  if (tokp != NULL) strcpy(name,tokp);

  strcpy(paramfile,name);
  strcat(paramfile,".param\0");

  fp_param=fopen(paramfile, "r");

  /**************************************************************************
  * read parameter file                                                     *
  **************************************************************************/

  param->name = name;

  if (!fp_param) {
       fprintf(stderr,"Param file %s cannot be opened: %s\n",paramfile,
	       strerror(errno));
    exit(1);
  }
  
  fprintf(fp_log," Reading parameter file:  <%s> ... \n", paramfile);
  if (read_param(param, fp_param, fp_log)==0) fprintf(fp_log," Starting calculation...\n");
  else{ fprintf(fp_log,"   failed! -> exit\n"); exit(1);}
  fclose(fp_param);
  
  fclose(fp_log);
}
  
static int read_stringline(FILE *fp, char *a){
  /* read whole line as string */
  int i = 0;    /* index */
  int c;
  a[0] = '\0';  /* erase string */
  while(((c=getc(fp)) != EOF) && (c!='\n')){
    a[i++]=c;
  }
  a[i]='\0';
  return 0;
}

static void do_command(param_t *param, char *paramfile, char *logfile, 
  char *command){

#define MIN(a, b)   (((a) < (b)) ? (a):(b))  

  FILE *fp;
  FILE *fp_param;
  char filename[160];
  int c;
  int job=-67;

  if (streq(command, "v")) {
    verbosity = !verbosity; 
  } else if (streq(command, "ae")){ 
    do_ae(param, logfile);
    if (verbosity) do_ae_report(stdout);
  } else if (streq(command, "ps")){
    do_ps(param, logfile);
    if (verbosity) do_ps_report(stdout);
  } else if (streq(command, "nl")){
    do_nl(param, logfile);
    if (verbosity) do_nl_report(stdout);
  } else if (streq(command, "tc")){
    do_tc(param, logfile, job);
    if (verbosity) do_tc_report(stdout);
  } else if (streq(command, "pwf")){
    fp = fopen(paramfile, "r"); 
    do_pwf(param, fp, logfile); 
    fclose(fp);
  } else if (streq(command, "recpot")){
    fp = fopen(paramfile, "r"); 
    do_recpot(param, fp, logfile); 
    fclose(fp);
  } else if (streq(command, "upf")) {
    do_upf(param, logfile);
  } else if (streq(command, "fhi")) {
    fp = fopen(paramfile, "r"); 
    do_fhi(param, fp, logfile);
    fclose(fp);
  } else if (streq(command, "ncpp")) {
    fp = fopen(paramfile, "r"); 
    do_ncpp(param, fp,logfile);
    fclose(fp);
  } else if (streq(command, "rpt")) {
    fp = fopen(logfile, "a");
    fprintf(fp,"<<<do_rpt>>>\n");
    fclose(fp);
    sprintf(filename, "%s.rpt", param->name);
    fp = fopen(filename, "w");
    fprintf(fp, "##########################################################\n");
    fprintf(fp, "#    Opium Report File                                   #\n");
    fprintf(fp, "##########################################################\n");
    fprintf(fp, "     Opium version: %10s\n", VERSION);
    fprintf(fp, "\n### copy of the parameter file #######################\n\n");
    fp_param = fopen(paramfile, "r");
    while((c=fgetc(fp_param)) != EOF) fputc(c, fp);
    fclose(fp_param);
    fprintf(fp, "\n### AE report ########################################\n\n");
    do_ae_report(fp);
    fprintf(fp, "\n### PS report ########################################\n\n");
    do_ps_report(fp);
    fprintf(fp, "\n### NL report ########################################\n\n");
    do_nl_report(fp);
    fprintf(fp, "\n### TC report ########################################\n\n");
    do_tc_report(fp);
    fclose(fp);
  }else if (streq(command, "help"))
    do_help();
  else if (streq(command, "plot"))
    do_phelp();
  else if (streq(command, "comm"))
    do_chelp();
  else if (streq(command, "keys"))
    do_khelp();

  else if (streq(command, "all")){
    do_ae(param, logfile); if (verbosity) do_ae_report(stdout);
    do_ps(param, logfile); if (verbosity) do_ps_report(stdout);
    do_nl(param, logfile); if (verbosity) do_nl_report(stdout);
    do_tc(param, logfile, job); if (verbosity) do_tc_report(stdout);
  }else{
    fp = fopen(logfile, "a");
    fprintf(fp,"option [%s] not supported!\n", command);
    fclose(fp);
    printf(">> option [%s] not supported!\n  Type \"opium help\" for help\n\n", command);
  }
}


static void do_help(){
  printf("\t===========================Opium help============================\n");
  printf("\t              version: %10s\n", VERSION);
  printf("\t=================================================================\n");
  printf("\n\t usage:\n");
  printf("\t opium <parameterfile> <logfile> <command line options>\n\n");
  printf("\t Opium can be executed interactively or non-interactively\n");
  printf("\t To enter an interactive session just enter: \n\n");
  printf("\t\t opium <parameterfile> <logfile> \n\n");
  printf("\t In an interactive session... \n");
  printf("\t enter \"comm\" for command line help\n");
  printf("\t enter \"plot\" for plotting help\n");
  printf("\t enter \"keys\" for keyblock help\n\n");

  printf("\t To run Opium non-interactively enter: \n\n");
  printf("\t\t opium <parameterfile> <logfile> <commands>\n\n\n");
  printf("\t In a non-interactive session... \n");
  printf("\t enter \"opium -c\" or \"opium comm\" for command line help\n");
  printf("\t enter \"opium -p\" or \"opium plot\" for plotting help\n");
  printf("\t enter \"opium -k\" or \"opium keys\" for keyblock help\n");

  printf("\t=================================================================\n\n\n");
}

static void do_chelp(){
  printf("\t====================Opium command line help======================\n");
  printf("\t              version: %10s\n", VERSION);
  printf("\t=================================================================\n");
  printf("\n\tatomic solves and pseudopotential construction \n");
  printf("\tae                  - all electron solve of the atom\n");
  printf("\tps                  - generate optimized pseudopotential\n");
  printf("\tnl                  - pseudopotential solve of the \n"
         "\t                        atomic ref. configuration\n");
  printf("\ttc                  - test additional configurations\n");
  printf("\tall                 - run through the complete cycle (ae ps nl tc)\n");
  printf("\n\tpseudo file output style \n");
  /*  printf("\tupf      - generate *.upf  output\n"); */
  printf("\tpwf                 - generate *.pwf  output\n");
  printf("\trecpot              - generate *.recpot  output\n");
  printf("\tfhi                 - generate *.fhi  output\n");
  printf("\tncpp                - generate *.ncpp output\n");
  printf("\n\tmiscellaneous options \n");
  printf("\trpt                 - generate report file\n");
  printf("\tplot [plot_type]    - make a plot of type [plot_type]\n" 
         "\t                        (do \"opium -p\" for more information) \n");
  printf("\tv                   - toggle verbosity flag (interactive mode only)\n");
  printf("\t=================================================================\n\n\n");
}

static void do_phelp(){
  printf("\t==================Opium plotting commands help===================\n");
  printf("\t              version: %10s\n", VERSION);
  printf("\t=================================================================\n");
  printf("\n\tusage examples: \n");
  printf("\n\t\t >> opium c log ae ps nl plot wp \n");
  printf("\n\t\t >> opium o log ae plot wa \n");
  printf("\n\t\t >> opium o log ae ps nl plot vi rpt \n");
  printf("\n\n\tCurrent plot types: \n");
  printf("\n\n\t wa   - all-electron wavefunctions \n");
  printf("\t wp   - pseudo & all-electron wavefunctions \n");
  printf("\t pcc  - core, valence density and/or partial core density\n");
  printf("\t den  - same as pcc \n");
  printf("\t vs   - screened pseudopotentials \n");
  printf("\t vi   - ionic psuedopotentials (descreened) \n");
  printf("\t qp   - q-space pseudowavefunctions and potentials \n");
  printf("\t logd - log. deriv. for state indentified in [Loginfo] keyblock \n");
  printf("\t=================================================================\n\n\n");
}

static void do_khelp(){
  printf("\t=======================Opium keyblock help=======================\n");
  printf("\t              version: %10s\n", VERSION);
  printf("\t=================================================================\n");
  printf("\n\tkeyblock help: \n");
  printf("\n\tThe following keyblocks are avaliable in the version of Opium...\n\n");
  printf("\t[Atom]                                                          \n");
  printf("\t  atom symbol(1 or 2 characters)                                  \n");
  printf("\t  number of reference orbitals(int)                             \n");
  printf("\t  nlm(int), occupation(float), eig. guess(float or - ) for orb 1 \n");
  printf("\t  .                                 .                            . \n");
  printf("\t  .                                 .                            . \n");
  printf("\t  nlm(int), occupation(float), eig. guess(float or - ) for orb n \n\n ");
  printf("\t  NOTE: To indicate a virtual/unbound/Hamann type orbital, make occupation < 0.0 \n ");
  printf("\t  and specify an eigenvalue in the guess column. \n\n\n");
  printf("\t[Pseudo]\n");
  printf("\t  number of orbitals in pseudopotetntial(int) \n");                         
  printf("\t  cut-off radius for pseudo orbital 1 (float) \n");
  printf("\t  .                   .             .         \n");
  printf("\t  cut-off radius for pseudo orbital n (float) \n");
  printf("\t  (o)ptimized or (k)erker--pseud. method (only 1st letter needed) \n\n\n");
  printf("\t[Optinfo]                                                        \n");
  printf("\t  cut-off wavevector(float), # bessel fxns(int) for pseudo orb 1 \n");  
  printf("\t  .                            .                                . \n");
  printf("\t  cut-off wavevector(float), # bessel fxns(int) for pseudo orb n \n\n\n");  
  printf("\t[XC]                                                                 \n");
  printf("\t  pzlda or pwlda or gga  --xc funct. (default is pzlda, lda=pzlda) \n\n\n");
  printf("\t[Pcc]                                               \n");
  printf("\t partial core radius(float) (default = no pcc)  \n");
  printf("\t partial core method (character) lfc or fuchs  \n");
  printf("\t lfc=Louie, Froyen & Cohen / fuchs=Fuchs & Scheffler  (default=lfc)\n\n\n");
  printf("\t[Relativity] \n");
  printf("\t  nrl or srl -- non-rel or scalar-rel solve (default = nrl) \n\n\n");
  printf("\t[Grid]                                                 \n");
  printf("\t  max grid points(int), a(float)(see documentation), b log. spacing (float) \n");
  printf("\t  default = 1201 0.0001 0.013 \n");
  printf("\t    note: this grid is used for all non-rel solves and psp construction \n\n\n");
  printf("\t[Relgrid]                                                         \n");
  printf("\t  same as [Grid] but for the non-relativisitc solve\n");
  printf("\t  default = 2201 0.0001 0.0143 \n\n\n");
  printf("\t[Tol]                                                             \n");
  printf("\t  tol. for eig.(float), tol for V(r) \n");
  printf("\t  default 1e-8 and 1e-6 \n");
  printf("\t  used for AE and NL solve stopping conditions \n\n\n");
  printf("\t[Configs]                                                         \n");
  printf("\t  number of configurations to test(int) \n");
  printf("\t  nlm(int), occupation(float), eig. guess(float or - ) for orb 1 config 1 \n");
  printf("\t  .                                 .                                   . \n");
  printf("\t  nlm(int), occupation(float), eig. guess(float or - ) for orb n config 1 \n\n\n");
  printf("\t  nlm(int), occupation(float), eig. guess(float or - ) for orb 1 config 2 \n");
  printf("\t  .                                 .                                   . \n");
  printf("\t  .                                 .                                   . \n");
  printf("\t  nlm(int), occupation(float), eig. guess(float or - ) for orb n config n \n\n");
  printf("\t  NOTE:  valence orbitals must be in same order as the reference state \n\n\n");

  printf("\t[KBdesign]                                                        \n");
  printf("\t  angular momentum of the local potential used for KB construction \n");
  printf("\t  specified by letter or l-value (i.e. 0 or s, 1 or p etc.) (default=0)\n\n");
  printf("\t  number of boxes used in designed non-local procedure(int) (default=0, no boxes)\n");
  printf("\t  units used for 1st box placement, should be either gp (grid points) or au (bohr) \n");
  printf("\t  left edge of 1st box(int or float), right edge of 1st box(int or float) \n");
  printf("\t  depth of 1st box(float) in rydberg units \n");
  printf("\t  .                                       .                                      . \n");
  printf("\t  .                                       .                                      . \n");
  printf("\t  units used for n\'th box placement, should be either gp (grid points) or au (bohr) \n");
  printf("\t  left edge of n\'th box(int or float), right edge of n\'th box(int or float) \n");
  printf("\t  depth of n'th box(float) in rydberg units \n\n\n");


  printf("\t[Loginfo]                                            \n");
  printf("\t  configuration number used in log. deriv. plot(int) (reference=0) \n");
  printf("\t  radius for log. deriv.(float), min and max energy for log deriv.(2 floats) \n\n");

  printf("\t  See $OPIUM/docs/keyblocks or http://opium.sourceforge.net for more help \n\n");
  printf("\t=================================================================\n\n\n");
}

