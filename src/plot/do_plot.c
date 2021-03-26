/*
 * $Id: do_plot.c,v 1.7 2004/07/06 15:29:51 ewalter Exp $
 */

#include <stdio.h>
#include <string.h>
#include <math.h>

#include "parameter.h"        /* defines structure: 'param_t' */
#include "nlm.h"              /* nlm_label call */
#include "fortparam.h"        /* fortran code parameters */
#include "do_plot.h"            /* the module's own header */
#include "common_blocks.h"    /* fortran common blocks */
#include "energy.h"           /* this is for saving the energies */
#include "do_wplot.h"
#include "do_vplot.h"
#include "do_pccplt.h"
#include "do_logplt.h"
#include "do_qplot.h"

#define streq(a,b) (*a==*b && !strcmp(a+1,b+1))

int do_plot(param_t *param, char *logfile, char *plot){
  
  if (streq(plot, "wa")) {
    do_wplot(param, logfile,"a");
  }else if (streq(plot, "wp")) {
    do_wplot(param, logfile,"n"); 
  }else if (streq(plot, "vs")) {
    do_vplot(param, logfile,"s"); 
  }else if (streq(plot, "vi")) {
    do_vplot(param, logfile,"i"); 
  }else if ((streq(plot, "pcc"))||(streq(plot, "den"))) {
    do_pccplt(param, logfile); 
  }else if (streq(plot, "logd")) {
    do_logplt(param, logfile); 
  }else if (streq(plot, "qp")) {
    do_qplot(param, logfile);
  }else{
    fprintf(stderr,"Wrong plot type %s\n", plot);
    return -1;
  }

  return 0;

}


/* vim: cindent sw=2 showmatch
 */
