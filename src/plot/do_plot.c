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
 * $Id: do_plot.c,v 1.9 2004/10/02 18:34:49 ewalter Exp $
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

#define streq(a,b) (!strcasecmp(a,b))

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