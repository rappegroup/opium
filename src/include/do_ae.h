/*
 * $Id: do_ae.h,v 1.2 2004/06/16 20:46:17 mbarnes Exp $
 */

/****************************************************************************
  Copyright (C) 2002 Andrew M Rappe group, University of Pennsylvania
  This file is distributed under the terms of the GNU General Public
  License as described in the file 'License' in the current directory.
****************************************************************************/

#include "parameter.h"

int do_ae(param_t *param, char *logfile);
void do_ae_report(FILE *fp);
