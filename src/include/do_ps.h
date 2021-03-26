/*
 * $Id: do_ps.h,v 1.3 2004/06/16 20:46:17 mbarnes Exp $
 */

#include "parameter.h"

int do_ps(param_t *param, char *logfile);
void do_ps_report(FILE *fp);
void readAE(param_t *param);

