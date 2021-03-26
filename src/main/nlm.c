/*
 * $Id: nlm.c,v 1.3 2004/06/16 20:46:17 mbarnes Exp $
 */

#include "nlm.h"


struct nlm nlm_label(int nlm){
  /* determine the labels for nlm */
static  struct nlm labels;
  labels.n = nlm/100; /* main quantum numbe */
  nlm -= labels.n * 100;
  labels.l = nlm/10;  /* orbital quantum number */
  nlm -= labels.l * 10;
  labels.m = nlm;     /* magnetic quantum number */
  return labels;
}
