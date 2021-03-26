/*
 * $Id: nwrite.c,v 1.3 2004/06/16 20:46:17 mbarnes Exp $
 */

/****************************************************************************
* This file is part of the NEW pseudo code                                  *
*****************************************************************************
*                                                                           *
* Simulate Fortran output of implicite loops within write statements        *
* -> This functionality is necessary in order to create BH comlient PWFs    *
****************************************************************************/


#include <stdio.h>
#include <string.h>


static int counter;

void nwrite_(FILE *fp, char *linetag, int *n, char *format, double *x){
  
  if (counter==0) fprintf(fp, "%s", linetag);
  
  fprintf(fp, format, *x);
  
  if (counter<*n-1)
    ++counter;
  else{
    fprintf(fp, "\n");
    counter=0;
  }
}

void nclear_(FILE *fp){
  if (counter != 0){
    fprintf(fp, "\n");
    counter = 0;
  }
}  
