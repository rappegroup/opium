/*
 * Copyright (c) 1998-2012 The OPIUM Group
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
  
  if (counter==0) fprintf(fp, "%s", "E");

  fprintf(fp, format, *x);
  
  if (counter<*n-1)
    ++counter;
  else{
    fprintf(fp, "\n");
    counter=0;
  }
}
void nwrite2_(FILE *fp, int *n, char *format, double *x){
  
  fprintf(fp, format, *x);
  
  if (counter<*n-1)
    ++counter;
  else{
    fprintf(fp, "\n");
    counter=0;
  }
}

void iwrite_(FILE *fp, char *linetag, int *n, char *format, int *x){
  
  if (counter==0) fprintf(fp, "%s", "");
  
  fprintf(fp, format, *x);
  
  if (counter<*n-1)
    ++counter;
  else{
    fprintf(fp, "\n");
    counter=0;
  }
}
void iwrite2_(FILE *fp, int *n, char *format, int *x){
  
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
