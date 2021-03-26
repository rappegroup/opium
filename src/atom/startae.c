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

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>

#include "parameter.h"
#include "cdim.h"
#include "common_blocks.h"   

void startae(param_t *param, int norb) {
  
  double y=0.0;
  int i,j;


  if ((param->z - adat_.xion - 1.0) > 0.0) 
    y=pow((param->z - adat_.xion - 1.0),0.4);
  
  for (i=0; i<param->ngrid+1; i++) {
    double t=grid_.r[i]/0.675;
    if (t>100) 
      t=100.0;
    t=2.0*(param->z - adat_.xion - 1.0) * (1.0 - 1.0/(0.675*y*(exp(t)-1.0)+1.0))-2*param->z;
    for (j=0; j<norb; j++){
      totpot_.rvps[j][i] = t;
      totpot_.rvcore[j][i] = -2.0* param->z; 
    }
    totpot_.rvcoul[i] =  totpot_.rvps[0][i]-totpot_.rvcore[0][i];
  }
}

