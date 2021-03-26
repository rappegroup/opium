/*
 */

/****************************************************************************
* Part of uniPP Library providing access to universal pseudopotential files.*
*                                                                           *
*                                                                           *
* last touched: 13.12.00 gjt                                                *
****************************************************************************/


/****************************************************************************
*  'uniPPlib' - Universal PseudoPotential Library                           *
*  Copyright (C) 2000  University of California, Santa Barbara              *
*                                                                           *
*  This program is free software; you can redistribute it and/or modify     *
*  it under the terms of the GNU General Public License as published by     *
*  the Free Software Foundation; either version 2 of the License, or        *
*  (at your option) any later version.                                      *
*                                                                           *
*  This program is distributed in the hope that it will be useful,          *
*  but WITHOUT ANY WARRANTY; without even the implied warranty of           *
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            *
*  GNU General Public License for more details.                             *
*                                                                           *
*  You should have received a copy of the GNU General Public License        *
*  along with this program; if not, write to the Free Software              *
*  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA  *
****************************************************************************/


#include <stdio.h>
#include <stdlib.h>

#include "xmcomplex.h"

#include "uniPPlib.h"


void uniPP_free(uniPP *unipp){

  int i;

  for (i=0; i<unipp->l_max; i++){
    free(unipp->u_ps[i]);
    free(unipp->v_ps[i]);
  }
  
  free(unipp->r_m);
  free(unipp->u_ps);
  free(unipp->v_ps);
  if (unipp->nlcc)
    free(unipp->n_pc);
}
  
