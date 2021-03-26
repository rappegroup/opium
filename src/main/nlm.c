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
