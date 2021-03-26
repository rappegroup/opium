/*
 * Copyright (c) 1998-2005 The OPIUM Group
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
 * $Id: parameter.h,v 1.7 2004/10/02 18:34:49 ewalter Exp $
 */

#ifndef __INCLUDE_PARAMETER_H
#define __INCLUDE_PARAMETER_H

/****************************************************************************
* Read the parameter settings from 'fp' using 'FlexiLib' and set default    *
* values.                                                                   *
*                                                                           *
****************************************************************************/


/* define param type structure */

typedef struct param_t{

  /* [Element] */
  char
    *name,        /* name */
    *symbol,      /* symbol */
    *longname;

  /* [Relativity] */
  char
    *reltype;     /* relativistic type */

  /* [Grid] */
  int  
    ngrid;        /* number of grid points */
  double  
    a,b;          /* grid parameters: r[i] = a*exp(b*i) */

  /* [RelGrid] */
  int  
    ngrid2;        /* number of grid points */
  double  
    a2,b2;          /* grid parameters: r[i] = a*exp(b*i) */
    
  /* [Atom] */
  int  
    norb,         /* number of orbitals */
    *nlm,         /* orbital index */
    *ibound;
  double
    z,            /* atomic number Z */  
    *wnl,         /* orbital occupation */
    *en,
    *ensave;          /* orbital energy guess */

  /* [Configs] */
  int  
    nconfigs,     /* number of test configurations */
    **nlm_conf;   /* orbital index */
  double
    **wnl_conf,   /* orbital occupation */
    **en_conf;    /* orbital energy guess */

  /* [Pseudo] */
  int 
    nval,         /* number of valence orbitals */
    *nb;          /* Nb [?] */
  double  
    *rc,          /* rc [A] */
    *qc,          /* qc [?] */
    rpcc;         /* rppc [A] */
  char psmeth;
  char optmeth;

  /* [LogOpt] */
  int  
    *numen;       /* number of points */
  double
    *wnorm,       /* weighting factor */
    *wke,         /* rapid weighting factor */
    **escat,      /* matching energy */
    **rscat,      /* matching radius */
    **wscat;      /* matching weighting factor */

  /* [LogInfo] */
  int  
    ilogder;      /* ? */
  double  
    rphas,        /* ? */
    emin,         /* ? */
    emax;         /* ? */

  /* [KBdesign] */
  int  
    local,localind,        /* local orbital [0...< nval[ */
    nboxes,       /* number of design boxes */
    *box_start,
    *box_end;
  
  double
    *box_height;

  /* [XC] */
  char *xcparam;
  int  ixc;     /* xc parametrization */
  double rxccut; /* origin behavior */

  /* [Misc] */
  int  
  ioutput,      /* output flag */
    irelvirt,     /* relativ. virt. crystal stuff ???? */
    ist1,         /* local pot only flag */
    ncut;         /* number of grid points in *.pwf file arrays */
  
  /* [Conmax] */
  int
    switch1,
    switch2,
    switch3,
    limsm,
    nstep;
  double
  encsm,
    tolcon;
  
  
/* the following block(s) are temporary and will be removed in later versions */
  
  /* [AERelMethod] */
  int  
    aerelmethod,  /* method flag */
    aereltransf;  /* transfer flag */
  
  
  /* auxiliary information NOT directly read from parameter file */
  int
  *ipot,
    nll;          /* number of different angular momenta */
    
}param_t;


/* read the parameter settings from 'fp' */

int read_param(param_t *param, FILE *fp, FILE *fp_log);

#endif /* __INCLUDE_PARAMETER_H */
