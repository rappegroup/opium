/*
 * $Id: parameter.h,v 1.4 2004/06/16 20:46:17 mbarnes Exp $
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
    **en_conf,    /* orbital energy guess */
    **ensave_conf;    /* orbital energy guess */

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
    local,        /* local orbital [0...< nval[ */
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
