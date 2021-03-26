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
 * $Id: common_blocks.h,v 1.7 2004/10/02 18:34:48 ewalter Exp $
 */

/***************************************************************************
* these are the common blocks from the fortran code                        *
****************************************************************************/
#ifndef COMMON_BLOCKS_H
#define COMMON_BLOCKS_H


/* 37  */
extern struct{
  double h, r1, z, r[NPDM];
  int np;
}grid_;  

/* 23  */
extern struct{
  double xion,rnl[N0][NPDM];
  int nlm[N0];
  double wnl[N0],en[N0];
  int norb;
}atomic_;  

/* 18  */
extern struct{
  double rvcore[N0][NPDM],rvps[N0][NPDM],rvcoul[NPDM];
}totpot_;

/* 19  */
extern struct{
  int nmax[NVALE0];
  int maxim;
}nmax_;

/* 13  */
extern struct{
  char file_log[80];
}filenames_;

/* 11  */
extern struct{
  int ncores, nvales;
}np_;  

/* 10  */
extern struct{
  int ibd[N0];
}ibound_;  

/* 10  */
extern struct{
  int inl,indrc[N0],IB,ID,IG;
}nlpot2_;  

/* 10  */
extern struct{
  double rcall[N0];
}rcall_;

/* 8  */
extern struct{
  double rpcc,rpccz;
}rpcc_;  

/* 8  */
extern struct{
  double rscore[NPDM],rdd[NPDM],rddd[NPDM],rscoretot[NPDM];
}rscore_;  

/* 7  */
extern struct{
  int ncores, nvales, iskip, maxip;
}npm_;  

/* 7  */
extern struct{
  double rvloc[NPDM];
}nlcore_;      

/* 6  */
extern struct{
  int    iboxstart[10],iboxend[10];
  double boxheight[10];
  int    numbox;
}box_;  

/* 6  */
extern struct{
  double rphas,elogmax,elogmin,dlwf[N0][NPL0];
}logarith_;  

/* 5  */
extern struct{
  double rnorm[30],etot;
  int lghost[30];
}results_;  

/* 5  */
extern struct{
  int  norb, ncore, nval;
  int  no[N0], lo[N0];
}reli_;  

/* 5  */
extern struct{
  double zcore,zval;
  double zo[N0],so[N0];
}reld_;

/* 5  */
extern struct{
  double a,b,r[NPDM],rab[NPDM];
}rgrid_;

/* 5  */
extern struct{
  int nr;
}nrgrid_;

/* 5  */
extern struct{
  int  ilogder;
}ilogder_;  

/* 4  */
extern struct{
  int iloc, idesign;
}local_;  

/* 3  */
extern struct{
  double etol, vtol;
  int maxit, isoft;
}consts_;  

/* 3  */
extern struct{
  double rsval[NPDM];
}valden_;

/* 3  */
extern struct{
  int nll;
}nll_;  

/* 3  */
extern struct{
  int     nnr,numorb;
  double  rr[NPDM],aar[N0][NPDM],bbr[N0][NPDM],ccdc[NPDM],eev[N0];
  int     llo[N0];
  double  sso[N0];
  int     nno[N0];
}atmwave_;

/* 3  */
extern struct{
  double zeff, rpcc;
  double spacing, first;
  int    nval, igrid, inl, ist1, ncut, ilocal,ilocalind;
}pwf_;      

/* 2  */
extern struct{
  int meth;
}opt_;  

/* 2  */
extern struct{
  double rcrel[30],relnorm[30];
}rcrel_;

/* 2  */
extern struct{
  int npots;
}npp_;  

/* 2  */
extern struct{
  double wfu[NVALE0][NPDM],wfd[NVALE0][NPDM];
}relwf_;  

/* 2  */
extern struct{
  double qcl[10];
  int    nbl[10];
  /*  int    ifnll[10];
  int    numenl[10];
  double wnorml[10],wkel[10];
  double escatl[10][NUMEN0],rscatl[10][NUMEN0],wscatl[10][NUMEN0];*/
}lparam_;  

/* 2  */
extern struct{
  double enc, tolc;
  int    lim,nstp;
  int    isw1,isw2,isw3;
}scrconmax_;

/* 2  */
extern struct{
  double wsumt[NVALE0],sumc[NVALE0];
  int lghost[NVALE0];
}convrpt_;  

/* 1  */
/*extern struct{
 double rxccut;
 }rxccut_;  */

/* 1  */
extern struct{
  int ipccmeth;
}ipccmeth_;  

/* 1  */
extern struct{
  int nlmp[NVALE0]; 
  double wnlp[NVALE0],ev[NVALE0*2],rnor[NVALE0*2];
}frlwf_;  

/* 1  */
extern struct{
  double ensave[N0];
}ensave_;  





#endif
/* COMMON_BLOCKS_H */
