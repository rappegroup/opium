/*
 * Copyright (c) 1998-2004 The OPIUM Group
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

/* thfs COMMON blocks */

extern struct{
  double etol, vtol;
  int maxit, isoft;
}consts_;  

extern struct{
  double rxccut;
}rxccut_;  

extern struct{
  double rsval[NPDM];
}valden_;

extern struct{
  int meth;
}opt_;  

extern struct{
  int ibd[N0];
}ibound_;  

extern struct{
  double h, r1, z, r[NPDM];
  int np;
}grid_;  

extern struct{
  int igga;
}igga_;  

extern struct{
  int  ilogder;
}ilogder_;  

extern struct{
  int    iboxstart[10],iboxend[10];
  double boxheight[10];
  int    numbox;
}box_;  

extern struct{
  int iloc, idesign;
}local_;  

extern struct{
  int ipccmeth;
}ipccmeth_;  

extern struct{
  double rphas,elogmax,elogmin,dlwf[N0][NPL0];
}logarith_;  

extern struct{
  int inl,indrc[N0],IB,ID,IG;
}nlpot2_;  

extern struct{
  int ncores, nvales, iskip, maxip;
}npm_;  

extern struct{
  int ioutput, irelvirt;
}output_;  

extern struct{
  double rcall[N0];
}rcall_;

extern struct{
  double rcrel[30],relnorm[30];
}rcrel_;

extern struct{
  double rpcc,rpccz;
}rpcc_;  

extern struct{
  double rnorm[30],etot;
  int lghost[30];
}results_;  

extern struct{
  int nll;
}nll_;  

extern struct{
  double rscore[NPDM],rdd[NPDM],rddd[NPDM],rscoretot[NPDM];
}rscore_;  

extern struct{
  double rvloc[NPDM];
}nlcore_;      

/* atm COMMON blocks */

extern struct{
  int  norb, ncore, nval;
  int  no[40], lo[40];
}reli_;  

extern struct{
  double zcore,zval;
  double zo[40],so[40];
}reld_;

extern struct{
  double a,b,r[NRMAX],rab[NRMAX];
}rgrid_;

extern struct{
  int nr;
}nrgrid_;

extern struct{
  int     nnr,numorb;
  double  rr[NRMAX],aar[30][NRMAX],bbr[30][NRMAX],ccdc[NRMAX],eev[NORBP];
  int     llo[NORBP];
  double  sso[NORBP];
  int     nno[NORBP];
}atmwave_;

/* new COMMON blocks */

extern struct{
  int iwriteae;
}flags_;

extern struct{
  char file_log[80];
}filenames_;

extern struct{
  char filev[80];
}temp_;

extern struct{
  int iwritepcc;
}iwrpcc_;

/* interp COMMON blocks */

extern struct{
  double title,c0,rlam,rcutof,rphas,emin,emax;
  int    nnt,igo,idate,icore,ispin,ipratt,ifprt,iphas,iuserc;
}atom1_;  
          
extern struct{
  int ncores, nvales;
}np_;  

extern struct{
  int npots;
}npp_;  

extern struct{
  double wfu[NVALE0][NPDM],wfd[NVALE0][NPDM];
}relwf_;  

extern struct{
  int nlmp[NVALE0]; 
  double wnlp[NVALE0],ev[2*NVALE0],rnor[2*NVALE0];
}frlwf_;  


/* pseudo COMMON blocks */

extern struct{
  double xion,rnl[N0][NPDM];
  int nlm[N0];
  double wnl[N0],en[N0];
  int norb;
}atomic_;  

extern struct{
  double ensave[N0];
}ensave_;  

extern struct{
  double qcl[10];
  int    nbl[10],ifnll[10];
  int    numenl[10];
  double wnorml[10],wkel[10];
  double escatl[10][NUMEN0],rscatl[10][NUMEN0],wscatl[10][NUMEN0];
}lparam_;  

extern struct{
  double rvps[NVALE0+1][NPDM],rvcore[NVALE0+1][NPDM];
}rvps_;  

extern struct{
  double enc, tolc;
  int    lim,nstp;
  int    isw1,isw2,isw3;
}scrconmax_;

extern struct{
  double wsumt[NVALE0],sumc[NVALE0];
  int lghost[NVALE0];
}convrpt_;  
          
/* pwf COMMON blocks */

extern struct{
  double zeff, rpcc;
  double spacing, first;
  int    nval, igrid, inl, ist1, ncut, ilocal,ilocalind;
}pwf_;      


/* pot COMMONS */

extern struct{
  double rvcore[N0][NPDM],rvps[N0][NPDM],rvcoul[NPDM];
}totpot_;

extern struct{
  int nmax[NVALE0];
  int maxim;
}nmax_;



#endif
/* COMMON_BLOCKS_H */
