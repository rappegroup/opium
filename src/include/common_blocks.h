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
/***************************************************************************
* these are the common blocks from the fortran code                        *
****************************************************************************/
#ifndef COMMON_BLOCKS_H
#define COMMON_BLOCKS_H

extern struct{
  char file_log[80];
}filenames_;

extern struct{
  double etol, vtol;
  int maxit, isoft;
}consts_;  

extern struct{
  double h, r1, z, r[NPDM];
  int np;
}grid_;  

extern struct{
  double rnl[N0][NPDM];
}wfn_;  

extern struct{
  double rnla[N0][NPDM],rnlb[N0][NPDM];
}wfnrel_;  

extern struct{
  double rscore[NPDM],rdd[NPDM],rddd[NPDM],rscoretot[NPDM],rsval[NPDM];
}rscore_;  

extern struct{
  double rvcore[N0][NPDM],rvps[N0][NPDM],rvcoul[NPDM];
}totpot_;

extern struct{
  int ncore,nval,norb,nlm[N0],no[N0],lo[N0],nmax[N0],maxim;
}aorb_;  

extern struct{
  double wnl[N0],en[N0],so[N0],xion;
}adat_;  

extern struct{
  double rcall[N0],rvap[N0],rnorm[N0];
  int ibd[N0];
  double etot;
}aval_;

extern struct{
  double etrial[N0];
}etrial_;  

extern struct{
  int nll,meth;
}psdat_;

extern struct{
  double qcl[N0];
  int    nbl[N0];
}optparam_;  

extern struct{
  double enc, tolc;
  int    lim,nstp;
  int    isw1,isw2,isw3;
}scrconmax_;

extern struct{
  double wsumt[N0],sumc[N0];
  int npsghost[N0];
}psout_;  


/* DNL vars */
extern struct{
  double rvloc[NPDM];
}nlcore_;      

extern struct{
  int    iboxstart[N0],iboxend[N0];
  double boxheight[N0];
  int    numbox;
}box_;  

extern struct{
  int inl,indrc[N0],IB,ID,IG;
}nlpot2_;  

extern struct{
  int nlghost[N0];
  int iloc, idesign;
}local_;  


/* pcc */
extern struct{
  double rpcc,rpccz;
}rpcc_;  

extern struct{
  int ipccmeth;
}ipccmeth_;  

extern struct{
  int iavgtype;
}iavgtype_;

/* logder */
extern struct{
  double rphas,elogmax,elogmin,dlwf[N0][NPL0];
}logarith_;  

extern struct{
  int  ilogder;
}ilogder_;  

extern struct{
  double zeff, rpcc;
  double spacing, first;
  int    nval, igrid, inl, ist1, ncut, ilocal,ilocalind;
}pwf_;      

extern struct{
  double a,b,r[NPDM],rab[NPDM];
}rgrid_;

extern struct{
  int nr;
}nrgrid_;




#endif
/* COMMON_BLOCKS_H */


