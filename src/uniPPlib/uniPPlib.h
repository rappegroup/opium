/* uniPPlib.h */

/****************************************************************************
* uniPP library provides access to universal pseudopotential files.         *
*                                                                           *
* gjt                                                                       *
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

#ifndef XMCOMPLEX
#include "xmcomplex.h"
#endif

/****************************************************************************
* This is the uniPP data type definition:                                   *
****************************************************************************/

typedef struct uniPP{

  char name[80];    /* name of the pseudopotential */

  double z_ion;     /* number of valence electrons */
  int l_max;        /* number of pseudopotential components */
  int l_loc;        /* index of local pseudopotential component */
  
  int rel;          /* flag: relativistic potentials */
  int nlcc;         /* flag: non-linear core correction */
    
  int m_mesh;       /* number of radial mesh points */
  double a_mesh;    /* mesh increment r_(m+1)/r_m */
  
  double *r_m;      /* [m_mesh] radial coordinate [a.u.] */
  
  double *v_loc;    /* [m_mesh] radial local pseudopotential [Hartree] */
  
  double ***v_ps;   /* [l_max][j][m_mesh] radial pseudopotential [Hartree]*/
  double ***u_ps;   /* [l_max][j][m_mesh] rad.pseudo wavefunction * r_m 
                                          normalized --> [1./sqrt(a.u.)]*/
  
  double *n_pc;     /* [m_mesh] partial core density [1./(a.u.)^3] */

  /* EJW added these */ 
  double *n_pc1;     /* [m_mesh] partial core density [1./(a.u.)^3] deriv*/
  double *n_pc2;     /* [m_mesh] partial core density [1./(a.u.)^3] 2nd deriv*/
  
}uniPP;

/****************************************************************************
* Method prototypes follow:                                                 *
****************************************************************************/

int uniPP_read(uniPP *unipp, FILE *fp);
int uniPP_write(uniPP *unipp, FILE *fp);

int uniPP_writefhi(uniPP *unipp, FILE *fp);

void uniPP_free(uniPP *unipp);                                                        
double uniPP_alpha(uniPP *unipp);

void uniPP_ffloc(uniPP *unipp, double *ffloc, double *g, int n);

void uniPP_ffnlcc(uniPP *unipp, double *ffnlcc, double *g, int n);
  
void uniPP_kbden(uniPP *unipp, double ***kbden, double nproj);

void uniPP_kbfnl(uniPP *unipp, dcomplex *****kbfnl, int nproj, 
                 double *gk_abs, double **gk, int nmax, int flag);
