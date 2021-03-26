/*
 */

/****************************************************************************
* Part of uniPP Library providing access to universal pseudopotential files.*
*                                                                           *
* if flag==0:                                                               *
*                                                                           *
* - uniPP_kbfnl() calculates the planewave projection factors onto the      *
*   pseudo wavefunctions. Returned units: [sqrt((a.u.)^3)]                  *
*                                                                           *
* if flag==1:                                                               *
*                                                                           *
* - uniPP_kbfnl() calculates the KB factors of the non-local part of the    *
*   pseudopotential. Returned units: [Hartree * sqrt((a.u.)^3)]             *
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
#include <math.h>

#include "xmcomplex.h"
#include "xmintegral.h"

#include "uniPPlib.h"


#ifdef HAVE_SIMPSON

void uniPP_kbfnl(uniPP *unipp, dcomplex *****kbfnl, int nproj, 
                 double *gk_abs, double **gk, int nmax, int flag){
  
  int i, j, j_max, ig, l, l_loc, proj;
  double r, u, dv, q, z, gka;
  double cosx, cosy, cosz;
  double eta, lambda_eta;
  double *g = (double *)malloc(unipp->m_mesh * sizeof(double));
  double *f = (double *)malloc(unipp->m_mesh * sizeof(double));
  
  l_loc = unipp->l_loc;
    
  for (l=0; l<unipp->l_max; l++){
    if ((l!=l_loc) | (flag==0)){
            
      if (unipp->rel && l>0) 
        j_max=2;
      else 
        j_max=1;
            
      for (j=0; j<j_max; j++){
        for (proj=0; proj<nproj; proj++){
                
          /******************************************************************
          * Prepare the two arrays: f[] and g[]                             *
          ******************************************************************/
 
          for (i=0; i<unipp->m_mesh; i++){
            r = unipp->r_m[i];
            
            u = unipp->u_ps[l][j][i];
            
            if (flag==1)
              dv = unipp->v_ps[l][j][i] - unipp->v_loc[i]; /* [Hartree] */
            else  
              dv = 1.;                                     /* [unitless] */
 
            
            f[i] = r*u*u*dv;
            g[i] = r*r*u*dv;
             
          }
           
          /********************************************************************
          * Handle different angular momenta:                                 *
          ********************************************************************/
 
          if (l==0){
 
            /******************************************************************
            * s-part                                                          *
            ******************************************************************/
 
            if (proj==0){
 
              for (ig=0; ig<nmax; ig++){
 
                gka = gk_abs[ig] * 0.529177249;  /* 1/A -> 1/au */

                if (gka==0.)
                  q = simpson(unipp->m_mesh, g, log(unipp->a_mesh)) + .5*g[0];
                else{
                  for (i=0; i<unipp->m_mesh; i++){
                    r = unipp->r_m[i];
                    z = gka*r;
                    f[i] = (sin(z)/z) * g[i];
                  }
                  q = simpson(unipp->m_mesh, f, log(unipp->a_mesh)) + .5*f[0];
 
 
                }
                kbfnl[0][j][0][0][ig].r = q;
                kbfnl[0][j][0][0][ig].i = 0.;
              }
 
            }else if (proj==1){
 
              eta = simpson(unipp->m_mesh, f, log(unipp->a_mesh)) + .5*f[0];
 
              for (i=0; i<unipp->m_mesh; i++){
                r = unipp->r_m[i];
                g[i] *= r*r;
                f[i] *= r*r;
              }
 
              lambda_eta = eta/
                (simpson(unipp->m_mesh, f, log(unipp->a_mesh)) + .5*f[0]);
 
              for (ig=0; ig<nmax; ig++){
 
                gka = gk_abs[ig] * 0.529177249;  /* 1/A -> 1/au */

                if (gka==0.)
                  q = simpson(unipp->m_mesh, g, log(unipp->a_mesh)) + .5*g[0];
                else{
                  for (i=0; i<unipp->m_mesh; i++){
                    r = unipp->r_m[i];
                    z = gka*r;
                    f[i] = (sin(z)/z) * g[i];
                  }
                  q = simpson(unipp->m_mesh, f, log(unipp->a_mesh)) + .5*f[0];
 
 
                }
                kbfnl[0][j][1][0][ig].r = q -lambda_eta*kbfnl[0][j][0][0][ig].r;
                kbfnl[0][j][1][0][ig].i = 0.;
              }
 
            }else{
              printf("uniPP_kbfnl: requested projector of order"
                     " higher than two -> exit(0)\n");
              exit(0);
            }
 
 
          }else if (l==1){
 
            /******************************************************************
            * p-part                                                          *
            ******************************************************************/
 
            if (proj==0){
 
              for (ig=0; ig<nmax; ig++){
 
                gka = gk_abs[ig] * 0.529177249;  /* 1/A -> 1/au */
 
                if (gka==0.){
                  kbfnl[1][j][0][0][ig].r = kbfnl[1][j][0][0][ig].i = 0.;
                  kbfnl[1][j][0][1][ig].r = kbfnl[1][j][0][1][ig].i = 0.;
                  kbfnl[1][j][0][2][ig].r = kbfnl[1][j][0][2][ig].i = 0.;
                }else{
                  for (i=0; i<unipp->m_mesh; i++){
                    r = unipp->r_m[i];
                    z = gka*r;
                    f[i] = ((sin(z)/z-cos(z))/z) * g[i];
                  }
                  q = simpson(unipp->m_mesh, f, log(unipp->a_mesh)) + .5*f[0];
 
 
                  cosx  = gk[ig][0]/gk_abs[ig];
                  cosy  = gk[ig][1]/gk_abs[ig];
                  cosz  = gk[ig][2]/gk_abs[ig];

                  kbfnl[1][j][0][0][ig].r =   q * cosx / sqrt(2.);
                  kbfnl[1][j][0][0][ig].i = - q * cosy / sqrt(2.);
                  kbfnl[1][j][0][1][ig].r =   q * cosz;
                  kbfnl[1][j][0][1][ig].i =   0.;
                  kbfnl[1][j][0][2][ig].r = - q * cosx / sqrt(2.);
                  kbfnl[1][j][0][2][ig].i = - q * cosy / sqrt(2.);
                  
                }
              }
 
            }else if (proj==1){
 
              eta = simpson(unipp->m_mesh, f, log(unipp->a_mesh)) + .5*f[0];
 
              /* printf("eta-p = %g\n", eta); */
 
              for (i=0; i<unipp->m_mesh; i++){
                r = unipp->r_m[i];
                g[i] *= r*r;
                f[i] *= r*r;
              }
              
 
              lambda_eta =
                (simpson(unipp->m_mesh, f, log(unipp->a_mesh)) + .5*f[0])/eta;
 
              for (ig=0; ig<nmax; ig++){
 
                gka = gk_abs[ig] * 0.529177249;  /* 1/A -> 1/au */

                if (gka==0.){
                  kbfnl[1][j][1][0][ig].r = kbfnl[1][j][1][0][ig].i = 0.;
                  kbfnl[1][j][1][1][ig].r = kbfnl[1][j][1][1][ig].i = 0.;
                  kbfnl[1][j][1][2][ig].r = kbfnl[1][j][1][2][ig].i = 0.;
                }else{
 
                  for (i=0; i<unipp->m_mesh; i++){
                    r = unipp->r_m[i];
                    z = gka*r;
                    f[i] = ((sin(z)/z-cos(z))/z) * g[i];
                  }
 
                  q = simpson(unipp->m_mesh, f, log(unipp->a_mesh)) + .5*f[0];
 
                  cosx  = gk[ig][0]/gk_abs[ig];
                  cosy  = gk[ig][1]/gk_abs[ig];
                  cosz  = gk[ig][2]/gk_abs[ig];
 
                  kbfnl[1][j][1][0][ig].r =   q * cosx / sqrt(2.)
                                         - lambda_eta * kbfnl[1][j][0][0][ig].r;
                  kbfnl[1][j][1][0][ig].i = - q * cosy / sqrt(2.)
                                         - lambda_eta * kbfnl[1][j][0][0][ig].i;
                  kbfnl[1][j][1][1][ig].r =   q * cosz
                                         - lambda_eta * kbfnl[1][j][0][1][ig].r;
                  kbfnl[1][j][1][1][ig].i =   0.;
                  kbfnl[1][j][1][2][ig].r = - q * cosx / sqrt(2.)
                                         - lambda_eta * kbfnl[1][j][0][2][ig].r;
                  kbfnl[1][j][1][2][ig].i = - q * cosy / sqrt(2.)
                                         - lambda_eta * kbfnl[1][j][0][2][ig].i;
                }
              }
 
            }else{
              printf("uniPP_kbfnl: requested projector of order"
                     " higher than two -> exit(0)\n");
              exit(0);
            }
 
 
          }else if (l==2){
 
            /******************************************************************
            * d-part                                                          *
            ******************************************************************/
 
            if (proj==0){
 
              for (ig=0; ig<nmax; ig++){
 
                gka = gk_abs[ig] * 0.529177249;  /* 1/A -> 1/au */

                if (gka==0.){
                  kbfnl[2][j][0][0][ig].r = kbfnl[2][j][0][0][ig].i = 0.;
                  kbfnl[2][j][0][1][ig].r = kbfnl[2][j][0][1][ig].i = 0.;
                  kbfnl[2][j][0][2][ig].r = kbfnl[2][j][0][2][ig].i = 0.;
                  kbfnl[2][j][0][3][ig].r = kbfnl[2][j][0][3][ig].i = 0.;
                  kbfnl[2][j][0][4][ig].r = kbfnl[2][j][0][4][ig].i = 0.;
                }else{
                  for (i=0; i<unipp->m_mesh; i++){
                    r = unipp->r_m[i];
                    z = gka*r;
                    f[i] = (((3./(z*z)-1)*sin(z) - (3./z)*cos(z))/z) * g[i];
                  }
                  q = simpson(unipp->m_mesh, f, log(unipp->a_mesh)) + .5*f[0];
  
                  cosx  = gk[ig][0]/gk_abs[ig];
                  cosy  = gk[ig][1]/gk_abs[ig];
                  cosz  = gk[ig][2]/gk_abs[ig];
 
                  kbfnl[2][j][0][0][ig].r = q*(cosx*cosx-cosy*cosy)*sqrt(3./8.);
                  kbfnl[2][j][0][0][ig].i = -q * (2.*cosx*cosy) * sqrt(3./8.);
                  kbfnl[2][j][0][1][ig].r =  q * (cosx*cosz) * sqrt(3./2.);
                  kbfnl[2][j][0][1][ig].i = -q * (cosy*cosz) * sqrt(3./2.);
                  kbfnl[2][j][0][2][ig].r =  q * (3.*cosz*cosz-1.) * .5;
                  kbfnl[2][j][0][2][ig].i =  0.;
                  kbfnl[2][j][0][3][ig].r = -q * (cosx*cosz) * sqrt(3./2.);
                  kbfnl[2][j][0][3][ig].i = -q * (cosy*cosz) * sqrt(3./2.);
                  kbfnl[2][j][0][4][ig].r = q*(cosx*cosx-cosy*cosy)*sqrt(3./8.);
                  kbfnl[2][j][0][4][ig].i =  q * (2.*cosx*cosy) * sqrt(3./8.);
 
                }
              }
 
            }else if (proj==1){
 
              eta = simpson(unipp->m_mesh, f, log(unipp->a_mesh)) + .5*f[0];
 
              for (i=0; i<unipp->m_mesh; i++){
                r = unipp->r_m[i];
                g[i] *= r*r;
                f[i] *= r*r;
              }
 
              lambda_eta =
                (simpson(unipp->m_mesh, f, log(unipp->a_mesh)) + .5*f[0])/eta;
 
              for (ig=0; ig<nmax; ig++){
 
                gka = gk_abs[ig] * 0.529177249;  /* 1/A -> 1/au */

                if (gka==0.){
                  kbfnl[2][j][1][0][ig].r = kbfnl[2][j][1][0][ig].i = 0.;
                  kbfnl[2][j][1][1][ig].r = kbfnl[2][j][1][1][ig].i = 0.;
                  kbfnl[2][j][1][2][ig].r = kbfnl[2][j][1][2][ig].i = 0.;
                  kbfnl[2][j][1][3][ig].r = kbfnl[2][j][1][3][ig].i = 0.;
                  kbfnl[2][j][1][4][ig].r = kbfnl[2][j][1][4][ig].i = 0.;
                }else{
                  for (i=0; i<unipp->m_mesh; i++){
                    r = unipp->r_m[i];
                    z = gka*r;
                    f[i] = (((3./(z*z)-1)*sin(z) - (3./z)*cos(z))/z) * g[i];
                  }
                  q = simpson(unipp->m_mesh, f, log(unipp->a_mesh)) + .5*f[0];
  
                  cosx  = gk[ig][0]/gk_abs[ig];
                  cosy  = gk[ig][1]/gk_abs[ig];
                  cosz  = gk[ig][2]/gk_abs[ig];
 
 
                  kbfnl[2][j][1][0][ig].r =q*(cosx*cosx-cosy*cosy) * sqrt(3./8.)
                                         - lambda_eta * kbfnl[2][j][0][0][ig].r;
                  kbfnl[2][j][1][0][ig].i = -q * (2.*cosx*cosy) * sqrt(3./8.)
                                         - lambda_eta * kbfnl[2][j][0][0][ig].i;
                  kbfnl[2][j][1][1][ig].r =  q * (cosx*cosz) * sqrt(3./2.)
                                         - lambda_eta * kbfnl[2][j][0][1][ig].r;
                  kbfnl[2][j][1][1][ig].i = -q * (cosy*cosz) * sqrt(3./2.)
                                         - lambda_eta * kbfnl[2][j][0][1][ig].i;
                  kbfnl[2][j][1][2][ig].r =  q * (3.*cosz*cosz-1.) * .5
                                         - lambda_eta * kbfnl[2][j][0][2][ig].r;
                  kbfnl[2][j][1][2][ig].i =  0.;
                  kbfnl[2][j][1][3][ig].r = -q * (cosx*cosz) * sqrt(3./2.)
                                         - lambda_eta * kbfnl[2][j][0][3][ig].r;
                  kbfnl[2][j][1][3][ig].i = -q * (cosy*cosz) * sqrt(3./2.)
                                         - lambda_eta * kbfnl[2][j][0][3][ig].i;
                  kbfnl[2][j][1][4][ig].r =q*(cosx*cosx-cosy*cosy) * sqrt(3./8.)
                                         - lambda_eta * kbfnl[2][j][0][4][ig].r;
                  kbfnl[2][j][1][4][ig].i =  q * (2.*cosx*cosy) * sqrt(3./8.)
                                         - lambda_eta * kbfnl[2][j][0][4][ig].i;
                }
              }
 
            }else{
              printf("uniPP_kbfnl: requested projector of order"
                     " higher than two -> exit(0)\n");
              exit(0);
            }
 
 
 
 
          }else{
            printf("uniPP_kbfnl.c: l=%d > 2 is not supported!\n", l);
            exit(0);
          }
        }
      }
    }    
  }
  free(f);
  free(g);
 
}

#endif /* HAVE_SIMPSON */


