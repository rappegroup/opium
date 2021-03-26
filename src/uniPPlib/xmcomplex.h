/* xmcomplex.h */
/****************************************************************************
* Part of the eXtended Math Library (Xmath):                                *
*    - Complex numbers and arithmetic                                       *
*                                                                           *
* References:                                                               *
*    (1) "numerical recipes in c", ISBN 0-521-43108-5                       *
*  =======================================================================  *
*                                                                           *
*  last touched: 19.07.00 gjt                                               *
****************************************************************************/

/****************************************************************************
*  'Xmath' - The eXtended Math Library                                      *
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

#define XMCOMPLEX

/******************** complex numbers ****************************************/

typedef struct {float r,i;} fcomplex;    /* float complex as in (1) */
typedef struct {double r,i;} dcomplex;   /* double complex */


/******************** complex arithmetic *************************************/

/* MACROS have advantages:
   - up to 3 times faster than functions
   - accept any type consistent set of variables
   and disadvantages:
   - cannot be used within another expression (except for CABS2)
   - expressions on the "righthand side" of the macro will be evaluated 
     multiple times (in other words: not a good idea to gain speed!) */

/* c = a + ib */
#define COMPLEX(c, a, b) {c.r = a; c.i = b;} 

/* c = a + b */
#define CADD(c, a, b) {c.r = a.r + b.r; c.i = a.i + b.i;}

/* c = a - b */
#define CSUB(c, a, b) {c.r = a.r - b.r; c.i = a.i - b.i;}

/* c = a * b */
#define CMUL(c, a, b) {c.r=a.r*b.r-a.i*b.i; c.i=a.i*b.r+a.r*b.i;}

/* c = z* */
#define CONJG(c, z) {c.r = z.r; c.i = -z.i;}

/* (z * z*) */
#define CABS2(z) (z.r * z.r + z.i * z.i)

/* c = a(real) * b(complex) */
#define RCMUL(c, a, b) {c.r=a*b.r; c.i=a*b.i;}


/* float complex functions */

/* r + j*i */
fcomplex fComplex(float r, float i);

/* a + b */
fcomplex fCadd(fcomplex a, fcomplex b);

/* a - b */
fcomplex fCsub(fcomplex a, fcomplex b);

/* a * b */
fcomplex fCmul(fcomplex a, fcomplex b);

/* a/b */
fcomplex fCdiv(fcomplex a, fcomplex b);

/* z* */
fcomplex fConjg(fcomplex z);

/* sqrt(a * a*) */
float fCabs(fcomplex z);

/* a * a* */
float fCabs2(fcomplex z);

/* sqrt(z) */
fcomplex fCsqrt(fcomplex z);

/* x(real) * a(complex) */
fcomplex fRCmul(float x, fcomplex a);


/* double complex functions */

/* r + j*i */
dcomplex dComplex(double r, double i);

/* a + b */
dcomplex dCadd(dcomplex a, dcomplex b);

/* a - b */
dcomplex dCsub(dcomplex a, dcomplex b);

/* a * b */
dcomplex dCmul(dcomplex a, dcomplex b);

/* Re{a * b} */
double dCmul_r(dcomplex a, dcomplex b);

/* a/b */
dcomplex dCdiv(dcomplex a, dcomplex b);

/* 1/b */
dcomplex dC1div(dcomplex b);

/* z* */
dcomplex dConjg(dcomplex z);

/* sqrt(a * a*) */
double dCabs(dcomplex z);

/* a * a* */
double dCabs2(dcomplex z);

/* sqrt(z) */
dcomplex dCsqrt(dcomplex z);

/* -z */
dcomplex dCneg(dcomplex z);

/* exp(ix) */
dcomplex dCexp(double x);

/* x(real) * a(complex) */
dcomplex dRCmul(double x, dcomplex a);

#endif
