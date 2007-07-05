/* $Id: Tersoff_common_defines_surf.h,v 1.1 2007-07-05 20:37:58 hspark Exp $ */
#ifndef TERSOFF_COMMON_DEFINES_SURF_H
#define TERSOFF_COMMON_DEFINES_SURF_H

/* Sequence of parameters is:
 * A
 * B
 * Mass
 * lambda
 * mu
 * beta
 * n
 * c
 * d
 * h
 * chi
 * R
 * S
 */

/* common variable defintion/mappings for auto-generated C code */

/* potential parameters */
double A      = params[ 0];
double B      = params[ 1];
double Mass   = params[ 2];
double lam    = params[ 3];
double mu     = params[ 4];
double beta   = params[ 5];
double n      = params[ 6];
double c      = params[ 7];
double d      = params[ 8];
double h      = params[ 9];
double chi    = params[10];
double R      = params[11];
double S      = params[12];

/* degrees of freedom */
/* ARE THESE ACTUALLY IN XS1,XS2,YS1,YS2,ZS1,ZS2 ORDER? */
double Xs1 = Xsi[0];
double Ys1 = Xsi[1];
double Zs1 = Xsi[2];
double Xs2 = Xsi[3];
double Ys2 = Xsi[4];
double Zs2 = Xsi[5];

/* atom coordinates (reference) */
double X1 = Xa[0];	
double X2 = Xa[1];	
double X3 = Xa[2];	
double X4 = Xa[3];	
double X5 = Xa[4];	
double X6 = Xa[5];	
double X7 = Xa[6];	
double X8 = Xa[7];	
double X9 = Xa[8];	

double Y1 = Ya[0];	
double Y2 = Ya[1];	
double Y3 = Ya[2];	
double Y4 = Ya[3];	
double Y5 = Ya[4];	
double Y6 = Ya[5];	
double Y7 = Ya[6];	
double Y8 = Ya[7];	
double Y9 = Ya[8];	

double Z1 = Za[0];	
double Z2 = Za[1];	
double Z3 = Za[2];	
double Z4 = Za[3];	
double Z5 = Za[4];	
double Z6 = Za[5];	
double Z7 = Za[6];	
double Z8 = Za[7];	
double Z9 = Za[8];	

/* deformation */
double C11 = Cmat[0];
double C12 = Cmat[1];
double C13 = Cmat[2];
double C21 = Cmat[3];
double C22 = Cmat[4];
double C23 = Cmat[5];
double C31 = Cmat[6];
double C32 = Cmat[7];
double C33 = Cmat[8];

#endif /* TERSOFF_COMMON_DEFINES_SURF_H */