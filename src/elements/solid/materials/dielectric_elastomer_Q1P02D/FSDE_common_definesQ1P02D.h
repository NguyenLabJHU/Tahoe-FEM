#ifndef FSDE_COMMON_DEFINESQ1P02D_H
#define FSDE_COMMON_DEFINESQ1P02D_H

/* Sequence of parameters is:
 * mu
 * lambda
 * epsilon
 * Nrig
 */

/* common variable defintion/mappings for auto-generated C code */

/* potential parameters */
double mu   		= params[ 0];
double lambda		= params[ 1];
double epsilon      = params[ 2];
double Nrig	  	 	= params[ 3];

/* E-field */
double ex = Xsi[0];
double ey = Xsi[1];
double ez = Xsi[2];

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

/* more deformation */
double F11 = Fmat[0];
double F12 = Fmat[1];
double F13 = Fmat[2];
double F21 = Fmat[3];
double F22 = Fmat[4];
double F23 = Fmat[5];
double F31 = Fmat[6];
double F32 = Fmat[7];
double F33 = Fmat[8];

#endif /* FSDE_COMMON_DEFINESQ1P02D_H */