#ifndef FSDE_COMMON_DEFINES_H
#define FSDE_COMMON_DEFINES_H

/* Sequence of parameters is:
 * epsilon
 * mu
 * Nrig
 * Lambda
 */

/* common variable defintion/mappings for auto-generated C code */

/* potential parameters */
double epsilon      = params[ 0];
double mu   		= params[ 1];
double Nrig	  	 	= params[ 2];
double lambda		= params[ 3];

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

#endif /* FSDE_COMMON_DEFINES_H */