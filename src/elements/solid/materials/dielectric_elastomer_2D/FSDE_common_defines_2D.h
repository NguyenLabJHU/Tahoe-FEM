#ifndef FSDE_COMMON_DEFINES_2D_H
#define FSDE_COMMON_DEFINES_2D_H

/* Sequence of parameters is:
 * epsilon
 * mu
 * Nrig
 * lambda
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

/* deformation */
double C11 = Cmat[0];
double C12 = Cmat[1];
double C21 = Cmat[2];
double C22 = Cmat[3];

#endif /* FSDE_COMMON_DEFINES_2D_H */