/* $Id: Tersoff_inc.h,v 1.2 2007-06-17 04:02:34 paklein Exp $ */
#ifndef TERSOFF_INC_H
#define TERSOFF_INC_H

#ifdef __cplusplus
extern "C" {
#endif

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

/* function to compute derivatives of the potential function wrt to the
 * internal degrees of freedom */
void get_dXsi(double* params, double *Xsi, double *Xa, double *Ya, double *Za, double* Cmat, double* dXsi, double* ddXsi); 

/* function to compute derivatives of the potential function wrt to the
 * internal degrees of freedom */
void get_dUdC(double* params, double *Xsi, double *Xa, double *Ya, double *Za, double* Cmat, double* dUdC); 

#ifdef __cplusplus
}
#endif

#endif /* TERSOFF_INC_H */
