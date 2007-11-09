/* $Id: TersoffDimer_inc.h,v 1.1 2007-11-09 16:53:55 hspark Exp $ */
#ifndef TERSOFFDIMER_INC_H
#define TERSOFFDIMER_INC_H

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
void get_dXsi(const double* params, const double *Xsi, const double *Xa, const double *Ya, const double *Za, const double* Cmat, double* dXsi, double* ddXsi); 

/* function to compute derivatives of the potential function wrt to the stretch tensor */
void get_dUdC(const double* params, const double *Xsi, const double *Xa, const double *Ya, const double *Za, const double* Cmat, double* dUdC); 

/* function to compute all second order derivatives of the potential function needed for the modulus */
void get_ddC(const double* params, const double *Xsi, const double *Xa, const double *Ya, const double *Za, const double* Cmat, double* dCdC, double* dCdXsi);

#ifdef __cplusplus
}
#endif

#endif /* TERSOFFDIMER_INC_H */
