/* $Id: Wurtzite_inc.h,v 1.2 2007-11-09 21:31:36 hspark Exp $ */
#ifndef WURTZITE_INC_H
#define WURTZITE_INC_H

#ifdef __cplusplus
extern "C" {
#endif

/* Sequence of parameters is:
 * C-axis
 * A-axis
 * Mass
 * D0
 * S0
 * r0
 * beta
 * gamma
 * c
 * d
 * h
 * R
 * D
 */

/* function to compute derivatives of the potential function wrt to the
 * internal degrees of freedom */
void WZget_dXsi(const double* params, const double *Xsi, const double *Xa, const double *Ya, const double *Za, const double* Cmat, double* dXsi, double* ddXsi); 

/* function to compute derivatives of the potential function wrt to the stretch tensor */
void WZget_dUdC(const double* params, const double *Xsi, const double *Xa, const double *Ya, const double *Za, const double* Cmat, double* dUdC); 

/* function to compute all second order derivatives of the potential function needed for the modulus */
void WZget_ddC(const double* params, const double *Xsi, const double *Xa, const double *Ya, const double *Za, const double* Cmat, double* dCdC, double* dCdXsi);

#ifdef __cplusplus
}
#endif

#endif /* WURTZITE_INC_H */
