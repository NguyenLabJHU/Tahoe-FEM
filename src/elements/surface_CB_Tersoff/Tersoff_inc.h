/* $Id: Tersoff_inc.h,v 1.1 2007-06-17 03:27:30 paklein Exp $ */
#ifndef TERSOFF_INC_H
#define TERSOFF_INC_H

#ifdef __cplusplus
extern "C" {
#endif

/* function to compute derivatives of the potential function wrt to the
 * internal degrees of freedom 
 *
 * Sequence of parameters is:
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
void get_dXsi(double* params, double *Xsi, double *Xa, double *Ya, double *Za, double* Cmat, double* dXsi, double* ddXsi); 

#ifdef __cplusplus
}
#endif

#endif /* TERSOFF_INC_H */
