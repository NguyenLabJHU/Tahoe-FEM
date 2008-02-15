/* $Id: Tersoff_inc_surf.h,v 1.2 2008-02-15 05:15:49 hspark Exp $ */
#ifndef TERSOFF_INC_SURF_H
#define TERSOFF_INC_SURF_H

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
void get_dXsi_surf(const double* params, const double *Xsi, const double *Xa, const double *Ya, const double *Za, const double* Cmat, double* dXsi, double* ddXsi); 

/* function to compute derivatives of the potential function wrt to the stretch tensor */
void get_dUdC_surf(const double* params, const double *Xsi, const double *Xa, const double *Ya, const double *Za, const double* Cmat, double* dUdC); 

/* function to compute all second order derivatives of the potential function needed for the modulus */
void get_ddC_surf(const double* params, const double *Xsi, const double *Xa, const double *Ya, const double *Za, const double* Cmat, double* dCdC, double* dCdXsi);

/* function to compute the surface strain energy density */
double get_energy_surf(const double* params, const double *Xsi, const double *Xa, const double *Ya, const double *Za, const double* Cmat); 

#ifdef __cplusplus
}
#endif

#endif /* TERSOFF_INC_H */
