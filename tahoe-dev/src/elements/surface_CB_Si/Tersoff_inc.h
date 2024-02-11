/* $Id: Tersoff_inc.h,v 1.3 2008-08-08 19:25:40 hspark Exp $ */
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
void get_dXsi(const double* params, const double *Xsi, const double *Xa, const double *Ya, const double *Za, const double* Cmat, double* dXsi, double* ddXsi); 

/* function to compute derivatives of the potential function wrt to the stretch tensor */
void get_dUdC(const double* params, const double *Xsi, const double *Xa, const double *Ya, const double *Za, const double* Cmat, double* dUdC); 

/* function to compute all second order derivatives of the potential function needed for the modulus */
//void get_ddC(const double* params, const double *Xsi, const double *Xa, const double *Ya, const double *Za, const double* Cmat, double* dCdC, double* dCdXsi);
void ddC_driver(const double* params,const double* Xsi,const double* Xa,const double* Ya,const double* Za,const double* Cmat,double* dCdC,double* dCdXsi);

/* function to get the bulk strain energy density */
double get_energy(const double* params, const double *Xsi, const double *Xa, const double *Ya, const double *Za, const double* Cmat); 

#ifdef __cplusplus
}
#endif

#endif /* TERSOFF_INC_H */
