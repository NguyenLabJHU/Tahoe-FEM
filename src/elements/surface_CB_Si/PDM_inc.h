/* $Id: PDM_inc.h,v 1.2 2010-09-29 15:18:06 hspark Exp $ */
#ifndef PDM_INC_H
#define PDM_INC_H

#ifdef __cplusplus
extern "C" {
#endif

/* Sequence of parameters is:
 * Ex
 * Ey
 * Ez
 * econv
 * alphatot
 * alpha1
 */

/* function to compute derivatives of the potential function wrt to the
 * internal degrees of freedom */
void get_dXsi_pdm(const double* params, const double *Xsi, const double *Xa, const double *Ya, const double *Za, const double* Cmat, double* dXsi, double* ddXsi); 

/* function to compute derivatives of the potential function wrt to the stretch tensor */
void get_dUdC_pdm(const double* params, const double *Xsi, const double *Xa, const double *Ya, const double *Za, const double* Cmat, double* dUdC); 

/* function to compute all second order derivatives of the potential function needed for the modulus */
void get_ddC_pdm(const double* params,const double* Xsi,const double* Xa,const double* Ya,const double* Za,const double* Cmat,double* dCdC,double* dCdXsi);

/* function to get the bulk strain energy density */
double get_energy_pdm(const double* params, const double *Xsi, const double *Xa, const double *Ya, const double *Za, const double* Cmat); 

#ifdef __cplusplus
}
#endif

#endif /* PDM_INC_H */
