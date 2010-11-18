#ifndef FSDE_INC_H
#define FSDE_INC_H

#ifdef __cplusplus
extern "C" {
#endif

/* Sequence of parameters is:
 * epsilon
 * mu
 * Nrig
 * lambda
 */

/* Gets first derivative of free energy with respect to E (electric displacement) */
void get_dUdE(const double* params, const double *Xsi, const double* Cmat, double* dXsi); 

/* Gets second derivative of free energy with respect to E */
void get_ddE(const double* params, const double *Xsi, const double* Cmat, double* ddXsi); 

/* function to compute first derivative of free energy wrt to the stretch tensor C */
void get_dUdCmech(const double* params, const double *Xsi, const double* Cmat, double* dUdC); 

/* function to compute second derivative of the potential function for mixed electromechanical modulus */
void get_ddCE(const double* params, const double *Xsi, const double* Cmat, double* dCdXsi);

/* function to compute second derivative of the potential function for both purely mechanical modulus */
void get_ddCmech(const double* params, const double *Xsi, const double* Cmat, double* dCdC);

/* function to compute mechanical and electromechanical tangent modulus */
void get_ddCmech_elec(const double* params, const double *Xsi, const double* Cmat, double* dCdC, double* dCdXsi);

/* function to compute electric displacement and tangent modulus */
void get_ddC_sc_elec(const double* params, const double *Xsi, const double* Cmat, double* dE, double* dEdE);

/* function to get the bulk strain energy density - not implemented yet */
//double get_energy(const double* params, const double *Xsi, const double* Cmat); 

#ifdef __cplusplus
}
#endif

#endif /* FSDE_INC_H */
