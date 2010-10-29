#ifndef FSDE_INC_H
#define FSDE_INC_H

#ifdef __cplusplus
extern "C" {
#endif

/* Sequence of parameters is:
 * Nrig
 * Epsilon
 * Mu
 */

/* Gets first and second derivatives of free energy with respect to E */
void get_dXsi(const double* params, const double *Xsi, const double* Cmat, double* dXsi, double* ddXsi); 

/* function to compute first derivative of free energy wrt to the stretch tensor C */
void get_dUdC(const double* params, const double *Xsi, const double* Cmat, double* dUdC); 

/* function to compute second derivative of the potential function for both purely mechanical and mixed electromechanical modulus */
void get_ddC(const double* params, const double *Xsi, const double* Cmat, double* dCdC, double* dCdXsi);

/* function to get the bulk strain energy density - not implemented yet */
//double get_energy(const double* params, const double *Xsi, const double* Cmat); 

#ifdef __cplusplus
}
#endif

#endif /* FSDE_INC_H */
