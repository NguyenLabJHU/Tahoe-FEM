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

/* function to compute first derivative of free energy wrt to the stretch tensor C */
void mech_pk2_ab(const double* params, const double* Xsi, const double* Cmat, double J, double I1, double* dUdCmech); 
void me_pk2_ab(const double* params, const double* Xsi, const double* Cmat, double J, double* dUdCmechelec); 

/* function to compute second derivative of the potential function for mixed electromechanical modulus */
void me_mixedmodulus_ab(const double* params, const double *Xsi, const double* Cmat, double J, double* dCdXsi);

/* function to compute second derivative of the potential function for both purely mechanical modulus */
void mech_tanmod_ab(const double* params, const double* Xsi, const double* Cmat, double J, double I1, double* ddCmech);
void me_tanmod_ab(const double* params, const double* Xsi, const double* Cmat, double J, double* ddCmechelec);

/* function to get the bulk strain energy density - not implemented yet */
//double get_energy(const double* params, const double *Xsi, const double* Cmat); 

#ifdef __cplusplus
}
#endif

#endif /* FSDE_INC_H */
