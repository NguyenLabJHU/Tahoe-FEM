#ifndef FSDEVISCO_INC_H
#define FSDEVISCO_INC_H

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
void mech_pk2_ab_visco(const double* params, const double* Xsi, const double* Cmat, double J, double I1, double* dUdCmech); 
void me_pk2_ab_visco(const double* params, const double* Xsi, const double* Cmat, double J, double* dUdCmechelec); 

/* function to compute second derivative of the potential function for mixed electromechanical modulus */
void me_mixedmodulus_ab_visco(const double* params, const double *Xsi, const double* Cmat, double J, double* dCdXsi);

/* function to compute second derivative of the potential function for both purely mechanical modulus */
void mech_tanmod_ab_visco(const double* params, const double* Xsi, const double* Cmat, double J, double I1, double* ddCmech);
void me_tanmod_ab_visco(const double* params, const double* Xsi, const double* Cmat, double J, double* ddCmechelec);

/* function to get electric displacement */
void elec_pk2_visco(const double* params, const double *Xsi, const double* Cmat, double J, double* dUdCelec); 

#ifdef __cplusplus
}
#endif

#endif /* FSDEVISCO_INC_H */
