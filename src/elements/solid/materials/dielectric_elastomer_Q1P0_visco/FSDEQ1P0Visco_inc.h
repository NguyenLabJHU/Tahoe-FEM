#ifndef FSDEQ1P0VISCO_INC_H
#define FSDEQ1P0VISCO_INC_H

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
void mech_pk2_q1p0visco(const double* params, const double* Xsi, const double* Cmat, const double* Fmat, double J, double I1, double* dUdCmech); 
void me_pk2_q1p0visco(const double* params, const double* Xsi, const double* Cmat, const double* Fmat, double J, double* dUdCmechelec); 

/* function to compute second derivative of the potential function for mixed electromechanical modulus */
void me_mixedmodulus_q1p0visco(const double* params, const double *Xsi, const double* Cmat, const double* Fmat, double J, double* dCdXsi);

/* function to compute second derivative of the potential function for mixed electromechanical modulus */
void me_mixedmodulus_q1p0viscospatial(const double* params, const double *Xsi, const double* Cmat, const double* Fmat, double J, double* Kempf);

/* function to compute second derivative of the potential function for both purely mechanical modulus */
void mech_tanmod_q1p0visco(const double* params, const double* Xsi, const double* Cmat, const double* Fmat, double J, double I1, double* ddCmech);
void me_tanmod_q1p0visco(const double* params, const double* Xsi, const double* Cmat, const double* Fmat, double J, double* ddCmechelec);

/* function to get electric displacement */
void elec_pk2_q1p0visco(const double* params, const double *Xsi, const double* Cmat, const double* Fmat, double J, double* dUdCelec); 

/* function to get electromechanical energy density */
double me_energy_q1p0visco(const double* params, const double *Xsi, const double* Cmat, const double* Fmat, double J, double I1); 

#ifdef __cplusplus
}
#endif

#endif /* FSDEQ1P0VISCO_INC_H */
