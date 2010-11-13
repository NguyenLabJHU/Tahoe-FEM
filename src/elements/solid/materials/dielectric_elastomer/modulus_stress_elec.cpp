#include "FSDE_inc.h"

#include <math.h>

static double z[26];

/* function to compute electrical displacement */
void get_ddC_sc_elec(const double* params, const double *Xsi, const double* Cmat, double* dE, double* dEdE) { 

/* common definitions */
#include "FSDE_common_defines.h"
	
	/* Stress code */
	z[1] = -C12*C21;
	z[2] = C13*C21;
	z[3] = C11*C22;
	z[4] = -C13*C22;
	z[5] = -C11*C23;
	z[6] = C12*C23;
	z[7] = C12*C31;
	z[8] = -C13*C31;
	z[9] = -C22*C31;
	z[10] = C23*C31;
	z[11] = -C11*C32;
	z[12] = C13*C32;
	z[13] = C21*C32;
	z[14] = -C23*C32;
	z[15] = C11*C33;
	z[16] = -C12*C33;
	z[17] = -C21*C33;
	z[18] = C22*C33;
	z[19] = C33*z[1];
	z[20] = C32*z[2];
	z[21] = C33*z[3];
	z[22] = C31*z[4];
	z[23] = C32*z[5];
	z[24] = C31*z[6];
	z[12] = z[12] + z[16];
	z[10] = z[10] + z[17];
	z[14] = z[14] + z[18];
	z[16] = z[19] + z[20] + z[21] + z[22] + z[23] + z[24];
	z[1] = z[1] + z[3];
	z[2] = z[2] + z[5];
	z[3] = z[4] + z[6];
	z[4] = z[11] + z[7];
	z[5] = z[15] + z[8];
	z[6] = z[13] + z[9];
	z[7] = 1./z[16];
	z[8] = 1./sqrt(z[16]);
	z[9] = sqrt(z[16]);
	z[11] = z[12]*z[7];
	z[10] = z[10]*z[7];
	z[12] = 2.*ex*z[14]*z[7];
	z[13] = 2.*ez*z[1]*z[7];
	z[2] = z[2]*z[7];
	z[3] = z[3]*z[7];
	z[4] = z[4]*z[7];
	z[15] = 2.*ey*z[5]*z[7];
	z[6] = z[6]*z[7];
	z[7] = -epsilon*z[14]*z[8];
	z[1] = -epsilon*z[1]*z[8];
	z[5] = -epsilon*z[5]*z[8];
	z[8] = ex*z[11];
	z[14] = ey*z[11];
	z[16] = ex*z[10];
	z[17] = ey*z[10];
	z[18] = ey*z[2];
	z[19] = ez*z[2];
	z[20] = ex*z[3];
	z[21] = ez*z[3];
	z[22] = ey*z[4];
	z[23] = ez*z[4];
	z[24] = ex*z[6];
	z[25] = ez*z[6];
	z[10] = z[10] + z[11];
	z[2] = z[2] + z[4];
	z[3] = z[3] + z[6];
	z[4] = z[15] + z[16] + z[19] + z[23] + z[8];
	z[6] = z[13] + z[18] + z[20] + z[22] + z[24];
	z[8] = z[12] + z[14] + z[17] + z[21] + z[25];
	z[9] = -0.5*epsilon*z[9];
	z[10] = z[10]*z[9];
	z[2] = z[2]*z[9];
	z[3] = z[3]*z[9];
	z[4] = z[4]*z[9];
	z[6] = z[6]*z[9];
	z[8] = z[8]*z[9];


	/* Output stress */
	dE[0] = z[8];
	dE[1] = z[6];
	dE[2] = z[4];

	/* Output stiffness */
	dEdE[0] = z[7];
	dEdE[1] = z[10];
	dEdE[2] = z[3];
	dEdE[3] = z[10];
	dEdE[4] = z[5];
	dEdE[5] = z[2];
	dEdE[6] = z[3];
	dEdE[7] = z[2];
	dEdE[8] = z[1];
}
