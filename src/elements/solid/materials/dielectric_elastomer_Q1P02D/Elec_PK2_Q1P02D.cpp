#include "FSDE_incQ1P02D.h"

#include <math.h>

static double z[25];

/* function to compute electric displacement */
void elec_pk2_q1p02D(const double* params, const double* Xsi, const double* Cmat, const double* Fmat, double J, double* dUdCelec) { 

/* common definitions */
#include "FSDE_common_definesQ1P02D.h"

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
	z[7] = epsilon*J*z[7];
	z[8] = ex*z[7];
	z[9] = ey*z[7];
	z[7] = ez*z[7];
	z[10] = z[10]*z[8];
	z[11] = z[14]*z[8];
	z[6] = z[6]*z[8];
	z[8] = z[12]*z[9];
	z[4] = z[4]*z[9];
	z[5] = z[5]*z[9];
	z[1] = z[1]*z[7];
	z[2] = z[2]*z[7];
	z[3] = z[3]*z[7];
	z[1] = z[1] + z[4] + z[6];
	z[2] = z[10] + z[2] + z[5];
	z[3] = z[11] + z[3] + z[8];

	/* return values */
	dUdCelec[0] = z[3];
	dUdCelec[1] = z[2];
	dUdCelec[2] = z[1];
}