#include "FSDEQ1P0Visco_inc.h"

#include <math.h>

static double z[33];

/* function to compute mechanical part of PK2 stress */
void mech_pk2_q1p0visco(const double* params, const double* Xsi, const double* Cmat, const double* Fmat, double J, double I1, double* dUdCmech) { 

/* common definitions */
#include "FSDEQ1P0Visco_common_defines.h"

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
	z[19] = I1*I1;
	z[20] = pow(I1,3.);
	z[21] = pow(I1,4.);
	z[22] = pow(Nrig,-4.);
	z[23] = pow(Nrig,-3.);
	z[24] = 1./(Nrig*Nrig);
	z[25] = 1./Nrig;
	z[26] = C33*z[1];
	z[27] = C32*z[2];
	z[28] = C33*z[3];
	z[29] = C31*z[4];
	z[30] = C32*z[5];
	z[31] = C31*z[6];
	z[32] = log(J);
	z[12] = z[12] + z[16];
	z[10] = z[10] + z[17];
	z[14] = z[14] + z[18];
	z[16] = 0.3119777365491651*z[22];
	z[17] = 0.0038515769944341373*z[21]*z[22];
	z[18] = 0.29314285714285715*z[23];
	z[20] = 0.010857142857142857*z[20]*z[23];
	z[21] = 0.28285714285714286*z[24];
	z[19] = 0.03142857142857143*z[19]*z[24];
	z[22] = 0.3*z[25];
	z[23] = 0.1*I1*z[25];
	z[1] = z[1] + z[3];
	z[3] = z[26] + z[27] + z[28] + z[29] + z[30] + z[31];
	z[24] = 0.5*lambda*z[32];
	z[2] = z[2] + z[5];
	z[4] = z[4] + z[6];
	z[5] = z[11] + z[7];
	z[6] = z[15] + z[8];
	z[7] = z[13] + z[9];
	z[8] = 0.5 + z[16] + z[18] + z[21] + z[22];
	z[9] = 0.5 + z[17] + z[19] + z[20] + z[23];
	z[3] = 1./z[3];
	z[8] = -mu*z[8];
	z[9] = 2.*mu*z[9];
	z[8] = z[24] + z[8];
	z[3] = 2.*z[3]*z[8];
	z[8] = z[12]*z[3];
	z[10] = z[10]*z[3];
	z[11] = z[14]*z[3];
	z[1] = z[1]*z[3];
	z[2] = z[2]*z[3];
	z[4] = z[3]*z[4];
	z[5] = z[3]*z[5];
	z[6] = z[3]*z[6];
	z[3] = z[3]*z[7];
	z[7] = z[11] + z[9];
	z[1] = z[1] + z[9];
	z[6] = z[6] + z[9];

	/* return values */
	dUdCmech[0] = z[7];
	dUdCmech[1] = z[10];
	dUdCmech[2] = z[3];
	dUdCmech[3] = z[8];
	dUdCmech[4] = z[6];
	dUdCmech[5] = z[5];
	dUdCmech[6] = z[4];
	dUdCmech[7] = z[2];
	dUdCmech[8] = z[1];	
}