#include "FSDEQ1P0Visco_inc.h"

#include <math.h>

static double z[37];

/* function to compute electrical displacement */
double me_energy_q1p0visco(const double* params, const double *Xsi, const double* Cmat, const double* Fmat, double J, double I1) { 

/* common definitions */
#include "FSDEQ1P0Visco_common_defines.h"

	/* EM energy code */
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
	z[19] = ex*ex;
	z[20] = ey*ey;
	z[21] = ez*ez;
	z[22] = -3. + I1;
	z[23] = I1*I1;
	z[24] = pow(I1,3.);
	z[25] = pow(I1,4.);
	z[26] = pow(I1,5.);
	z[27] = pow(Nrig,-4.);
	z[28] = pow(Nrig,-3.);
	z[29] = 1./(Nrig*Nrig);
	z[30] = 1./Nrig;
	z[31] = C33*z[1];
	z[32] = C32*z[2];
	z[33] = C33*z[3];
	z[34] = C31*z[4];
	z[35] = C32*z[5];
	z[36] = C31*z[6];
	z[12] = z[12] + z[16];
	z[10] = z[10] + z[17];
	z[14] = z[14] + z[18];
	z[16] = 0.5*z[22];
	z[17] = -9. + z[23];
	z[18] = -27. + z[24];
	z[22] = -81. + z[25];
	z[23] = -243. + z[26];
	z[24] = 0.3119777365491651*z[27];
	z[25] = 0.29314285714285715*z[28];
	z[26] = 0.28285714285714286*z[29];
	z[1] = z[1] + z[3];
	z[3] = 0.3*z[30];
	z[31] = z[31] + z[32] + z[33] + z[34] + z[35] + z[36];
	z[2] = z[2] + z[5];
	z[4] = z[4] + z[6];
	z[5] = z[11] + z[7];
	z[6] = z[15] + z[8];
	z[7] = z[13] + z[9];
	z[8] = 0.05*z[17]*z[30];
	z[9] = 0.010476190476190476*z[18]*z[29];
	z[11] = 0.0027142857142857142*z[22]*z[28];
	z[13] = 0.0007703153988868275*z[23]*z[27];
	z[3] = 0.5 + z[24] + z[25] + z[26] + z[3];
	z[15] = 1./z[31];
	z[17] = sqrt(z[31]);
	z[8] = z[11] + z[13] + z[16] + z[8] + z[9];
	z[9] = ex*ey*z[15];
	z[11] = z[14]*z[15]*z[19];
	z[1] = z[1]*z[15]*z[21];
	z[2] = ey*ez*z[15]*z[2];
	z[4] = ex*ez*z[15]*z[4];
	z[5] = ey*ez*z[15]*z[5];
	z[6] = z[15]*z[20]*z[6];
	z[7] = ex*ez*z[15]*z[7];
	z[12] = z[12]*z[9];
	z[9] = z[10]*z[9];
	z[10] = log(z[17]);
	z[8] = mu*z[8];
	z[1] = z[1] + z[11] + z[12] + z[2] + z[4] + z[5] + z[6] + z[7] + z[9];
	z[2] = -2.*mu*z[10]*z[3];
	z[3] = z[10]*z[10];
	z[1] = -0.5*epsilon*z[1]*z[17];
	z[3] = 0.5*lambda*z[3];
	z[1] = z[1] + z[2] + z[3] + z[8];

	/* Output energy */
	return z[1];
}