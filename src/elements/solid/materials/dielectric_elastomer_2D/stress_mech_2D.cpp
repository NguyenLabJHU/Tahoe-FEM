#include "FSDE_inc_2D.h"
#include <math.h>
static double z[22];

/* function to compute mechanical stress */
void get_dUdCmech_2D(const double* params, const double *Xsi, const double* Cmat, double* dUdC) { 

/* common definitions */
#include "FSDE_common_defines_2D.h"
	
	/* Stress code */
	z[1] = C11*C11;
	z[2] = C12*C12;
	z[3] = -C12*C21;
	z[4] = C21*C21;
	z[5] = C11*C22;
	z[6] = C22*C22;
	z[7] = C11 + C22;
	z[8] = ex*ex;
	z[9] = ey*ey;
	z[10] = 1./(Nrig*Nrig);
	z[11] = 1./Nrig;
	z[12] = z[3] + z[5];
	z[11] = 0.1*z[11]*z[7];
	z[7] = z[7]*z[7];
	z[13] = 1./(z[12]*z[12]);
	z[12] = 1./z[12];
	z[7] = 0.03142857142857143*z[10]*z[7];
	z[10] = ex*ey;
	z[7] = 0.5 + z[11] + z[7];
	z[11] = z[10]*z[13];
	z[10] = -z[10]*z[12];
	z[14] = C11*z[11];
	z[15] = C12*z[11];
	z[16] = C21*C22*z[11];
	z[2] = -z[11]*z[2];
	z[17] = C12*z[14];
	z[14] = C21*z[14];
	z[15] = C22*z[15];
	z[3] = z[11]*z[3];
	z[4] = -z[11]*z[4];
	z[11] = C12*C22*z[13]*z[8];
	z[18] = C21*C22*z[13]*z[8];
	z[19] = z[12]*z[8];
	z[20] = -z[13]*z[5]*z[8];
	z[6] = -z[13]*z[6]*z[8];
	z[8] = C11*C12*z[13]*z[9];
	z[21] = C11*C21*z[13]*z[9];
	z[1] = -z[1]*z[13]*z[9];
	z[12] = z[12]*z[9];
	z[5] = -z[13]*z[5]*z[9];
	z[7] = mu*z[7];
	z[3] = z[10] + z[3];
	z[1] = z[1] + z[14] + z[17] + z[19] + z[20];
	z[5] = z[12] + z[15] + z[16] + z[5] + z[6];
	z[2] = z[11] + z[2] + z[3] + z[8];
	z[3] = z[18] + z[21] + z[3] + z[4];
	z[4] = -0.5*epsilon;
	z[1] = z[1]*z[4];
	z[5] = z[4]*z[5];
	z[2] = z[2]*z[4];
	z[3] = z[3]*z[4];
	z[1] = z[1] + z[7];
	z[4] = z[5] + z[7];

	/* return values */
	dUdC[0] = z[4];
	dUdC[1] = z[2];
	dUdC[2] = z[3];
	dUdC[3] = z[1];
}