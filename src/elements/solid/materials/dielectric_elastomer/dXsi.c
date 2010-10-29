#include "FSDE_inc.h"

#include <math.h>

static double z[25];

/* function to compute first and second derivatives of free energy wrt E-field */
/* Gives electric displacement D and electrical tangent modulus */
void get_dXsi(const double* params, const double *Xsi, const double* Cmat, double* dXsi, double* ddXsi) { 

/* common definitions */
#include "FSDE_common_defines.h"

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
	z[8] = z[12]*z[7];
	z[9] = z[10]*z[7];
	z[10] = 1.*epsilon*z[14]*z[7];
	z[11] = 2.*ex*z[14]*z[7];
	z[12] = 1.*epsilon*z[1]*z[7];
	z[1] = 2.*ez*z[1]*z[7];
	z[2] = z[2]*z[7];
	z[3] = z[3]*z[7];
	z[4] = z[4]*z[7];
	z[13] = 1.*epsilon*z[5]*z[7];
	z[5] = 2.*ey*z[5]*z[7];
	z[6] = z[6]*z[7];
	z[7] = ex*z[8];
	z[14] = ey*z[8];
	z[15] = ex*z[9];
	z[16] = ey*z[9];
	z[17] = ey*z[2];
	z[18] = ez*z[2];
	z[19] = ex*z[3];
	z[20] = ez*z[3];
	z[21] = ey*z[4];
	z[22] = ez*z[4];
	z[23] = ex*z[6];
	z[24] = ez*z[6];
	z[8] = z[8] + z[9];
	z[2] = z[2] + z[4];
	z[3] = z[3] + z[6];
	z[4] = z[15] + z[18] + z[22] + z[5] + z[7];
	z[1] = z[1] + z[17] + z[19] + z[21] + z[23];
	z[5] = z[11] + z[14] + z[16] + z[20] + z[24];
	z[6] = 0.5*epsilon;
	z[7] = z[6]*z[8];
	z[2] = z[2]*z[6];
	z[3] = z[3]*z[6];
	z[4] = z[4]*z[6];
	z[1] = z[1]*z[6];
	z[5] = z[5]*z[6];

	/* output */
	/* dXsi  = {z5, z4, z1} */
	/* ddXsi = {{z10, z7, z3},
	            {z7, z13, z2},
				{z3, z2, z12}} */

	/* return values */
	dXsi[0] = z[5];
	dXsi[1] = z[4];
	dXsi[2] = z[1];
	
	ddXsi[0] = z[10];
	ddXsi[1] = z[7];
	ddXsi[2] = z[3];
	ddXsi[3] = z[7];
	ddXsi[4] = z[13];
	ddXsi[5] = z[2];
	ddXsi[6] = z[3];
	ddXsi[7] = z[2];
	ddXsi[8] = z[12];
}