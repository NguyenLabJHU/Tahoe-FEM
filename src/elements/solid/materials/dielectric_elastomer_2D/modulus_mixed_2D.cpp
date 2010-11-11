#include "FSDE_inc_2D.h"

#include <math.h>

static double z[21];

/* function to compute mixed electromechanical modulus */
void get_ddCE_2D(const double* params, const double *Xsi, const double* Cmat, double* dCdXsi) { 

/* common definitions */
#include "FSDE_common_defines_2D.h"
	
	/* Stress code */
	z[1] = C11*C11;
	z[2] = C12*C12;
	z[3] = -C12*C21;
	z[4] = C21*C21;
	z[5] = C11*C22;
	z[6] = C22*C22;
	z[7] = z[3] + z[5];
	z[8] = 1./(z[7]*z[7]);
	z[7] = 1./z[7];
	z[9] = ex*z[8];
	z[10] = C11*z[9];
	z[11] = C12*z[9];
	z[12] = C21*C22*z[9];
	z[13] = C12*z[10];
	z[10] = C21*z[10];
	z[11] = C22*z[11];
	z[14] = 2.*z[12];
	z[2] = -z[2]*z[9];
	z[15] = z[3]*z[9];
	z[16] = -z[4]*z[9];
	z[6] = -2.*z[6]*z[9];
	z[9] = 2.*C11*C12*ey*z[8];
	z[17] = 2.*C11*C21*ey*z[8];
	z[18] = C12*C22*ey*z[8];
	z[19] = C21*C22*ey*z[8];
	z[1] = -2.*ey*z[1]*z[8];
	z[3] = ey*z[3]*z[8];
	z[4] = -ey*z[4]*z[8];
	z[5] = -2.*ey*z[5]*z[8];
	z[8] = -ex*z[7];
	z[20] = -ey*z[7];
	z[7] = 2.*ey*z[7];
	z[6] = z[18] + z[19] + z[6];
	z[1] = z[1] + z[10] + z[13];
	z[2] = z[15] + z[2] + z[8] + z[9];
	z[8] = z[15] + z[16] + z[17] + z[8];
	z[3] = z[14] + z[20] + z[3] + z[4];
	z[4] = z[11] + z[12] + z[5] + z[7];
	z[5] = -0.5*epsilon;
	z[6] = z[5]*z[6];
	z[1] = z[1]*z[5];
	z[2] = z[2]*z[5];
	z[7] = z[5]*z[8];
	z[3] = z[3]*z[5];
	z[4] = z[4]*z[5];

	/* Output stress */
	/* dCdE:  3 x 2 */
	dCdXsi[0] = z[6];
	dCdXsi[1] = z[2];
	dCdXsi[2] = z[4];
	dCdXsi[3] = z[3];
	dCdXsi[4] = z[1];
	dCdXsi[5] = z[7];
}