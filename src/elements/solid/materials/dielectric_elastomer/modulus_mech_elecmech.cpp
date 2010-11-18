// Arruda-Boyce 3D
#include "FSDE_inc.h"

#include <math.h>

static double z[8];

/* function to compute mechanical modulus */
void get_ddCmech_elec(const double* params, const double *Xsi, const double* Cmat, double* dCdC, double* dCdXsi) { 

/* common definitions */
#include "FSDE_common_defines.h"
	
	/* Stress code */
	z[1] = C11 + C22 + C33;
	z[2] = pow(Nrig,-4.);
	z[3] = pow(Nrig,-3.);
	z[4] = 1./(Nrig*Nrig);
	z[5] = 1./Nrig;
	z[6] = z[1]*z[1];
	z[7] = pow(z[1],3.);
	z[1] = 0.06285714285714286*z[1]*z[4];
	z[4] = 0.1*z[5];
	z[3] = 0.03257142857142857*z[3]*z[6];
	z[2] = 0.015406307977736549*z[2]*z[7];
	z[1] = z[1] + z[2] + z[3] + z[4];
	z[1] = mu*z[1];

	/* dCdE:  6 x 3 */
	dCdXsi[0] = 0.0;
	dCdXsi[1] = 0.0;
	dCdXsi[2] = 0.0;
	dCdXsi[3] = 0.0;
	dCdXsi[4] = 0.0;
	dCdXsi[5] = 0.0;
	
	dCdXsi[6] = 0.0;
	dCdXsi[7] = 0.0;
	dCdXsi[8] = 0.0;
	dCdXsi[9] = 0.0;
	dCdXsi[10] = 0.0;
	dCdXsi[11] = 0.0;
	
	dCdXsi[12] = 0.0;
	dCdXsi[13] = 0.0;
	dCdXsi[14] = 0.0;
	dCdXsi[15] = 0.0;
	dCdXsi[16] = 0.0;
	dCdXsi[17] = 0.0;	

	/* dCdC: 6 x 6 */
	dCdC[ 0] = z[1];
	dCdC[ 1] = z[1];
	dCdC[ 2] = z[1];
	dCdC[ 3] = 0.0;
	dCdC[ 4] = 0.0;
	dCdC[ 5] = 0.0;

	dCdC[ 6] = z[1];
	dCdC[ 7] = z[1];
	dCdC[ 8] = z[1];
	dCdC[ 9] = 0.0;
	dCdC[10] = 0.0;
	dCdC[11] = 0.0;

	dCdC[12] = z[1];
	dCdC[13] = z[1];
	dCdC[14] = z[1];
	dCdC[15] = 0.0;
	dCdC[16] = 0.0;
	dCdC[17] = 0.0;

	dCdC[18] = 0.0;
	dCdC[19] = 0.0;
	dCdC[20] = 0.0;
	dCdC[21] = 0.0;
	dCdC[22] = 0.0;
	dCdC[23] = 0.0;

	dCdC[24] = 0.0;
	dCdC[25] = 0.0;
	dCdC[26] = 0.0;
	dCdC[27] = 0.0;
	dCdC[28] = 0.0;
	dCdC[29] = 0.0;

	dCdC[30] = 0.0;
	dCdC[31] = 0.0;
	dCdC[32] = 0.0;
	dCdC[33] = 0.0;
	dCdC[34] = 0.0;
	dCdC[35] = 0.0;
}