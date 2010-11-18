// Arruda-Boyce 3D
#include "FSDE_inc.h"
#include <math.h>
static double z[9];

/* function to compute mechanical stress */
void get_dUdCmech(const double* params, const double *Xsi, const double* Cmat, double* dUdC) { 

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
	z[8] = pow(z[1],4.);
	z[1] = 0.1*z[1]*z[5];
	z[4] = 0.03142857142857143*z[4]*z[6];
	z[3] = 0.010857142857142857*z[3]*z[7];
	z[2] = 0.0038515769944341373*z[2]*z[8];
	z[1] = 0.5 + z[1] + z[2] + z[3] + z[4];
	z[1] = mu*z[1];

	/* return values */
	dUdC[0] = z[1];
	dUdC[1] = 0.0;
	dUdC[2] = 0.0;
	dUdC[3] = 0.0;
	dUdC[4] = z[1];
	dUdC[5] = 0.0;
	dUdC[6] = 0.0;
	dUdC[7] = 0.0;
	dUdC[8] = z[1];
}