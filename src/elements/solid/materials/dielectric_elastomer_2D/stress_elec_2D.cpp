#include "FSDE_inc_2D.h"

#include <math.h>

static double z[7];

/* function to compute electrical displacement */
void get_dUdE_2D(const double* params, const double *Xsi, const double* Cmat, double* dXsi) { 

/* common definitions */
#include "FSDE_common_defines_2D.h"

	/* Stress code */
	z[1] = -C12*C21;
	z[2] = C11*C22;
	z[1] = z[1] + z[2];
	z[1] = 1./z[1];
	z[2] = -ex*z[1];
	z[3] = ex*z[1];
	z[4] = 2.*C11*ey*z[1];
	z[5] = -C12*ey*z[1];
	z[1] = -C21*ey*z[1];
	z[6] = C12*z[2];
	z[2] = C21*z[2];
	z[3] = 2.*C22*z[3];
	z[2] = z[2] + z[4] + z[6];
	z[1] = z[1] + z[3] + z[5];
	z[3] = -0.5*epsilon;
	z[2] = z[2]*z[3];
	z[1] = z[1]*z[3];

	/* return values */
	dXsi[0] = z[1];
	dXsi[1] = z[2];
}