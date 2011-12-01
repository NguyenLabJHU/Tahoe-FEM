#include "FSDE_inc_2D.h"

#include <cmath>

static double z[5];

/* function to compute electrical displacement */
void get_ddE_2D(const double* params, const double *Xsi, const double* Cmat, double* ddXsi) { 

/* common definitions */
#include "FSDE_common_defines_2D.h"
	
	/* Stress code */
	z[1] = -C12*C21;
	z[2] = C11*C22;
	z[1] = z[1] + z[2];
	z[1] = 1./z[1];
	z[2] = -z[1];
	z[3] = -C11*epsilon*z[1];
	z[1] = -C22*epsilon*z[1];
	z[4] = C12*z[2];
	z[2] = C21*z[2];
	z[2] = z[2] + z[4];
	z[2] = -0.5*epsilon*z[2];

	/* Output stiffness */
	ddXsi[0] = z[1];
	ddXsi[1] = z[2];
	ddXsi[2] = z[2];
	ddXsi[3] = z[3];
}
