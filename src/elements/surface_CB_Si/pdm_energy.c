/* $Id: pdm_energy.c,v 1.2 2010-09-29 15:19:38 hspark Exp $ */
#include "PDM_inc.h"
#include <math.h>
static double z[64];

/* function to compute bulk strain energy density for point dipole method (thole dipole interaction tensor) */
double get_energy_pdm(const double* params, const double *Xsi, const double *Xa, const double *Ya, const double *Za, const double* Cmat) { 
	
/* common definitions */
#include "PDM_common_defines.h"
	
	/* Energy code */
	z[1] = 1./alpha1;
	z[2] = pow(alphatot,-4.);
	z[3] = pow(alphatot,-3.);
	z[4] = -X2;
	z[5] = -X3;
	z[6] = -X4;
	z[7] = -X5;
	z[8] = -Y2;
	z[9] = -Y3;
	z[10] = -Y4;
	z[11] = -Y5;
	z[12] = -Z2;
	z[13] = -Z3;
	z[14] = -Z4;
	z[15] = -Z5;
	z[16] = Y1 + Ys1;
	z[1] = -z[1];
	z[10] = z[10] + z[16];
	z[11] = z[11] + z[16];
	z[3] = -16.*z[3];
	z[4] = X1 + Xs1 + z[4];
	z[5] = X1 + Xs1 + z[5];
	z[6] = X1 + Xs1 + z[6];
	z[7] = X1 + Xs1 + z[7];
	z[8] = z[16] + z[8];
	z[9] = z[16] + z[9];
	z[12] = Z1 + z[12] + Zs1;
	z[13] = Z1 + z[13] + Zs1;
	z[14] = Z1 + z[14] + Zs1;
	z[15] = Z1 + z[15] + Zs1;
	z[16] = C21*z[10];
	z[17] = C22*z[10];
	z[18] = C23*z[10];
	z[19] = z[10]*z[10];
	z[20] = C21*z[11];
	z[21] = C22*z[11];
	z[22] = C23*z[11];
	z[23] = z[11]*z[11];
	z[24] = C11*z[4];
	z[25] = C12*z[4];
	z[26] = C13*z[4];
	z[27] = z[4]*z[4];
	z[28] = C11*z[5];
	z[29] = C12*z[5];
	z[30] = C13*z[5];
	z[31] = z[5]*z[5];
	z[32] = C11*z[6];
	z[33] = C12*z[6];
	z[34] = C13*z[6];
	z[35] = z[6]*z[6];
	z[36] = C11*z[7];
	z[37] = C12*z[7];
	z[38] = C13*z[7];
	z[39] = z[7]*z[7];
	z[40] = C21*z[8];
	z[41] = C22*z[8];
	z[42] = C23*z[8];
	z[43] = z[8]*z[8];
	z[44] = C21*z[9];
	z[45] = C22*z[9];
	z[46] = C23*z[9];
	z[47] = z[9]*z[9];
	z[48] = C31*z[12];
	z[49] = C32*z[12];
	z[50] = C33*z[12];
	z[51] = z[12]*z[12];
	z[52] = C31*z[13];
	z[53] = C32*z[13];
	z[54] = C33*z[13];
	z[55] = z[13]*z[13];
	z[56] = C31*z[14];
	z[57] = C32*z[14];
	z[58] = C33*z[14];
	z[59] = z[14]*z[14];
	z[60] = C31*z[15];
	z[61] = C32*z[15];
	z[62] = C33*z[15];
	z[63] = z[15]*z[15];
	z[24] = z[24] + z[40] + z[48];
	z[25] = z[25] + z[41] + z[49];
	z[26] = z[26] + z[42] + z[50];
	z[28] = z[28] + z[44] + z[52];
	z[29] = z[29] + z[45] + z[53];
	z[30] = z[30] + z[46] + z[54];
	z[16] = z[16] + z[32] + z[56];
	z[17] = z[17] + z[33] + z[57];
	z[18] = z[18] + z[34] + z[58];
	z[20] = z[20] + z[36] + z[60];
	z[21] = z[21] + z[37] + z[61];
	z[22] = z[22] + z[38] + z[62];
	z[24] = z[24]*z[4];
	z[25] = z[25]*z[8];
	z[26] = z[12]*z[26];
	z[28] = z[28]*z[5];
	z[29] = z[29]*z[9];
	z[30] = z[13]*z[30];
	z[16] = z[16]*z[6];
	z[17] = z[10]*z[17];
	z[18] = z[14]*z[18];
	z[20] = z[20]*z[7];
	z[21] = z[11]*z[21];
	z[22] = z[15]*z[22];
	z[20] = z[20] + z[21] + z[22];
	z[21] = z[24] + z[25] + z[26];
	z[22] = z[28] + z[29] + z[30];
	z[16] = z[16] + z[17] + z[18];
	z[17] = 1./sqrt(z[20]);
	z[18] = sqrt(z[20]);
	z[20] = 1./sqrt(z[21]);
	z[21] = sqrt(z[21]);
	z[24] = 1./sqrt(z[22]);
	z[22] = sqrt(z[22]);
	z[25] = 1./sqrt(z[16]);
	z[16] = sqrt(z[16]);
	z[2] = 3.*z[2];
	z[18] = z[18]*z[2];
	z[21] = z[2]*z[21];
	z[22] = z[2]*z[22];
	z[16] = z[16]*z[2];
	z[26] = z[10]*z[2]*z[25]*z[6];
	z[28] = z[11]*z[17]*z[2]*z[7];
	z[29] = z[2]*z[20]*z[4]*z[8];
	z[30] = z[2]*z[24]*z[5]*z[9];
	z[4] = z[12]*z[2]*z[20]*z[4];
	z[8] = z[12]*z[2]*z[20]*z[8];
	z[5] = z[13]*z[2]*z[24]*z[5];
	z[9] = z[13]*z[2]*z[24]*z[9];
	z[10] = z[10]*z[14]*z[2]*z[25];
	z[6] = z[14]*z[2]*z[25]*z[6];
	z[11] = z[11]*z[15]*z[17]*z[2];
	z[7] = z[15]*z[17]*z[2]*z[7];
	z[12] = z[19]*z[2]*z[25];
	z[13] = z[17]*z[2]*z[23];
	z[14] = z[2]*z[20]*z[27];
	z[15] = z[2]*z[24]*z[31];
	z[19] = z[2]*z[25]*z[35];
	z[23] = z[17]*z[2]*z[39];
	z[27] = z[2]*z[20]*z[43];
	z[31] = z[2]*z[24]*z[47];
	z[20] = z[2]*z[20]*z[51];
	z[24] = z[2]*z[24]*z[55];
	z[25] = z[2]*z[25]*z[59];
	z[2] = z[17]*z[2]*z[63];
	z[17] = z[26] + z[28] + z[29] + z[30];
	z[8] = z[10] + z[11] + z[8] + z[9];
	z[4] = z[4] + z[5] + z[6] + z[7];
	z[5] = z[1] + z[14] + z[15] + z[16] + z[18] + z[19] + z[21] + z[22] + z[23] + z[3];
	z[6] = z[1] + z[12] + z[13] + z[16] + z[18] + z[21] + z[22] + z[27] + z[3] + z[31];
	z[1] = z[1] + z[16] + z[18] + z[2] + z[20] + z[21] + z[22] + z[24] + z[25] + z[3];
	z[2] = z[17]*z[17];
	z[3] = z[17]*z[8];
	z[7] = z[8]*z[8];
	z[9] = -z[17]*z[4];
	z[10] = z[17]*z[4];
	z[11] = z[4]*z[8];
	z[12] = z[4]*z[4];
	z[13] = -z[5]*z[8];
	z[14] = z[5]*z[8];
	z[15] = -z[4]*z[6];
	z[16] = z[5]*z[6];
	z[17] = -z[1]*z[17];
	z[5] = z[1]*z[5];
	z[6] = z[1]*z[6];
	z[2] = -z[2];
	z[7] = -z[7];
	z[12] = -z[12];
	z[10] = z[10] + z[13];
	z[9] = z[14] + z[9];
	z[3] = z[15] + z[3];
	z[11] = z[11] + z[17];
	z[2] = z[16] + z[2];
	z[6] = z[6] + z[7];
	z[5] = z[12] + z[5];
	z[7] = -z[8]*z[9];
	z[4] = z[3]*z[4];
	z[1] = z[1]*z[2];
	z[1] = z[1] + z[4] + z[7];
	z[1] = 1./z[1];
	z[4] = z[1]*z[10];
	z[7] = ex*z[1]*z[3];
	z[3] = ez*z[1]*z[3];
	z[8] = ex*z[1]*z[11];
	z[9] = ey*z[1]*z[11];
	z[2] = ez*z[1]*z[2];
	z[6] = ex*z[1]*z[6];
	z[1] = ey*z[1]*z[5];
	z[5] = ey*z[4];
	z[4] = ez*z[4];
	z[3] = z[3] + z[6] + z[9];
	z[2] = z[2] + z[5] + z[7];
	z[1] = z[1] + z[4] + z[8];
	z[3] = ex*z[3];
	z[2] = ez*z[2];
	z[1] = ey*z[1];
	z[1] = z[1] + z[2] + z[3];
	z[1] = -0.5*econv*z[1];

	/* Return value */
	return z[1];
}