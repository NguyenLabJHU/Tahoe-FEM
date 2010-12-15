#include "FSDE_inc.h"

#include <math.h>

static double z[76];

/* function to compute electrical displacement */
void me_mixedmodulus_ab(const double* params, const double *Xsi, const double* Cmat, double J, double* dCdXsi) { 

/* common definitions */
#include "FSDE_common_defines.h"
	
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
	z[7] = z[12]*z[12];
	z[8] = z[10]*z[10];
	z[9] = z[14]*z[14];
	z[11] = 1./(z[16]*z[16]);
	z[13] = z[1]*z[1];
	z[15] = z[2]*z[2];
	z[16] = z[3]*z[3];
	z[17] = z[4]*z[4];
	z[18] = z[5]*z[5];
	z[19] = z[6]*z[6];
	z[12] = z[11]*z[12];
	z[20] = -z[1]*z[10]*z[11];
	z[21] = z[1]*z[10]*z[11];
	z[22] = -z[1]*z[11]*z[14];
	z[23] = z[1]*z[11]*z[14];
	z[24] = -z[10]*z[11]*z[2];
	z[25] = z[10]*z[11]*z[2];
	z[26] = -z[11]*z[14]*z[2];
	z[27] = z[11]*z[14]*z[2];
	z[28] = -2.*z[1]*z[11]*z[2];
	z[29] = -z[10]*z[11]*z[3];
	z[30] = z[10]*z[11]*z[3];
	z[31] = epsilon*ex*J*z[11]*z[14]*z[3];
	z[32] = epsilon*ez*J*z[11]*z[14]*z[3];
	z[33] = -2.*z[1]*z[11]*z[3];
	z[34] = -z[11]*z[2]*z[3];
	z[35] = -z[10]*z[11]*z[4];
	z[36] = z[10]*z[11]*z[4];
	z[37] = -z[11]*z[14]*z[4];
	z[38] = z[11]*z[14]*z[4];
	z[39] = z[1]*z[11]*z[4];
	z[40] = -z[11]*z[2]*z[4];
	z[41] = z[11]*z[2]*z[4];
	z[42] = -z[11]*z[3]*z[4];
	z[43] = z[11]*z[3]*z[4];
	z[44] = z[10]*z[11]*z[5];
	z[45] = -z[11]*z[14]*z[5];
	z[46] = z[11]*z[14]*z[5];
	z[47] = -z[1]*z[11]*z[5];
	z[48] = z[1]*z[11]*z[5];
	z[49] = epsilon*ey*J*z[11]*z[2]*z[5];
	z[50] = epsilon*ez*J*z[11]*z[2]*z[5];
	z[51] = -z[11]*z[3]*z[5];
	z[52] = z[11]*z[3]*z[5];
	z[53] = -z[10]*z[11]*z[6];
	z[54] = z[10]*z[11]*z[6];
	z[55] = z[1]*z[11]*z[6];
	z[56] = -z[11]*z[2]*z[6];
	z[57] = z[11]*z[2]*z[6];
	z[58] = -z[11]*z[3]*z[6];
	z[59] = z[11]*z[3]*z[6];
	z[60] = -z[11]*z[4]*z[6];
	z[61] = z[11]*z[4]*z[6];
	z[62] = -z[11]*z[5]*z[6];
	z[63] = z[11]*z[5]*z[6];
	z[7] = -z[11]*z[7];
	z[64] = -z[11]*z[8];
	z[8] = z[11]*z[8];
	z[9] = epsilon*ex*J*z[11]*z[9];
	z[13] = epsilon*ez*J*z[11]*z[13];
	z[15] = -z[11]*z[15];
	z[16] = -z[11]*z[16];
	z[65] = -z[11]*z[17];
	z[17] = z[11]*z[17];
	z[18] = epsilon*ey*J*z[11]*z[18];
	z[66] = -z[11]*z[19];
	z[11] = z[11]*z[19];
	z[19] = -z[10]*z[12];
	z[10] = z[10]*z[12];
	z[67] = epsilon*ex*J*z[12]*z[14];
	z[14] = epsilon*ey*J*z[12]*z[14];
	z[68] = -z[1]*z[12];
	z[1] = z[1]*z[12];
	z[2] = -z[12]*z[2];
	z[3] = -z[12]*z[3];
	z[69] = -z[12]*z[4];
	z[4] = z[12]*z[4];
	z[5] = -2.*z[12]*z[5];
	z[70] = -z[12]*z[6];
	z[6] = z[12]*z[6];
	z[12] = epsilon*ey*J*z[39];
	z[71] = epsilon*ez*J*z[39];
	z[72] = epsilon*ex*J*z[44];
	z[73] = epsilon*ey*J*z[44];
	z[74] = epsilon*ex*J*z[55];
	z[75] = epsilon*ez*J*z[55];
	z[27] = z[27] + z[29] + z[3];
	z[30] = z[26] + z[3] + z[30];
	z[6] = z[3] + z[37] + z[6];
	z[3] = z[3] + z[38] + z[70];
	z[28] = z[28] + z[39];
	z[39] = z[21] + z[34] + z[42];
	z[43] = z[20] + z[34] + z[43];
	z[5] = z[44] + z[5];
	z[25] = z[2] + z[25] + z[51];
	z[4] = z[35] + z[4] + z[51];
	z[44] = z[36] + z[51] + z[69];
	z[51] = z[2] + z[24] + z[52];
	z[38] = z[29] + z[38] + z[53];
	z[26] = z[26] + z[54] + z[70];
	z[29] = z[29] + z[37] + z[54];
	z[33] = z[33] + z[55];
	z[34] = z[34] + z[57] + z[68];
	z[11] = z[11] + z[22] + z[58];
	z[37] = z[23] + z[58] + z[66];
	z[1] = z[1] + z[42] + z[60];
	z[21] = z[21] + z[56] + z[60];
	z[42] = z[42] + z[61] + z[68];
	z[20] = z[20] + z[56] + z[61];
	z[36] = z[24] + z[36] + z[62];
	z[2] = z[2] + z[63] + z[69];
	z[24] = z[24] + z[35] + z[63];
	z[10] = z[10] + z[45] + z[7];
	z[7] = z[19] + z[46] + z[7];
	z[35] = z[19] + z[46] + z[64];
	z[8] = z[19] + z[45] + z[8];
	z[9] = z[14] + z[32] + z[9];
	z[12] = z[12] + z[13] + z[74];
	z[13] = z[15] + z[41] + z[47];
	z[14] = z[15] + z[40] + z[48];
	z[15] = z[16] + z[23] + z[58];
	z[16] = z[16] + z[22] + z[59];
	z[19] = z[40] + z[48] + z[65];
	z[17] = z[17] + z[40] + z[47];
	z[18] = z[18] + z[50] + z[72];
	z[22] = -epsilon*J;
	z[23] = ex*z[22];
	z[32] = ey*z[22];
	z[22] = ez*z[22];
	z[38] = z[23]*z[38];
	z[26] = z[23]*z[26];
	z[29] = z[23]*z[29];
	z[11] = z[11]*z[23];
	z[37] = z[23]*z[37];
	z[21] = z[21]*z[23];
	z[20] = z[20]*z[23];
	z[36] = z[23]*z[36];
	z[24] = z[23]*z[24];
	z[35] = z[23]*z[35];
	z[8] = z[23]*z[8];
	z[6] = z[32]*z[6];
	z[3] = z[3]*z[32];
	z[5] = z[32]*z[5];
	z[4] = z[32]*z[4];
	z[23] = z[32]*z[44];
	z[1] = z[1]*z[32];
	z[40] = z[32]*z[42];
	z[2] = z[2]*z[32];
	z[10] = z[10]*z[32];
	z[7] = z[32]*z[7];
	z[19] = z[19]*z[32];
	z[17] = z[17]*z[32];
	z[27] = z[22]*z[27];
	z[30] = z[22]*z[30];
	z[28] = z[22]*z[28];
	z[32] = z[22]*z[39];
	z[39] = z[22]*z[43];
	z[25] = z[22]*z[25];
	z[41] = z[22]*z[51];
	z[33] = z[22]*z[33];
	z[34] = z[22]*z[34];
	z[13] = z[13]*z[22];
	z[14] = z[14]*z[22];
	z[15] = z[15]*z[22];
	z[16] = z[16]*z[22];
	z[1] = z[1] + z[37] + z[75];
	z[19] = z[19] + z[21] + z[71];
	z[7] = z[27] + z[67] + z[7];
	z[10] = z[10] + z[30] + z[67];
	z[17] = z[17] + z[20] + z[28];
	z[20] = z[23] + z[29] + z[32];
	z[4] = z[38] + z[39] + z[4];
	z[5] = z[25] + z[5] + z[8];
	z[8] = z[35] + z[41] + z[73];
	z[11] = z[11] + z[33] + z[40];
	z[2] = z[2] + z[26] + z[34];
	z[3] = z[15] + z[3] + z[31];
	z[6] = z[16] + z[31] + z[6];
	z[13] = z[13] + z[36] + z[49];
	z[14] = z[14] + z[24] + z[49];

	/* Output stiffness */
	dCdXsi[0] = z[9];
	dCdXsi[1] = z[8];
	dCdXsi[2] = z[1];
	dCdXsi[3] = z[4];
	dCdXsi[4] = z[6];
	dCdXsi[5] = z[10];
	
	dCdXsi[6] = z[7];
	dCdXsi[7] = z[18];
	dCdXsi[8] = z[19];
	dCdXsi[9] = z[13];
	dCdXsi[10] = z[2];
	dCdXsi[11] = z[5];
	
	dCdXsi[12] = z[3];
	dCdXsi[13] = z[14];
	dCdXsi[14] = z[12];
	dCdXsi[15] = z[17];
	dCdXsi[16] = z[11];
	dCdXsi[17] = z[20];	
}