#include "FSDE_inc.h"

#include <math.h>

static double z[109];

/* function to compute electrical displacement */
void me_pk2_ab(const double* params, const double *Xsi, const double* Cmat, double J, double* dUdCmechelec) { 

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
	z[19] = ex*ex;
	z[20] = ey*ey;
	z[21] = ez*ez;
	z[22] = C33*z[1];
	z[23] = C32*z[2];
	z[24] = C33*z[3];
	z[25] = C31*z[4];
	z[26] = C32*z[5];
	z[27] = C31*z[6];
	z[12] = z[12] + z[16];
	z[10] = z[10] + z[17];
	z[14] = z[14] + z[18];
	z[16] = z[22] + z[23] + z[24] + z[25] + z[26] + z[27];
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
	z[22] = z[6]*z[6];
	z[12] = z[11]*z[12];
	z[23] = -0.5*z[10]*z[11]*z[14];
	z[24] = z[10]*z[11]*z[14];
	z[25] = -0.5*z[1]*z[10]*z[11];
	z[26] = z[1]*z[10]*z[11];
	z[27] = -0.5*z[1]*z[11]*z[14];
	z[28] = z[1]*z[11]*z[14];
	z[29] = -0.5*z[10]*z[11]*z[2];
	z[30] = z[10]*z[11]*z[2];
	z[31] = -0.5*z[11]*z[14]*z[2];
	z[32] = z[11]*z[14]*z[2];
	z[33] = -0.5*z[1]*z[11]*z[2];
	z[34] = z[1]*z[11]*z[2];
	z[35] = -0.5*z[10]*z[11]*z[3];
	z[36] = z[10]*z[11]*z[3];
	z[37] = -0.5*z[11]*z[14]*z[3];
	z[38] = z[11]*z[14]*z[3];
	z[39] = -0.5*z[1]*z[11]*z[3];
	z[40] = z[1]*z[11]*z[3];
	z[41] = -0.5*z[11]*z[2]*z[3];
	z[42] = z[11]*z[2]*z[3];
	z[43] = -0.5*z[10]*z[11]*z[4];
	z[44] = z[10]*z[11]*z[4];
	z[45] = -0.5*z[11]*z[14]*z[4];
	z[46] = z[11]*z[14]*z[4];
	z[47] = -0.5*z[1]*z[11]*z[4];
	z[48] = z[1]*z[11]*z[4];
	z[49] = 0.5*epsilon*ey*ez*J*z[11]*z[2]*z[4];
	z[50] = -0.5*z[11]*z[3]*z[4];
	z[51] = z[11]*z[3]*z[4];
	z[52] = -0.5*z[10]*z[11]*z[5];
	z[53] = z[10]*z[11]*z[5];
	z[54] = -0.5*z[11]*z[14]*z[5];
	z[55] = z[11]*z[14]*z[5];
	z[56] = -0.5*z[1]*z[11]*z[5];
	z[57] = z[1]*z[11]*z[5];
	z[58] = -0.5*z[11]*z[2]*z[5];
	z[59] = z[11]*z[2]*z[5];
	z[60] = -0.5*z[11]*z[3]*z[5];
	z[61] = z[11]*z[3]*z[5];
	z[62] = -0.5*z[11]*z[4]*z[5];
	z[63] = z[11]*z[4]*z[5];
	z[64] = -0.5*z[10]*z[11]*z[6];
	z[65] = z[10]*z[11]*z[6];
	z[66] = -0.5*z[11]*z[14]*z[6];
	z[67] = z[11]*z[14]*z[6];
	z[68] = -0.5*z[1]*z[11]*z[6];
	z[69] = z[1]*z[11]*z[6];
	z[70] = -0.5*z[11]*z[2]*z[6];
	z[71] = z[11]*z[2]*z[6];
	z[72] = 0.5*epsilon*ex*ez*J*z[11]*z[3]*z[6];
	z[73] = -0.5*z[11]*z[4]*z[6];
	z[74] = z[11]*z[4]*z[6];
	z[75] = -0.5*z[11]*z[5]*z[6];
	z[76] = z[11]*z[5]*z[6];
	z[77] = -0.5*z[11]*z[7];
	z[7] = z[11]*z[7];
	z[78] = -0.5*z[11]*z[8];
	z[8] = z[11]*z[8];
	z[9] = 0.5*epsilon*J*z[11]*z[19]*z[9];
	z[13] = 0.5*epsilon*J*z[11]*z[13]*z[21];
	z[79] = -0.5*z[11]*z[15];
	z[15] = z[11]*z[15];
	z[80] = -0.5*z[11]*z[16];
	z[16] = z[11]*z[16];
	z[81] = -0.5*z[11]*z[17];
	z[17] = z[11]*z[17];
	z[18] = 0.5*epsilon*J*z[11]*z[18]*z[20];
	z[82] = -0.5*z[11]*z[22];
	z[11] = z[11]*z[22];
	z[10] = 0.5*epsilon*ex*ey*J*z[10]*z[12];
	z[22] = -0.5*z[12]*z[14];
	z[14] = z[12]*z[14];
	z[83] = 0.5*epsilon*J*z[14]*z[19];
	z[84] = -0.5*z[1]*z[12];
	z[1] = z[1]*z[12];
	z[85] = -0.5*z[12]*z[2];
	z[2] = z[12]*z[2];
	z[86] = 0.5*epsilon*ex*ey*J*z[2];
	z[87] = -0.5*z[12]*z[3];
	z[3] = z[12]*z[3];
	z[88] = -0.5*z[12]*z[4];
	z[4] = z[12]*z[4];
	z[89] = -0.5*z[12]*z[5];
	z[5] = z[12]*z[5];
	z[90] = 0.5*epsilon*ex*ey*J*z[5];
	z[91] = -0.5*z[12]*z[6];
	z[6] = z[12]*z[6];
	z[12] = 0.5*epsilon*ex*ez*J*z[6];
	z[92] = 0.5*epsilon*ex*ey*J*z[24];
	z[93] = 0.5*epsilon*ey*ez*J*z[34];
	z[94] = 0.5*epsilon*ex*ey*J*z[36];
	z[95] = 0.5*epsilon*J*z[19]*z[38];
	z[96] = 0.5*epsilon*ex*ez*J*z[40];
	z[97] = 0.5*epsilon*ey*ez*J*z[44];
	z[98] = 0.5*epsilon*J*z[21]*z[48];
	z[99] = 0.5*epsilon*ex*ez*J*z[51];
	z[100] = 0.5*epsilon*J*z[20]*z[53];
	z[101] = 0.5*epsilon*J*z[20]*z[59];
	z[102] = 0.5*epsilon*ey*ez*J*z[63];
	z[103] = 0.5*epsilon*ex*ez*J*z[67];
	z[104] = 0.5*epsilon*J*z[21]*z[69];
	z[105] = 0.5*epsilon*ey*ez*J*z[71];
	z[106] = z[74] + z[84];
	z[107] = z[76] + z[88];
	z[108] = z[4] + z[75];
	z[14] = z[14] + z[23];
	z[22] = z[22] + z[24];
	z[23] = z[25] + z[74];
	z[16] = z[16] + z[27];
	z[11] = z[11] + z[27];
	z[24] = z[28] + z[80];
	z[27] = z[28] + z[82];
	z[28] = z[29] + z[76];
	z[74] = z[30] + z[75];
	z[75] = z[3] + z[31];
	z[76] = z[32] + z[87];
	z[6] = z[35] + z[6];
	z[35] = z[36] + z[91];
	z[36] = z[1] + z[41];
	z[41] = z[26] + z[41];
	z[80] = z[42] + z[84];
	z[25] = z[25] + z[42];
	z[2] = z[2] + z[43];
	z[42] = z[44] + z[85];
	z[3] = z[3] + z[45];
	z[43] = z[46] + z[87];
	z[34] = z[34] + z[47];
	z[33] = z[33] + z[48];
	z[5] = z[5] + z[52];
	z[44] = z[53] + z[89];
	z[7] = z[54] + z[7];
	z[8] = z[54] + z[8];
	z[47] = z[55] + z[77];
	z[48] = z[55] + z[78];
	z[15] = z[15] + z[56];
	z[17] = z[17] + z[56];
	z[52] = z[57] + z[79];
	z[53] = z[57] + z[81];
	z[4] = z[4] + z[60];
	z[30] = z[30] + z[60];
	z[54] = z[61] + z[88];
	z[29] = z[29] + z[61];
	z[55] = z[59] + z[62];
	z[56] = z[58] + z[63];
	z[32] = z[32] + z[64];
	z[46] = z[46] + z[64];
	z[31] = z[31] + z[65];
	z[45] = z[45] + z[65];
	z[38] = z[38] + z[66];
	z[37] = z[37] + z[67];
	z[40] = z[40] + z[68];
	z[39] = z[39] + z[69];
	z[51] = z[51] + z[70];
	z[50] = z[50] + z[71];
	z[1] = z[1] + z[73];
	z[26] = z[26] + z[73];
	z[57] = epsilon*J;
	z[58] = ey*z[57];
	z[59] = ez*z[57];
	z[14] = z[14]*z[19]*z[57];
	z[16] = z[16]*z[19]*z[57];
	z[60] = z[19]*z[57]*z[75];
	z[3] = z[19]*z[3]*z[57];
	z[7] = z[19]*z[57]*z[7];
	z[19] = z[19]*z[38]*z[57];
	z[38] = z[20]*z[57]*z[74];
	z[44] = z[20]*z[44]*z[57];
	z[8] = z[20]*z[57]*z[8];
	z[15] = z[15]*z[20]*z[57];
	z[30] = z[20]*z[30]*z[57];
	z[20] = z[20]*z[55]*z[57];
	z[55] = ez*z[58];
	z[22] = ex*z[22]*z[58];
	z[61] = ex*z[58]*z[76];
	z[35] = ex*z[35]*z[58];
	z[62] = ex*z[58]*z[80];
	z[25] = ex*z[25]*z[58];
	z[2] = ex*z[2]*z[58];
	z[5] = ex*z[5]*z[58];
	z[47] = ex*z[47]*z[58];
	z[48] = ex*z[48]*z[58];
	z[54] = ex*z[54]*z[58];
	z[29] = ex*z[29]*z[58];
	z[32] = ex*z[32]*z[58];
	z[58] = ex*z[108]*z[59];
	z[24] = ex*z[24]*z[59];
	z[27] = ex*z[27]*z[59];
	z[6] = ex*z[59]*z[6];
	z[36] = ex*z[36]*z[59];
	z[43] = ex*z[43]*z[59];
	z[4] = ex*z[4]*z[59];
	z[46] = ex*z[46]*z[59];
	z[37] = ex*z[37]*z[59];
	z[40] = ex*z[40]*z[59];
	z[51] = ex*z[51]*z[59];
	z[1] = ex*z[1]*z[59];
	z[59] = z[106]*z[21]*z[57];
	z[23] = z[21]*z[23]*z[57];
	z[11] = z[11]*z[21]*z[57];
	z[33] = z[21]*z[33]*z[57];
	z[17] = z[17]*z[21]*z[57];
	z[21] = z[21]*z[39]*z[57];
	z[39] = z[107]*z[55];
	z[28] = z[28]*z[55];
	z[41] = z[41]*z[55];
	z[42] = z[42]*z[55];
	z[34] = z[34]*z[55];
	z[52] = z[52]*z[55];
	z[53] = z[53]*z[55];
	z[56] = z[55]*z[56];
	z[31] = z[31]*z[55];
	z[45] = z[45]*z[55];
	z[50] = z[50]*z[55];
	z[26] = z[26]*z[55];
	z[6] = z[10] + z[100] + z[14] + z[23] + z[28] + z[46] + z[48] + z[6] + z[97];
	z[10] = z[10] + z[12] + z[39] + z[42] + z[43] + z[44] + z[47] + z[59] + z[83];
	z[12] = z[13] + z[15] + z[16] + z[25] + z[34] + z[40] + z[62] + z[93] + z[96];
	z[4] = z[102] + z[17] + z[18] + z[4] + z[5] + z[56] + z[58] + z[7] + z[90];
	z[5] = z[103] + z[11] + z[22] + z[31] + z[37] + z[45] + z[8] + z[9] + z[92];
	z[7] = z[101] + z[29] + z[33] + z[36] + z[49] + z[51] + z[52] + z[60] + z[86];
	z[1] = z[1] + z[2] + z[20] + z[3] + z[49] + z[53] + z[54] + z[98] + z[99];
	z[2] = z[21] + z[24] + z[30] + z[41] + z[50] + z[61] + z[72] + z[94] + z[95];
	z[3] = z[104] + z[105] + z[19] + z[26] + z[27] + z[32] + z[35] + z[38] + z[72];

	/* return values */
	dUdCmechelec[0] = z[5];
	dUdCmechelec[1] = z[6];
	dUdCmechelec[2] = z[3];
	dUdCmechelec[3] = z[10];
	dUdCmechelec[4] = z[4];
	dUdCmechelec[5] = z[1];
	dUdCmechelec[6] = z[2];
	dUdCmechelec[7] = z[7];
	dUdCmechelec[8] = z[12];	
}