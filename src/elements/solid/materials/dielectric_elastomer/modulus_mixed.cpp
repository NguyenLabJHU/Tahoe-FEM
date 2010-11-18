#include "FSDE_inc.h"

#include <math.h>

static double z[120];

/* function to compute mixed electromechanical modulus */
void get_ddCE(const double* params, const double *Xsi, const double* Cmat, double* dCdXsi) { 

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
	z[13] = 1./z[16];
	z[15] = z[1]*z[1];
	z[16] = z[2]*z[2];
	z[17] = z[3]*z[3];
	z[18] = z[4]*z[4];
	z[19] = z[5]*z[5];
	z[20] = z[6]*z[6];
	z[21] = -z[11]*z[12];
	z[12] = z[11]*z[12];
	z[22] = -z[10]*z[11];
	z[23] = ex*z[10]*z[11];
	z[24] = -2.*ez*z[1]*z[10]*z[11];
	z[25] = -2.*ez*z[1]*z[11]*z[14];
	z[26] = -ey*z[11]*z[14]*z[2];
	z[27] = -ez*z[11]*z[14]*z[2];
	z[28] = -ey*z[1]*z[11]*z[2];
	z[29] = -2.*ez*z[1]*z[11]*z[2];
	z[30] = -ex*z[11]*z[14]*z[3];
	z[31] = -ez*z[11]*z[14]*z[3];
	z[32] = -ex*z[1]*z[11]*z[3];
	z[33] = -2.*ez*z[1]*z[11]*z[3];
	z[34] = -ex*z[11]*z[2]*z[3];
	z[35] = -ey*z[11]*z[2]*z[3];
	z[36] = -ey*z[11]*z[14]*z[4];
	z[37] = -ez*z[11]*z[14]*z[4];
	z[38] = -ey*z[1]*z[11]*z[4];
	z[39] = -2.*ez*z[1]*z[11]*z[4];
	z[40] = -ey*z[11]*z[2]*z[4];
	z[41] = -ez*z[11]*z[2]*z[4];
	z[42] = -ex*z[11]*z[3]*z[4];
	z[43] = -ey*z[11]*z[3]*z[4];
	z[44] = -2.*ey*z[10]*z[11]*z[5];
	z[45] = -2.*ey*z[11]*z[14]*z[5];
	z[46] = -2.*ez*z[1]*z[11]*z[5];
	z[47] = -ey*z[11]*z[2]*z[5];
	z[48] = -ez*z[11]*z[2]*z[5];
	z[49] = -ex*z[11]*z[3]*z[5];
	z[50] = -2.*ey*z[11]*z[4]*z[5];
	z[51] = -ey*z[11]*z[4]*z[5];
	z[52] = -ez*z[11]*z[4]*z[5];
	z[53] = -2.*ex*z[11]*z[14]*z[6];
	z[54] = -ex*z[11]*z[14]*z[6];
	z[55] = -ez*z[11]*z[14]*z[6];
	z[56] = -ex*z[1]*z[11]*z[6];
	z[57] = -2.*ez*z[1]*z[11]*z[6];
	z[58] = -ex*z[11]*z[2]*z[6];
	z[59] = -ey*z[11]*z[2]*z[6];
	z[60] = -ez*z[11]*z[2]*z[6];
	z[61] = -ex*z[11]*z[3]*z[6];
	z[62] = -ez*z[11]*z[3]*z[6];
	z[63] = -ex*z[11]*z[4]*z[6];
	z[64] = -ey*z[11]*z[4]*z[6];
	z[65] = -ez*z[11]*z[4]*z[6];
	z[66] = -ex*z[11]*z[5]*z[6];
	z[67] = -2.*ey*z[11]*z[5]*z[6];
	z[7] = -ex*z[11]*z[7];
	z[68] = -ex*z[11]*z[8];
	z[8] = -ey*z[11]*z[8];
	z[9] = -2.*ex*z[11]*z[9];
	z[69] = C12*ex*z[13];
	z[70] = -C13*ex*z[13];
	z[71] = C21*ex*z[13];
	z[72] = -C22*ex*z[13];
	z[73] = C23*ex*z[13];
	z[74] = -C31*ex*z[13];
	z[75] = C31*ex*z[13];
	z[76] = C32*ex*z[13];
	z[77] = -C33*ex*z[13];
	z[78] = -C11*ey*z[13];
	z[79] = C12*ey*z[13];
	z[80] = C13*ey*z[13];
	z[81] = C21*ey*z[13];
	z[82] = -C23*ey*z[13];
	z[83] = -2.*C31*ey*z[13];
	z[84] = C31*ey*z[13];
	z[85] = -C32*ey*z[13];
	z[86] = C32*ey*z[13];
	z[87] = -C33*ey*z[13];
	z[88] = 2.*C33*ey*z[13];
	z[89] = -C11*ez*z[13];
	z[90] = 2.*C11*ez*z[13];
	z[91] = -2.*C12*ez*z[13];
	z[92] = C13*ez*z[13];
	z[93] = -2.*C21*ez*z[13];
	z[94] = C21*ez*z[13];
	z[95] = -C22*ez*z[13];
	z[96] = 2.*C22*ez*z[13];
	z[97] = -C23*ez*z[13];
	z[98] = C23*ez*z[13];
	z[99] = C31*ez*z[13];
	z[13] = -C32*ez*z[13];
	z[15] = -2.*ez*z[11]*z[15];
	z[16] = -ey*z[11]*z[16];
	z[17] = -ex*z[11]*z[17];
	z[100] = -ey*z[11]*z[18];
	z[18] = -ez*z[11]*z[18];
	z[19] = -2.*ey*z[11]*z[19];
	z[101] = -ex*z[11]*z[20];
	z[11] = -ez*z[11]*z[20];
	z[20] = ex*z[21];
	z[102] = z[10]*z[20];
	z[103] = z[14]*z[20];
	z[104] = z[20]*z[3];
	z[105] = z[20]*z[4];
	z[106] = z[20]*z[5];
	z[20] = z[20]*z[6];
	z[10] = z[10]*z[21];
	z[10] = ey*z[10];
	z[107] = ey*z[14]*z[21];
	z[108] = ey*z[2]*z[21];
	z[109] = ez*z[2]*z[21];
	z[110] = ey*z[21]*z[4];
	z[111] = ez*z[21]*z[4];
	z[21] = ey*z[21]*z[6];
	z[1] = -2.*ez*z[1]*z[12];
	z[12] = -2.*ey*z[12]*z[5];
	z[112] = ex*z[22];
	z[113] = z[112]*z[14];
	z[114] = z[112]*z[3];
	z[115] = z[112]*z[4];
	z[5] = z[112]*z[5];
	z[112] = z[112]*z[6];
	z[116] = ey*z[14]*z[22];
	z[117] = ey*z[2]*z[22];
	z[2] = ez*z[2]*z[22];
	z[3] = ez*z[22]*z[3];
	z[118] = ey*z[22]*z[4];
	z[4] = ez*z[22]*z[4];
	z[119] = ey*z[22]*z[6];
	z[6] = ez*z[22]*z[6];
	z[14] = -2.*z[14]*z[23];
	z[1] = z[1] + z[104] + z[108] + z[110] + z[20] + z[76] + z[80] + z[91];
	z[22] = z[112] + z[114] + z[117] + z[118] + z[24] + z[73] + z[84] + z[93];
	z[13] = z[103] + z[113] + z[13] + z[27] + z[37] + z[45] + z[88] + z[97];
	z[18] = z[105] + z[115] + z[18] + z[41] + z[50] + z[75] + z[89];
	z[5] = z[106] + z[19] + z[48] + z[5] + z[52];
	z[19] = z[25] + z[26] + z[30] + z[36] + z[54] + z[82] + z[85] + z[96];
	z[15] = z[15] + z[28] + z[32] + z[38] + z[56];
	z[16] = z[16] + z[29] + z[34] + z[40] + z[58] + z[71] + z[78];
	z[17] = z[17] + z[33] + z[35] + z[43] + z[61] + z[72] + z[79];
	z[11] = z[11] + z[119] + z[21] + z[53] + z[62] + z[86] + z[95];
	z[21] = z[101] + z[57] + z[59] + z[61] + z[64] + z[72] + z[81];
	z[20] = z[112] + z[20] + z[60] + z[65] + z[67] + z[76] + z[83] + z[94];
	z[7] = z[102] + z[109] + z[111] + z[12] + z[7] + z[77] + z[92];
	z[2] = z[102] + z[2] + z[4] + z[44] + z[68] + z[77] + z[99];
	z[3] = z[10] + z[14] + z[3] + z[6] + z[8] + z[87] + z[98];
	z[4] = z[107] + z[116] + z[31] + z[55] + z[9];
	z[6] = z[100] + z[39] + z[40] + z[42] + z[63] + z[69] + z[78];
	z[8] = z[46] + z[47] + z[49] + z[51] + z[66] + z[70] + z[74] + z[90];
	z[9] = -0.5*epsilon;
	z[1] = z[1]*z[9];
	z[10] = z[22]*z[9];
	z[12] = z[13]*z[9];
	z[13] = z[18]*z[9];
	z[5] = z[5]*z[9];
	z[14] = z[19]*z[9];
	z[15] = z[15]*z[9];
	z[16] = z[16]*z[9];
	z[17] = z[17]*z[9];
	z[11] = z[11]*z[9];
	z[18] = z[21]*z[9];
	z[19] = z[20]*z[9];
	z[7] = z[7]*z[9];
	z[2] = z[2]*z[9];
	z[3] = z[3]*z[9];
	z[4] = z[4]*z[9];
	z[6] = z[6]*z[9];
	z[8] = z[8]*z[9];


	/* Output stress */
	/* dCdE:  6 x 3 */
	dCdXsi[0] = z[4];
	dCdXsi[1] = z[7];
	dCdXsi[2] = z[17];
	dCdXsi[3] = z[1];
	dCdXsi[4] = z[14];
	dCdXsi[5] = z[12];
	
	dCdXsi[6] = z[3];
	dCdXsi[7] = z[5];
	dCdXsi[8] = z[16];
	dCdXsi[9] = z[8];
	dCdXsi[10] = z[10];
	dCdXsi[11] = z[2];
	
	dCdXsi[12] = z[11];
	dCdXsi[13] = z[13];
	dCdXsi[14] = z[15];
	dCdXsi[15] = z[6];
	dCdXsi[16] = z[18];
	dCdXsi[17] = z[19];	
}