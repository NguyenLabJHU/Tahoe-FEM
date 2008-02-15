#include "Tersoff_inc.h"

#include <math.h>

static double z[197];

/* function to compute bulk strain energy density for silicon */
double get_energy_surf_bulk(const double* params, const double *Xsi, const double *Xa, const double *Ya, const double *Za, const double* Cmat) { 

/* common definitions */
#include "Tersoff_common_defines.h"

	z[1] = c*c;
	z[2] = 1./(d*d);
	z[3] = d*d;
	z[4] = 1./n;
	z[5] = -X2;
	z[6] = -X3;
	z[7] = -X4;
	z[8] = -X5;
	z[9] = -Y2;
	z[10] = -Y3;
	z[11] = -Y4;
	z[12] = -Y5;
	z[13] = -Z2;
	z[14] = -Z3;
	z[15] = -Z4;
	z[16] = -Z5;
	z[17] = Y2 + z[10];
	z[18] = Y4 + z[10];
	z[19] = Y5 + z[10];
	z[10] = Y1 + Ys1 + z[10];
	z[20] = Y2 + z[11];
	z[21] = Y3 + z[11];
	z[22] = Y5 + z[11];
	z[11] = Y1 + Ys1 + z[11];
	z[23] = Y2 + z[12];
	z[24] = Y3 + z[12];
	z[25] = Y4 + z[12];
	z[12] = Y1 + Ys1 + z[12];
	z[2] = z[1]*z[2];
	z[26] = z[14] + Z2;
	z[27] = z[15] + Z2;
	z[28] = z[16] + Z2;
	z[29] = z[13] + Z3;
	z[30] = z[15] + Z3;
	z[31] = z[16] + Z3;
	z[4] = -0.5*z[4];
	z[32] = z[13] + Z4;
	z[33] = z[14] + Z4;
	z[34] = z[16] + Z4;
	z[35] = X3 + z[5];
	z[36] = X4 + z[5];
	z[37] = X5 + z[5];
	z[5] = X1 + Xs1 + z[5];
	z[38] = z[13] + Z5;
	z[39] = z[14] + Z5;
	z[40] = z[15] + Z5;
	z[41] = X2 + z[6];
	z[42] = X4 + z[6];
	z[43] = X5 + z[6];
	z[6] = X1 + Xs1 + z[6];
	z[44] = X2 + z[7];
	z[45] = X3 + z[7];
	z[46] = X5 + z[7];
	z[7] = X1 + Xs1 + z[7];
	z[47] = X2 + z[8];
	z[48] = X3 + z[8];
	z[49] = X4 + z[8];
	z[8] = X1 + Xs1 + z[8];
	z[50] = Y3 + z[9];
	z[51] = Y4 + z[9];
	z[52] = Y5 + z[9];
	z[9] = Y1 + Ys1 + z[9];
	z[13] = Z1 + z[13] + Zs1;
	z[14] = Z1 + z[14] + Zs1;
	z[15] = Z1 + z[15] + Zs1;
	z[16] = Z1 + z[16] + Zs1;
	z[53] = C21*z[17];
	z[54] = C22*z[17];
	z[55] = C23*z[17];
	z[56] = C21*z[18];
	z[57] = C22*z[18];
	z[58] = C23*z[18];
	z[59] = C21*z[19];
	z[60] = C22*z[19];
	z[61] = C23*z[19];
	z[62] = C21*z[10];
	z[63] = C22*z[10];
	z[64] = C23*z[10];
	z[65] = C21*z[20];
	z[66] = C22*z[20];
	z[67] = C23*z[20];
	z[68] = C21*z[21];
	z[69] = C22*z[21];
	z[70] = C23*z[21];
	z[71] = C21*z[22];
	z[72] = C22*z[22];
	z[73] = C23*z[22];
	z[74] = C21*z[11];
	z[75] = C22*z[11];
	z[76] = C23*z[11];
	z[77] = C21*z[23];
	z[78] = C22*z[23];
	z[79] = C23*z[23];
	z[80] = C21*z[24];
	z[81] = C22*z[24];
	z[82] = C23*z[24];
	z[83] = C21*z[25];
	z[84] = C22*z[25];
	z[85] = C23*z[25];
	z[86] = C21*z[12];
	z[87] = C22*z[12];
	z[88] = C23*z[12];
	z[89] = C31*z[26];
	z[90] = C32*z[26];
	z[91] = C33*z[26];
	z[92] = C31*z[27];
	z[93] = C32*z[27];
	z[94] = C33*z[27];
	z[95] = C31*z[28];
	z[96] = C32*z[28];
	z[97] = C33*z[28];
	z[98] = C31*z[29];
	z[99] = C32*z[29];
	z[100] = C33*z[29];
	z[101] = C31*z[30];
	z[102] = C32*z[30];
	z[103] = C33*z[30];
	z[104] = C31*z[31];
	z[105] = C32*z[31];
	z[106] = C33*z[31];
	z[107] = C31*z[32];
	z[108] = C32*z[32];
	z[109] = C33*z[32];
	z[110] = C31*z[33];
	z[111] = C32*z[33];
	z[112] = C33*z[33];
	z[113] = C31*z[34];
	z[114] = C32*z[34];
	z[115] = C33*z[34];
	z[116] = C11*z[35];
	z[117] = C12*z[35];
	z[118] = C13*z[35];
	z[119] = C11*z[36];
	z[120] = C12*z[36];
	z[121] = C13*z[36];
	z[122] = C11*z[37];
	z[123] = C12*z[37];
	z[124] = C13*z[37];
	z[125] = C11*z[5];
	z[126] = C12*z[5];
	z[127] = C13*z[5];
	z[128] = C31*z[38];
	z[129] = C32*z[38];
	z[130] = C33*z[38];
	z[131] = C31*z[39];
	z[132] = C32*z[39];
	z[133] = C33*z[39];
	z[134] = C31*z[40];
	z[135] = C32*z[40];
	z[136] = C33*z[40];
	z[137] = C11*z[41];
	z[138] = C12*z[41];
	z[139] = C13*z[41];
	z[140] = C11*z[42];
	z[141] = C12*z[42];
	z[142] = C13*z[42];
	z[143] = C11*z[43];
	z[144] = C12*z[43];
	z[145] = C13*z[43];
	z[146] = C11*z[6];
	z[147] = C12*z[6];
	z[148] = C13*z[6];
	z[149] = C11*z[44];
	z[150] = C12*z[44];
	z[151] = C13*z[44];
	z[152] = C11*z[45];
	z[153] = C12*z[45];
	z[154] = C13*z[45];
	z[155] = C11*z[46];
	z[156] = C12*z[46];
	z[157] = C13*z[46];
	z[158] = C11*z[7];
	z[159] = C12*z[7];
	z[160] = C13*z[7];
	z[161] = C11*z[47];
	z[162] = C12*z[47];
	z[163] = C13*z[47];
	z[164] = C11*z[48];
	z[165] = C12*z[48];
	z[166] = C13*z[48];
	z[167] = C11*z[49];
	z[168] = C12*z[49];
	z[169] = C13*z[49];
	z[170] = C11*z[8];
	z[171] = C12*z[8];
	z[172] = C13*z[8];
	z[173] = C21*z[50];
	z[174] = C22*z[50];
	z[175] = C23*z[50];
	z[176] = C21*z[51];
	z[177] = C22*z[51];
	z[178] = C23*z[51];
	z[179] = C21*z[52];
	z[180] = C22*z[52];
	z[181] = C23*z[52];
	z[182] = C21*z[9];
	z[183] = C22*z[9];
	z[184] = C23*z[9];
	z[185] = C31*z[13];
	z[186] = C32*z[13];
	z[187] = C33*z[13];
	z[188] = C31*z[14];
	z[189] = C32*z[14];
	z[190] = C33*z[14];
	z[191] = C31*z[15];
	z[192] = C32*z[15];
	z[193] = C33*z[15];
	z[194] = C31*z[16];
	z[195] = C32*z[16];
	z[196] = C33*z[16];
	z[98] = z[116] + z[173] + z[98];
	z[99] = z[117] + z[174] + z[99];
	z[100] = z[100] + z[118] + z[175];
	z[107] = z[107] + z[119] + z[176];
	z[108] = z[108] + z[120] + z[177];
	z[109] = z[109] + z[121] + z[178];
	z[116] = z[122] + z[128] + z[179];
	z[117] = z[123] + z[129] + z[180];
	z[118] = z[124] + z[130] + z[181];
	z[119] = z[125] + z[182] + z[185];
	z[120] = z[126] + z[183] + z[186];
	z[121] = z[127] + z[184] + z[187];
	z[86] = z[170] + z[194] + z[86];
	z[87] = z[171] + z[195] + z[87];
	z[88] = z[172] + z[196] + z[88];
	z[53] = z[137] + z[53] + z[89];
	z[54] = z[138] + z[54] + z[90];
	z[55] = z[139] + z[55] + z[91];
	z[56] = z[110] + z[140] + z[56];
	z[57] = z[111] + z[141] + z[57];
	z[58] = z[112] + z[142] + z[58];
	z[59] = z[131] + z[143] + z[59];
	z[60] = z[132] + z[144] + z[60];
	z[61] = z[133] + z[145] + z[61];
	z[62] = z[146] + z[188] + z[62];
	z[63] = z[147] + z[189] + z[63];
	z[64] = z[148] + z[190] + z[64];
	z[65] = z[149] + z[65] + z[92];
	z[66] = z[150] + z[66] + z[93];
	z[67] = z[151] + z[67] + z[94];
	z[68] = z[101] + z[152] + z[68];
	z[69] = z[102] + z[153] + z[69];
	z[70] = z[103] + z[154] + z[70];
	z[71] = z[134] + z[155] + z[71];
	z[72] = z[135] + z[156] + z[72];
	z[73] = z[136] + z[157] + z[73];
	z[74] = z[158] + z[191] + z[74];
	z[75] = z[159] + z[192] + z[75];
	z[76] = z[160] + z[193] + z[76];
	z[77] = z[161] + z[77] + z[95];
	z[78] = z[162] + z[78] + z[96];
	z[79] = z[163] + z[79] + z[97];
	z[80] = z[104] + z[164] + z[80];
	z[81] = z[105] + z[165] + z[81];
	z[82] = z[106] + z[166] + z[82];
	z[83] = z[113] + z[167] + z[83];
	z[84] = z[114] + z[168] + z[84];
	z[85] = z[115] + z[169] + z[85];
	z[17] = -z[17]*z[54];
	z[18] = -z[18]*z[57];
	z[19] = -z[19]*z[60];
	z[10] = z[10]*z[63];
	z[20] = -z[20]*z[66];
	z[21] = -z[21]*z[69];
	z[22] = -z[22]*z[72];
	z[11] = z[11]*z[75];
	z[23] = -z[23]*z[78];
	z[24] = -z[24]*z[81];
	z[25] = -z[25]*z[84];
	z[12] = z[12]*z[87];
	z[26] = -z[26]*z[55];
	z[27] = -z[27]*z[67];
	z[28] = -z[28]*z[79];
	z[29] = -z[100]*z[29];
	z[30] = -z[30]*z[70];
	z[31] = -z[31]*z[82];
	z[32] = -z[109]*z[32];
	z[33] = -z[33]*z[58];
	z[34] = -z[34]*z[85];
	z[35] = -z[35]*z[98];
	z[36] = -z[107]*z[36];
	z[37] = -z[116]*z[37];
	z[5] = z[119]*z[5];
	z[38] = -z[118]*z[38];
	z[39] = -z[39]*z[61];
	z[40] = -z[40]*z[73];
	z[41] = -z[41]*z[53];
	z[42] = -z[42]*z[56];
	z[43] = -z[43]*z[59];
	z[6] = z[6]*z[62];
	z[44] = -z[44]*z[65];
	z[45] = -z[45]*z[68];
	z[46] = -z[46]*z[71];
	z[7] = z[7]*z[74];
	z[47] = -z[47]*z[77];
	z[48] = -z[48]*z[80];
	z[49] = -z[49]*z[83];
	z[8] = z[8]*z[86];
	z[50] = -z[50]*z[99];
	z[51] = -z[108]*z[51];
	z[52] = -z[117]*z[52];
	z[9] = z[120]*z[9];
	z[13] = z[121]*z[13];
	z[14] = z[14]*z[64];
	z[15] = z[15]*z[76];
	z[16] = z[16]*z[88];
	z[5] = z[13] + z[5] + z[9];
	z[6] = z[10] + z[14] + z[6];
	z[7] = z[11] + z[15] + z[7];
	z[8] = z[12] + z[16] + z[8];
	z[9] = z[17] + z[26] + z[41] + z[5] + z[6];
	z[10] = z[29] + z[35] + z[5] + z[50] + z[6];
	z[11] = z[20] + z[27] + z[44] + z[5] + z[7];
	z[12] = z[32] + z[36] + z[5] + z[51] + z[7];
	z[13] = z[18] + z[33] + z[42] + z[6] + z[7];
	z[14] = z[21] + z[30] + z[45] + z[6] + z[7];
	z[15] = z[23] + z[28] + z[47] + z[5] + z[8];
	z[16] = z[37] + z[38] + z[5] + z[52] + z[8];
	z[17] = z[19] + z[39] + z[43] + z[6] + z[8];
	z[18] = z[24] + z[31] + z[48] + z[6] + z[8];
	z[19] = z[22] + z[40] + z[46] + z[7] + z[8];
	z[20] = z[25] + z[34] + z[49] + z[7] + z[8];
	z[21] = 1./sqrt(z[5]);
	z[5] = sqrt(z[5]);
	z[22] = 1./sqrt(z[6]);
	z[6] = sqrt(z[6]);
	z[23] = 1./sqrt(z[7]);
	z[7] = sqrt(z[7]);
	z[24] = 1./sqrt(z[8]);
	z[8] = sqrt(z[8]);
	z[5] = -z[5];
	z[9] = -0.5*z[21]*z[22]*z[9];
	z[10] = -0.5*z[10]*z[21]*z[22];
	z[25] = -lam*z[6];
	z[6] = -mu*z[6];
	z[11] = -0.5*z[11]*z[21]*z[23];
	z[12] = -0.5*z[12]*z[21]*z[23];
	z[13] = -0.5*z[13]*z[22]*z[23];
	z[14] = -0.5*z[14]*z[22]*z[23];
	z[26] = -lam*z[7];
	z[7] = -mu*z[7];
	z[15] = -0.5*z[15]*z[21]*z[24];
	z[16] = -0.5*z[16]*z[21]*z[24];
	z[17] = -0.5*z[17]*z[22]*z[24];
	z[18] = -0.5*z[18]*z[22]*z[24];
	z[19] = -0.5*z[19]*z[23]*z[24];
	z[20] = -0.5*z[20]*z[23]*z[24];
	z[21] = -lam*z[8];
	z[8] = -mu*z[8];
	z[22] = lam*z[5];
	z[5] = mu*z[5];
	z[23] = exp(z[25]);
	z[6] = exp(z[6]);
	z[24] = exp(z[26]);
	z[7] = exp(z[7]);
	z[21] = exp(z[21]);
	z[8] = exp(z[8]);
	z[22] = exp(z[22]);
	z[5] = exp(z[5]);
	z[9] = h + z[9];
	z[10] = h + z[10];
	z[11] = h + z[11];
	z[12] = h + z[12];
	z[13] = h + z[13];
	z[14] = h + z[14];
	z[15] = h + z[15];
	z[16] = h + z[16];
	z[17] = h + z[17];
	z[18] = h + z[18];
	z[19] = h + z[19];
	z[20] = h + z[20];
	z[23] = A*z[23];
	z[24] = A*z[24];
	z[21] = A*z[21];
	z[22] = A*z[22];
	z[9] = z[9]*z[9];
	z[10] = z[10]*z[10];
	z[11] = z[11]*z[11];
	z[12] = z[12]*z[12];
	z[13] = z[13]*z[13];
	z[14] = z[14]*z[14];
	z[15] = z[15]*z[15];
	z[16] = z[16]*z[16];
	z[17] = z[17]*z[17];
	z[18] = z[18]*z[18];
	z[19] = z[19]*z[19];
	z[20] = z[20]*z[20];
	z[9] = z[3] + z[9];
	z[10] = z[10] + z[3];
	z[11] = z[11] + z[3];
	z[12] = z[12] + z[3];
	z[13] = z[13] + z[3];
	z[14] = z[14] + z[3];
	z[15] = z[15] + z[3];
	z[16] = z[16] + z[3];
	z[17] = z[17] + z[3];
	z[18] = z[18] + z[3];
	z[19] = z[19] + z[3];
	z[3] = z[20] + z[3];
	z[9] = 1./z[9];
	z[10] = 1./z[10];
	z[11] = 1./z[11];
	z[12] = 1./z[12];
	z[13] = 1./z[13];
	z[14] = 1./z[14];
	z[15] = 1./z[15];
	z[16] = 1./z[16];
	z[17] = 1./z[17];
	z[18] = 1./z[18];
	z[19] = 1./z[19];
	z[3] = 1./z[3];
	z[1] = -z[1];
	z[9] = z[1]*z[9];
	z[10] = z[1]*z[10];
	z[11] = z[1]*z[11];
	z[12] = z[1]*z[12];
	z[13] = z[1]*z[13];
	z[14] = z[1]*z[14];
	z[15] = z[1]*z[15];
	z[16] = z[1]*z[16];
	z[17] = z[1]*z[17];
	z[18] = z[1]*z[18];
	z[19] = z[1]*z[19];
	z[1] = z[1]*z[3];
	z[2] = 1. + z[2];
	z[3] = z[2] + z[9];
	z[9] = z[10] + z[2];
	z[10] = z[11] + z[2];
	z[11] = z[12] + z[2];
	z[12] = z[13] + z[2];
	z[13] = z[14] + z[2];
	z[14] = z[15] + z[2];
	z[15] = z[16] + z[2];
	z[16] = z[17] + z[2];
	z[17] = z[18] + z[2];
	z[18] = z[19] + z[2];
	z[1] = z[1] + z[2];
	z[2] = beta*z[3];
	z[3] = beta*z[9];
	z[9] = beta*z[10];
	z[10] = beta*z[11];
	z[11] = beta*z[12];
	z[12] = beta*z[13];
	z[13] = beta*z[14];
	z[14] = beta*z[15];
	z[15] = beta*z[16];
	z[16] = beta*z[17];
	z[17] = beta*z[18];
	z[1] = beta*z[1];
	z[2] = z[13] + z[2] + z[9];
	z[3] = z[12] + z[16] + z[3];
	z[9] = z[14] + z[15] + z[17];
	z[1] = z[1] + z[10] + z[11];
	z[2] = pow(z[2],n);
	z[3] = pow(z[3],n);
	z[9] = pow(z[9],n);
	z[1] = pow(z[1],n);
	z[2] = 1. + z[2];
	z[3] = 1. + z[3];
	z[9] = 1. + z[9];
	z[1] = 1. + z[1];
	z[2] = pow(z[2],z[4]);
	z[3] = pow(z[3],z[4]);
	z[9] = pow(z[9],z[4]);
	z[1] = pow(z[1],z[4]);
	z[4] = -B*chi;
	z[2] = z[2]*z[4]*z[5];
	z[3] = z[3]*z[4]*z[6];
	z[5] = z[4]*z[8]*z[9];
	z[1] = z[1]*z[4]*z[7];
	z[1] = z[1] + z[2] + z[21] + z[22] + z[23] + z[24] + z[3] + z[5];
	z[1] = 0.5*z[1];

	/* return values */
	return z[1];
}