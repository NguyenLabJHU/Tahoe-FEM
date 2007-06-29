/* $Id: dUdC.c,v 1.5 2007-06-29 02:55:07 hspark Exp $ */
#include "Tersoff_inc.h"

#include <math.h>

static double z[296];

/* function to compute derivatives of the potential function wrt to the
 * internal degrees of freedom */
void get_dUdC(const double* params, const double *Xsi, const double *Xa, const double *Ya, const double *Za, const double* Cmat, double* dUdC) {

/* common definitions */
#include "Tersoff_common_defines.h"

	z[1] = c*c;
	z[2] = 1./(d*d);
	z[3] = d*d;
	z[4] = -1. + n;
	z[5] = 1./n;
	z[6] = -X2;
	z[7] = -X3;
	z[8] = -X4;
	z[9] = -X5;
	z[10] = -Y2;
	z[11] = -Y3;
	z[12] = -Y4;
	z[13] = -Y5;
	z[14] = -Z2;
	z[15] = -Z3;
	z[16] = -Z4;
	z[17] = -Z5;
	z[18] = Y3 + z[10];
	z[19] = Y4 + z[10];
	z[20] = Y5 + z[10];
	z[10] = Y1 + Ys1 + z[10];
	z[21] = Y2 + z[11];
	z[22] = Y4 + z[11];
	z[23] = Y5 + z[11];
	z[11] = Y1 + Ys1 + z[11];
	z[24] = Y2 + z[12];
	z[25] = Y3 + z[12];
	z[26] = Y5 + z[12];
	z[12] = Y1 + Ys1 + z[12];
	z[27] = Y2 + z[13];
	z[28] = Y3 + z[13];
	z[29] = Y4 + z[13];
	z[13] = Y1 + Ys1 + z[13];
	z[2] = z[1]*z[2];
	z[30] = z[15] + Z2;
	z[31] = z[16] + Z2;
	z[32] = z[17] + Z2;
	z[33] = z[14] + Z3;
	z[34] = z[16] + Z3;
	z[35] = z[17] + Z3;
	z[36] = z[14] + Z4;
	z[37] = z[15] + Z4;
	z[38] = z[17] + Z4;
	z[5] = -0.5*z[5];
	z[39] = z[14] + Z5;
	z[40] = z[15] + Z5;
	z[41] = z[16] + Z5;
	z[42] = X3 + z[6];
	z[43] = X4 + z[6];
	z[44] = X5 + z[6];
	z[6] = X1 + Xs1 + z[6];
	z[45] = X2 + z[7];
	z[46] = X4 + z[7];
	z[47] = X5 + z[7];
	z[7] = X1 + Xs1 + z[7];
	z[48] = X2 + z[8];
	z[49] = X3 + z[8];
	z[50] = X5 + z[8];
	z[8] = X1 + Xs1 + z[8];
	z[51] = X2 + z[9];
	z[52] = X3 + z[9];
	z[53] = X4 + z[9];
	z[9] = X1 + Xs1 + z[9];
	z[14] = Z1 + z[14] + Zs1;
	z[15] = Z1 + z[15] + Zs1;
	z[16] = Z1 + z[16] + Zs1;
	z[17] = Z1 + z[17] + Zs1;
	z[54] = C21*z[18];
	z[55] = C22*z[18];
	z[56] = C23*z[18];
	z[57] = z[18]*z[18];
	z[58] = C21*z[19];
	z[59] = C22*z[19];
	z[60] = C23*z[19];
	z[61] = z[19]*z[19];
	z[62] = C21*z[20];
	z[63] = C22*z[20];
	z[64] = C23*z[20];
	z[65] = z[20]*z[20];
	z[66] = C21*z[10];
	z[67] = C22*z[10];
	z[68] = C23*z[10];
	z[69] = z[10]*z[10];
	z[70] = C21*z[21];
	z[71] = C22*z[21];
	z[72] = C23*z[21];
	z[73] = z[21]*z[21];
	z[74] = C21*z[22];
	z[75] = C22*z[22];
	z[76] = C23*z[22];
	z[77] = z[22]*z[22];
	z[78] = C21*z[23];
	z[79] = C22*z[23];
	z[80] = C23*z[23];
	z[81] = z[23]*z[23];
	z[82] = C21*z[11];
	z[83] = C22*z[11];
	z[84] = C23*z[11];
	z[85] = z[11]*z[11];
	z[86] = C21*z[24];
	z[87] = C22*z[24];
	z[88] = C23*z[24];
	z[89] = z[24]*z[24];
	z[90] = C21*z[25];
	z[91] = C22*z[25];
	z[92] = C23*z[25];
	z[93] = z[25]*z[25];
	z[94] = C21*z[26];
	z[95] = C22*z[26];
	z[96] = C23*z[26];
	z[97] = z[26]*z[26];
	z[98] = C21*z[12];
	z[99] = C22*z[12];
	z[100] = C23*z[12];
	z[101] = z[12]*z[12];
	z[102] = C21*z[27];
	z[103] = C22*z[27];
	z[104] = C23*z[27];
	z[105] = z[27]*z[27];
	z[106] = C21*z[28];
	z[107] = C22*z[28];
	z[108] = C23*z[28];
	z[109] = z[28]*z[28];
	z[110] = C21*z[29];
	z[111] = C22*z[29];
	z[112] = C23*z[29];
	z[113] = z[29]*z[29];
	z[114] = C21*z[13];
	z[115] = C22*z[13];
	z[116] = C23*z[13];
	z[117] = z[13]*z[13];
	z[118] = C31*z[30];
	z[119] = C32*z[30];
	z[120] = C33*z[30];
	z[121] = -z[21]*z[30];
	z[122] = z[30]*z[30];
	z[123] = C31*z[31];
	z[124] = C32*z[31];
	z[125] = C33*z[31];
	z[126] = -z[24]*z[31];
	z[127] = z[31]*z[31];
	z[128] = C31*z[32];
	z[129] = C32*z[32];
	z[130] = C33*z[32];
	z[131] = -z[27]*z[32];
	z[132] = z[32]*z[32];
	z[133] = C31*z[33];
	z[134] = C32*z[33];
	z[135] = C33*z[33];
	z[136] = -z[18]*z[33];
	z[137] = z[33]*z[33];
	z[138] = C31*z[34];
	z[139] = C32*z[34];
	z[140] = C33*z[34];
	z[141] = -z[25]*z[34];
	z[142] = z[34]*z[34];
	z[143] = C31*z[35];
	z[144] = C32*z[35];
	z[145] = C33*z[35];
	z[146] = -z[28]*z[35];
	z[147] = z[35]*z[35];
	z[148] = C31*z[36];
	z[149] = C32*z[36];
	z[150] = C33*z[36];
	z[151] = -z[19]*z[36];
	z[152] = z[36]*z[36];
	z[153] = C31*z[37];
	z[154] = C32*z[37];
	z[155] = C33*z[37];
	z[156] = -z[22]*z[37];
	z[157] = z[37]*z[37];
	z[158] = C31*z[38];
	z[159] = C32*z[38];
	z[160] = C33*z[38];
	z[161] = -z[29]*z[38];
	z[162] = z[38]*z[38];
	z[163] = -1. + z[5];
	z[164] = C31*z[39];
	z[165] = C32*z[39];
	z[166] = C33*z[39];
	z[167] = -z[20]*z[39];
	z[168] = z[39]*z[39];
	z[169] = C31*z[40];
	z[170] = C32*z[40];
	z[171] = C33*z[40];
	z[172] = -z[23]*z[40];
	z[173] = z[40]*z[40];
	z[174] = C31*z[41];
	z[175] = C32*z[41];
	z[176] = C33*z[41];
	z[177] = -z[26]*z[41];
	z[178] = z[41]*z[41];
	z[179] = C11*z[42];
	z[180] = C12*z[42];
	z[181] = C13*z[42];
	z[182] = -z[18]*z[42];
	z[183] = -z[33]*z[42];
	z[184] = z[42]*z[42];
	z[185] = C11*z[43];
	z[186] = C12*z[43];
	z[187] = C13*z[43];
	z[188] = -z[19]*z[43];
	z[189] = -z[36]*z[43];
	z[190] = z[43]*z[43];
	z[191] = C11*z[44];
	z[192] = C12*z[44];
	z[193] = C13*z[44];
	z[194] = -z[20]*z[44];
	z[195] = -z[39]*z[44];
	z[196] = z[44]*z[44];
	z[197] = C11*z[6];
	z[198] = C12*z[6];
	z[199] = C13*z[6];
	z[200] = z[10]*z[6];
	z[201] = z[6]*z[6];
	z[202] = C11*z[45];
	z[203] = C12*z[45];
	z[204] = C13*z[45];
	z[205] = -z[21]*z[45];
	z[206] = -z[30]*z[45];
	z[207] = z[45]*z[45];
	z[208] = C11*z[46];
	z[209] = C12*z[46];
	z[210] = C13*z[46];
	z[211] = -z[22]*z[46];
	z[212] = -z[37]*z[46];
	z[213] = z[46]*z[46];
	z[214] = C11*z[47];
	z[215] = C12*z[47];
	z[216] = C13*z[47];
	z[217] = -z[23]*z[47];
	z[218] = -z[40]*z[47];
	z[219] = z[47]*z[47];
	z[220] = C11*z[7];
	z[221] = C12*z[7];
	z[222] = C13*z[7];
	z[223] = z[11]*z[7];
	z[224] = z[7]*z[7];
	z[225] = C11*z[48];
	z[226] = C12*z[48];
	z[227] = C13*z[48];
	z[228] = -z[24]*z[48];
	z[229] = -z[31]*z[48];
	z[230] = z[48]*z[48];
	z[231] = C11*z[49];
	z[232] = C12*z[49];
	z[233] = C13*z[49];
	z[234] = -z[25]*z[49];
	z[235] = -z[34]*z[49];
	z[236] = z[49]*z[49];
	z[237] = C11*z[50];
	z[238] = C12*z[50];
	z[239] = C13*z[50];
	z[240] = -z[26]*z[50];
	z[241] = -z[41]*z[50];
	z[242] = z[50]*z[50];
	z[243] = C11*z[8];
	z[244] = C12*z[8];
	z[245] = C13*z[8];
	z[246] = z[12]*z[8];
	z[247] = z[8]*z[8];
	z[248] = C11*z[51];
	z[249] = C12*z[51];
	z[250] = C13*z[51];
	z[251] = -z[27]*z[51];
	z[252] = -z[32]*z[51];
	z[253] = z[51]*z[51];
	z[254] = C11*z[52];
	z[255] = C12*z[52];
	z[256] = C13*z[52];
	z[257] = -z[28]*z[52];
	z[258] = -z[35]*z[52];
	z[259] = z[52]*z[52];
	z[260] = C11*z[53];
	z[261] = C12*z[53];
	z[262] = C13*z[53];
	z[263] = -z[29]*z[53];
	z[264] = -z[38]*z[53];
	z[265] = z[53]*z[53];
	z[266] = C11*z[9];
	z[267] = C12*z[9];
	z[268] = C13*z[9];
	z[269] = z[13]*z[9];
	z[270] = z[9]*z[9];
	z[271] = C31*z[14];
	z[272] = C32*z[14];
	z[273] = C33*z[14];
	z[274] = z[10]*z[14];
	z[275] = z[14]*z[6];
	z[276] = z[14]*z[14];
	z[277] = C31*z[15];
	z[278] = C32*z[15];
	z[279] = C33*z[15];
	z[280] = z[11]*z[15];
	z[281] = z[15]*z[7];
	z[282] = z[15]*z[15];
	z[283] = C31*z[16];
	z[284] = C32*z[16];
	z[285] = C33*z[16];
	z[286] = z[12]*z[16];
	z[287] = z[16]*z[8];
	z[288] = z[16]*z[16];
	z[289] = C31*z[17];
	z[290] = C32*z[17];
	z[291] = C33*z[17];
	z[292] = z[13]*z[17];
	z[293] = z[17]*z[9];
	z[294] = z[17]*z[17];
	z[89] = -z[89];
	z[93] = -z[93];
	z[97] = -z[97];
	z[105] = -z[105];
	z[109] = -z[109];
	z[113] = -z[113];
	z[122] = -z[122];
	z[127] = -z[127];
	z[132] = -z[132];
	z[137] = -z[137];
	z[142] = -z[142];
	z[147] = -z[147];
	z[152] = -z[152];
	z[157] = -z[157];
	z[162] = -z[162];
	z[168] = -z[168];
	z[173] = -z[173];
	z[178] = -z[178];
	z[184] = -z[184];
	z[190] = -z[190];
	z[196] = -z[196];
	z[207] = -z[207];
	z[213] = -z[213];
	z[219] = -z[219];
	z[295] = z[200] + z[223];
	z[86] = z[123] + z[225] + z[86];
	z[87] = z[124] + z[226] + z[87];
	z[88] = z[125] + z[227] + z[88];
	z[123] = -z[230];
	z[90] = z[138] + z[231] + z[90];
	z[91] = z[139] + z[232] + z[91];
	z[92] = z[140] + z[233] + z[92];
	z[124] = -z[236];
	z[94] = z[174] + z[237] + z[94];
	z[95] = z[175] + z[238] + z[95];
	z[96] = z[176] + z[239] + z[96];
	z[125] = -z[242];
	z[138] = z[188] + z[200] + z[246];
	z[139] = z[211] + z[223] + z[246];
	z[140] = z[200] + z[228] + z[246];
	z[174] = z[223] + z[234] + z[246];
	z[102] = z[102] + z[128] + z[248];
	z[103] = z[103] + z[129] + z[249];
	z[104] = z[104] + z[130] + z[250];
	z[128] = -z[253];
	z[106] = z[106] + z[143] + z[254];
	z[107] = z[107] + z[144] + z[255];
	z[108] = z[108] + z[145] + z[256];
	z[129] = -z[259];
	z[110] = z[110] + z[158] + z[260];
	z[111] = z[111] + z[159] + z[261];
	z[112] = z[112] + z[160] + z[262];
	z[130] = -z[265];
	z[143] = z[194] + z[200] + z[269];
	z[144] = z[217] + z[223] + z[269];
	z[145] = z[240] + z[246] + z[269];
	z[158] = z[200] + z[251] + z[269];
	z[159] = z[223] + z[257] + z[269];
	z[160] = z[246] + z[263] + z[269];
	z[121] = z[121] + z[274] + z[280];
	z[136] = z[136] + z[274] + z[280];
	z[175] = z[183] + z[275] + z[281];
	z[176] = z[206] + z[275] + z[281];
	z[98] = z[243] + z[283] + z[98];
	z[99] = z[244] + z[284] + z[99];
	z[100] = z[100] + z[245] + z[285];
	z[126] = z[126] + z[274] + z[286];
	z[151] = z[151] + z[274] + z[286];
	z[141] = z[141] + z[280] + z[286];
	z[156] = z[156] + z[280] + z[286];
	z[183] = z[189] + z[275] + z[287];
	z[188] = z[229] + z[275] + z[287];
	z[189] = z[212] + z[281] + z[287];
	z[194] = z[235] + z[281] + z[287];
	z[114] = z[114] + z[266] + z[289];
	z[115] = z[115] + z[267] + z[290];
	z[116] = z[116] + z[268] + z[291];
	z[131] = z[131] + z[274] + z[292];
	z[167] = z[167] + z[274] + z[292];
	z[146] = z[146] + z[280] + z[292];
	z[172] = z[172] + z[280] + z[292];
	z[161] = z[161] + z[286] + z[292];
	z[177] = z[177] + z[286] + z[292];
	z[195] = z[195] + z[275] + z[293];
	z[206] = z[252] + z[275] + z[293];
	z[211] = z[218] + z[281] + z[293];
	z[212] = z[258] + z[281] + z[293];
	z[217] = z[241] + z[287] + z[293];
	z[218] = z[264] + z[287] + z[293];
	z[182] = z[182] + z[295];
	z[205] = z[205] + z[295];
	z[54] = z[133] + z[179] + z[54];
	z[55] = z[134] + z[180] + z[55];
	z[56] = z[135] + z[181] + z[56];
	z[57] = -z[57];
	z[58] = z[148] + z[185] + z[58];
	z[59] = z[149] + z[186] + z[59];
	z[60] = z[150] + z[187] + z[60];
	z[61] = -z[61];
	z[62] = z[164] + z[191] + z[62];
	z[63] = z[165] + z[192] + z[63];
	z[64] = z[166] + z[193] + z[64];
	z[65] = -z[65];
	z[66] = z[197] + z[271] + z[66];
	z[67] = z[198] + z[272] + z[67];
	z[68] = z[199] + z[273] + z[68];
	z[70] = z[118] + z[202] + z[70];
	z[71] = z[119] + z[203] + z[71];
	z[72] = z[120] + z[204] + z[72];
	z[73] = -z[73];
	z[74] = z[153] + z[208] + z[74];
	z[75] = z[154] + z[209] + z[75];
	z[76] = z[155] + z[210] + z[76];
	z[77] = -z[77];
	z[78] = z[169] + z[214] + z[78];
	z[79] = z[170] + z[215] + z[79];
	z[80] = z[171] + z[216] + z[80];
	z[81] = -z[81];
	z[82] = z[220] + z[277] + z[82];
	z[83] = z[221] + z[278] + z[83];
	z[84] = z[222] + z[279] + z[84];
	z[118] = z[101] + z[117];
	z[119] = z[122] + z[276] + z[282];
	z[120] = z[127] + z[276] + z[288];
	z[122] = z[132] + z[276] + z[294];
	z[127] = z[137] + z[276] + z[282];
	z[132] = z[142] + z[282] + z[288];
	z[133] = z[147] + z[282] + z[294];
	z[134] = z[152] + z[276] + z[288];
	z[135] = z[157] + z[282] + z[288];
	z[137] = z[162] + z[288] + z[294];
	z[142] = z[168] + z[276] + z[294];
	z[147] = z[173] + z[282] + z[294];
	z[148] = z[178] + z[288] + z[294];
	z[149] = z[184] + z[201] + z[224];
	z[150] = z[190] + z[201] + z[247];
	z[152] = z[196] + z[201] + z[270];
	z[153] = z[201] + z[207] + z[224];
	z[154] = z[213] + z[224] + z[247];
	z[155] = z[219] + z[224] + z[270];
	z[24] = -z[24]*z[87];
	z[87] = z[123] + z[201] + z[247];
	z[25] = -z[25]*z[91];
	z[91] = z[124] + z[224] + z[247];
	z[26] = -z[26]*z[95];
	z[95] = z[125] + z[247] + z[270];
	z[27] = -z[103]*z[27];
	z[103] = z[128] + z[201] + z[270];
	z[28] = -z[107]*z[28];
	z[107] = z[129] + z[224] + z[270];
	z[29] = -z[111]*z[29];
	z[31] = -z[31]*z[88];
	z[88] = z[130] + z[247] + z[270];
	z[32] = -z[104]*z[32];
	z[12] = z[12]*z[99];
	z[13] = z[115]*z[13];
	z[34] = -z[34]*z[92];
	z[35] = -z[108]*z[35];
	z[18] = -z[18]*z[55];
	z[33] = -z[33]*z[56];
	z[19] = -z[19]*z[59];
	z[20] = -z[20]*z[63];
	z[36] = -z[36]*z[60];
	z[10] = z[10]*z[67];
	z[21] = -z[21]*z[71];
	z[30] = -z[30]*z[72];
	z[22] = -z[22]*z[75];
	z[37] = -z[37]*z[76];
	z[23] = -z[23]*z[79];
	z[11] = z[11]*z[83];
	z[38] = -z[112]*z[38];
	z[55] = z[118] + z[97];
	z[56] = z[113] + z[118];
	z[39] = -z[39]*z[64];
	z[40] = -z[40]*z[80];
	z[41] = -z[41]*z[96];
	z[42] = -z[42]*z[54];
	z[43] = -z[43]*z[58];
	z[44] = -z[44]*z[62];
	z[6] = z[6]*z[66];
	z[45] = -z[45]*z[70];
	z[46] = -z[46]*z[74];
	z[47] = -z[47]*z[78];
	z[7] = z[7]*z[82];
	z[48] = -z[48]*z[86];
	z[49] = -z[49]*z[90];
	z[50] = -z[50]*z[94];
	z[8] = z[8]*z[98];
	z[51] = -z[102]*z[51];
	z[52] = -z[106]*z[52];
	z[53] = -z[110]*z[53];
	z[9] = z[114]*z[9];
	z[14] = z[14]*z[68];
	z[15] = z[15]*z[84];
	z[16] = z[100]*z[16];
	z[17] = z[116]*z[17];
	z[54] = z[101] + z[69] + z[89];
	z[58] = z[105] + z[117] + z[69];
	z[59] = z[101] + z[61] + z[69];
	z[60] = z[117] + z[65] + z[69];
	z[61] = z[101] + z[85] + z[93];
	z[62] = z[109] + z[117] + z[85];
	z[63] = z[101] + z[77] + z[85];
	z[64] = z[117] + z[81] + z[85];
	z[57] = z[57] + z[69] + z[85];
	z[65] = z[69] + z[73] + z[85];
	z[6] = z[10] + z[14] + z[6];
	z[7] = z[11] + z[15] + z[7];
	z[8] = z[12] + z[16] + z[8];
	z[9] = z[13] + z[17] + z[9];
	z[10] = z[18] + z[33] + z[42] + z[6] + z[7];
	z[11] = z[21] + z[30] + z[45] + z[6] + z[7];
	z[12] = z[19] + z[36] + z[43] + z[6] + z[8];
	z[13] = z[24] + z[31] + z[48] + z[6] + z[8];
	z[14] = z[22] + z[37] + z[46] + z[7] + z[8];
	z[15] = z[25] + z[34] + z[49] + z[7] + z[8];
	z[16] = z[20] + z[39] + z[44] + z[6] + z[9];
	z[17] = z[27] + z[32] + z[51] + z[6] + z[9];
	z[18] = z[23] + z[40] + z[47] + z[7] + z[9];
	z[19] = z[28] + z[35] + z[52] + z[7] + z[9];
	z[20] = z[26] + z[41] + z[50] + z[8] + z[9];
	z[21] = z[29] + z[38] + z[53] + z[8] + z[9];
	z[22] = pow(z[6],-1.5);
	z[23] = 1./sqrt(z[6]);
	z[6] = sqrt(z[6]);
	z[24] = pow(z[7],-1.5);
	z[25] = 1./sqrt(z[7]);
	z[7] = sqrt(z[7]);
	z[26] = pow(z[8],-1.5);
	z[27] = 1./sqrt(z[8]);
	z[8] = sqrt(z[8]);
	z[28] = pow(z[9],-1.5);
	z[29] = 1./sqrt(z[9]);
	z[9] = sqrt(z[9]);
	z[6] = -z[6];
	z[30] = 0.25*z[10]*z[223]*z[23]*z[24];
	z[31] = 0.25*z[10]*z[224]*z[23]*z[24];
	z[32] = 0.25*z[10]*z[23]*z[24]*z[280];
	z[33] = 0.25*z[10]*z[23]*z[24]*z[281];
	z[34] = 0.25*z[10]*z[23]*z[24]*z[282];
	z[35] = 0.25*z[11]*z[223]*z[23]*z[24];
	z[36] = 0.25*z[11]*z[224]*z[23]*z[24];
	z[37] = 0.25*z[11]*z[23]*z[24]*z[280];
	z[38] = 0.25*z[11]*z[23]*z[24]*z[281];
	z[39] = 0.25*z[11]*z[23]*z[24]*z[282];
	z[40] = 0.25*z[10]*z[200]*z[22]*z[25];
	z[41] = 0.25*z[10]*z[201]*z[22]*z[25];
	z[42] = 0.25*z[10]*z[22]*z[25]*z[274];
	z[43] = 0.25*z[10]*z[22]*z[25]*z[275];
	z[44] = 0.25*z[10]*z[22]*z[25]*z[276];
	z[45] = 0.25*z[11]*z[200]*z[22]*z[25];
	z[46] = 0.25*z[11]*z[201]*z[22]*z[25];
	z[47] = 0.25*z[11]*z[22]*z[25]*z[274];
	z[48] = 0.25*z[11]*z[22]*z[25]*z[275];
	z[49] = 0.25*z[11]*z[22]*z[25]*z[276];
	z[50] = -0.5*z[121]*z[23]*z[25];
	z[51] = -0.5*z[136]*z[23]*z[25];
	z[52] = -0.5*z[175]*z[23]*z[25];
	z[53] = -0.5*z[176]*z[23]*z[25];
	z[66] = -0.5*z[182]*z[23]*z[25];
	z[67] = -0.5*z[205]*z[23]*z[25];
	z[68] = -0.5*z[119]*z[23]*z[25];
	z[70] = -0.5*z[127]*z[23]*z[25];
	z[71] = -0.5*z[149]*z[23]*z[25];
	z[72] = -0.5*z[153]*z[23]*z[25];
	z[57] = -0.5*z[23]*z[25]*z[57];
	z[65] = -0.5*z[23]*z[25]*z[65];
	z[73] = -0.5*z[10]*z[23]*z[25];
	z[74] = -0.5*z[11]*z[23]*z[25];
	z[75] = -lam*z[7];
	z[7] = -mu*z[7];
	z[76] = 0.25*z[101]*z[12]*z[23]*z[26];
	z[77] = 0.25*z[12]*z[23]*z[246]*z[26];
	z[78] = 0.25*z[12]*z[23]*z[247]*z[26];
	z[79] = 0.25*z[12]*z[23]*z[26]*z[286];
	z[80] = 0.25*z[12]*z[23]*z[26]*z[287];
	z[81] = 0.25*z[12]*z[23]*z[26]*z[288];
	z[82] = 0.25*z[101]*z[13]*z[23]*z[26];
	z[83] = 0.25*z[13]*z[23]*z[246]*z[26];
	z[84] = 0.25*z[13]*z[23]*z[247]*z[26];
	z[86] = 0.25*z[13]*z[23]*z[26]*z[286];
	z[89] = 0.25*z[13]*z[23]*z[26]*z[287];
	z[90] = 0.25*z[13]*z[23]*z[26]*z[288];
	z[92] = 0.25*z[101]*z[14]*z[25]*z[26];
	z[93] = 0.25*z[14]*z[246]*z[25]*z[26];
	z[94] = 0.25*z[14]*z[247]*z[25]*z[26];
	z[96] = 0.25*z[14]*z[25]*z[26]*z[286];
	z[97] = 0.25*z[14]*z[25]*z[26]*z[287];
	z[98] = 0.25*z[14]*z[25]*z[26]*z[288];
	z[99] = 0.25*z[101]*z[15]*z[25]*z[26];
	z[100] = 0.25*z[15]*z[246]*z[25]*z[26];
	z[102] = 0.25*z[15]*z[247]*z[25]*z[26];
	z[104] = 0.25*z[15]*z[25]*z[26]*z[286];
	z[105] = 0.25*z[15]*z[25]*z[26]*z[287];
	z[106] = 0.25*z[15]*z[25]*z[26]*z[288];
	z[108] = 0.25*z[12]*z[200]*z[22]*z[27];
	z[109] = 0.25*z[12]*z[201]*z[22]*z[27];
	z[110] = 0.25*z[12]*z[22]*z[27]*z[274];
	z[111] = 0.25*z[12]*z[22]*z[27]*z[275];
	z[112] = 0.25*z[12]*z[22]*z[27]*z[276];
	z[113] = 0.25*z[13]*z[200]*z[22]*z[27];
	z[114] = 0.25*z[13]*z[201]*z[22]*z[27];
	z[115] = 0.25*z[13]*z[22]*z[27]*z[274];
	z[116] = 0.25*z[13]*z[22]*z[27]*z[275];
	z[118] = 0.25*z[13]*z[22]*z[27]*z[276];
	z[119] = -0.5*z[138]*z[23]*z[27];
	z[121] = -0.5*z[140]*z[23]*z[27];
	z[123] = -0.5*z[126]*z[23]*z[27];
	z[124] = -0.5*z[151]*z[23]*z[27];
	z[125] = -0.5*z[183]*z[23]*z[27];
	z[126] = -0.5*z[188]*z[23]*z[27];
	z[120] = -0.5*z[120]*z[23]*z[27];
	z[127] = -0.5*z[134]*z[23]*z[27];
	z[128] = -0.5*z[150]*z[23]*z[27];
	z[87] = -0.5*z[23]*z[27]*z[87];
	z[54] = -0.5*z[23]*z[27]*z[54];
	z[59] = -0.5*z[23]*z[27]*z[59];
	z[129] = -0.5*z[12]*z[23]*z[27];
	z[130] = -0.5*z[13]*z[23]*z[27];
	z[134] = 0.25*z[14]*z[223]*z[24]*z[27];
	z[136] = 0.25*z[14]*z[224]*z[24]*z[27];
	z[138] = 0.25*z[14]*z[24]*z[27]*z[280];
	z[140] = 0.25*z[14]*z[24]*z[27]*z[281];
	z[149] = 0.25*z[14]*z[24]*z[27]*z[282];
	z[150] = 0.25*z[15]*z[223]*z[24]*z[27];
	z[151] = 0.25*z[15]*z[224]*z[24]*z[27];
	z[153] = 0.25*z[15]*z[24]*z[27]*z[280];
	z[157] = 0.25*z[15]*z[24]*z[27]*z[281];
	z[162] = 0.25*z[15]*z[24]*z[27]*z[282];
	z[139] = -0.5*z[139]*z[25]*z[27];
	z[164] = -0.5*z[174]*z[25]*z[27];
	z[141] = -0.5*z[141]*z[25]*z[27];
	z[156] = -0.5*z[156]*z[25]*z[27];
	z[165] = -0.5*z[189]*z[25]*z[27];
	z[166] = -0.5*z[194]*z[25]*z[27];
	z[132] = -0.5*z[132]*z[25]*z[27];
	z[135] = -0.5*z[135]*z[25]*z[27];
	z[154] = -0.5*z[154]*z[25]*z[27];
	z[91] = -0.5*z[25]*z[27]*z[91];
	z[61] = -0.5*z[25]*z[27]*z[61];
	z[63] = -0.5*z[25]*z[27]*z[63];
	z[168] = -0.5*z[14]*z[25]*z[27];
	z[169] = -0.5*z[15]*z[25]*z[27];
	z[170] = -lam*z[8];
	z[8] = -mu*z[8];
	z[171] = 0.25*z[117]*z[16]*z[23]*z[28];
	z[173] = 0.25*z[16]*z[23]*z[269]*z[28];
	z[174] = 0.25*z[16]*z[23]*z[270]*z[28];
	z[175] = 0.25*z[16]*z[23]*z[28]*z[292];
	z[176] = 0.25*z[16]*z[23]*z[28]*z[293];
	z[178] = 0.25*z[16]*z[23]*z[28]*z[294];
	z[179] = 0.25*z[117]*z[17]*z[23]*z[28];
	z[180] = 0.25*z[17]*z[23]*z[269]*z[28];
	z[181] = 0.25*z[17]*z[23]*z[270]*z[28];
	z[182] = 0.25*z[17]*z[23]*z[28]*z[292];
	z[183] = 0.25*z[17]*z[23]*z[28]*z[293];
	z[184] = 0.25*z[17]*z[23]*z[28]*z[294];
	z[185] = 0.25*z[117]*z[18]*z[25]*z[28];
	z[186] = 0.25*z[18]*z[25]*z[269]*z[28];
	z[187] = 0.25*z[18]*z[25]*z[270]*z[28];
	z[188] = 0.25*z[18]*z[25]*z[28]*z[292];
	z[189] = 0.25*z[18]*z[25]*z[28]*z[293];
	z[190] = 0.25*z[18]*z[25]*z[28]*z[294];
	z[191] = 0.25*z[117]*z[19]*z[25]*z[28];
	z[192] = 0.25*z[19]*z[25]*z[269]*z[28];
	z[193] = 0.25*z[19]*z[25]*z[270]*z[28];
	z[194] = 0.25*z[19]*z[25]*z[28]*z[292];
	z[196] = 0.25*z[19]*z[25]*z[28]*z[293];
	z[197] = 0.25*z[19]*z[25]*z[28]*z[294];
	z[198] = 0.25*z[117]*z[20]*z[27]*z[28];
	z[199] = 0.25*z[20]*z[269]*z[27]*z[28];
	z[202] = 0.25*z[20]*z[27]*z[270]*z[28];
	z[203] = 0.25*z[20]*z[27]*z[28]*z[292];
	z[204] = 0.25*z[20]*z[27]*z[28]*z[293];
	z[205] = 0.25*z[20]*z[27]*z[28]*z[294];
	z[207] = 0.25*z[117]*z[21]*z[27]*z[28];
	z[208] = 0.25*z[21]*z[269]*z[27]*z[28];
	z[209] = 0.25*z[21]*z[27]*z[270]*z[28];
	z[210] = 0.25*z[21]*z[27]*z[28]*z[292];
	z[213] = 0.25*z[21]*z[27]*z[28]*z[293];
	z[28] = 0.25*z[21]*z[27]*z[28]*z[294];
	z[214] = 0.25*z[16]*z[200]*z[22]*z[29];
	z[215] = 0.25*z[16]*z[201]*z[22]*z[29];
	z[216] = 0.25*z[16]*z[22]*z[274]*z[29];
	z[219] = 0.25*z[16]*z[22]*z[275]*z[29];
	z[220] = 0.25*z[16]*z[22]*z[276]*z[29];
	z[221] = 0.25*z[17]*z[200]*z[22]*z[29];
	z[222] = 0.25*z[17]*z[201]*z[22]*z[29];
	z[225] = 0.25*z[17]*z[22]*z[274]*z[29];
	z[226] = 0.25*z[17]*z[22]*z[275]*z[29];
	z[227] = 0.25*z[17]*z[22]*z[276]*z[29];
	z[143] = -0.5*z[143]*z[23]*z[29];
	z[158] = -0.5*z[158]*z[23]*z[29];
	z[131] = -0.5*z[131]*z[23]*z[29];
	z[167] = -0.5*z[167]*z[23]*z[29];
	z[195] = -0.5*z[195]*z[23]*z[29];
	z[206] = -0.5*z[206]*z[23]*z[29];
	z[122] = -0.5*z[122]*z[23]*z[29];
	z[142] = -0.5*z[142]*z[23]*z[29];
	z[152] = -0.5*z[152]*z[23]*z[29];
	z[103] = -0.5*z[103]*z[23]*z[29];
	z[58] = -0.5*z[23]*z[29]*z[58];
	z[60] = -0.5*z[23]*z[29]*z[60];
	z[228] = -0.5*z[16]*z[23]*z[29];
	z[229] = -0.5*z[17]*z[23]*z[29];
	z[230] = 0.25*z[18]*z[223]*z[24]*z[29];
	z[231] = 0.25*z[18]*z[224]*z[24]*z[29];
	z[232] = 0.25*z[18]*z[24]*z[280]*z[29];
	z[233] = 0.25*z[18]*z[24]*z[281]*z[29];
	z[234] = 0.25*z[18]*z[24]*z[282]*z[29];
	z[235] = 0.25*z[19]*z[223]*z[24]*z[29];
	z[236] = 0.25*z[19]*z[224]*z[24]*z[29];
	z[237] = 0.25*z[19]*z[24]*z[280]*z[29];
	z[238] = 0.25*z[19]*z[24]*z[281]*z[29];
	z[239] = 0.25*z[19]*z[24]*z[282]*z[29];
	z[144] = -0.5*z[144]*z[25]*z[29];
	z[159] = -0.5*z[159]*z[25]*z[29];
	z[146] = -0.5*z[146]*z[25]*z[29];
	z[172] = -0.5*z[172]*z[25]*z[29];
	z[211] = -0.5*z[211]*z[25]*z[29];
	z[212] = -0.5*z[212]*z[25]*z[29];
	z[133] = -0.5*z[133]*z[25]*z[29];
	z[147] = -0.5*z[147]*z[25]*z[29];
	z[155] = -0.5*z[155]*z[25]*z[29];
	z[107] = -0.5*z[107]*z[25]*z[29];
	z[62] = -0.5*z[25]*z[29]*z[62];
	z[64] = -0.5*z[25]*z[29]*z[64];
	z[240] = -0.5*z[18]*z[25]*z[29];
	z[241] = -0.5*z[19]*z[25]*z[29];
	z[242] = 0.25*z[101]*z[20]*z[26]*z[29];
	z[243] = 0.25*z[20]*z[246]*z[26]*z[29];
	z[244] = 0.25*z[20]*z[247]*z[26]*z[29];
	z[245] = 0.25*z[20]*z[26]*z[286]*z[29];
	z[248] = 0.25*z[20]*z[26]*z[287]*z[29];
	z[249] = 0.25*z[20]*z[26]*z[288]*z[29];
	z[250] = 0.25*z[101]*z[21]*z[26]*z[29];
	z[251] = 0.25*z[21]*z[246]*z[26]*z[29];
	z[252] = 0.25*z[21]*z[247]*z[26]*z[29];
	z[253] = 0.25*z[21]*z[26]*z[286]*z[29];
	z[254] = 0.25*z[21]*z[26]*z[287]*z[29];
	z[26] = 0.25*z[21]*z[26]*z[288]*z[29];
	z[145] = -0.5*z[145]*z[27]*z[29];
	z[160] = -0.5*z[160]*z[27]*z[29];
	z[161] = -0.5*z[161]*z[27]*z[29];
	z[177] = -0.5*z[177]*z[27]*z[29];
	z[217] = -0.5*z[217]*z[27]*z[29];
	z[218] = -0.5*z[218]*z[27]*z[29];
	z[137] = -0.5*z[137]*z[27]*z[29];
	z[148] = -0.5*z[148]*z[27]*z[29];
	z[95] = -0.5*z[27]*z[29]*z[95];
	z[88] = -0.5*z[27]*z[29]*z[88];
	z[55] = -0.5*z[27]*z[29]*z[55];
	z[56] = -0.5*z[27]*z[29]*z[56];
	z[20] = -0.5*z[20]*z[27]*z[29];
	z[21] = -0.5*z[21]*z[27]*z[29];
	z[255] = -lam*z[9];
	z[9] = -mu*z[9];
	z[256] = lam*z[6];
	z[6] = mu*z[6];
	z[257] = 0.25*z[10]*z[22]*z[25]*z[69];
	z[258] = 0.25*z[11]*z[22]*z[25]*z[69];
	z[12] = 0.25*z[12]*z[22]*z[27]*z[69];
	z[13] = 0.25*z[13]*z[22]*z[27]*z[69];
	z[16] = 0.25*z[16]*z[22]*z[29]*z[69];
	z[17] = 0.25*z[17]*z[22]*z[29]*z[69];
	z[10] = 0.25*z[10]*z[23]*z[24]*z[85];
	z[11] = 0.25*z[11]*z[23]*z[24]*z[85];
	z[14] = 0.25*z[14]*z[24]*z[27]*z[85];
	z[15] = 0.25*z[15]*z[24]*z[27]*z[85];
	z[18] = 0.25*z[18]*z[24]*z[29]*z[85];
	z[19] = 0.25*z[19]*z[24]*z[29]*z[85];
	z[22] = exp(z[75]);
	z[7] = exp(z[7]);
	z[24] = exp(z[170]);
	z[8] = exp(z[8]);
	z[75] = exp(z[255]);
	z[9] = exp(z[9]);
	z[170] = exp(z[256]);
	z[6] = exp(z[6]);
	z[37] = z[37] + z[47] + z[50];
	z[32] = z[32] + z[42] + z[51];
	z[33] = z[33] + z[43] + z[52];
	z[38] = z[38] + z[48] + z[53];
	z[30] = z[30] + z[40] + z[66];
	z[35] = z[35] + z[45] + z[67];
	z[39] = z[39] + z[49] + z[68];
	z[34] = z[34] + z[44] + z[70];
	z[31] = z[31] + z[41] + z[71];
	z[36] = z[36] + z[46] + z[72];
	z[40] = h + z[73];
	z[41] = h + z[74];
	z[42] = z[108] + z[119] + z[77];
	z[43] = z[113] + z[121] + z[83];
	z[44] = z[115] + z[123] + z[86];
	z[45] = z[110] + z[124] + z[79];
	z[46] = z[111] + z[125] + z[80];
	z[47] = z[116] + z[126] + z[89];
	z[48] = z[118] + z[120] + z[90];
	z[49] = z[112] + z[127] + z[81];
	z[50] = z[109] + z[128] + z[78];
	z[51] = z[114] + z[84] + z[87];
	z[52] = h + z[129];
	z[53] = h + z[130];
	z[66] = z[134] + z[139] + z[93];
	z[67] = z[100] + z[150] + z[164];
	z[68] = z[104] + z[141] + z[153];
	z[70] = z[138] + z[156] + z[96];
	z[71] = z[140] + z[165] + z[97];
	z[72] = z[105] + z[157] + z[166];
	z[73] = z[106] + z[132] + z[162];
	z[74] = z[135] + z[149] + z[98];
	z[77] = z[136] + z[154] + z[94];
	z[78] = z[102] + z[151] + z[91];
	z[79] = h + z[168];
	z[80] = h + z[169];
	z[81] = z[143] + z[173] + z[214];
	z[83] = z[158] + z[180] + z[221];
	z[84] = z[131] + z[182] + z[225];
	z[86] = z[167] + z[175] + z[216];
	z[87] = z[176] + z[195] + z[219];
	z[89] = z[183] + z[206] + z[226];
	z[90] = z[122] + z[184] + z[227];
	z[91] = z[142] + z[178] + z[220];
	z[93] = z[152] + z[174] + z[215];
	z[94] = z[103] + z[181] + z[222];
	z[96] = h + z[228];
	z[97] = h + z[229];
	z[98] = z[144] + z[186] + z[230];
	z[100] = z[159] + z[192] + z[235];
	z[102] = z[146] + z[194] + z[237];
	z[103] = z[172] + z[188] + z[232];
	z[104] = z[189] + z[211] + z[233];
	z[105] = z[196] + z[212] + z[238];
	z[106] = z[133] + z[197] + z[239];
	z[108] = z[147] + z[190] + z[234];
	z[109] = z[155] + z[187] + z[231];
	z[107] = z[107] + z[193] + z[236];
	z[110] = h + z[240];
	z[111] = h + z[241];
	z[112] = z[145] + z[199] + z[243];
	z[113] = z[160] + z[208] + z[251];
	z[114] = z[161] + z[210] + z[253];
	z[115] = z[177] + z[203] + z[245];
	z[116] = z[204] + z[217] + z[248];
	z[118] = z[213] + z[218] + z[254];
	z[26] = z[137] + z[26] + z[28];
	z[28] = z[148] + z[205] + z[249];
	z[95] = z[202] + z[244] + z[95];
	z[88] = z[209] + z[252] + z[88];
	z[55] = z[198] + z[242] + z[55];
	z[56] = z[207] + z[250] + z[56];
	z[20] = h + z[20];
	z[21] = h + z[21];
	z[12] = z[12] + z[59] + z[76];
	z[13] = z[13] + z[54] + z[82];
	z[16] = z[16] + z[171] + z[60];
	z[17] = z[17] + z[179] + z[58];
	z[10] = z[10] + z[257] + z[57];
	z[11] = z[11] + z[258] + z[65];
	z[14] = z[14] + z[63] + z[92];
	z[15] = z[15] + z[61] + z[99];
	z[18] = z[18] + z[185] + z[64];
	z[19] = z[19] + z[191] + z[62];
	z[54] = -0.5*A*lam;
	z[57] = z[40]*z[40];
	z[58] = z[41]*z[41];
	z[59] = z[52]*z[52];
	z[60] = z[53]*z[53];
	z[61] = z[79]*z[79];
	z[62] = z[80]*z[80];
	z[63] = z[96]*z[96];
	z[64] = z[97]*z[97];
	z[65] = z[110]*z[110];
	z[76] = z[111]*z[111];
	z[82] = z[20]*z[20];
	z[92] = z[21]*z[21];
	z[22] = z[22]*z[25]*z[54];
	z[24] = z[24]*z[27]*z[54];
	z[75] = z[29]*z[54]*z[75];
	z[54] = z[170]*z[23]*z[54];
	z[99] = z[22]*z[223];
	z[119] = z[22]*z[224];
	z[120] = z[22]*z[280];
	z[121] = z[22]*z[281];
	z[122] = z[22]*z[282];
	z[123] = z[101]*z[24];
	z[124] = z[24]*z[246];
	z[125] = z[24]*z[247];
	z[126] = z[24]*z[286];
	z[127] = z[24]*z[287];
	z[24] = z[24]*z[288];
	z[128] = z[117]*z[75];
	z[129] = z[269]*z[75];
	z[130] = z[270]*z[75];
	z[131] = z[292]*z[75];
	z[132] = z[293]*z[75];
	z[75] = z[294]*z[75];
	z[133] = z[200]*z[54];
	z[134] = z[201]*z[54];
	z[135] = z[274]*z[54];
	z[136] = z[275]*z[54];
	z[137] = z[276]*z[54];
	z[54] = z[54]*z[69];
	z[22] = z[22]*z[85];
	z[57] = z[3] + z[57];
	z[58] = z[3] + z[58];
	z[59] = z[3] + z[59];
	z[60] = z[3] + z[60];
	z[61] = z[3] + z[61];
	z[62] = z[3] + z[62];
	z[63] = z[3] + z[63];
	z[64] = z[3] + z[64];
	z[65] = z[3] + z[65];
	z[76] = z[3] + z[76];
	z[82] = z[3] + z[82];
	z[3] = z[3] + z[92];
	z[92] = 1./(z[57]*z[57]);
	z[57] = 1./z[57];
	z[138] = 1./(z[58]*z[58]);
	z[58] = 1./z[58];
	z[139] = 1./(z[59]*z[59]);
	z[59] = 1./z[59];
	z[140] = 1./(z[60]*z[60]);
	z[60] = 1./z[60];
	z[141] = 1./(z[61]*z[61]);
	z[61] = 1./z[61];
	z[142] = 1./(z[62]*z[62]);
	z[62] = 1./z[62];
	z[143] = 1./(z[63]*z[63]);
	z[63] = 1./z[63];
	z[144] = 1./(z[64]*z[64]);
	z[64] = 1./z[64];
	z[145] = 1./(z[65]*z[65]);
	z[65] = 1./z[65];
	z[146] = 1./(z[76]*z[76]);
	z[76] = 1./z[76];
	z[147] = 1./(z[82]*z[82]);
	z[82] = 1./z[82];
	z[148] = 1./(z[3]*z[3]);
	z[3] = 1./z[3];
	z[40] = 2.*beta*z[1]*z[40]*z[92];
	z[57] = -z[1]*z[57];
	z[37] = 2.*beta*z[1]*z[138]*z[37]*z[41];
	z[38] = 2.*beta*z[1]*z[138]*z[38]*z[41];
	z[35] = 2.*beta*z[1]*z[138]*z[35]*z[41];
	z[39] = 2.*beta*z[1]*z[138]*z[39]*z[41];
	z[36] = 2.*beta*z[1]*z[138]*z[36]*z[41];
	z[11] = 2.*beta*z[1]*z[11]*z[138]*z[41];
	z[41] = -z[1]*z[58];
	z[42] = 2.*beta*z[1]*z[139]*z[42]*z[52];
	z[45] = 2.*beta*z[1]*z[139]*z[45]*z[52];
	z[46] = 2.*beta*z[1]*z[139]*z[46]*z[52];
	z[49] = 2.*beta*z[1]*z[139]*z[49]*z[52];
	z[50] = 2.*beta*z[1]*z[139]*z[50]*z[52];
	z[12] = 2.*beta*z[1]*z[12]*z[139]*z[52];
	z[52] = -z[1]*z[59];
	z[43] = 2.*beta*z[1]*z[140]*z[43]*z[53];
	z[44] = 2.*beta*z[1]*z[140]*z[44]*z[53];
	z[47] = 2.*beta*z[1]*z[140]*z[47]*z[53];
	z[48] = 2.*beta*z[1]*z[140]*z[48]*z[53];
	z[51] = 2.*beta*z[1]*z[140]*z[51]*z[53];
	z[13] = 2.*beta*z[1]*z[13]*z[140]*z[53];
	z[53] = -z[1]*z[60];
	z[58] = 2.*beta*z[1]*z[141]*z[66]*z[79];
	z[59] = 2.*beta*z[1]*z[141]*z[70]*z[79];
	z[60] = 2.*beta*z[1]*z[141]*z[71]*z[79];
	z[66] = 2.*beta*z[1]*z[141]*z[74]*z[79];
	z[70] = 2.*beta*z[1]*z[141]*z[77]*z[79];
	z[14] = 2.*beta*z[1]*z[14]*z[141]*z[79];
	z[61] = -z[1]*z[61];
	z[67] = 2.*beta*z[1]*z[142]*z[67]*z[80];
	z[68] = 2.*beta*z[1]*z[142]*z[68]*z[80];
	z[71] = 2.*beta*z[1]*z[142]*z[72]*z[80];
	z[72] = 2.*beta*z[1]*z[142]*z[73]*z[80];
	z[73] = 2.*beta*z[1]*z[142]*z[78]*z[80];
	z[15] = 2.*beta*z[1]*z[142]*z[15]*z[80];
	z[62] = -z[1]*z[62];
	z[74] = 2.*beta*z[1]*z[143]*z[81]*z[96];
	z[77] = 2.*beta*z[1]*z[143]*z[86]*z[96];
	z[78] = 2.*beta*z[1]*z[143]*z[87]*z[96];
	z[79] = 2.*beta*z[1]*z[143]*z[91]*z[96];
	z[80] = 2.*beta*z[1]*z[143]*z[93]*z[96];
	z[16] = 2.*beta*z[1]*z[143]*z[16]*z[96];
	z[63] = -z[1]*z[63];
	z[81] = 2.*beta*z[1]*z[144]*z[83]*z[97];
	z[83] = 2.*beta*z[1]*z[144]*z[84]*z[97];
	z[84] = 2.*beta*z[1]*z[144]*z[89]*z[97];
	z[86] = 2.*beta*z[1]*z[144]*z[90]*z[97];
	z[87] = 2.*beta*z[1]*z[144]*z[94]*z[97];
	z[17] = 2.*beta*z[1]*z[144]*z[17]*z[97];
	z[64] = -z[1]*z[64];
	z[89] = 2.*beta*z[1]*z[110]*z[145]*z[98];
	z[90] = 2.*beta*z[1]*z[103]*z[110]*z[145];
	z[91] = 2.*beta*z[1]*z[104]*z[110]*z[145];
	z[92] = 2.*beta*z[1]*z[108]*z[110]*z[145];
	z[93] = 2.*beta*z[1]*z[109]*z[110]*z[145];
	z[18] = 2.*beta*z[1]*z[110]*z[145]*z[18];
	z[65] = -z[1]*z[65];
	z[94] = 2.*beta*z[1]*z[100]*z[111]*z[146];
	z[96] = 2.*beta*z[1]*z[102]*z[111]*z[146];
	z[97] = 2.*beta*z[1]*z[105]*z[111]*z[146];
	z[98] = 2.*beta*z[1]*z[106]*z[111]*z[146];
	z[100] = 2.*beta*z[1]*z[107]*z[111]*z[146];
	z[19] = 2.*beta*z[1]*z[111]*z[146]*z[19];
	z[76] = -z[1]*z[76];
	z[102] = 2.*beta*z[1]*z[112]*z[147]*z[20];
	z[103] = 2.*beta*z[1]*z[115]*z[147]*z[20];
	z[104] = 2.*beta*z[1]*z[116]*z[147]*z[20];
	z[28] = 2.*beta*z[1]*z[147]*z[20]*z[28];
	z[95] = 2.*beta*z[1]*z[147]*z[20]*z[95];
	z[20] = 2.*beta*z[1]*z[147]*z[20]*z[55];
	z[55] = -z[1]*z[82];
	z[82] = 2.*beta*z[1]*z[113]*z[148]*z[21];
	z[105] = 2.*beta*z[1]*z[114]*z[148]*z[21];
	z[106] = 2.*beta*z[1]*z[118]*z[148]*z[21];
	z[26] = 2.*beta*z[1]*z[148]*z[21]*z[26];
	z[88] = 2.*beta*z[1]*z[148]*z[21]*z[88];
	z[21] = 2.*beta*z[1]*z[148]*z[21]*z[56];
	z[1] = -z[1]*z[3];
	z[3] = z[32]*z[40];
	z[32] = z[33]*z[40];
	z[30] = z[30]*z[40];
	z[33] = z[34]*z[40];
	z[31] = z[31]*z[40];
	z[10] = z[10]*z[40];
	z[2] = 1. + z[2];
	z[34] = z[2] + z[65];
	z[40] = z[2] + z[76];
	z[55] = z[2] + z[55];
	z[1] = z[1] + z[2];
	z[56] = z[2] + z[57];
	z[41] = z[2] + z[41];
	z[52] = z[2] + z[52];
	z[53] = z[2] + z[53];
	z[42] = z[42] + z[58] + z[82];
	z[45] = z[105] + z[45] + z[59];
	z[46] = z[106] + z[46] + z[60];
	z[26] = z[26] + z[49] + z[66];
	z[49] = z[50] + z[70] + z[88];
	z[12] = z[12] + z[14] + z[21];
	z[14] = z[2] + z[61];
	z[21] = z[30] + z[67] + z[94];
	z[3] = z[3] + z[68] + z[96];
	z[30] = z[32] + z[71] + z[97];
	z[32] = z[33] + z[72] + z[98];
	z[31] = z[100] + z[31] + z[73];
	z[10] = z[10] + z[15] + z[19];
	z[15] = z[2] + z[62];
	z[19] = z[102] + z[74] + z[89];
	z[33] = z[103] + z[77] + z[90];
	z[50] = z[104] + z[78] + z[91];
	z[28] = z[28] + z[79] + z[92];
	z[57] = z[80] + z[93] + z[95];
	z[16] = z[16] + z[18] + z[20];
	z[18] = z[2] + z[63];
	z[20] = z[35] + z[43] + z[81];
	z[35] = z[37] + z[44] + z[83];
	z[37] = z[38] + z[47] + z[84];
	z[38] = z[39] + z[48] + z[86];
	z[36] = z[36] + z[51] + z[87];
	z[11] = z[11] + z[13] + z[17];
	z[2] = z[2] + z[64];
	z[13] = beta*z[34];
	z[17] = beta*z[40];
	z[34] = beta*z[55];
	z[1] = beta*z[1];
	z[39] = beta*z[56];
	z[40] = beta*z[41];
	z[41] = beta*z[52];
	z[43] = beta*z[53];
	z[14] = beta*z[14];
	z[15] = beta*z[15];
	z[18] = beta*z[18];
	z[2] = beta*z[2];
	z[1] = z[1] + z[14] + z[41];
	z[14] = z[15] + z[17] + z[39];
	z[13] = z[13] + z[18] + z[34];
	z[2] = z[2] + z[40] + z[43];
	z[15] = pow(z[1],n);
	z[1] = pow(z[1],z[4]);
	z[17] = pow(z[14],n);
	z[14] = pow(z[14],z[4]);
	z[18] = pow(z[13],n);
	z[13] = pow(z[13],z[4]);
	z[34] = pow(z[2],n);
	z[2] = pow(z[2],z[4]);
	z[4] = 1. + z[15];
	z[15] = 1. + z[17];
	z[17] = 1. + z[18];
	z[18] = 1. + z[34];
	z[34] = pow(z[4],z[163]);
	z[4] = pow(z[4],z[5]);
	z[39] = pow(z[15],z[163]);
	z[15] = pow(z[15],z[5]);
	z[40] = pow(z[17],z[163]);
	z[17] = pow(z[17],z[5]);
	z[41] = pow(z[18],z[163]);
	z[5] = pow(z[18],z[5]);
	z[18] = 0.5*B*chi;
	z[43] = mu*z[18];
	z[44] = z[18]*z[7];
	z[14] = z[14]*z[39]*z[44];
	z[21] = z[14]*z[21];
	z[3] = z[14]*z[3];
	z[30] = z[14]*z[30];
	z[32] = z[14]*z[32];
	z[31] = z[14]*z[31];
	z[10] = z[10]*z[14];
	z[7] = z[15]*z[25]*z[43]*z[7];
	z[14] = z[223]*z[7];
	z[15] = z[224]*z[7];
	z[25] = z[280]*z[7];
	z[39] = z[281]*z[7];
	z[44] = z[282]*z[7];
	z[47] = z[18]*z[8];
	z[1] = z[1]*z[34]*z[47];
	z[34] = z[1]*z[42];
	z[42] = z[1]*z[45];
	z[45] = z[1]*z[46];
	z[26] = z[1]*z[26];
	z[46] = z[1]*z[49];
	z[1] = z[1]*z[12];
	z[4] = z[27]*z[4]*z[43]*z[8];
	z[8] = z[101]*z[4];
	z[12] = z[246]*z[4];
	z[27] = z[247]*z[4];
	z[47] = z[286]*z[4];
	z[48] = z[287]*z[4];
	z[4] = z[288]*z[4];
	z[49] = z[18]*z[9];
	z[13] = z[13]*z[40]*z[49];
	z[19] = z[13]*z[19];
	z[33] = z[13]*z[33];
	z[40] = z[13]*z[50];
	z[28] = z[13]*z[28];
	z[49] = z[13]*z[57];
	z[13] = z[13]*z[16];
	z[9] = z[17]*z[29]*z[43]*z[9];
	z[16] = z[117]*z[9];
	z[17] = z[269]*z[9];
	z[29] = z[270]*z[9];
	z[50] = z[292]*z[9];
	z[51] = z[293]*z[9];
	z[9] = z[294]*z[9];
	z[18] = z[18]*z[6];
	z[2] = z[18]*z[2]*z[41];
	z[18] = z[2]*z[20];
	z[20] = z[2]*z[35];
	z[35] = z[2]*z[37];
	z[37] = z[2]*z[38];
	z[36] = z[2]*z[36];
	z[2] = z[11]*z[2];
	z[5] = z[23]*z[43]*z[5]*z[6];
	z[6] = z[200]*z[5];
	z[11] = z[201]*z[5];
	z[23] = z[274]*z[5];
	z[38] = z[275]*z[5];
	z[41] = z[276]*z[5];
	z[5] = z[5]*z[69];
	z[7] = z[7]*z[85];
	z[6] = z[12] + z[124] + z[14] + z[17] + z[18] + z[19] + z[21] + z[34] + z[6] + z[99];
	z[6] = z[129] + z[133] + z[6];
	z[11] = z[11] + z[119] + z[125] + z[15] + z[27] + z[29] + z[31] + z[36] + z[46] + z[49];
	z[11] = z[11] + z[130] + z[134];
	z[3] = z[120] + z[126] + z[20] + z[23] + z[25] + z[3] + z[33] + z[42] + z[47] + z[50];
	z[3] = z[131] + z[135] + z[3];
	z[12] = z[121] + z[127] + z[30] + z[35] + z[38] + z[39] + z[40] + z[45] + z[48] + z[51];
	z[12] = z[12] + z[132] + z[136];
	z[4] = z[122] + z[24] + z[26] + z[28] + z[32] + z[37] + z[4] + z[41] + z[44] + z[9];
	z[4] = z[137] + z[4] + z[75];
	z[1] = z[1] + z[10] + z[123] + z[128] + z[13] + z[16] + z[2] + z[5] + z[7] + z[8];
	z[1] = z[1] + z[22] + z[54];
	z[2] = 0.5*z[6];
	z[5] = 0.5*z[11];
	z[3] = 0.5*z[3];
	z[6] = 0.5*z[12];
	z[4] = 0.5*z[4];
	z[1] = 0.5*z[1];
	
	/* output */
	/* {{z5, z2, z6},
	 *  {z2, z1, z3},
	 *  {z6, z3, z4}}
	 */
	 
	/* return values */
	dUdC[0] = z[5];
	dUdC[1] = z[2];
	dUdC[2] = z[6];
	dUdC[3] = z[2];
	dUdC[4] = z[1];
	dUdC[5] = z[3];
	dUdC[6] = z[6];
	dUdC[7] = z[3];
	dUdC[8] = z[4];
}
