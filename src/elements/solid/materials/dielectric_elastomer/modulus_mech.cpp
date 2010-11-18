#include "FSDE_inc.h"

#include <math.h>

static double z[394];

/* function to compute mechanical modulus */
void get_ddCmech(const double* params, const double *Xsi, const double* Cmat, double* dCdC) { 

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
	z[19] = C11 + C22 + C33;
	z[20] = ex*ex;
	z[21] = ey*ey;
	z[22] = ez*ez;
	z[23] = 1./(Nrig*Nrig);
	z[24] = 1./Nrig;
	z[25] = C33*z[1];
	z[26] = C32*z[2];
	z[27] = C33*z[3];
	z[28] = C31*z[4];
	z[29] = C32*z[5];
	z[30] = C31*z[6];
	z[12] = z[12] + z[16];
	z[10] = z[10] + z[17];
	z[14] = z[14] + z[18];
	z[16] = 0.06285714285714286*z[19]*z[23];
	z[17] = 0.1*z[24];
	z[1] = z[1] + z[3];
	z[3] = z[25] + z[26] + z[27] + z[28] + z[29] + z[30];
	z[2] = z[2] + z[5];
	z[4] = z[4] + z[6];
	z[5] = z[11] + z[7];
	z[6] = z[15] + z[8];
	z[7] = z[13] + z[9];
	z[8] = z[12]*z[12];
	z[9] = pow(z[12],3.);
	z[11] = z[14]*z[14];
	z[13] = pow(z[14],3.);
	z[15] = z[16] + z[17];
	z[16] = z[1]*z[1];
	z[17] = pow(z[1],3.);
	z[18] = pow(z[3],-3.);
	z[19] = 1./(z[3]*z[3]);
	z[3] = 1./z[3];
	z[23] = z[2]*z[2];
	z[24] = pow(z[2],3.);
	z[25] = z[4]*z[4];
	z[26] = pow(z[4],3.);
	z[27] = z[6]*z[6];
	z[28] = pow(z[6],3.);
	z[15] = mu*z[15];
	z[29] = 2.*ex*ey*z[10]*z[12]*z[14]*z[18];
	z[30] = 2.*ex*ey*z[1]*z[10]*z[12]*z[18];
	z[31] = 2.*ex*ey*z[1]*z[12]*z[14]*z[18];
	z[32] = 2.*z[1]*z[12]*z[14]*z[18]*z[20];
	z[33] = 2.*z[1]*z[12]*z[14]*z[18]*z[22];
	z[34] = 2.*ex*ey*z[1]*z[10]*z[14]*z[18];
	z[35] = 2.*ex*ey*z[10]*z[12]*z[18]*z[2];
	z[36] = 2.*ex*ey*z[12]*z[14]*z[18]*z[2];
	z[37] = 2.*ey*ez*z[12]*z[14]*z[18]*z[2];
	z[38] = 2.*z[12]*z[14]*z[18]*z[2]*z[20];
	z[39] = 2.*ex*ey*z[10]*z[14]*z[18]*z[2];
	z[40] = 2.*ex*ey*z[1]*z[12]*z[18]*z[2];
	z[41] = 2.*ey*ez*z[1]*z[12]*z[18]*z[2];
	z[42] = 2.*z[1]*z[12]*z[18]*z[2]*z[22];
	z[43] = 2.*ex*ey*z[1]*z[10]*z[18]*z[2];
	z[44] = 2.*ey*ez*z[1]*z[14]*z[18]*z[2];
	z[45] = 2.*z[1]*z[14]*z[18]*z[2]*z[20];
	z[46] = 2.*z[1]*z[14]*z[18]*z[2]*z[22];
	z[47] = 2.*ex*ey*z[10]*z[12]*z[18]*z[4];
	z[48] = 2.*ex*ey*z[12]*z[14]*z[18]*z[4];
	z[49] = 2.*ex*ez*z[12]*z[14]*z[18]*z[4];
	z[50] = 2.*z[12]*z[14]*z[18]*z[20]*z[4];
	z[51] = 2.*ex*ey*z[10]*z[14]*z[18]*z[4];
	z[52] = 2.*ex*ey*z[1]*z[12]*z[18]*z[4];
	z[53] = 2.*ex*ez*z[1]*z[12]*z[18]*z[4];
	z[54] = 2.*z[1]*z[12]*z[18]*z[22]*z[4];
	z[55] = 2.*ex*ey*z[1]*z[10]*z[18]*z[4];
	z[56] = 2.*ex*ez*z[1]*z[14]*z[18]*z[4];
	z[57] = 2.*z[1]*z[14]*z[18]*z[20]*z[4];
	z[58] = 2.*z[1]*z[14]*z[18]*z[22]*z[4];
	z[59] = 2.*ex*ey*z[12]*z[18]*z[2]*z[4];
	z[60] = 2.*ex*ez*z[12]*z[18]*z[2]*z[4];
	z[61] = 2.*ey*ez*z[12]*z[18]*z[2]*z[4];
	z[62] = 2.*ex*ey*z[10]*z[18]*z[2]*z[4];
	z[63] = 2.*ex*ez*z[14]*z[18]*z[2]*z[4];
	z[64] = 2.*ey*ez*z[14]*z[18]*z[2]*z[4];
	z[65] = 2.*z[14]*z[18]*z[2]*z[20]*z[4];
	z[66] = 2.*ex*ez*z[1]*z[18]*z[2]*z[4];
	z[67] = 2.*ey*ez*z[1]*z[18]*z[2]*z[4];
	z[68] = 2.*z[1]*z[18]*z[2]*z[22]*z[4];
	z[69] = 2.*ey*ez*z[12]*z[14]*z[18]*z[5];
	z[70] = 2.*ey*ez*z[1]*z[12]*z[18]*z[5];
	z[71] = 2.*ey*ez*z[1]*z[14]*z[18]*z[5];
	z[72] = 2.*ey*ez*z[12]*z[18]*z[2]*z[5];
	z[73] = 2.*ey*ez*z[14]*z[18]*z[2]*z[5];
	z[74] = 2.*ey*ez*z[1]*z[18]*z[2]*z[5];
	z[75] = 2.*ey*ez*z[12]*z[18]*z[4]*z[5];
	z[76] = 2.*ey*ez*z[14]*z[18]*z[4]*z[5];
	z[77] = 2.*ey*ez*z[1]*z[18]*z[4]*z[5];
	z[78] = 2.*ey*ez*z[18]*z[2]*z[4]*z[5];
	z[79] = 2.*ex*ey*z[10]*z[12]*z[18]*z[6];
	z[80] = 2.*ex*ey*z[12]*z[14]*z[18]*z[6];
	z[81] = 2.*z[12]*z[14]*z[18]*z[20]*z[6];
	z[82] = 2.*z[12]*z[14]*z[18]*z[21]*z[6];
	z[83] = 2.*ex*ey*z[10]*z[14]*z[18]*z[6];
	z[84] = 2.*ex*ey*z[1]*z[12]*z[18]*z[6];
	z[85] = 2.*z[1]*z[12]*z[18]*z[21]*z[6];
	z[86] = 2.*z[1]*z[12]*z[18]*z[22]*z[6];
	z[87] = 2.*ex*ey*z[1]*z[10]*z[18]*z[6];
	z[88] = 2.*z[1]*z[14]*z[18]*z[20]*z[6];
	z[89] = 2.*z[1]*z[14]*z[18]*z[21]*z[6];
	z[90] = 2.*z[1]*z[14]*z[18]*z[22]*z[6];
	z[91] = 2.*ex*ey*z[12]*z[18]*z[2]*z[6];
	z[92] = 2.*ey*ez*z[12]*z[18]*z[2]*z[6];
	z[93] = 2.*z[12]*z[18]*z[2]*z[21]*z[6];
	z[94] = 2.*ex*ey*z[10]*z[18]*z[2]*z[6];
	z[95] = 2.*ey*ez*z[14]*z[18]*z[2]*z[6];
	z[96] = 2.*z[14]*z[18]*z[2]*z[20]*z[6];
	z[97] = 2.*z[14]*z[18]*z[2]*z[21]*z[6];
	z[98] = 2.*ey*ez*z[1]*z[18]*z[2]*z[6];
	z[99] = 2.*z[1]*z[18]*z[2]*z[21]*z[6];
	z[100] = 2.*z[1]*z[18]*z[2]*z[22]*z[6];
	z[101] = 2.*ex*ey*z[12]*z[18]*z[4]*z[6];
	z[102] = 2.*ex*ez*z[12]*z[18]*z[4]*z[6];
	z[103] = 2.*z[12]*z[18]*z[21]*z[4]*z[6];
	z[104] = 2.*ex*ey*z[10]*z[18]*z[4]*z[6];
	z[105] = 2.*ex*ez*z[14]*z[18]*z[4]*z[6];
	z[106] = 2.*z[14]*z[18]*z[20]*z[4]*z[6];
	z[107] = 2.*z[14]*z[18]*z[21]*z[4]*z[6];
	z[108] = 2.*ex*ez*z[1]*z[18]*z[4]*z[6];
	z[109] = 2.*z[1]*z[18]*z[21]*z[4]*z[6];
	z[110] = 2.*z[1]*z[18]*z[22]*z[4]*z[6];
	z[111] = 2.*ex*ez*z[18]*z[2]*z[4]*z[6];
	z[112] = 2.*ey*ez*z[18]*z[2]*z[4]*z[6];
	z[113] = 2.*z[18]*z[2]*z[21]*z[4]*z[6];
	z[114] = 2.*ey*ez*z[12]*z[18]*z[5]*z[6];
	z[115] = 2.*ey*ez*z[14]*z[18]*z[5]*z[6];
	z[116] = 2.*ey*ez*z[1]*z[18]*z[5]*z[6];
	z[117] = 2.*ey*ez*z[18]*z[2]*z[5]*z[6];
	z[118] = 2.*ey*ez*z[18]*z[4]*z[5]*z[6];
	z[119] = 2.*ex*ez*z[12]*z[14]*z[18]*z[7];
	z[120] = 2.*ex*ez*z[1]*z[12]*z[18]*z[7];
	z[121] = 2.*ex*ez*z[1]*z[14]*z[18]*z[7];
	z[122] = 2.*ex*ez*z[12]*z[18]*z[2]*z[7];
	z[123] = 2.*ex*ez*z[14]*z[18]*z[2]*z[7];
	z[124] = 2.*ex*ez*z[1]*z[18]*z[2]*z[7];
	z[125] = 2.*ex*ez*z[12]*z[18]*z[4]*z[7];
	z[126] = 2.*ex*ez*z[14]*z[18]*z[4]*z[7];
	z[127] = 2.*ex*ez*z[1]*z[18]*z[4]*z[7];
	z[128] = 2.*ex*ez*z[18]*z[2]*z[4]*z[7];
	z[129] = 2.*ex*ez*z[12]*z[18]*z[6]*z[7];
	z[130] = 2.*ex*ez*z[14]*z[18]*z[6]*z[7];
	z[131] = 2.*ex*ez*z[1]*z[18]*z[6]*z[7];
	z[132] = 2.*ex*ez*z[18]*z[2]*z[6]*z[7];
	z[133] = 2.*ex*ez*z[18]*z[4]*z[6]*z[7];
	z[134] = 2.*ex*ey*z[10]*z[18]*z[8];
	z[135] = 2.*ex*ey*z[14]*z[18]*z[8];
	z[136] = 2.*z[14]*z[18]*z[20]*z[8];
	z[137] = 2.*ex*ey*z[1]*z[18]*z[8];
	z[138] = 2.*z[1]*z[18]*z[22]*z[8];
	z[139] = 2.*ex*ey*z[18]*z[2]*z[8];
	z[140] = 2.*ey*ez*z[18]*z[2]*z[8];
	z[141] = 2.*ex*ey*z[18]*z[4]*z[8];
	z[142] = 2.*ex*ez*z[18]*z[4]*z[8];
	z[143] = 2.*ey*ez*z[18]*z[5]*z[8];
	z[144] = 2.*ex*ey*z[18]*z[6]*z[8];
	z[145] = 2.*z[18]*z[21]*z[6]*z[8];
	z[8] = 2.*ex*ez*z[18]*z[7]*z[8];
	z[9] = 2.*ex*ey*z[18]*z[9];
	z[146] = 2.*ex*ey*z[11]*z[12]*z[18];
	z[147] = 2.*z[11]*z[12]*z[18]*z[20];
	z[148] = 2.*ex*ey*z[10]*z[11]*z[18];
	z[149] = 2.*z[1]*z[11]*z[18]*z[20];
	z[150] = 2.*z[1]*z[11]*z[18]*z[22];
	z[151] = 2.*ey*ez*z[11]*z[18]*z[2];
	z[152] = 2.*z[11]*z[18]*z[2]*z[20];
	z[153] = 2.*ex*ez*z[11]*z[18]*z[4];
	z[154] = 2.*z[11]*z[18]*z[20]*z[4];
	z[155] = 2.*ey*ez*z[11]*z[18]*z[5];
	z[156] = 2.*z[11]*z[18]*z[20]*z[6];
	z[157] = 2.*z[11]*z[18]*z[21]*z[6];
	z[11] = 2.*ex*ez*z[11]*z[18]*z[7];
	z[13] = 2.*z[13]*z[18]*z[20];
	z[158] = 2.*ex*ey*z[12]*z[16]*z[18];
	z[159] = 2.*z[12]*z[16]*z[18]*z[22];
	z[160] = 2.*ex*ey*z[10]*z[16]*z[18];
	z[161] = 2.*z[14]*z[16]*z[18]*z[20];
	z[162] = 2.*z[14]*z[16]*z[18]*z[22];
	z[163] = 2.*ey*ez*z[16]*z[18]*z[2];
	z[164] = 2.*z[16]*z[18]*z[2]*z[22];
	z[165] = 2.*ex*ez*z[16]*z[18]*z[4];
	z[166] = 2.*z[16]*z[18]*z[22]*z[4];
	z[167] = 2.*ey*ez*z[16]*z[18]*z[5];
	z[168] = 2.*z[16]*z[18]*z[21]*z[6];
	z[169] = 2.*z[16]*z[18]*z[22]*z[6];
	z[16] = 2.*ex*ez*z[16]*z[18]*z[7];
	z[17] = 2.*z[17]*z[18]*z[22];
	z[170] = -C11*ex*ey*z[12]*z[19];
	z[171] = 2.*C12*ex*ey*z[12]*z[19];
	z[172] = -2.*C13*ex*ey*z[12]*z[19];
	z[173] = C13*ex*ey*z[12]*z[19];
	z[174] = C21*ex*ey*z[12]*z[19];
	z[175] = -C22*ex*ey*z[12]*z[19];
	z[176] = -C23*ex*ey*z[12]*z[19];
	z[177] = C23*ex*ey*z[12]*z[19];
	z[178] = -C33*ex*ey*z[12]*z[19];
	z[179] = 2.*C33*ex*ey*z[12]*z[19];
	z[180] = C13*ex*ez*z[12]*z[19];
	z[181] = -C21*ex*ez*z[12]*z[19];
	z[182] = C22*ex*ez*z[12]*z[19];
	z[183] = C31*ex*ez*z[12]*z[19];
	z[184] = -2.*C32*ex*ez*z[12]*z[19];
	z[185] = C11*ey*ez*z[12]*z[19];
	z[186] = -C12*ey*ez*z[12]*z[19];
	z[187] = -2.*C13*ey*ez*z[12]*z[19];
	z[188] = C23*ey*ez*z[12]*z[19];
	z[189] = C32*ey*ez*z[12]*z[19];
	z[190] = -C22*z[12]*z[19]*z[20];
	z[191] = C23*z[12]*z[19]*z[20];
	z[192] = -C33*z[12]*z[19]*z[20];
	z[193] = -C11*z[12]*z[19]*z[21];
	z[194] = C13*z[12]*z[19]*z[21];
	z[195] = -C33*z[12]*z[19]*z[21];
	z[196] = -C11*z[12]*z[19]*z[22];
	z[197] = 2.*C12*z[12]*z[19]*z[22];
	z[198] = -C22*z[12]*z[19]*z[22];
	z[199] = -C11*ex*ey*z[10]*z[19];
	z[200] = C12*ex*ey*z[10]*z[19];
	z[201] = -C13*ex*ey*z[10]*z[19];
	z[202] = C13*ex*ey*z[10]*z[19];
	z[203] = -C22*ex*ey*z[10]*z[19];
	z[204] = C23*ex*ey*z[10]*z[19];
	z[205] = -C33*ex*ey*z[10]*z[19];
	z[206] = C12*ex*ey*z[14]*z[19];
	z[207] = -C13*ex*ey*z[14]*z[19];
	z[208] = C21*ex*ey*z[14]*z[19];
	z[209] = -C23*ex*ey*z[14]*z[19];
	z[210] = C33*ex*ey*z[14]*z[19];
	z[211] = C13*ex*ez*z[14]*z[19];
	z[212] = -C21*ex*ez*z[14]*z[19];
	z[213] = C22*ex*ez*z[14]*z[19];
	z[214] = C31*ex*ez*z[14]*z[19];
	z[215] = -C32*ex*ez*z[14]*z[19];
	z[216] = C11*ey*ez*z[14]*z[19];
	z[217] = -C12*ey*ez*z[14]*z[19];
	z[218] = -C13*ey*ez*z[14]*z[19];
	z[219] = 2.*C23*ey*ez*z[14]*z[19];
	z[220] = 2.*C32*ey*ez*z[14]*z[19];
	z[221] = -C11*z[14]*z[19]*z[20];
	z[222] = C12*z[14]*z[19]*z[20];
	z[223] = -C13*z[14]*z[19]*z[20];
	z[224] = C13*z[14]*z[19]*z[20];
	z[225] = -2.*C22*z[14]*z[19]*z[20];
	z[226] = 2.*C23*z[14]*z[19]*z[20];
	z[227] = -2.*C33*z[14]*z[19]*z[20];
	z[228] = -C11*z[14]*z[19]*z[21];
	z[229] = C13*z[14]*z[19]*z[21];
	z[230] = -2.*C33*z[14]*z[19]*z[21];
	z[231] = -C11*z[14]*z[19]*z[22];
	z[232] = C12*z[14]*z[19]*z[22];
	z[233] = -2.*C22*z[14]*z[19]*z[22];
	z[234] = 2.*C12*ex*ey*z[1]*z[19];
	z[235] = -C13*ex*ey*z[1]*z[19];
	z[236] = 2.*C21*ex*ey*z[1]*z[19];
	z[237] = -C23*ex*ey*z[1]*z[19];
	z[238] = C33*ex*ey*z[1]*z[19];
	z[239] = C13*ex*ez*z[1]*z[19];
	z[240] = -C21*ex*ez*z[1]*z[19];
	z[241] = C22*ex*ez*z[1]*z[19];
	z[242] = C31*ex*ez*z[1]*z[19];
	z[243] = -C32*ex*ez*z[1]*z[19];
	z[244] = C11*ey*ez*z[1]*z[19];
	z[245] = -C12*ey*ez*z[1]*z[19];
	z[246] = -C13*ey*ez*z[1]*z[19];
	z[247] = C23*ey*ez*z[1]*z[19];
	z[248] = C32*ey*ez*z[1]*z[19];
	z[249] = -2.*C22*z[1]*z[19]*z[20];
	z[250] = C23*z[1]*z[19]*z[20];
	z[251] = -C33*z[1]*z[19]*z[20];
	z[252] = -2.*C11*z[1]*z[19]*z[21];
	z[253] = C13*z[1]*z[19]*z[21];
	z[254] = -C33*z[1]*z[19]*z[21];
	z[255] = -2.*C11*z[1]*z[19]*z[22];
	z[256] = 2.*C12*z[1]*z[19]*z[22];
	z[257] = -C13*z[1]*z[19]*z[22];
	z[258] = C13*z[1]*z[19]*z[22];
	z[259] = -2.*C22*z[1]*z[19]*z[22];
	z[260] = C23*z[1]*z[19]*z[22];
	z[261] = -C33*z[1]*z[19]*z[22];
	z[262] = C12*ex*ey*z[19]*z[2];
	z[263] = -2.*C13*ex*ey*z[19]*z[2];
	z[264] = C21*ex*ey*z[19]*z[2];
	z[265] = -C23*ex*ey*z[19]*z[2];
	z[266] = C33*ex*ey*z[19]*z[2];
	z[267] = C13*ex*ez*z[19]*z[2];
	z[268] = -2.*C21*ex*ez*z[19]*z[2];
	z[269] = C22*ex*ez*z[19]*z[2];
	z[270] = C31*ex*ez*z[19]*z[2];
	z[271] = -C32*ex*ez*z[19]*z[2];
	z[272] = -C11*ey*ez*z[19]*z[2];
	z[273] = 2.*C11*ey*ez*z[19]*z[2];
	z[274] = -C12*ey*ez*z[19]*z[2];
	z[275] = C12*ey*ez*z[19]*z[2];
	z[276] = -2.*C13*ey*ez*z[19]*z[2];
	z[277] = C13*ey*ez*z[19]*z[2];
	z[278] = -C22*ey*ez*z[19]*z[2];
	z[279] = 2.*C23*ey*ez*z[19]*z[2];
	z[280] = C32*ey*ez*z[19]*z[2];
	z[281] = -C33*ey*ez*z[19]*z[2];
	z[282] = -C22*z[19]*z[2]*z[20];
	z[283] = 2.*C23*z[19]*z[2]*z[20];
	z[284] = -C33*z[19]*z[2]*z[20];
	z[285] = -C11*z[19]*z[2]*z[21];
	z[286] = C13*z[19]*z[2]*z[21];
	z[287] = -C33*z[19]*z[2]*z[21];
	z[288] = -C11*z[19]*z[2]*z[22];
	z[289] = C12*z[19]*z[2]*z[22];
	z[290] = -C22*z[19]*z[2]*z[22];
	z[291] = C12*ex*ey*z[19]*z[4];
	z[292] = -C13*ex*ey*z[19]*z[4];
	z[293] = C21*ex*ey*z[19]*z[4];
	z[294] = -2.*C23*ex*ey*z[19]*z[4];
	z[295] = C33*ex*ey*z[19]*z[4];
	z[296] = -C11*ex*ez*z[19]*z[4];
	z[297] = C12*ex*ez*z[19]*z[4];
	z[298] = -C13*ex*ez*z[19]*z[4];
	z[299] = 2.*C13*ex*ez*z[19]*z[4];
	z[300] = -C21*ex*ez*z[19]*z[4];
	z[301] = -C22*ex*ez*z[19]*z[4];
	z[302] = 2.*C22*ex*ez*z[19]*z[4];
	z[303] = C23*ex*ez*z[19]*z[4];
	z[304] = C31*ex*ez*z[19]*z[4];
	z[305] = -C32*ex*ez*z[19]*z[4];
	z[306] = -C33*ex*ez*z[19]*z[4];
	z[307] = C11*ey*ez*z[19]*z[4];
	z[308] = -2.*C12*ey*ez*z[19]*z[4];
	z[309] = -C13*ey*ez*z[19]*z[4];
	z[310] = C23*ey*ez*z[19]*z[4];
	z[311] = C32*ey*ez*z[19]*z[4];
	z[312] = -C22*z[19]*z[20]*z[4];
	z[313] = C23*z[19]*z[20]*z[4];
	z[314] = -C33*z[19]*z[20]*z[4];
	z[315] = -C11*z[19]*z[21]*z[4];
	z[316] = 2.*C13*z[19]*z[21]*z[4];
	z[317] = -C33*z[19]*z[21]*z[4];
	z[318] = -C11*z[19]*z[22]*z[4];
	z[319] = C12*z[19]*z[22]*z[4];
	z[320] = -C22*z[19]*z[22]*z[4];
	z[321] = -C11*ey*ez*z[19]*z[5];
	z[322] = C12*ey*ez*z[19]*z[5];
	z[323] = -C13*ey*ez*z[19]*z[5];
	z[324] = C13*ey*ez*z[19]*z[5];
	z[325] = -C22*ey*ez*z[19]*z[5];
	z[326] = C23*ey*ez*z[19]*z[5];
	z[327] = -C33*ey*ez*z[19]*z[5];
	z[328] = C12*ex*ey*z[19]*z[6];
	z[329] = -C13*ex*ey*z[19]*z[6];
	z[330] = C21*ex*ey*z[19]*z[6];
	z[331] = -C23*ex*ey*z[19]*z[6];
	z[332] = C33*ex*ey*z[19]*z[6];
	z[333] = 2.*C13*ex*ez*z[19]*z[6];
	z[334] = -C21*ex*ez*z[19]*z[6];
	z[335] = C22*ex*ez*z[19]*z[6];
	z[336] = 2.*C31*ex*ez*z[19]*z[6];
	z[337] = -C32*ex*ez*z[19]*z[6];
	z[338] = C11*ey*ez*z[19]*z[6];
	z[339] = -C12*ey*ez*z[19]*z[6];
	z[340] = -C13*ey*ez*z[19]*z[6];
	z[341] = C23*ey*ez*z[19]*z[6];
	z[342] = C32*ey*ez*z[19]*z[6];
	z[343] = -C22*z[19]*z[20]*z[6];
	z[344] = C23*z[19]*z[20]*z[6];
	z[345] = -2.*C33*z[19]*z[20]*z[6];
	z[346] = -2.*C11*z[19]*z[21]*z[6];
	z[347] = C12*z[19]*z[21]*z[6];
	z[348] = -C13*z[19]*z[21]*z[6];
	z[349] = 2.*C13*z[19]*z[21]*z[6];
	z[350] = -C22*z[19]*z[21]*z[6];
	z[351] = C23*z[19]*z[21]*z[6];
	z[352] = -2.*C33*z[19]*z[21]*z[6];
	z[353] = -2.*C11*z[19]*z[22]*z[6];
	z[354] = C12*z[19]*z[22]*z[6];
	z[355] = -C22*z[19]*z[22]*z[6];
	z[356] = -C11*ex*ez*z[19]*z[7];
	z[357] = C12*ex*ez*z[19]*z[7];
	z[358] = -C13*ex*ez*z[19]*z[7];
	z[359] = C13*ex*ez*z[19]*z[7];
	z[360] = -C22*ex*ez*z[19]*z[7];
	z[361] = C23*ex*ez*z[19]*z[7];
	z[19] = -C33*ex*ez*z[19]*z[7];
	z[362] = -ex*ey*z[3];
	z[363] = -ex*ez*z[3];
	z[364] = ex*ez*z[3];
	z[365] = -ey*ez*z[3];
	z[366] = z[20]*z[3];
	z[367] = z[21]*z[3];
	z[3] = z[22]*z[3];
	z[368] = 2.*ex*ey*z[12]*z[18]*z[23];
	z[369] = 2.*ey*ez*z[12]*z[18]*z[23];
	z[370] = 2.*ex*ey*z[10]*z[18]*z[23];
	z[371] = 2.*ey*ez*z[14]*z[18]*z[23];
	z[372] = 2.*z[14]*z[18]*z[20]*z[23];
	z[373] = 2.*ey*ez*z[1]*z[18]*z[23];
	z[374] = 2.*z[1]*z[18]*z[22]*z[23];
	z[375] = 2.*ex*ez*z[18]*z[23]*z[4];
	z[376] = 2.*ey*ez*z[18]*z[23]*z[4];
	z[377] = 2.*ey*ez*z[18]*z[23]*z[5];
	z[378] = 2.*ey*ez*z[18]*z[23]*z[6];
	z[379] = 2.*z[18]*z[21]*z[23]*z[6];
	z[23] = 2.*ex*ez*z[18]*z[23]*z[7];
	z[24] = 2.*ey*ez*z[18]*z[24];
	z[380] = 2.*ex*ey*z[12]*z[18]*z[25];
	z[381] = 2.*ex*ez*z[12]*z[18]*z[25];
	z[382] = 2.*ex*ey*z[10]*z[18]*z[25];
	z[383] = 2.*ex*ez*z[14]*z[18]*z[25];
	z[384] = 2.*z[14]*z[18]*z[20]*z[25];
	z[385] = 2.*ex*ez*z[1]*z[18]*z[25];
	z[386] = 2.*z[1]*z[18]*z[22]*z[25];
	z[387] = 2.*ex*ez*z[18]*z[2]*z[25];
	z[388] = 2.*ey*ez*z[18]*z[2]*z[25];
	z[389] = 2.*ey*ez*z[18]*z[25]*z[5];
	z[390] = 2.*ex*ez*z[18]*z[25]*z[6];
	z[6] = 2.*z[18]*z[21]*z[25]*z[6];
	z[25] = 2.*ex*ez*z[18]*z[25]*z[7];
	z[26] = 2.*ex*ez*z[18]*z[26];
	z[391] = 2.*ex*ey*z[12]*z[18]*z[27];
	z[12] = 2.*z[12]*z[18]*z[21]*z[27];
	z[10] = 2.*ex*ey*z[10]*z[18]*z[27];
	z[20] = 2.*z[14]*z[18]*z[20]*z[27];
	z[14] = 2.*z[14]*z[18]*z[21]*z[27];
	z[392] = 2.*z[1]*z[18]*z[21]*z[27];
	z[1] = 2.*z[1]*z[18]*z[22]*z[27];
	z[22] = 2.*ey*ez*z[18]*z[2]*z[27];
	z[2] = 2.*z[18]*z[2]*z[21]*z[27];
	z[393] = 2.*ex*ez*z[18]*z[27]*z[4];
	z[4] = 2.*z[18]*z[21]*z[27]*z[4];
	z[5] = 2.*ey*ez*z[18]*z[27]*z[5];
	z[7] = 2.*ex*ez*z[18]*z[27]*z[7];
	z[18] = 2.*z[18]*z[21]*z[28];
	z[8] = z[134] + z[136] + z[138] + z[140] + z[142] + z[143] + z[145] + z[179] + z[8] + z[9];
	z[9] = z[11] + z[13] + z[146] + z[148] + z[150] + z[151] + z[153] + z[155] + z[157] + z[219];
	z[11] = z[158] + z[16] + z[160] + z[161] + z[163] + z[165] + z[167] + z[168] + z[17] + z[234];
	z[13] = z[263] + z[268] + z[273] + z[283] + z[368] + z[370] + z[372] + z[374] + z[375] + z[377];
	z[16] = z[294] + z[302] + z[308] + z[316] + z[380] + z[382] + z[384] + z[386] + z[388] + z[389];
	z[17] = z[102] + z[114] + z[129] + z[144] + z[180] + z[183] + z[79] + z[81] + z[86] + z[92];
	z[21] = z[105] + z[115] + z[130] + z[156] + z[178] + z[205] + z[80] + z[83] + z[90] + z[95];
	z[27] = z[108] + z[116] + z[131] + z[169] + z[170] + z[199] + z[84] + z[87] + z[88] + z[98];
	z[28] = z[100] + z[111] + z[117] + z[132] + z[267] + z[270] + z[284] + z[91] + z[94] + z[96];
	z[79] = z[101] + z[104] + z[106] + z[110] + z[112] + z[118] + z[133] + z[173] + z[202] + z[224];
	z[1] = z[1] + z[10] + z[20] + z[22] + z[333] + z[336] + z[345] + z[353] + z[391] + z[393];
	z[8] = z[184] + z[187] + z[197] + z[8];
	z[9] = z[220] + z[230] + z[233] + z[9];
	z[10] = z[11] + z[236] + z[249] + z[252];
	z[11] = z[13] + z[23] + z[24] + z[379];
	z[6] = z[16] + z[25] + z[26] + z[6];
	z[12] = z[12] + z[17] + z[192] + z[196] + z[332] + z[337] + z[340] + z[354];
	z[13] = z[211] + z[214] + z[227] + z[231] + z[261] + z[281] + z[306] + z[327] + z[341] + z[342];
	z[16] = z[221] + z[239] + z[242] + z[251] + z[255] + z[272] + z[296] + z[321] + z[328] + z[330];
	z[2] = z[2] + z[28] + z[288] + z[329] + z[334] + z[338] + z[344] + z[378];
	z[17] = z[258] + z[277] + z[299] + z[304] + z[314] + z[318] + z[324] + z[331] + z[335] + z[339];
	z[1] = z[1] + z[18] + z[5] + z[7];
	z[3] = z[13] + z[14] + z[19] + z[21] + z[3] + z[352] + z[355];
	z[5] = z[16] + z[27] + z[343] + z[346] + z[356] + z[366] + z[392];
	z[4] = z[17] + z[349] + z[359] + z[363] + z[390] + z[4] + z[79];
	z[7] = z[119] + z[135] + z[147] + z[188] + z[189] + z[195] + z[198] + z[210] + z[69] + z[82];
	z[7] = z[215] + z[218] + z[232] + z[29] + z[33] + z[37] + z[49] + z[7];
	z[13] = z[120] + z[137] + z[159] + z[171] + z[174] + z[190] + z[193] + z[200] + z[70] + z[85];
	z[14] = z[222] + z[238] + z[243] + z[246] + z[256] + z[275] + z[297] + z[322] + z[347] + z[357];
	z[13] = z[13] + z[14] + z[30] + z[32] + z[362] + z[41] + z[53];
	z[14] = z[121] + z[149] + z[162] + z[175] + z[203] + z[206] + z[208] + z[225] + z[71] + z[89];
	z[16] = z[228] + z[247] + z[248] + z[254] + z[259] + z[278] + z[301] + z[325] + z[350] + z[360];
	z[14] = z[14] + z[16] + z[31] + z[34] + z[367] + z[44] + z[56];
	z[16] = z[122] + z[139] + z[172] + z[181] + z[185] + z[191] + z[201] + z[223] + z[72] + z[93];
	z[17] = z[257] + z[266] + z[271] + z[276] + z[289] + z[298] + z[323] + z[348] + z[358] + z[364];
	z[16] = z[16] + z[17] + z[35] + z[369] + z[38] + z[42] + z[60];
	z[17] = z[103] + z[125] + z[141] + z[176] + z[182] + z[186] + z[194] + z[295] + z[305] + z[75];
	z[17] = z[17] + z[309] + z[319] + z[381] + z[47] + z[50] + z[54] + z[61];
	z[18] = z[123] + z[152] + z[177] + z[204] + z[207] + z[212] + z[216] + z[226] + z[73] + z[97];
	z[19] = z[260] + z[279] + z[280] + z[287] + z[290] + z[303] + z[326] + z[351] + z[361] + z[365];
	z[18] = z[18] + z[19] + z[36] + z[371] + z[39] + z[46] + z[63];
	z[19] = z[107] + z[126] + z[154] + z[209] + z[213] + z[217] + z[229] + z[310] + z[311] + z[76];
	z[19] = z[19] + z[317] + z[320] + z[383] + z[48] + z[51] + z[58] + z[64];
	z[20] = z[124] + z[164] + z[235] + z[240] + z[244] + z[250] + z[262] + z[264] + z[74] + z[99];
	z[20] = z[20] + z[282] + z[285] + z[373] + z[40] + z[43] + z[45] + z[66];
	z[21] = z[109] + z[127] + z[166] + z[237] + z[241] + z[245] + z[253] + z[291] + z[293] + z[77];
	z[21] = z[21] + z[312] + z[315] + z[385] + z[52] + z[55] + z[57] + z[67];
	z[22] = z[113] + z[128] + z[265] + z[269] + z[274] + z[286] + z[292] + z[300] + z[307] + z[78];
	z[22] = z[22] + z[313] + z[376] + z[387] + z[59] + z[62] + z[65] + z[68];
	z[23] = -0.5*epsilon;
	z[8] = z[23]*z[8];
	z[9] = z[23]*z[9];
	z[10] = z[10]*z[23];
	z[11] = z[11]*z[23];
	z[6] = z[23]*z[6];
	z[12] = z[12]*z[23];
	z[2] = z[2]*z[23];
	z[1] = z[1]*z[23];
	z[3] = z[23]*z[3];
	z[5] = z[23]*z[5];
	z[4] = z[23]*z[4];
	z[7] = z[23]*z[7];
	z[13] = z[13]*z[23];
	z[14] = z[14]*z[23];
	z[16] = z[16]*z[23];
	z[17] = z[17]*z[23];
	z[18] = z[18]*z[23];
	z[19] = z[19]*z[23];
	z[20] = z[20]*z[23];
	z[21] = z[21]*z[23];
	z[22] = z[22]*z[23];
	z[9] = z[15] + z[9];
	z[10] = z[10] + z[15];
	z[1] = z[1] + z[15];
	z[3] = z[15] + z[3];
	z[5] = z[15] + z[5];
	z[14] = z[14] + z[15];

	/* dCdC: 6 x 6 */
	dCdC[ 0] = z[9];
	dCdC[ 1] = z[3];
	dCdC[ 2] = z[14];
	dCdC[ 3] = z[18];
	dCdC[ 4] = z[19];
	dCdC[ 5] = z[7];

	dCdC[ 6] = z[3];
	dCdC[ 7] = z[1];
	dCdC[ 8] = z[5];
	dCdC[ 9] = z[2];
	dCdC[10] = z[4];
	dCdC[11] = z[12];

	dCdC[12] = z[14];
	dCdC[13] = z[5];
	dCdC[14] = z[10];
	dCdC[15] = z[20];
	dCdC[16] = z[21];
	dCdC[17] = z[13];

	dCdC[18] = z[18];
	dCdC[19] = z[2];
	dCdC[20] = z[20];
	dCdC[21] = z[11];
	dCdC[22] = z[22];
	dCdC[23] = z[16];

	dCdC[24] = z[19];
	dCdC[25] = z[4];
	dCdC[26] = z[21];
	dCdC[27] = z[22];
	dCdC[28] = z[6];
	dCdC[29] = z[17];

	dCdC[30] = z[7];
	dCdC[31] = z[12];
	dCdC[32] = z[13];
	dCdC[33] = z[16];
	dCdC[34] = z[17];
	dCdC[35] = z[8];
}