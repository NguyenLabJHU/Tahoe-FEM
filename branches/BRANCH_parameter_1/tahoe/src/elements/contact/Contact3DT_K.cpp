/* $Id: Contact3DT_K.cpp,v 1.2 2003-03-02 19:35:07 paklein Exp $ */
#include "Contact3DT.h"

void Contact3DT::DDg_tri_facet(
	double* X1, double* X2, double* X3, double* pX,
	double* d1, double* d2, double* d3, double* pd,
	dMatrixT& K) const
{
	/* dimension checks */
	if (K.Rows() < 12 || K.Cols() < 12) 
		ExceptionT::SizeMismatch("Contact3DT::DDg_tri_facet");

	/* translation */
	double X11 = X1[0];
	double X12 = X1[1];
	double X13 = X1[2];

	double X21 = X2[0];
	double X22 = X2[1];
	double X23 = X2[2];

	double X31 = X3[0];
	double X32 = X3[1];
	double X33 = X3[2];

	double pX1 = pX[0];
	double pX2 = pX[1];
	double pX3 = pX[2];

	double d11 = d1[0];
	double d12 = d1[1];
	double d13 = d1[2];

	double d21 = d2[0];
	double d22 = d2[1];
	double d23 = d2[2];

	double d31 = d3[0];
	double d32 = d3[1];
	double d33 = d3[2];

	double pd1 = pd[0];
	double pd2 = pd[1];
	double pd3 = pd[2];

	/* work space */
	double z[544];

	/* the Z's - max[Z] = 543 */
	z[1] = -d11;
	z[2] = -d12;
	z[3] = -d13;
	z[4] = -d21;
	z[5] = d13*d21;
	z[6] = -d22;
	z[7] = d11*d22;
	z[8] = -d23;
	z[9] = d12*d23;
	z[10] = -d31;
	z[11] = d12*d31;
	z[12] = d23*d31;
	z[13] = -d32;
	z[14] = d13*d32;
	z[15] = d21*d32;
	z[16] = -d33;
	z[17] = d11*d33;
	z[18] = d22*d33;
	z[19] = -X11;
	z[20] = d22*X11;
	z[21] = d33*X11;
	z[22] = -X12;
	z[23] = d23*X12;
	z[24] = d31*X12;
	z[25] = -X13;
	z[26] = d21*X13;
	z[27] = d32*X13;
	z[28] = -X21;
	z[29] = d13*X21;
	z[30] = d32*X21;
	z[31] = X13*X21;
	z[32] = -X22;
	z[33] = d11*X22;
	z[34] = d33*X22;
	z[35] = X11*X22;
	z[36] = -X23;
	z[37] = d12*X23;
	z[38] = d31*X23;
	z[39] = X12*X23;
	z[40] = -X31;
	z[41] = d12*X31;
	z[42] = d23*X31;
	z[43] = X12*X31;
	z[44] = X23*X31;
	z[45] = -X32;
	z[46] = d13*X32;
	z[47] = d21*X32;
	z[48] = X13*X32;
	z[49] = X21*X32;
	z[50] = -X33;
	z[51] = d11*X33;
	z[52] = d22*X33;
	z[53] = X11*X33;
	z[54] = X22*X33;
	z[55] = d23*z[1];
	z[56] = d32*z[1];
	z[57] = X23*z[1];
	z[58] = X32*z[1];
	z[59] = X13*z[10];
	z[60] = X22*z[10];
	z[61] = X11*z[13];
	z[62] = X23*z[13];
	z[63] = X12*z[16];
	z[64] = X21*z[16];
	z[65] = X23*z[19];
	z[66] = X32*z[19];
	z[67] = d21*z[2];
	z[68] = d33*z[2];
	z[69] = X21*z[2];
	z[70] = X33*z[2];
	z[71] = X21*z[22];
	z[72] = X33*z[22];
	z[73] = X22*z[25];
	z[74] = X31*z[25];
	z[75] = X33*z[28];
	z[76] = d22*z[3];
	z[77] = d31*z[3];
	z[78] = X22*z[3];
	z[79] = X31*z[3];
	z[80] = X31*z[32];
	z[81] = X32*z[36];
	z[82] = d33*z[4];
	z[83] = X12*z[4];
	z[84] = X33*z[4];
	z[85] = d31*z[6];
	z[86] = X13*z[6];
	z[87] = X31*z[6];
	z[88] = d32*z[8];
	z[89] = X11*z[8];
	z[90] = X32*z[8];
	z[1] = z[1] + z[19];
	z[19] = d22 + X22 + z[2] + z[22];
	z[91] = d32 + X32 + z[2] + z[22];
	z[92] = d23 + X23 + z[25] + z[3];
	z[93] = d33 + X33 + z[25] + z[3];
	z[94] = d11 + X11 + z[28] + z[4];
	z[95] = d31 + X31 + z[28] + z[4];
	z[96] = d11 + X11 + z[10] + z[40];
	z[97] = d21 + X21 + z[10] + z[40];
	z[98] = d12 + X12 + z[13] + z[45];
	z[99] = d22 + X22 + z[13] + z[45];
	z[100] = d13 + X13 + z[16] + z[50];
	z[101] = d23 + X23 + z[16] + z[50];
	z[102] = d12 + X12 + z[32] + z[6];
	z[103] = d32 + X32 + z[32] + z[6];
	z[2] = z[13] + z[2] + z[22] + z[32] + z[45] + z[6];
	z[6] = d13 + X13 + z[36] + z[8];
	z[13] = d33 + X33 + z[36] + z[8];
	z[3] = z[16] + z[25] + z[3] + z[36] + z[50] + z[8];
	z[8] = z[11] + z[15] + z[20] + z[24] + z[30] + z[33] + z[35] + z[41] + z[43] + z[47];
	z[7] = z[49] + z[56] + z[58] + z[60] + z[61] + z[66] + z[67] + z[69] + z[7] + z[8];
	z[7] = z[7] + z[71] + z[80] + z[83] + z[85] + z[87];
	z[5] = z[12] + z[17] + z[21] + z[26] + z[29] + z[31] + z[38] + z[42] + z[44] + z[5];
	z[5] = z[5] + z[51] + z[53] + z[55] + z[57] + z[59] + z[64] + z[65] + z[74] + z[75];
	z[5] = z[5] + z[77] + z[79] + z[82] + z[84] + z[89];
	z[8] = z[14] + z[18] + z[23] + z[27] + z[34] + z[37] + z[39] + z[46] + z[48] + z[52];
	z[8] = z[54] + z[62] + z[63] + z[68] + z[70] + z[72] + z[73] + z[76] + z[78] + z[8];
	z[8] = z[8] + z[81] + z[86] + z[88] + z[9] + z[90];
	z[9] = d21 + X21 + z[1];
	z[11] = d31 + X31 + z[1];
	z[1] = z[1] + z[10] + z[28] + z[4] + z[40];
	z[4] = 2.*z[98];
	z[10] = z[98]*z[98];
	z[12] = z[99]*z[99];
	z[14] = z[100]*z[100];
	z[15] = z[101]*z[101];
	z[16] = 2.*z[102]*z[99];
	z[17] = z[102]*z[102];
	z[18] = z[103]*z[103];
	z[2] = 0.3333333333333333*z[2];
	z[20] = 2.*z[101]*z[6];
	z[21] = 2.*z[103]*z[6];
	z[22] = z[6]*z[6];
	z[23] = 2.*z[100]*z[13];
	z[24] = z[13]*z[13];
	z[3] = 0.3333333333333333*z[3];
	z[25] = -2.*z[7];
	z[26] = 2.*z[7];
	z[27] = z[7]*z[7];
	z[28] = -2.*z[5];
	z[29] = 2.*z[5];
	z[30] = z[5]*z[5];
	z[31] = -2.*z[8];
	z[32] = 2.*z[8];
	z[33] = z[8]*z[8];
	z[34] = 2.*z[9]*z[99];
	z[35] = z[9]*z[9];
	z[36] = 2.*z[11]*z[13];
	z[37] = z[11]*z[11];
	z[1] = 0.3333333333333333*z[1];
	z[38] = z[101]*z[4];
	z[39] = z[103]*z[4];
	z[40] = z[4]*z[6];
	z[41] = z[4]*z[8];
	z[42] = z[26]*z[99];
	z[43] = z[102]*z[26];
	z[44] = z[26]*z[9];
	z[45] = z[100]*z[29];
	z[46] = z[13]*z[29];
	z[47] = z[11]*z[29];
	z[48] = z[101]*z[32];
	z[49] = z[103]*z[32];
	z[50] = z[32]*z[6];
	z[51] = 2.*z[101]*z[19];
	z[52] = 2.*z[103]*z[19];
	z[4] = z[19]*z[4];
	z[53] = z[19]*z[32];
	z[54] = z[19]*z[19];
	z[55] = 2.*z[91]*z[99];
	z[56] = 2.*z[102]*z[91];
	z[57] = 2.*z[9]*z[91];
	z[58] = z[26]*z[91];
	z[59] = z[91]*z[91];
	z[60] = 2.*z[100]*z[92];
	z[61] = 2.*z[13]*z[92];
	z[62] = 2.*z[11]*z[92];
	z[63] = z[29]*z[92];
	z[64] = z[92]*z[92];
	z[65] = 2.*z[101]*z[93];
	z[66] = 2.*z[103]*z[93];
	z[67] = 2.*z[6]*z[93];
	z[68] = z[32]*z[93];
	z[69] = 2.*z[19]*z[93];
	z[70] = z[93]*z[93];
	z[71] = 2.*z[100]*z[94];
	z[72] = 2.*z[13]*z[94];
	z[73] = 2.*z[11]*z[94];
	z[74] = z[29]*z[94];
	z[75] = z[94]*z[94];
	z[76] = 2.*z[102]*z[95];
	z[77] = 2.*z[9]*z[95];
	z[78] = z[26]*z[95];
	z[79] = 2.*z[91]*z[95];
	z[80] = z[95]*z[95];
	z[81] = 2.*z[96]*z[99];
	z[82] = 2.*z[102]*z[96];
	z[83] = 2.*z[9]*z[96];
	z[84] = z[26]*z[96];
	z[85] = 2.*z[95]*z[96];
	z[86] = z[96]*z[96];
	z[87] = 2.*z[100]*z[97];
	z[88] = 2.*z[11]*z[97];
	z[89] = z[29]*z[97];
	z[90] = 2.*z[92]*z[97];
	z[104] = 2.*z[94]*z[97];
	z[105] = z[97]*z[97];
	z[10] = 2.*z[10];
	z[12] = 2.*z[12];
	z[14] = 2.*z[14];
	z[15] = 2.*z[15];
	z[17] = 2.*z[17];
	z[18] = 2.*z[18];
	z[2] = pd2 + pX2 + z[2];
	z[22] = 2.*z[22];
	z[24] = 2.*z[24];
	z[3] = pd3 + pX3 + z[3];
	z[21] = z[21] + z[32];
	z[27] = z[27] + z[30] + z[33];
	z[30] = z[25] + z[34];
	z[33] = 2.*z[35];
	z[34] = z[28] + z[36];
	z[35] = 2.*z[37];
	z[1] = pd1 + pX1 + z[1];
	z[36] = z[32] + z[38];
	z[37] = z[31] + z[40];
	z[38] = z[42] + z[46];
	z[40] = z[41] + z[47];
	z[41] = z[44] + z[50];
	z[42] = z[31] + z[51];
	z[44] = 2.*z[54];
	z[23] = z[23] + z[55];
	z[46] = z[26] + z[57];
	z[45] = z[45] + z[58];
	z[47] = 2.*z[59];
	z[50] = z[56] + z[60];
	z[16] = z[16] + z[61];
	z[51] = z[29] + z[62];
	z[43] = z[43] + z[63];
	z[54] = 2.*z[64];
	z[31] = z[31] + z[66];
	z[32] = z[32] + z[69];
	z[55] = 2.*z[70];
	z[56] = z[28] + z[71];
	z[57] = z[29] + z[72];
	z[4] = z[4] + z[73];
	z[53] = z[53] + z[74];
	z[58] = 2.*z[75];
	z[59] = z[26] + z[76];
	z[20] = z[20] + z[77];
	z[48] = z[48] + z[78];
	z[60] = z[25] + z[79];
	z[61] = 2.*z[80];
	z[26] = z[26] + z[81];
	z[25] = z[25] + z[82];
	z[62] = z[67] + z[83];
	z[63] = z[68] + z[84];
	z[64] = z[65] + z[85];
	z[65] = 2.*z[86];
	z[29] = z[29] + z[87];
	z[39] = z[39] + z[88];
	z[49] = z[49] + z[89];
	z[28] = z[28] + z[90];
	z[52] = z[104] + z[52];
	z[66] = 2.*z[105];
	z[12] = z[12] + z[24];
	z[24] = pow(z[27],-2.5);
	z[67] = pow(z[27],-1.5);
	z[27] = 1./sqrt(z[27]);
	z[22] = z[22] + z[33];
	z[10] = z[10] + z[35];
	z[33] = z[38]*z[38];
	z[35] = z[40]*z[40];
	z[68] = z[41]*z[41];
	z[69] = z[45]*z[45];
	z[14] = z[14] + z[47];
	z[47] = z[43]*z[43];
	z[17] = z[17] + z[54];
	z[54] = z[53]*z[53];
	z[44] = z[44] + z[58];
	z[58] = z[48]*z[48];
	z[15] = z[15] + z[61];
	z[61] = z[63]*z[63];
	z[55] = z[55] + z[65];
	z[65] = z[49]*z[49];
	z[18] = z[18] + z[66];
	z[24] = 0.75*z[24];
	z[66] = -z[101]*z[103]*z[2]*z[5]*z[67];
	z[70] = -z[102]*z[2]*z[5]*z[67]*z[9];
	z[71] = -z[100]*z[11]*z[2]*z[5]*z[67];
	z[72] = -z[101]*z[103]*z[3]*z[67]*z[7];
	z[73] = -z[102]*z[3]*z[67]*z[7]*z[9];
	z[74] = -z[100]*z[11]*z[3]*z[67]*z[7];
	z[75] = -0.5*z[2]*z[21]*z[5]*z[67];
	z[76] = -0.5*z[21]*z[3]*z[67]*z[7];
	z[77] = -0.5*z[2]*z[30]*z[5]*z[67];
	z[78] = -0.5*z[3]*z[30]*z[67]*z[7];
	z[79] = -0.5*z[2]*z[34]*z[5]*z[67];
	z[80] = -0.5*z[3]*z[34]*z[67]*z[7];
	z[81] = -z[1]*z[101]*z[103]*z[67]*z[8];
	z[82] = -z[1]*z[102]*z[67]*z[8]*z[9];
	z[83] = -z[1]*z[100]*z[11]*z[67]*z[8];
	z[21] = -0.5*z[1]*z[21]*z[67]*z[8];
	z[30] = -0.5*z[1]*z[30]*z[67]*z[8];
	z[34] = -0.5*z[1]*z[34]*z[67]*z[8];
	z[84] = -0.5*z[2]*z[36]*z[5]*z[67];
	z[85] = -0.5*z[3]*z[36]*z[67]*z[7];
	z[36] = -0.5*z[1]*z[36]*z[67]*z[8];
	z[86] = -0.5*z[2]*z[37]*z[5]*z[67];
	z[87] = -0.5*z[3]*z[37]*z[67]*z[7];
	z[37] = -0.5*z[1]*z[37]*z[67]*z[8];
	z[88] = -0.5*z[38]*z[67]*z[7];
	z[89] = 0.16666666666666666*z[38]*z[67]*z[7];
	z[90] = -0.5*z[38]*z[5]*z[67];
	z[104] = 0.16666666666666666*z[38]*z[5]*z[67];
	z[105] = -0.5*z[38]*z[67]*z[8];
	z[106] = 0.16666666666666666*z[38]*z[67]*z[8];
	z[107] = 0.3333333333333333*z[38]*z[67]*z[8];
	z[108] = -0.5*z[100]*z[2]*z[38]*z[67];
	z[109] = -z[13]*z[2]*z[38]*z[67];
	z[110] = -0.5*z[11]*z[2]*z[38]*z[67];
	z[111] = -z[3]*z[38]*z[67]*z[99];
	z[112] = -0.5*z[102]*z[3]*z[38]*z[67];
	z[113] = -0.5*z[3]*z[38]*z[67]*z[9];
	z[114] = -0.5*z[1]*z[38]*z[67]*z[98];
	z[115] = -0.5*z[1]*z[101]*z[38]*z[67];
	z[116] = -0.5*z[1]*z[103]*z[38]*z[67];
	z[117] = -0.5*z[1]*z[38]*z[6]*z[67];
	z[118] = -0.5*z[40]*z[67]*z[7];
	z[119] = 0.16666666666666666*z[40]*z[67]*z[7];
	z[120] = 0.3333333333333333*z[40]*z[67]*z[7];
	z[121] = -0.5*z[40]*z[5]*z[67];
	z[122] = 0.16666666666666666*z[40]*z[5]*z[67];
	z[123] = -0.5*z[40]*z[67]*z[8];
	z[124] = 0.16666666666666666*z[40]*z[67]*z[8];
	z[125] = -0.5*z[100]*z[2]*z[40]*z[67];
	z[126] = -0.5*z[13]*z[2]*z[40]*z[67];
	z[127] = -z[11]*z[2]*z[40]*z[67];
	z[128] = -0.5*z[3]*z[40]*z[67]*z[99];
	z[129] = -0.5*z[102]*z[3]*z[40]*z[67];
	z[130] = -0.5*z[3]*z[40]*z[67]*z[9];
	z[131] = -z[1]*z[40]*z[67]*z[98];
	z[132] = -0.5*z[1]*z[101]*z[40]*z[67];
	z[133] = -0.5*z[1]*z[103]*z[40]*z[67];
	z[134] = -0.5*z[1]*z[40]*z[6]*z[67];
	z[135] = -0.5*z[41]*z[67]*z[7];
	z[136] = 0.16666666666666666*z[41]*z[67]*z[7];
	z[137] = -0.5*z[41]*z[5]*z[67];
	z[138] = 0.16666666666666666*z[41]*z[5]*z[67];
	z[139] = 0.3333333333333333*z[41]*z[5]*z[67];
	z[140] = -0.5*z[41]*z[67]*z[8];
	z[141] = 0.16666666666666666*z[41]*z[67]*z[8];
	z[142] = -0.5*z[100]*z[2]*z[41]*z[67];
	z[143] = -0.5*z[13]*z[2]*z[41]*z[67];
	z[144] = -0.5*z[11]*z[2]*z[41]*z[67];
	z[145] = -0.5*z[3]*z[41]*z[67]*z[99];
	z[146] = -0.5*z[102]*z[3]*z[41]*z[67];
	z[147] = -z[3]*z[41]*z[67]*z[9];
	z[148] = -0.5*z[1]*z[41]*z[67]*z[98];
	z[149] = -0.5*z[1]*z[101]*z[41]*z[67];
	z[150] = -0.5*z[1]*z[103]*z[41]*z[67];
	z[151] = -z[1]*z[41]*z[6]*z[67];
	z[152] = -0.5*z[2]*z[42]*z[5]*z[67];
	z[153] = -0.5*z[3]*z[42]*z[67]*z[7];
	z[42] = -0.5*z[1]*z[42]*z[67]*z[8];
	z[154] = -0.5*z[2]*z[23]*z[5]*z[67];
	z[155] = -0.5*z[23]*z[3]*z[67]*z[7];
	z[23] = -0.5*z[1]*z[23]*z[67]*z[8];
	z[156] = -0.5*z[2]*z[46]*z[5]*z[67];
	z[157] = -0.5*z[3]*z[46]*z[67]*z[7];
	z[46] = -0.5*z[1]*z[46]*z[67]*z[8];
	z[158] = -0.5*z[45]*z[67]*z[7];
	z[159] = 0.16666666666666666*z[45]*z[67]*z[7];
	z[160] = -0.5*z[45]*z[5]*z[67];
	z[161] = 0.16666666666666666*z[45]*z[5]*z[67];
	z[162] = -0.5*z[45]*z[67]*z[8];
	z[163] = 0.16666666666666666*z[45]*z[67]*z[8];
	z[164] = 0.3333333333333333*z[45]*z[67]*z[8];
	z[165] = -z[100]*z[2]*z[45]*z[67];
	z[166] = -0.5*z[13]*z[2]*z[45]*z[67];
	z[167] = -0.5*z[11]*z[2]*z[45]*z[67];
	z[168] = -0.5*z[3]*z[45]*z[67]*z[99];
	z[169] = -0.5*z[102]*z[3]*z[45]*z[67];
	z[170] = -0.5*z[3]*z[45]*z[67]*z[9];
	z[171] = -0.5*z[1]*z[45]*z[67]*z[98];
	z[172] = -0.5*z[1]*z[101]*z[45]*z[67];
	z[173] = -0.5*z[1]*z[103]*z[45]*z[67];
	z[174] = -0.5*z[1]*z[45]*z[6]*z[67];
	z[175] = -0.5*z[2]*z[5]*z[50]*z[67];
	z[176] = -0.5*z[3]*z[50]*z[67]*z[7];
	z[50] = -0.5*z[1]*z[50]*z[67]*z[8];
	z[177] = -0.5*z[16]*z[2]*z[5]*z[67];
	z[178] = -0.5*z[16]*z[3]*z[67]*z[7];
	z[16] = -0.5*z[1]*z[16]*z[67]*z[8];
	z[179] = -0.5*z[2]*z[5]*z[51]*z[67];
	z[180] = -0.5*z[3]*z[51]*z[67]*z[7];
	z[51] = -0.5*z[1]*z[51]*z[67]*z[8];
	z[181] = -0.5*z[43]*z[67]*z[7];
	z[182] = 0.16666666666666666*z[43]*z[67]*z[7];
	z[183] = -0.5*z[43]*z[5]*z[67];
	z[184] = 0.16666666666666666*z[43]*z[5]*z[67];
	z[185] = -0.5*z[43]*z[67]*z[8];
	z[186] = 0.16666666666666666*z[43]*z[67]*z[8];
	z[187] = 0.3333333333333333*z[43]*z[67]*z[8];
	z[188] = -0.5*z[100]*z[2]*z[43]*z[67];
	z[189] = -0.5*z[13]*z[2]*z[43]*z[67];
	z[190] = -0.5*z[11]*z[2]*z[43]*z[67];
	z[191] = -0.5*z[3]*z[43]*z[67]*z[99];
	z[192] = -z[102]*z[3]*z[43]*z[67];
	z[193] = -0.5*z[3]*z[43]*z[67]*z[9];
	z[194] = -0.5*z[1]*z[43]*z[67]*z[98];
	z[195] = -0.5*z[1]*z[101]*z[43]*z[67];
	z[196] = -0.5*z[1]*z[103]*z[43]*z[67];
	z[197] = -0.5*z[1]*z[43]*z[6]*z[67];
	z[198] = -0.5*z[2]*z[31]*z[5]*z[67];
	z[199] = -0.5*z[3]*z[31]*z[67]*z[7];
	z[31] = -0.5*z[1]*z[31]*z[67]*z[8];
	z[200] = -0.5*z[2]*z[32]*z[5]*z[67];
	z[201] = -0.5*z[3]*z[32]*z[67]*z[7];
	z[32] = -0.5*z[1]*z[32]*z[67]*z[8];
	z[202] = -0.5*z[2]*z[5]*z[56]*z[67];
	z[203] = -0.5*z[3]*z[56]*z[67]*z[7];
	z[56] = -0.5*z[1]*z[56]*z[67]*z[8];
	z[204] = -0.5*z[2]*z[5]*z[57]*z[67];
	z[205] = -0.5*z[3]*z[57]*z[67]*z[7];
	z[57] = -0.5*z[1]*z[57]*z[67]*z[8];
	z[206] = -0.5*z[2]*z[4]*z[5]*z[67];
	z[207] = -0.5*z[3]*z[4]*z[67]*z[7];
	z[4] = -0.5*z[1]*z[4]*z[67]*z[8];
	z[208] = -0.5*z[53]*z[67]*z[7];
	z[209] = 0.16666666666666666*z[53]*z[67]*z[7];
	z[210] = 0.3333333333333333*z[53]*z[67]*z[7];
	z[211] = -0.5*z[5]*z[53]*z[67];
	z[212] = 0.16666666666666666*z[5]*z[53]*z[67];
	z[213] = -0.5*z[53]*z[67]*z[8];
	z[214] = 0.16666666666666666*z[53]*z[67]*z[8];
	z[215] = -0.5*z[100]*z[2]*z[53]*z[67];
	z[216] = -0.5*z[13]*z[2]*z[53]*z[67];
	z[217] = -0.5*z[11]*z[2]*z[53]*z[67];
	z[218] = -0.5*z[3]*z[53]*z[67]*z[99];
	z[219] = -0.5*z[102]*z[3]*z[53]*z[67];
	z[220] = -0.5*z[3]*z[53]*z[67]*z[9];
	z[221] = -0.5*z[1]*z[53]*z[67]*z[98];
	z[222] = -0.5*z[1]*z[101]*z[53]*z[67];
	z[223] = -0.5*z[1]*z[103]*z[53]*z[67];
	z[224] = -0.5*z[1]*z[53]*z[6]*z[67];
	z[225] = -0.5*z[2]*z[5]*z[59]*z[67];
	z[226] = -0.5*z[3]*z[59]*z[67]*z[7];
	z[59] = -0.5*z[1]*z[59]*z[67]*z[8];
	z[227] = -0.5*z[2]*z[20]*z[5]*z[67];
	z[228] = -0.5*z[20]*z[3]*z[67]*z[7];
	z[20] = -0.5*z[1]*z[20]*z[67]*z[8];
	z[229] = -0.5*z[48]*z[67]*z[7];
	z[230] = 0.16666666666666666*z[48]*z[67]*z[7];
	z[231] = -0.5*z[48]*z[5]*z[67];
	z[232] = 0.16666666666666666*z[48]*z[5]*z[67];
	z[233] = 0.3333333333333333*z[48]*z[5]*z[67];
	z[234] = -0.5*z[48]*z[67]*z[8];
	z[235] = 0.16666666666666666*z[48]*z[67]*z[8];
	z[236] = -0.5*z[100]*z[2]*z[48]*z[67];
	z[237] = -0.5*z[13]*z[2]*z[48]*z[67];
	z[238] = -0.5*z[11]*z[2]*z[48]*z[67];
	z[239] = -0.5*z[3]*z[48]*z[67]*z[99];
	z[240] = -0.5*z[102]*z[3]*z[48]*z[67];
	z[241] = -0.5*z[3]*z[48]*z[67]*z[9];
	z[242] = -0.5*z[1]*z[48]*z[67]*z[98];
	z[243] = -z[1]*z[101]*z[48]*z[67];
	z[244] = -0.5*z[1]*z[103]*z[48]*z[67];
	z[245] = -0.5*z[1]*z[48]*z[6]*z[67];
	z[246] = -0.5*z[2]*z[5]*z[60]*z[67];
	z[247] = -0.5*z[3]*z[60]*z[67]*z[7];
	z[60] = -0.5*z[1]*z[60]*z[67]*z[8];
	z[248] = -0.5*z[2]*z[26]*z[5]*z[67];
	z[249] = -0.5*z[26]*z[3]*z[67]*z[7];
	z[26] = -0.5*z[1]*z[26]*z[67]*z[8];
	z[250] = -0.5*z[2]*z[25]*z[5]*z[67];
	z[251] = -0.5*z[25]*z[3]*z[67]*z[7];
	z[25] = -0.5*z[1]*z[25]*z[67]*z[8];
	z[252] = -0.5*z[2]*z[5]*z[62]*z[67];
	z[253] = -0.5*z[3]*z[62]*z[67]*z[7];
	z[62] = -0.5*z[1]*z[62]*z[67]*z[8];
	z[254] = -0.5*z[63]*z[67]*z[7];
	z[255] = 0.16666666666666666*z[63]*z[67]*z[7];
	z[256] = -0.5*z[5]*z[63]*z[67];
	z[257] = 0.16666666666666666*z[5]*z[63]*z[67];
	z[258] = 0.3333333333333333*z[5]*z[63]*z[67];
	z[259] = -0.5*z[63]*z[67]*z[8];
	z[260] = 0.16666666666666666*z[63]*z[67]*z[8];
	z[261] = -0.5*z[100]*z[2]*z[63]*z[67];
	z[262] = -0.5*z[13]*z[2]*z[63]*z[67];
	z[263] = -0.5*z[11]*z[2]*z[63]*z[67];
	z[264] = -0.5*z[3]*z[63]*z[67]*z[99];
	z[265] = -0.5*z[102]*z[3]*z[63]*z[67];
	z[266] = -0.5*z[3]*z[63]*z[67]*z[9];
	z[267] = -0.5*z[1]*z[63]*z[67]*z[98];
	z[268] = -0.5*z[1]*z[101]*z[63]*z[67];
	z[269] = -0.5*z[1]*z[103]*z[63]*z[67];
	z[270] = -0.5*z[1]*z[6]*z[63]*z[67];
	z[271] = -0.5*z[2]*z[5]*z[64]*z[67];
	z[272] = -0.5*z[3]*z[64]*z[67]*z[7];
	z[64] = -0.5*z[1]*z[64]*z[67]*z[8];
	z[273] = -0.5*z[2]*z[29]*z[5]*z[67];
	z[274] = -0.5*z[29]*z[3]*z[67]*z[7];
	z[29] = -0.5*z[1]*z[29]*z[67]*z[8];
	z[275] = -0.5*z[2]*z[39]*z[5]*z[67];
	z[276] = -0.5*z[3]*z[39]*z[67]*z[7];
	z[39] = -0.5*z[1]*z[39]*z[67]*z[8];
	z[277] = -0.5*z[49]*z[67]*z[7];
	z[278] = 0.16666666666666666*z[49]*z[67]*z[7];
	z[279] = 0.3333333333333333*z[49]*z[67]*z[7];
	z[280] = -0.5*z[49]*z[5]*z[67];
	z[281] = 0.16666666666666666*z[49]*z[5]*z[67];
	z[282] = -0.5*z[49]*z[67]*z[8];
	z[283] = 0.16666666666666666*z[49]*z[67]*z[8];
	z[284] = -0.5*z[100]*z[2]*z[49]*z[67];
	z[285] = -0.5*z[13]*z[2]*z[49]*z[67];
	z[286] = -0.5*z[11]*z[2]*z[49]*z[67];
	z[287] = -0.5*z[3]*z[49]*z[67]*z[99];
	z[288] = -0.5*z[102]*z[3]*z[49]*z[67];
	z[289] = -0.5*z[3]*z[49]*z[67]*z[9];
	z[290] = -0.5*z[1]*z[49]*z[67]*z[98];
	z[291] = -0.5*z[1]*z[101]*z[49]*z[67];
	z[292] = -z[1]*z[103]*z[49]*z[67];
	z[293] = -0.5*z[1]*z[49]*z[6]*z[67];
	z[294] = -0.5*z[2]*z[28]*z[5]*z[67];
	z[295] = -0.5*z[28]*z[3]*z[67]*z[7];
	z[28] = -0.5*z[1]*z[28]*z[67]*z[8];
	z[296] = -0.5*z[2]*z[5]*z[52]*z[67];
	z[297] = -0.5*z[3]*z[52]*z[67]*z[7];
	z[52] = -0.5*z[1]*z[52]*z[67]*z[8];
	z[298] = -0.5*z[12]*z[2]*z[5]*z[67];
	z[299] = -0.5*z[12]*z[3]*z[67]*z[7];
	z[12] = -0.5*z[1]*z[12]*z[67]*z[8];
	z[300] = -0.3333333333333333*z[27]*z[98];
	z[301] = z[27]*z[98];
	z[302] = -0.3333333333333333*z[27]*z[99];
	z[303] = z[27]*z[99];
	z[304] = -0.3333333333333333*z[100]*z[27];
	z[100] = z[100]*z[27];
	z[305] = -0.3333333333333333*z[101]*z[27];
	z[101] = z[101]*z[27];
	z[306] = -0.3333333333333333*z[102]*z[27];
	z[102] = z[102]*z[27];
	z[307] = -0.3333333333333333*z[103]*z[27];
	z[103] = z[103]*z[27];
	z[308] = -0.3333333333333333*z[27]*z[6];
	z[309] = z[27]*z[6];
	z[310] = -0.3333333333333333*z[13]*z[27];
	z[311] = z[13]*z[27];
	z[312] = -0.3333333333333333*z[27]*z[9];
	z[9] = z[27]*z[9];
	z[313] = -0.3333333333333333*z[11]*z[27];
	z[11] = z[11]*z[27];
	z[314] = -z[2]*z[27];
	z[315] = z[2]*z[27];
	z[316] = -z[27]*z[3];
	z[317] = z[27]*z[3];
	z[318] = -z[1]*z[27];
	z[319] = z[1]*z[27];
	z[320] = -0.5*z[2]*z[22]*z[5]*z[67];
	z[321] = -0.5*z[22]*z[3]*z[67]*z[7];
	z[22] = -0.5*z[1]*z[22]*z[67]*z[8];
	z[322] = -0.5*z[10]*z[2]*z[5]*z[67];
	z[323] = -0.5*z[10]*z[3]*z[67]*z[7];
	z[10] = -0.5*z[1]*z[10]*z[67]*z[8];
	z[324] = -0.5*z[14]*z[2]*z[5]*z[67];
	z[325] = -0.5*z[14]*z[3]*z[67]*z[7];
	z[14] = -0.5*z[1]*z[14]*z[67]*z[8];
	z[326] = -0.5*z[17]*z[2]*z[5]*z[67];
	z[327] = -0.5*z[17]*z[3]*z[67]*z[7];
	z[17] = -0.5*z[1]*z[17]*z[67]*z[8];
	z[328] = -0.5*z[2]*z[44]*z[5]*z[67];
	z[329] = -0.5*z[3]*z[44]*z[67]*z[7];
	z[44] = -0.5*z[1]*z[44]*z[67]*z[8];
	z[330] = -0.5*z[15]*z[2]*z[5]*z[67];
	z[331] = -0.5*z[15]*z[3]*z[67]*z[7];
	z[15] = -0.5*z[1]*z[15]*z[67]*z[8];
	z[332] = -0.5*z[2]*z[5]*z[55]*z[67];
	z[333] = -0.5*z[3]*z[55]*z[67]*z[7];
	z[55] = -0.5*z[1]*z[55]*z[67]*z[8];
	z[334] = -0.5*z[18]*z[2]*z[5]*z[67];
	z[335] = -0.5*z[18]*z[3]*z[67]*z[7];
	z[18] = -0.5*z[1]*z[18]*z[67]*z[8];
	z[336] = z[2]*z[24]*z[5];
	z[337] = z[24]*z[38];
	z[338] = z[24]*z[40];
	z[339] = z[24]*z[3]*z[41]*z[45]*z[7];
	z[340] = z[1]*z[24]*z[41]*z[45]*z[8];
	z[341] = z[24]*z[3]*z[41]*z[43]*z[7];
	z[342] = z[1]*z[24]*z[41]*z[43]*z[8];
	z[343] = z[24]*z[3]*z[43]*z[45]*z[7];
	z[344] = z[1]*z[24]*z[43]*z[45]*z[8];
	z[345] = z[24]*z[3]*z[41]*z[53]*z[7];
	z[346] = z[1]*z[24]*z[41]*z[53]*z[8];
	z[347] = z[24]*z[3]*z[45]*z[53]*z[7];
	z[348] = z[1]*z[24]*z[45]*z[53]*z[8];
	z[349] = z[24]*z[3]*z[43]*z[53]*z[7];
	z[350] = z[1]*z[24]*z[43]*z[53]*z[8];
	z[351] = z[24]*z[3]*z[41]*z[48]*z[7];
	z[352] = z[1]*z[24]*z[41]*z[48]*z[8];
	z[353] = z[24]*z[3]*z[45]*z[48]*z[7];
	z[354] = z[1]*z[24]*z[45]*z[48]*z[8];
	z[355] = z[24]*z[3]*z[43]*z[48]*z[7];
	z[356] = z[1]*z[24]*z[43]*z[48]*z[8];
	z[357] = z[24]*z[3]*z[48]*z[53]*z[7];
	z[358] = z[1]*z[24]*z[48]*z[53]*z[8];
	z[359] = z[24]*z[3]*z[41]*z[63]*z[7];
	z[360] = z[1]*z[24]*z[41]*z[63]*z[8];
	z[361] = z[24]*z[3]*z[45]*z[63]*z[7];
	z[362] = z[1]*z[24]*z[45]*z[63]*z[8];
	z[363] = z[24]*z[3]*z[43]*z[63]*z[7];
	z[364] = z[1]*z[24]*z[43]*z[63]*z[8];
	z[365] = z[24]*z[3]*z[53]*z[63]*z[7];
	z[366] = z[1]*z[24]*z[53]*z[63]*z[8];
	z[367] = z[24]*z[3]*z[48]*z[63]*z[7];
	z[368] = z[1]*z[24]*z[48]*z[63]*z[8];
	z[369] = z[24]*z[3]*z[41]*z[49]*z[7];
	z[370] = z[1]*z[24]*z[41]*z[49]*z[8];
	z[371] = z[24]*z[3]*z[45]*z[49]*z[7];
	z[372] = z[1]*z[24]*z[45]*z[49]*z[8];
	z[373] = z[24]*z[3]*z[43]*z[49]*z[7];
	z[374] = z[1]*z[24]*z[43]*z[49]*z[8];
	z[375] = z[24]*z[3]*z[49]*z[53]*z[7];
	z[376] = z[1]*z[24]*z[49]*z[53]*z[8];
	z[377] = z[24]*z[3]*z[48]*z[49]*z[7];
	z[378] = z[1]*z[24]*z[48]*z[49]*z[8];
	z[379] = z[24]*z[3]*z[49]*z[63]*z[7];
	z[380] = z[1]*z[24]*z[49]*z[63]*z[8];
	z[381] = z[24]*z[3]*z[33]*z[7];
	z[382] = z[1]*z[24]*z[33]*z[8];
	z[383] = z[24]*z[3]*z[35]*z[7];
	z[384] = z[1]*z[24]*z[35]*z[8];
	z[385] = z[24]*z[3]*z[68]*z[7];
	z[386] = z[1]*z[24]*z[68]*z[8];
	z[387] = z[24]*z[3]*z[69]*z[7];
	z[388] = z[1]*z[24]*z[69]*z[8];
	z[389] = z[24]*z[3]*z[47]*z[7];
	z[390] = z[1]*z[24]*z[47]*z[8];
	z[391] = z[24]*z[3]*z[54]*z[7];
	z[392] = z[1]*z[24]*z[54]*z[8];
	z[393] = z[24]*z[3]*z[58]*z[7];
	z[394] = z[1]*z[24]*z[58]*z[8];
	z[395] = z[24]*z[3]*z[61]*z[7];
	z[396] = z[1]*z[24]*z[61]*z[8];
	z[397] = z[24]*z[3]*z[65]*z[7];
	z[24] = z[1]*z[24]*z[65]*z[8];
	z[398] = z[336]*z[38];
	z[399] = z[336]*z[40];
	z[400] = z[336]*z[41]*z[45];
	z[401] = z[336]*z[41]*z[43];
	z[402] = z[336]*z[43]*z[45];
	z[403] = z[336]*z[41]*z[53];
	z[404] = z[336]*z[45]*z[53];
	z[405] = z[336]*z[43]*z[53];
	z[406] = z[336]*z[41]*z[48];
	z[407] = z[336]*z[45]*z[48];
	z[408] = z[336]*z[43]*z[48];
	z[409] = z[336]*z[48]*z[53];
	z[410] = z[336]*z[41]*z[63];
	z[411] = z[336]*z[45]*z[63];
	z[412] = z[336]*z[43]*z[63];
	z[413] = z[336]*z[53]*z[63];
	z[414] = z[336]*z[48]*z[63];
	z[415] = z[336]*z[41]*z[49];
	z[416] = z[336]*z[45]*z[49];
	z[417] = z[336]*z[43]*z[49];
	z[418] = z[336]*z[49]*z[53];
	z[419] = z[336]*z[48]*z[49];
	z[420] = z[336]*z[49]*z[63];
	z[33] = z[33]*z[336];
	z[35] = z[336]*z[35];
	z[68] = z[336]*z[68];
	z[69] = z[336]*z[69];
	z[47] = z[336]*z[47];
	z[54] = z[336]*z[54];
	z[58] = z[336]*z[58];
	z[61] = z[336]*z[61];
	z[65] = z[336]*z[65];
	z[336] = z[337]*z[40];
	z[421] = z[3]*z[337]*z[41]*z[7];
	z[422] = z[1]*z[337]*z[41]*z[8];
	z[423] = z[3]*z[337]*z[45]*z[7];
	z[424] = z[1]*z[337]*z[45]*z[8];
	z[425] = z[3]*z[337]*z[43]*z[7];
	z[426] = z[1]*z[337]*z[43]*z[8];
	z[427] = z[3]*z[337]*z[53]*z[7];
	z[428] = z[1]*z[337]*z[53]*z[8];
	z[429] = z[3]*z[337]*z[48]*z[7];
	z[430] = z[1]*z[337]*z[48]*z[8];
	z[431] = z[3]*z[337]*z[63]*z[7];
	z[432] = z[1]*z[337]*z[63]*z[8];
	z[433] = z[3]*z[337]*z[49]*z[7];
	z[337] = z[1]*z[337]*z[49]*z[8];
	z[434] = z[3]*z[338]*z[41]*z[7];
	z[435] = z[1]*z[338]*z[41]*z[8];
	z[436] = z[3]*z[338]*z[45]*z[7];
	z[437] = z[1]*z[338]*z[45]*z[8];
	z[438] = z[3]*z[338]*z[43]*z[7];
	z[439] = z[1]*z[338]*z[43]*z[8];
	z[440] = z[3]*z[338]*z[53]*z[7];
	z[441] = z[1]*z[338]*z[53]*z[8];
	z[442] = z[3]*z[338]*z[48]*z[7];
	z[443] = z[1]*z[338]*z[48]*z[8];
	z[444] = z[3]*z[338]*z[63]*z[7];
	z[445] = z[1]*z[338]*z[63]*z[8];
	z[446] = z[3]*z[338]*z[49]*z[7];
	z[338] = z[1]*z[338]*z[49]*z[8];
	z[447] = z[398]*z[40];
	z[448] = z[398]*z[41];
	z[449] = z[398]*z[45];
	z[450] = z[398]*z[43];
	z[451] = z[398]*z[53];
	z[452] = z[398]*z[48];
	z[453] = z[398]*z[63];
	z[398] = z[398]*z[49];
	z[454] = z[399]*z[41];
	z[455] = z[399]*z[45];
	z[456] = z[399]*z[43];
	z[457] = z[399]*z[53];
	z[458] = z[399]*z[48];
	z[459] = z[399]*z[63];
	z[399] = z[399]*z[49];
	z[460] = z[3]*z[336]*z[7];
	z[336] = z[1]*z[336]*z[8];
	z[461] = -z[19]*z[2]*z[5]*z[6]*z[67];
	z[462] = -z[19]*z[3]*z[6]*z[67]*z[7];
	z[6] = -z[1]*z[19]*z[6]*z[67]*z[8];
	z[463] = -0.5*z[1]*z[19]*z[38]*z[67];
	z[464] = -0.5*z[1]*z[19]*z[40]*z[67];
	z[465] = -0.5*z[1]*z[19]*z[41]*z[67];
	z[466] = -0.5*z[1]*z[19]*z[45]*z[67];
	z[467] = -0.5*z[1]*z[19]*z[43]*z[67];
	z[468] = -z[1]*z[19]*z[53]*z[67];
	z[469] = -0.5*z[1]*z[19]*z[48]*z[67];
	z[470] = -0.5*z[1]*z[19]*z[63]*z[67];
	z[471] = -0.5*z[1]*z[19]*z[49]*z[67];
	z[472] = -0.3333333333333333*z[19]*z[27];
	z[19] = z[19]*z[27];
	z[473] = -0.5*z[3]*z[38]*z[67]*z[91];
	z[474] = -0.5*z[3]*z[40]*z[67]*z[91];
	z[475] = -0.5*z[3]*z[41]*z[67]*z[91];
	z[476] = -z[3]*z[45]*z[67]*z[91];
	z[477] = -0.5*z[3]*z[43]*z[67]*z[91];
	z[478] = -0.5*z[3]*z[53]*z[67]*z[91];
	z[479] = -0.5*z[3]*z[48]*z[67]*z[91];
	z[480] = -0.5*z[3]*z[63]*z[67]*z[91];
	z[481] = -0.5*z[3]*z[49]*z[67]*z[91];
	z[482] = -0.3333333333333333*z[27]*z[91];
	z[483] = z[27]*z[91];
	z[484] = -0.5*z[2]*z[38]*z[67]*z[92];
	z[485] = -0.5*z[2]*z[40]*z[67]*z[92];
	z[486] = -0.5*z[2]*z[41]*z[67]*z[92];
	z[487] = -0.5*z[2]*z[45]*z[67]*z[92];
	z[488] = -z[2]*z[43]*z[67]*z[92];
	z[489] = -0.5*z[2]*z[53]*z[67]*z[92];
	z[490] = -0.5*z[2]*z[48]*z[67]*z[92];
	z[491] = -0.5*z[2]*z[63]*z[67]*z[92];
	z[492] = -0.5*z[2]*z[49]*z[67]*z[92];
	z[493] = -0.3333333333333333*z[27]*z[92];
	z[494] = z[27]*z[92];
	z[495] = -z[2]*z[5]*z[67]*z[93]*z[98];
	z[496] = -z[3]*z[67]*z[7]*z[93]*z[98];
	z[98] = -z[1]*z[67]*z[8]*z[93]*z[98];
	z[497] = -0.5*z[1]*z[38]*z[67]*z[93];
	z[498] = -0.5*z[1]*z[40]*z[67]*z[93];
	z[499] = -0.5*z[1]*z[41]*z[67]*z[93];
	z[500] = -0.5*z[1]*z[45]*z[67]*z[93];
	z[501] = -0.5*z[1]*z[43]*z[67]*z[93];
	z[502] = -0.5*z[1]*z[53]*z[67]*z[93];
	z[503] = -0.5*z[1]*z[48]*z[67]*z[93];
	z[504] = -z[1]*z[63]*z[67]*z[93];
	z[505] = -0.5*z[1]*z[49]*z[67]*z[93];
	z[506] = -0.3333333333333333*z[27]*z[93];
	z[93] = z[27]*z[93];
	z[507] = -0.5*z[2]*z[38]*z[67]*z[94];
	z[508] = -0.5*z[2]*z[40]*z[67]*z[94];
	z[509] = -0.5*z[2]*z[41]*z[67]*z[94];
	z[510] = -0.5*z[2]*z[45]*z[67]*z[94];
	z[511] = -0.5*z[2]*z[43]*z[67]*z[94];
	z[512] = -z[2]*z[53]*z[67]*z[94];
	z[513] = -0.5*z[2]*z[48]*z[67]*z[94];
	z[514] = -0.5*z[2]*z[63]*z[67]*z[94];
	z[515] = -0.5*z[2]*z[49]*z[67]*z[94];
	z[516] = -0.3333333333333333*z[27]*z[94];
	z[517] = z[27]*z[94];
	z[518] = -z[2]*z[5]*z[67]*z[92]*z[94];
	z[519] = -z[3]*z[67]*z[7]*z[92]*z[94];
	z[92] = -z[1]*z[67]*z[8]*z[92]*z[94];
	z[94] = -z[2]*z[5]*z[67]*z[95]*z[99];
	z[520] = -z[3]*z[67]*z[7]*z[95]*z[99];
	z[99] = -z[1]*z[67]*z[8]*z[95]*z[99];
	z[521] = -0.5*z[3]*z[38]*z[67]*z[95];
	z[522] = -0.5*z[3]*z[40]*z[67]*z[95];
	z[523] = -0.5*z[3]*z[41]*z[67]*z[95];
	z[524] = -0.5*z[3]*z[45]*z[67]*z[95];
	z[525] = -0.5*z[3]*z[43]*z[67]*z[95];
	z[526] = -0.5*z[3]*z[53]*z[67]*z[95];
	z[527] = -z[3]*z[48]*z[67]*z[95];
	z[528] = -0.5*z[3]*z[63]*z[67]*z[95];
	z[529] = -0.5*z[3]*z[49]*z[67]*z[95];
	z[530] = -0.3333333333333333*z[27]*z[95];
	z[95] = z[27]*z[95];
	z[531] = -0.5*z[3]*z[38]*z[67]*z[96];
	z[532] = -0.5*z[3]*z[40]*z[67]*z[96];
	z[533] = -0.5*z[3]*z[41]*z[67]*z[96];
	z[534] = -0.5*z[3]*z[45]*z[67]*z[96];
	z[535] = -0.5*z[3]*z[43]*z[67]*z[96];
	z[536] = -0.5*z[3]*z[53]*z[67]*z[96];
	z[537] = -0.5*z[3]*z[48]*z[67]*z[96];
	z[538] = -z[3]*z[63]*z[67]*z[96];
	z[539] = -0.5*z[3]*z[49]*z[67]*z[96];
	z[540] = -0.3333333333333333*z[27]*z[96];
	z[541] = z[27]*z[96];
	z[542] = -z[2]*z[5]*z[67]*z[91]*z[96];
	z[543] = -z[3]*z[67]*z[7]*z[91]*z[96];
	z[91] = -z[1]*z[67]*z[8]*z[91]*z[96];
	z[5] = -z[13]*z[2]*z[5]*z[67]*z[97];
	z[3] = -z[13]*z[3]*z[67]*z[7]*z[97];
	z[1] = -z[1]*z[13]*z[67]*z[8]*z[97];
	z[7] = -0.5*z[2]*z[38]*z[67]*z[97];
	z[8] = -0.5*z[2]*z[40]*z[67]*z[97];
	z[13] = -0.5*z[2]*z[41]*z[67]*z[97];
	z[38] = -0.5*z[2]*z[45]*z[67]*z[97];
	z[40] = -0.5*z[2]*z[43]*z[67]*z[97];
	z[41] = -0.5*z[2]*z[53]*z[67]*z[97];
	z[43] = -0.5*z[2]*z[48]*z[67]*z[97];
	z[45] = -0.5*z[2]*z[63]*z[67]*z[97];
	z[2] = -z[2]*z[49]*z[67]*z[97];
	z[48] = -0.3333333333333333*z[27]*z[97];
	z[27] = z[27]*z[97];
	z[49] = z[123] + z[301];
	z[53] = z[303] + z[88];
	z[63] = z[100] + z[160];
	z[67] = z[101] + z[234];
	z[88] = z[102] + z[181];
	z[96] = z[103] + z[282];
	z[97] = z[141] + z[308];
	z[100] = z[140] + z[309];
	z[90] = z[311] + z[90];
	z[9] = z[135] + z[9];
	z[11] = z[11] + z[121];
	z[12] = z[107] + z[109] + z[111] + z[12] + z[298] + z[299] + z[33] + z[381] + z[382];
	z[10] = z[10] + z[120] + z[127] + z[131] + z[322] + z[323] + z[35] + z[383] + z[384];
	z[22] = z[139] + z[147] + z[151] + z[22] + z[320] + z[321] + z[385] + z[386] + z[68];
	z[33] = z[122] + z[130] + z[134] + z[136] + z[144] + z[148] + z[312] + z[37] + z[86] + z[87];
	z[34] = z[110] + z[114] + z[124] + z[126] + z[128] + z[300] + z[34] + z[79] + z[80] + z[89];
	z[19] = z[19] + z[213];
	z[23] = z[106] + z[108] + z[154] + z[155] + z[163] + z[166] + z[168] + z[23] + z[423] + z[424];
	z[14] = z[14] + z[164] + z[165] + z[324] + z[325] + z[387] + z[388] + z[476] + z[69];
	z[35] = z[124] + z[125] + z[159] + z[167] + z[171] + z[300] + z[436] + z[71] + z[74] + z[83];
	z[37] = z[158] + z[483];
	z[16] = z[106] + z[112] + z[16] + z[177] + z[178] + z[186] + z[189] + z[191] + z[425] + z[426];
	z[51] = z[124] + z[129] + z[179] + z[180] + z[182] + z[190] + z[194] + z[300] + z[306] + z[51];
	z[50] = z[163] + z[169] + z[175] + z[176] + z[186] + z[188] + z[343] + z[344] + z[402] + z[50];
	z[17] = z[17] + z[187] + z[192] + z[326] + z[327] + z[389] + z[390] + z[47] + z[488];
	z[47] = z[183] + z[494];
	z[68] = z[259] + z[93];
	z[57] = z[204] + z[205] + z[214] + z[216] + z[218] + z[302] + z[315] + z[427] + z[57] + z[89];
	z[4] = z[119] + z[206] + z[207] + z[209] + z[217] + z[221] + z[4] + z[440] + z[441] + z[457];
	z[56] = z[159] + z[202] + z[203] + z[214] + z[215] + z[314] + z[347] + z[348] + z[404] + z[56];
	z[44] = z[210] + z[328] + z[329] + z[391] + z[392] + z[44] + z[468] + z[512] + z[54];
	z[54] = z[136] + z[212] + z[220] + z[224] + z[312] + z[345] + z[346] + z[403] + z[461] + z[462];
	z[69] = z[211] + z[517];
	z[71] = z[182] + z[214] + z[219] + z[306] + z[349] + z[350] + z[405] + z[467] + z[472] + z[489];
	z[74] = z[104] + z[115] + z[235] + z[237] + z[239] + z[305] + z[310] + z[429] + z[430] + z[452];
	z[20] = z[138] + z[149] + z[20] + z[227] + z[228] + z[232] + z[241] + z[245] + z[351] + z[352];
	z[60] = z[161] + z[172] + z[235] + z[236] + z[246] + z[247] + z[304] + z[305] + z[316] + z[60];
	z[59] = z[184] + z[195] + z[225] + z[226] + z[235] + z[240] + z[305] + z[317] + z[355] + z[59];
	z[15] = z[15] + z[233] + z[243] + z[330] + z[331] + z[393] + z[394] + z[527] + z[58];
	z[36] = z[122] + z[132] + z[230] + z[238] + z[242] + z[313] + z[319] + z[36] + z[84] + z[85];
	z[42] = z[152] + z[153] + z[212] + z[222] + z[230] + z[318] + z[357] + z[358] + z[409] + z[42];
	z[58] = z[229] + z[95];
	z[26] = z[104] + z[248] + z[249] + z[26] + z[260] + z[262] + z[264] + z[310] + z[317] + z[431];
	z[62] = z[138] + z[252] + z[253] + z[257] + z[266] + z[270] + z[359] + z[360] + z[410] + z[62];
	z[25] = z[184] + z[25] + z[250] + z[251] + z[260] + z[265] + z[316] + z[363] + z[364] + z[412];
	z[64] = z[232] + z[257] + z[268] + z[271] + z[272] + z[367] + z[368] + z[414] + z[503] + z[64];
	z[55] = z[258] + z[332] + z[333] + z[395] + z[396] + z[504] + z[538] + z[55] + z[61];
	z[61] = z[122] + z[255] + z[263] + z[267] + z[313] + z[444] + z[445] + z[459] + z[495] + z[496];
	z[32] = z[200] + z[201] + z[212] + z[255] + z[319] + z[32] + z[365] + z[366] + z[413] + z[470];
	z[79] = z[254] + z[541];
	z[80] = z[161] + z[260] + z[261] + z[304] + z[361] + z[362] + z[411] + z[480] + z[500] + z[506];
	z[83] = z[116] + z[283] + z[285] + z[287] + z[302] + z[307] + z[337] + z[398] + z[433] + z[89];
	z[39] = z[119] + z[133] + z[275] + z[276] + z[278] + z[286] + z[290] + z[338] + z[39] + z[446];
	z[29] = z[159] + z[173] + z[273] + z[274] + z[283] + z[284] + z[29] + z[307] + z[315] + z[371];
	z[28] = z[182] + z[196] + z[28] + z[283] + z[288] + z[294] + z[295] + z[306] + z[307] + z[314];
	z[52] = z[209] + z[223] + z[278] + z[296] + z[297] + z[375] + z[376] + z[418] + z[471] + z[52];
	z[2] = z[18] + z[2] + z[24] + z[279] + z[292] + z[334] + z[335] + z[397] + z[65];
	z[18] = z[136] + z[150] + z[21] + z[281] + z[289] + z[293] + z[312] + z[319] + z[75] + z[76];
	z[21] = z[230] + z[244] + z[281] + z[291] + z[377] + z[378] + z[419] + z[66] + z[72] + z[81];
	z[24] = z[198] + z[199] + z[255] + z[269] + z[281] + z[31] + z[318] + z[379] + z[380] + z[420];
	z[27] = z[27] + z[280];
	z[30] = z[104] + z[113] + z[117] + z[143] + z[145] + z[30] + z[310] + z[316] + z[77] + z[78];
	z[31] = z[142] + z[156] + z[157] + z[161] + z[170] + z[174] + z[304] + z[317] + z[339] + z[46];
	z[46] = z[146] + z[184] + z[193] + z[197] + z[341] + z[342] + z[401] + z[70] + z[73] + z[82];
	z[33] = z[313] + z[318] + z[33] + z[434] + z[435] + z[454];
	z[34] = z[302] + z[314] + z[336] + z[34] + z[447] + z[460];
	z[23] = z[23] + z[449] + z[473];
	z[35] = z[35] + z[437] + z[455] + z[474] + z[482];
	z[16] = z[16] + z[450] + z[484];
	z[51] = z[315] + z[438] + z[439] + z[456] + z[485] + z[51];
	z[50] = z[477] + z[487] + z[50];
	z[57] = z[428] + z[451] + z[463] + z[472] + z[507] + z[57];
	z[4] = z[4] + z[464] + z[508];
	z[56] = z[466] + z[472] + z[478] + z[482] + z[510] + z[56];
	z[6] = z[465] + z[509] + z[516] + z[54] + z[6];
	z[54] = z[511] + z[518] + z[519] + z[71] + z[92];
	z[65] = z[520] + z[521] + z[74] + z[94] + z[99];
	z[20] = z[20] + z[406] + z[523];
	z[60] = z[353] + z[354] + z[407] + z[479] + z[524] + z[60];
	z[59] = z[356] + z[408] + z[490] + z[493] + z[525] + z[59];
	z[36] = z[36] + z[442] + z[443] + z[458] + z[522] + z[530];
	z[42] = z[42] + z[469] + z[513] + z[516] + z[526] + z[530];
	z[26] = z[26] + z[432] + z[453] + z[497] + z[506] + z[531];
	z[62] = z[499] + z[533] + z[62];
	z[25] = z[25] + z[491] + z[493] + z[501] + z[506] + z[535];
	z[64] = z[528] + z[537] + z[64];
	z[61] = z[498] + z[532] + z[540] + z[61] + z[98];
	z[32] = z[32] + z[502] + z[514] + z[516] + z[536] + z[540];
	z[66] = z[534] + z[542] + z[543] + z[80] + z[91];
	z[1] = z[1] + z[3] + z[5] + z[7] + z[83];
	z[3] = z[39] + z[399] + z[8];
	z[5] = z[29] + z[372] + z[38] + z[416] + z[481] + z[482];
	z[7] = z[28] + z[373] + z[374] + z[40] + z[417] + z[492];
	z[8] = z[41] + z[515] + z[52];
	z[13] = z[13] + z[18] + z[369] + z[370] + z[415] + z[48];
	z[18] = z[21] + z[43] + z[48] + z[529] + z[530];
	z[21] = z[24] + z[45] + z[48] + z[505] + z[539] + z[540];
	z[24] = z[30] + z[421] + z[422] + z[448] + z[97];
	z[28] = z[31] + z[340] + z[400] + z[475] + z[97];
	z[29] = z[46] + z[486] + z[493] + z[97];
	
	/* write into K */	
	K(0,0) = z[12];
	K(1,0) = z[65];
	K(2,0) = z[1];
	K(3,0) = z[23];
	K(4,0) = z[26];
	K(5,0) = z[34];
	K(6,0) = z[16];
	K(7,0) = z[24];
	K(8,0) = z[57];
	K(9,0) = z[105];
	K(10,0) = z[90];
	K(11,0) = z[53];

	K(0,1) = z[65];
	K(1,1) = z[15];
	K(2,1) = z[18];
	K(3,1) = z[60];
	K(4,1) = z[64];
	K(5,1) = z[36];
	K(6,1) = z[59];
	K(7,1) = z[20];
	K(8,1) = z[42];
	K(9,1) = z[67];
	K(10,1) = z[231];
	K(11,1) = z[58];

	K(0,2) = z[1];
	K(1,2) = z[18];
	K(2,2) = z[2];
	K(3,2) = z[5];
	K(4,2) = z[21];
	K(5,2) = z[3];
	K(6,2) = z[7];
	K(7,2) = z[13];
	K(8,2) = z[8];
	K(9,2) = z[96];
	K(10,2) = z[27];
	K(11,2) = z[277];

	K(0,3) = z[23];
	K(1,3) = z[60];
	K(2,3) = z[5];
	K(3,3) = z[14];
	K(4,3) = z[66];
	K(5,3) = z[35];
	K(6,3) = z[50];
	K(7,3) = z[28];
	K(8,3) = z[56];
	K(9,3) = z[162];
	K(10,3) = z[63];
	K(11,3) = z[37];

	K(0,4) = z[26];
	K(1,4) = z[64];
	K(2,4) = z[21];
	K(3,4) = z[66];
	K(4,4) = z[55];
	K(5,4) = z[61];
	K(6,4) = z[25];
	K(7,4) = z[62];
	K(8,4) = z[32];
	K(9,4) = z[68];
	K(10,4) = z[256];
	K(11,4) = z[79];

	K(0,5) = z[34];
	K(1,5) = z[36];
	K(2,5) = z[3];
	K(3,5) = z[35];
	K(4,5) = z[61];
	K(5,5) = z[10];
	K(6,5) = z[51];
	K(7,5) = z[33];
	K(8,5) = z[4];
	K(9,5) = z[49];
	K(10,5) = z[11];
	K(11,5) = z[118];

	K(0,6) = z[16];
	K(1,6) = z[59];
	K(2,6) = z[7];
	K(3,6) = z[50];
	K(4,6) = z[25];
	K(5,6) = z[51];
	K(6,6) = z[17];
	K(7,6) = z[29];
	K(8,6) = z[54];
	K(9,6) = z[185];
	K(10,6) = z[47];
	K(11,6) = z[88];

	K(0,7) = z[24];
	K(1,7) = z[20];
	K(2,7) = z[13];
	K(3,7) = z[28];
	K(4,7) = z[62];
	K(5,7) = z[33];
	K(6,7) = z[29];
	K(7,7) = z[22];
	K(8,7) = z[6];
	K(9,7) = z[100];
	K(10,7) = z[137];
	K(11,7) = z[9];

	K(0,8) = z[57];
	K(1,8) = z[42];
	K(2,8) = z[8];
	K(3,8) = z[56];
	K(4,8) = z[32];
	K(5,8) = z[4];
	K(6,8) = z[54];
	K(7,8) = z[6];
	K(8,8) = z[44];
	K(9,8) = z[19];
	K(10,8) = z[69];
	K(11,8) = z[208];

	K(0,9) = z[105];
	K(1,9) = z[67];
	K(2,9) = z[96];
	K(3,9) = z[162];
	K(4,9) = z[68];
	K(5,9) = z[49];
	K(6,9) = z[185];
	K(7,9) = z[100];
	K(8,9) = z[19];
	K(9,9) = 0;
	K(10,9) = 0;
	K(11,9) = 0;

	K(0,10) = z[90];
	K(1,10) = z[231];
	K(2,10) = z[27];
	K(3,10) = z[63];
	K(4,10) = z[256];
	K(5,10) = z[11];
	K(6,10) = z[47];
	K(7,10) = z[137];
	K(8,10) = z[69];
	K(9,10) = 0;
	K(10,10) = 0;
	K(11,10) = 0;

	K(0,11) = z[53];
	K(1,11) = z[58];
	K(2,11) = z[277];
	K(3,11) = z[37];
	K(4,11) = z[79];
	K(5,11) = z[118];
	K(6,11) = z[88];
	K(7,11) = z[9];
	K(8,11) = z[208];
	K(9,11) = 0;
	K(10,11) = 0;
	K(11,11) = 0;
}

