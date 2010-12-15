#include "FSDE_inc.h"

#include <math.h>

static double z[62];

/* function to compute mechanical part of tangent modulus */
void mech_tanmod_ab(const double* params, const double* Xsi, const double* Cmat, double J, double I1, double* ddCmech) { 

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
	z[19] = I1*I1;
	z[20] = pow(I1,3.);
	z[21] = pow(Nrig,-4.);
	z[22] = pow(Nrig,-3.);
	z[23] = 1./(Nrig*Nrig);
	z[24] = 1./Nrig;
	z[25] = C33*z[1];
	z[26] = C32*z[2];
	z[27] = C33*z[3];
	z[28] = C31*z[4];
	z[29] = C32*z[5];
	z[30] = C31*z[6];
	z[31] = log(J);
	z[12] = z[12] + z[16];
	z[10] = z[10] + z[17];
	z[14] = z[14] + z[18];
	z[16] = 0.3119777365491651*z[21];
	z[17] = 0.015406307977736549*z[20]*z[21];
	z[18] = 0.29314285714285715*z[22];
	z[19] = 0.03257142857142857*z[19]*z[22];
	z[20] = 0.28285714285714286*z[23];
	z[21] = 0.06285714285714286*I1*z[23];
	z[22] = 0.1*z[24];
	z[23] = 0.3*z[24];
	z[1] = z[1] + z[3];
	z[3] = z[25] + z[26] + z[27] + z[28] + z[29] + z[30];
	z[24] = -lambda*z[31];
	z[2] = z[2] + z[5];
	z[4] = z[4] + z[6];
	z[5] = z[11] + z[7];
	z[6] = z[15] + z[8];
	z[7] = z[13] + z[9];
	z[8] = z[12]*z[12];
	z[9] = z[10]*z[10];
	z[11] = z[14]*z[14];
	z[13] = z[17] + z[19] + z[21] + z[22];
	z[15] = 0.5 + z[16] + z[18] + z[20] + z[23];
	z[16] = z[1]*z[1];
	z[3] = 1./(z[3]*z[3]);
	z[17] = z[2]*z[2];
	z[18] = z[4]*z[4];
	z[19] = z[5]*z[5];
	z[20] = z[6]*z[6];
	z[21] = z[7]*z[7];
	z[13] = 4.*mu*z[13];
	z[15] = 2.*mu*z[15];
	z[22] = z[10]*z[12]*z[3];
	z[23] = lambda*z[10]*z[14]*z[3];
	z[25] = z[1]*z[12]*z[3];
	z[26] = z[1]*z[10]*z[3];
	z[27] = z[1]*z[14]*z[3];
	z[28] = z[12]*z[2]*z[3];
	z[29] = z[14]*z[2]*z[3];
	z[30] = z[10]*z[3]*z[4];
	z[31] = z[10]*z[3]*z[5];
	z[32] = z[14]*z[3]*z[5];
	z[33] = lambda*z[1]*z[3]*z[5];
	z[34] = z[2]*z[3]*z[5];
	z[35] = z[3]*z[4]*z[5];
	z[36] = lambda*z[10]*z[3]*z[6];
	z[37] = z[14]*z[3]*z[6];
	z[38] = z[1]*z[3]*z[6];
	z[39] = z[3]*z[4]*z[6];
	z[40] = lambda*z[3]*z[5]*z[6];
	z[41] = z[12]*z[3]*z[7];
	z[42] = lambda*z[10]*z[3]*z[7];
	z[43] = lambda*z[14]*z[3]*z[7];
	z[44] = lambda*z[1]*z[3]*z[7];
	z[45] = z[2]*z[3]*z[7];
	z[46] = z[3]*z[4]*z[7];
	z[47] = lambda*z[3]*z[5]*z[7];
	z[48] = z[3]*z[6]*z[7];
	z[49] = lambda*z[3]*z[9];
	z[50] = lambda*z[11]*z[3];
	z[51] = lambda*z[16]*z[3];
	z[52] = lambda*z[19]*z[3];
	z[53] = lambda*z[20]*z[3];
	z[54] = lambda*z[21]*z[3];
	z[55] = lambda*z[26];
	z[56] = lambda*z[27];
	z[57] = lambda*z[31];
	z[58] = lambda*z[32];
	z[59] = lambda*z[37];
	z[60] = lambda*z[38];
	z[61] = lambda*z[48];
	z[15] = z[15] + z[24];
	z[24] = z[29] + z[30];
	z[25] = z[25] + z[35];
	z[22] = z[22] + z[37];
	z[29] = z[34] + z[38];
	z[28] = z[28] + z[39];
	z[30] = z[32] + z[41];
	z[26] = z[26] + z[45];
	z[27] = z[27] + z[46];
	z[31] = z[31] + z[48];
	z[24] = z[15]*z[24];
	z[25] = z[15]*z[25];
	z[22] = z[15]*z[22];
	z[29] = z[15]*z[29];
	z[28] = z[15]*z[28];
	z[30] = z[15]*z[30];
	z[26] = z[15]*z[26];
	z[27] = z[15]*z[27];
	z[31] = z[15]*z[31];
	z[32] = 2.*z[12]*z[14]*z[15]*z[3];
	z[34] = 2.*z[10]*z[14]*z[15]*z[3];
	z[35] = 2.*z[10]*z[15]*z[2]*z[3];
	z[37] = 2.*z[1]*z[15]*z[2]*z[3];
	z[38] = 2.*z[12]*z[15]*z[3]*z[4];
	z[39] = 2.*z[14]*z[15]*z[3]*z[4];
	z[41] = 2.*z[1]*z[15]*z[3]*z[4];
	z[4] = 2.*z[15]*z[2]*z[3]*z[4];
	z[45] = 2.*z[12]*z[15]*z[3]*z[5];
	z[46] = 2.*z[1]*z[15]*z[3]*z[5];
	z[12] = 2.*z[12]*z[15]*z[3]*z[6];
	z[48] = 2.*z[10]*z[15]*z[3]*z[6];
	z[2] = 2.*z[15]*z[2]*z[3]*z[6];
	z[6] = 2.*z[15]*z[3]*z[5]*z[6];
	z[10] = 2.*z[10]*z[15]*z[3]*z[7];
	z[14] = 2.*z[14]*z[15]*z[3]*z[7];
	z[1] = 2.*z[1]*z[15]*z[3]*z[7];
	z[5] = 2.*z[15]*z[3]*z[5]*z[7];
	z[7] = 2.*z[15]*z[3]*z[8];
	z[8] = 2.*z[15]*z[3]*z[9];
	z[9] = 2.*z[11]*z[15]*z[3];
	z[11] = 2.*z[15]*z[16]*z[3];
	z[16] = 2.*z[15]*z[17]*z[3];
	z[17] = 2.*z[15]*z[18]*z[3];
	z[18] = 2.*z[15]*z[19]*z[3];
	z[19] = 2.*z[15]*z[20]*z[3];
	z[3] = 2.*z[15]*z[21]*z[3];
	z[15] = z[35] + z[61];
	z[20] = z[38] + z[58];
	z[21] = z[45] + z[61];
	z[10] = z[10] + z[58];
	z[7] = z[13] + z[59] + z[7];
	z[8] = z[13] + z[59] + z[8];
	z[16] = z[13] + z[16] + z[60];
	z[18] = z[13] + z[18] + z[60];
	z[32] = z[23] + z[32];
	z[23] = z[23] + z[34];
	z[34] = z[33] + z[37];
	z[33] = z[33] + z[46];
	z[12] = z[12] + z[36];
	z[35] = z[36] + z[48];
	z[2] = z[2] + z[40];
	z[6] = z[40] + z[6];
	z[24] = z[24] + z[42];
	z[30] = z[30] + z[42];
	z[36] = z[39] + z[43];
	z[14] = z[14] + z[43];
	z[37] = z[41] + z[44];
	z[1] = z[1] + z[44];
	z[25] = z[25] + z[47];
	z[26] = z[26] + z[47];
	z[22] = z[22] + z[49];
	z[9] = z[13] + z[50] + z[9];
	z[11] = z[11] + z[13] + z[51];
	z[29] = z[29] + z[52];
	z[19] = z[13] + z[19] + z[53];
	z[27] = z[27] + z[54];
	z[4] = z[4] + z[55];
	z[5] = z[5] + z[55];
	z[17] = z[13] + z[17] + z[56];
	z[3] = z[13] + z[3] + z[56];
	z[13] = z[28] + z[57];
	z[28] = z[31] + z[57];

	/* Output stress */
	/* dCdE:  6 x 6 */
	ddCmech[0] = z[9];
	ddCmech[1] = z[8];
	ddCmech[2] = z[3];
	ddCmech[3] = z[10];
	ddCmech[4] = z[14];
	ddCmech[5] = z[23];
	
	ddCmech[6] = z[7];
	ddCmech[7] = z[19];
	ddCmech[8] = z[18];
	ddCmech[9] = z[6];
	ddCmech[10] = z[21];
	ddCmech[11] = z[12];
	
	ddCmech[12] = z[17];
	ddCmech[13] = z[16];
	ddCmech[14] = z[11];
	ddCmech[15] = z[34];
	ddCmech[16] = z[37];
	ddCmech[17] = z[4];	
	
	ddCmech[18] = z[20];
	ddCmech[19] = z[2];
	ddCmech[20] = z[33];
	ddCmech[21] = z[29];
	ddCmech[22] = z[25];
	ddCmech[23] = z[13];	
	
	ddCmech[24] = z[36];
	ddCmech[25] = z[15];
	ddCmech[26] = z[1];
	ddCmech[27] = z[26];
	ddCmech[28] = z[27];
	ddCmech[29] = z[24];	
	
	ddCmech[30] = z[32];
	ddCmech[31] = z[35];
	ddCmech[32] = z[5];
	ddCmech[33] = z[28];
	ddCmech[34] = z[30];
	ddCmech[35] = z[22];		
}