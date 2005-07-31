/* $Id: LJFCC111.cpp,v 1.1.1.1 2001-01-29 08:20:23 paklein Exp $ */
/* created: paklein (07/31/1996)                                          */

#include "LJFCC111.h"

#include <math.h>
#include <iostream.h>

#include "fstreamT.h"

const double sqrt2 = sqrt(2.0);
const double sqrt3 = sqrt(3.0);

/* constructor */
LJFCC111::LJFCC111(ifstreamT& in, const ElasticT& element):
	NL_E_RotMat2DT(in, element, kPlaneStrain)
{
	in >> fScale;	if (fScale < 0.0) throw eBadInputValue;
}

/* I/O functions  */
void LJFCC111::Print(ostream& out) const
{
	/* inherited */
	NL_E_RotMat2DT::Print(out);

	out << " Scaling constant. . . . . . . . . . . . . . . . = " << fScale << '\n';
}

/*************************************************************************
* Protected
*************************************************************************/

void LJFCC111::PrintName(ostream& out) const
{
	NL_E_RotMat2DT::PrintName(out);

	out << "    LJ FCC <111> 2D\n";
}

/* compute the symetric Cij reduced index matrix */
void LJFCC111::ComputeModuli(const dSymMatrixT& E, dMatrixT& moduli)
{
	double z1, z2, z3, z4, z5, z6, z7, z8, z9, z10, z11, z12, z13, z14, z15;
	double z16, z17, z18, z19, z20, z21, z22, z23, z24, z25, z26, z27, z28;
	double z29, z30, z31, z32, z33, z34, z35, z36, z37, z38, z39, z40, z41;
	double z42, z43, z44, z45, z46;

	double E11 = E[0];
	double E22 = E[1];
	double E12 = E[2];
	
	double a = 1.0 + ThermalElongation();

	z1 = pow(a,6);
	z2 = pow(a,12);
	z3 = 2.*E11;
	z4 = 3.*E11;
	z5 = -3.464101615137754*E12;
	z6 = 3.464101615137754*E12;
	z7 = 2.*E22;
	z8 = 3.*E22;
	z3 = 1. + z3;
	z4 = 6. + E22 + z4;
	z9 = z4 + z5;
	z4 = z4 + z6;
	z7 = 3. + z7;
	z8 = 2. + E11 + z8;
	z5 = z5 + z8;
	z6 = z6 + z8;
	z8 = pow(z9,-8);
	z9 = pow(z9,-5);
	z10 = pow(z4,-8);
	z4 = pow(z4,-5);
	z11 = pow(z7,-8);
	z7 = pow(z7,-5);
	z12 = pow(z5,-8);
	z5 = pow(z5,-5);
	z13 = pow(z6,-8);
	z6 = pow(z6,-5);
	z14 = pow(z3,-8);
	z3 = pow(z3,-5);
	z15 = -23328.*z1*z9;
	z16 = -7776.000000000001*z1*z9;
	z17 = -2592.*z1*z9;
	z9 = z1*z9;
	z18 = -5.091065436109812e6*z2*z8;
	z19 = -1.697021812036604e6*z2*z8;
	z8 = z2*z8;
	z10 = z10*z2;
	z20 = -23328.*z1*z4;
	z21 = -13468.42707965559*z1*z4;
	z22 = -7776.000000000001*z1*z4;
	z23 = -4489.47569321853*z1*z4;
	z24 = -2592.*z1*z4;
	z4 = z1*z4;
	z11 = 61236.00000000001*z11*z2;
	z7 = -1296.*z1*z7;
	z25 = -6983.628856117712*z12*z2;
	z26 = -2327.876285372571*z12*z2;
	z12 = z12*z2;
	z27 = -864.*z1*z5;
	z28 = -288.*z1*z5;
	z29 = -96.*z1*z5;
	z5 = z1*z5;
	z13 = z13*z2;
	z30 = -864.*z1*z6;
	z31 = -498.8306325798366*z1*z6;
	z32 = -288.*z1*z6;
	z33 = -166.2768775266122*z1*z6;
	z34 = -96.*z1*z6;
	z6 = z1*z6;
	z2 = 84.*z14*z2;
	z1 = -48.*z1*z3;
	z3 = 4489.47569321853*z9;
	z9 = 13468.42707965559*z9;
	z14 = 979776.*z8;
	z35 = 2.939328e6*z8;
	z8 = 8.817984e6*z8;
	z36 = 979776.*z10;
	z37 = 1.697021812036604e6*z10;
	z38 = 2.939328e6*z10;
	z39 = 5.091065436109812e6*z10;
	z10 = 8.817984e6*z10;
	z40 = 1344.*z12;
	z41 = 4032.*z12;
	z12 = 12096.*z12;
	z42 = 166.2768775266122*z5;
	z5 = 498.8306325798366*z5;
	z43 = 1344.*z13;
	z44 = 2327.876285372571*z13;
	z45 = 4032.*z13;
	z46 = 6983.628856117712*z13;
	z13 = 12096.*z13;
	z7 = z11 + z7;
	z1 = z1 + z2;
	z2 = z19 + z3;
	z3 = z18 + z9;
	z9 = z14 + z17;
	z11 = z16 + z35;
	z8 = z15 + z8;
	z14 = z24 + z36;
	z15 = z23 + z37;
	z16 = z22 + z38;
	z17 = z21 + z39;
	z10 = z10 + z20;
	z18 = z29 + z40;
	z19 = z28 + z41;
	z12 = z12 + z27;
	z20 = z26 + z42;
	z5 = z25 + z5;
	z21 = z34 + z43;
	z22 = z33 + z44;
	z23 = z32 + z45;
	z24 = z31 + z46;
	z13 = z13 + z30;
	z25 = 3.*fScale;
	z23 = z23*z25;
	z24 = z24*z25;
	z13 = z13*z25;
	z7 = z25*z7;
	z1 = z1*z25;
	z2 = z2*z25;
	z3 = z25*z3;
	z9 = z25*z9;
	z11 = z11*z25;
	z8 = z25*z8;
	z14 = z14*z25;
	z15 = z15*z25;
	z16 = z16*z25;
	z17 = z17*z25;
	z10 = z10*z25;
	z18 = z18*z25;
	z19 = z19*z25;
	z12 = z12*z25;
	z20 = z20*z25;
	z5 = z25*z5;
	z21 = z21*z25;
	z22 = z22*z25;
	z11 = z11 + z16 + z19 + z23;
	z7 = z12 + z13 + z14 + z7 + z9;
	z2 = z15 + z2 + z24 + z5;
	z1 = z1 + z10 + z18 + z21 + z8;
	z3 = z17 + z20 + z22 + z3;
	z5 = 0.3535533905932737*z11;
	z7 = 0.3535533905932737*z7;
	z2 = 0.3535533905932737*z2;
	z1 = 0.3535533905932737*z1;
	z3 = 0.3535533905932737*z3;
			
	/* returns reduced index 3 x 3 */
	/* Also C_33 = C_12 for LJFCC111 */
	moduli(0,0) = z1;
	moduli(1,1) = z7;
	moduli(2,2) = z5;
	moduli(1,2) = z2;
	moduli(0,2) = z3;
	moduli(0,1) = z5;
	
	/* symmetric */
	moduli.CopySymmetric();
}

/* symetric 2nd Piola-Kirchhoff reduced index vector */
void LJFCC111::ComputePK2(const dSymMatrixT& E, dSymMatrixT& PK2)
{
	double z1, z2, z3, z4, z5, z6, z7, z8, z9, z10, z11, z12, z13, z14;
	double z15, z16, z17, z18, z19, z20, z21, z22, z23, z24, z25, z26, z27;
	double z28, z29, z30;

	double	E11 = E[0];
	double	E22 = E[1];
	double	E12 = E[2];
	
	double a = 1.0 + ThermalElongation();

	z1 = 1./sqrt2;
	z2 = sqrt3;
	z3 = pow(a,6);
	z4 = pow(a,12);
	z5 = 2*E11;
	z6 = 3*E11;
	z7 = 2*E22;
	z8 = 3*E22;
	z9 = -2*E12*z2;
	z10 = E12*z2;
	z10 = 2*z10;
	z5 = 1 + z5;
	z7 = 3 + z7;
	z11 = pow(z5,-7);
	z5 = pow(z5,-4);
	z12 = pow(z7,-7);
	z7 = pow(z7,-4);
	z13 = 6 + E22 + z10 + z6;
	z10 = 2 + E11 + z10 + z8;
	z6 = 6 + E22 + z6 + z9;
	z8 = 2 + E11 + z8 + z9;
	z9 = pow(z13,-7);
	z13 = pow(z13,-4);
	z14 = pow(z10,-7);
	z10 = pow(z10,-4);
	z15 = pow(z6,-7);
	z6 = pow(z6,-4);
	z16 = pow(z8,-7);
	z8 = pow(z8,-4);
	z5 = 6*z3*z5;
	z7 = 162*z3*z7;
	z11 = -6*z11*z4;
	z12 = -4374*z12*z4;
	z13 = z13*z3;
	z10 = z10*z3;
	z6 = z3*z6;
	z3 = z3*z8;
	z5 = z11 + z5;
	z7 = z12 + z7;
	z8 = 648*z13;
	z11 = 1944*z13;
	z12 = 24*z10;
	z10 = 72*z10;
	z13 = 648*z6;
	z17 = 1944*z6;
	z6 = -648*z2*z6;
	z18 = 24*z3;
	z19 = 72*z3;
	z3 = -24*z2*z3;
	z20 = -419904*z4*z9;
	z21 = -139968*z4*z9;
	z9 = z4*z9;
	z22 = -576*z14*z4;
	z23 = -192*z14*z4;
	z14 = z14*z4;
	z24 = -419904*z15*z4;
	z25 = -139968*z15*z4;
	z15 = z15*z4;
	z26 = -576*z16*z4;
	z27 = -192*z16*z4;
	z4 = z16*z4;
	z16 = z2*z8;
	z28 = z12*z2;
	z29 = z2*z21;
	z30 = z2*z23;
	z15 = 139968*z15*z2;
	z2 = 192*z2*z4;
	z4 = 3*fScale;
	z11 = z11 + z20;
	z8 = z21 + z8;
	z10 = z10 + z22;
	z12 = z12 + z23;
	z17 = z17 + z24;
	z13 = z13 + z25;
	z19 = z19 + z26;
	z18 = z18 + z27;
	z16 = z16 + z29;
	z20 = z28 + z30;
	z6 = z15 + z6;
	z2 = z2 + z3;
	z3 = z4*z5;
	z5 = z4*z7;
	z7 = z11*z4;
	z8 = z4*z8;
	z10 = z10*z4;
	z11 = z12*z4;
	z12 = z17*z4;
	z13 = z13*z4;
	z15 = z19*z4;
	z17 = z18*z4;
	z16 = z16*z4;
	z18 = z20*z4;
	z6 = z4*z6;
	z2 = z2*z4;
	z4 = z10 + z13 + z15 + z5 + z8;
	z3 = z11 + z12 + z17 + z3 + z7;
	z2 = z16 + z18 + z2 + z6;
	z1 = z1/2;
	z4 = z1*z4;
	z3 = z1*z3;
	z1 = z1*z2;

	PK2[0] = z3;
	PK2[1] = z4;
	PK2[2] = z1;
}

/* strain energy density */
double LJFCC111::ComputeEnergyDensity(const dSymMatrixT& E)
{
	double z1, z2, z3, z4, z5, z6, z7, z8, z9, z10, z11, z12, z13, z14, z15;	

	double E11 = E[0];
	double E22 = E[1];
	double E12 = E[2];

	double a = 1.0 + ThermalElongation();

	z1 = 1./sqrt2;
	z2 = sqrt3;
	z3 = pow(a,6);
	z4 = pow(a,12);
	z5 = 2*E11;
	z6 = 3*E11;
	z7 = 2*E22;
	z8 = 3*E22;
	z9 = -2*E12*z2;
	z2 = E12*z2;
	z2 = 2*z2;
	z5 = 1 + z5;
	z7 = 3 + z7;
	z10 = pow(z5,-6);
	z5 = pow(z5,-3);
	z11 = pow(z7,-6);
	z7 = pow(z7,-3);
	z12 = 6 + E22 + z2 + z6;
	z2 = 2 + E11 + z2 + z8;
	z6 = 6 + E22 + z6 + z9;
	z8 = 2 + E11 + z8 + z9;
	z9 = pow(z12,-6);
	z12 = pow(z12,-3);
	z13 = pow(z2,-6);
	z2 = pow(z2,-3);
	z14 = pow(z6,-6);
	z6 = pow(z6,-3);
	z15 = pow(z8,-6);
	z8 = pow(z8,-3);
	z5 = -(z3*z5);
	z7 = -27*z3*z7;
	z10 = z10*z4/2;
	z11 = 729*z11*z4/2;
	z12 = -216*z12*z3;
	z2 = -8*z2*z3;
	z6 = -216*z3*z6;
	z3 = -8*z3*z8;
	z5 = z10 + z5;
	z7 = z11 + z7;
	z8 = 23328*z4*z9;
	z9 = 32*z13*z4;
	z10 = 23328*z14*z4;
	z4 = 32*z15*z4;
	z11 = 3*fScale;
	z8 = z12 + z8;
	z2 = z2 + z9;
	z6 = z10 + z6;
	z3 = z3 + z4;
	z4 = z11*z5;
	z5 = z11*z7;
	z7 = z11*z8;
	z2 = z11*z2;
	z6 = z11*z6;
	z3 = z11*z3;
	z2 = z2 + z3 + z4 + z5 + z6 + z7;
	z1 = z1*z2/2;

	return z1;
}
