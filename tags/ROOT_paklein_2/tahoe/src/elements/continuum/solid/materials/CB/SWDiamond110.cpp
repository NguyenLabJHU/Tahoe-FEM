/* $Id: SWDiamond110.cpp,v 1.3 2002-07-02 19:55:34 cjkimme Exp $ */
/* created: paklein (08/25/1996)                                          */
/* (11/9/1996) : includes cut off terms                                   */

#include "SWDiamond110.h"
#include <math.h>
#include <iostream.h>

/* constructor */

using namespace Tahoe;

SWDiamond110::SWDiamond110(ifstreamT& in, const FiniteStrainT& element):
	SWMaterial2D(in, element)
{

}

/*************************************************************************
* Protected
*************************************************************************/

/* print name */
void SWDiamond110::PrintName(ostream& out) const
{
	/* inherited */
	SWMaterial2D::PrintName(out);

	out << "    SW Diamond <110> 2D\n";
}

/* symmetric Cij reduced index matrix */
void SWDiamond110::ComputeModuli(const dSymMatrixT& E, dMatrixT& moduli)
{
	double z1, z2, z3, z4, z5, z6, z7, z8, z9, z10, z11, z12;
	double z13, z14, z15, z16, z17, z18, z19, z20, z21, z22, z23, z24;
	double z25, z26, z27, z28, z29, z30, z31, z32, z33, z34, z35, z36;
	double z37, z38, z39, z40, z41, z42, z43, z44, z45, z46, z47, z48;
	double z49, z50, z51, z52, z53, z54, z55, z56, z57, z58, z59, z60;
	double z61, z62, z63, z64, z65, z66, z67, z68, z69, z70, z71, z72;
	double z73, z74, z75, z76, z77, z78, z79, z80, z81, z82, z83, z84;
	double z85, z86, z87, z88, z89, z90, z91, z92, z93, z94, z95, z96;
	double z97, z98, z99, z100, z101, z102, z103, z104, z105, z106;
	double z107, z108;//, z109;	
	
	double E11 = E[0];
	double E22 = E[1];
	double E12 = E[2];

	double fa = 1.0 + ThermalElongation();

	z1 = pow(2.,-0.6666666666666666);
	z2 = pow(2.,0.3333333333333333);
	z3 = pow(2.,0.6666666666666666);
	z4 = pow(3.,-0.5);
	z5 = pow(fa,-2.);
	z6 = pow(fa,-1.);
	z7 = pow(fa,4.);
	z8 = pow(fdelta,2.);
	z9 = -4.*E11;
	z10 = 4.*E11;
	z11 = -4.*E12;
	z12 = 4.*E12;
	z13 = -4.*E22;
	z14 = 4.*E22;
	z15 = pow(fgamma,2.);
	z16 = -1.*frcut;
	z10 = 3. + z10;
	z11 = -1. + z11;
	z17 = -1. + z12;
	z12 = 1. + z12;
	z13 = -1. + z13;
	z14 = 3. + z14;
	z9 = -1. + z9;
	z18 = pow(z10,-4.);
	z19 = pow(z10,-3.);
	z20 = pow(z10,-2.);
	z21 = pow(z10,-1.);
	z22 = pow(z11,2.);
	z23 = pow(z17,2.);
	z10 = z10*z2;
	z24 = pow(z12,2.);
	z25 = pow(z14,-4.);
	z26 = pow(z14,-3.);
	z27 = pow(z14,-2.);
	z28 = pow(z14,-1.);
	z14 = z14*z2;
	z29 = 64.00000000000001*z19*z9;
	z30 = 64.00000000000001*z20;
	z31 = -8.*z20*z9;
	z32 = -8.*z21;
	z9 = z21*z9;
	z33 = pow(z10,-2.5);
	z34 = pow(z10,-1.5);
	z35 = pow(z10,-0.5);
	z10 = pow(z10,0.5);
	z36 = 64.00000000000001*z13*z26;
	z37 = 64.00000000000001*z27;
	z38 = -8.*z13*z27;
	z39 = -8.*z28;
	z13 = z13*z28;
	z40 = pow(z14,-2.5);
	z41 = pow(z14,-1.5);
	z42 = pow(z14,-0.5);
	z14 = pow(z14,0.5);
	z1 = 9.*fB*z1*z7;
	z43 = z1*z20;
	z1 = z1*z27;
	z29 = z29 + z30;
	z30 = z31 + z32;
	z9 = 0.3333333333333333 + z9;
	z10 = z10*z4;
	z31 = z36 + z37;
	z36 = z38 + z39;
	z13 = 0.3333333333333333 + z13;
	z37 = z2*z35*z42;
	z14 = z14*z4;
	z38 = -1. + z43;
	z1 = -1. + z1;
	z43 = z4*z6;
	z44 = z10*z6;
	z45 = z11*z37;
	z46 = z17*z37;
	z47 = -1.*z12*z37;
	z6 = z14*z6;
	z48 = pow(z30,2.);
	z49 = pow(z9,2.);
	z50 = pow(z36,2.);
	z51 = pow(z13,2.);
	z44 = z16 + z44;
	z45 = 0.3333333333333333 + z45;
	z46 = 0.3333333333333333 + z46;
	z47 = 0.3333333333333333 + z47;
	z6 = z16 + z6;
	z16 = pow(z44,-4.);
	z52 = pow(z44,-3.);
	z53 = pow(z44,-2.);
	z44 = pow(z44,-1.);
	z54 = pow(z45,2.);
	z55 = pow(z46,2.);
	z56 = pow(z47,2.);
	z57 = pow(z6,-4.);
	z58 = pow(z6,-3.);
	z59 = pow(z6,-2.);
	z6 = pow(z6,-1.);
	z60 = 2.666666666666666*z2*z5;
	z61 = 4.*z3*z43;
	z62 = fdelta*z44;
	z44 = fgamma*z44;
	z63 = fdelta*z6;
	z6 = fgamma*z6;
	z52 = z21*z52*z60;
	z58 = z28*z58*z60;
	z60 = z34*z53*z61;
	z61 = z41*z59*z61;
	z6 = exp(z6);
	z62 = exp(z62);
	z44 = exp(z44);
	z63 = exp(z63);
	z52 = z52 + z60;
	z58 = z58 + z61;
	z60 = (z10 < frcut) ? z62 : 0;
	z61 = (z10 < frcut) ? z44 : 0;
	z64 = (z14 < frcut) ? z6 : 0;
	z65 = (z14 < frcut) ? z63 : 0;
	z66 = fdelta*z52*z62;
	z52 = fgamma*z44*z52;
	z67 = fgamma*z58*z6;
	z58 = fdelta*z58*z63;
	z68 = pow(z61,2.);
	z69 = -16.*flambda*z61*z64;
	z70 = flambda*z61*z64;
	z71 = pow(z64,2.);
	z72 = 8.*z70;
	z73 = 16.*z70;
	z74 = z11*z70;
	z75 = z11*z72;
	z76 = z17*z69;
	z77 = z17*z70;
	z78 = z12*z70;
	z79 = z20*z70;
	z79 = z20*z72;
	z80 = z21*z70;
	z81 = z21*z72;
	z82 = z22*z72;
	z23 = z23*z73;
	z83 = z3*z70;
	z83 = z3*z72;
	z84 = z24*z72;
	z81 = z26*z81;
	z85 = z22*z81;
	z81 = z24*z81;
	z86 = z21*z23*z26;
	z79 = z27*z79;
	z22 = z22*z79;
	z24 = z24*z79;
	z79 = z20*z23*z27;
	z87 = z21*z27*z75;
	z88 = z21*z27*z76;
	z89 = z28*z70;
	z72 = z28*z72;
	z72 = z19*z28*z82;
	z23 = z19*z23*z28;
	z82 = z19*z28*z84;
	z75 = z20*z28*z75;
	z76 = z20*z28*z76;
	z80 = 5461.333333333334*z28*z4*z80;
	z27 = z27*z32*z78;
	z32 = z34*z69;
	z84 = z34*z70;
	z84 = z35*z70;
	z20 = z20*z39*z78;
	z5 = 1.333333333333333*z2*z5;
	z15 = z15*z5;
	z39 = z41*z69;
	z69 = z41*z70;
	z69 = z35*z41*z83;
	z70 = z42*z70;
	z70 = z34*z42*z83;
	z83 = 864.*fA*fB*z2*z7;
	z18 = z18*z60*z83;
	z25 = z25*z65*z83;
	z60 = z5*z8;
	z60 = z45*z74;
	z60 = 48.*z60;
	z65 = z35*z40*z60;
	z60 = z33*z42*z60;
	z74 = z45*z69;
	z83 = z45*z70;
	z73 = z11*z34*z41*z45*z73;
	z77 = z46*z77;
	z84 = 96.*z77;
	z89 = z35*z40*z84;
	z77 = 32.*z34*z41*z77;
	z84 = z33*z42*z84;
	z39 = z3*z35*z39*z46;
	z90 = z3*z32*z42*z46;
	z78 = -48.*z47*z78;
	z40 = z35*z40*z78;
	z33 = z33*z42*z78;
	z69 = z47*z69;
	z70 = z47*z70;
	z32 = z12*z32*z41*z47;
	z78 = (z10 < frcut) ? -2.*fdelta*z2*z35*z43*z53*z62 : 0;
	z53 = (z10 < frcut) ? -2.*fgamma*z2*z35*z43*z44*z53 : 0;
	z91 = (z14 < frcut) ? -2.*fgamma*z2*z42*z43*z59*z6 : 0;
	z43 = (z14 < frcut) ? -2.*fdelta*z2*z42*z43*z59*z63 : 0;
	z59 = -4.*flambda;
	z92 = 0.5*flambda;
	z93 = 4.*flambda;
	z94 = 8.*flambda;
	z95 = flambda*z61;
	z96 = flambda*z64;
	z97 = flambda*z68;
	z98 = flambda*z71;
	z99 = flambda*z53;
	z100 = pow(z53,2.);
	z101 = pow(z91,2.);
	z102 = z53*z59*z64;
	z103 = z59*z61*z91;
	z104 = z61*z93;
	z93 = z64*z93;
	z64 = z53*z64*z94;
	z61 = z61*z91*z94;
	z94 = z91*z95;
	z105 = z53*z96;
	z99 = z91*z99;
	z106 = flambda*z3;
	z107 = -8.*z3*z94;
	z94 = z3*z94;
	z108 = -8.*z105*z3;
	z105 = z105*z3;
	// z109 = z106*z35*z41; not used again
	z94 = z35*z41*z94;
	z106 = z106*z34*z42;
	z105 = z105*z34*z42;
	z29 = z29*z9*z97;
	z9 = z104*z30*z53*z9;
	z30 = z13*z31*z98;
	z13 = z13*z36*z91*z93;
	z31 = z37*z59;
	z31 = z102*z37;
	z36 = z103*z37;
	z2 = -144.*fA*fB*z2*z7;
	z7 = z19*z2*z78;
	z2 = z2*z26*z43;
	z19 = z48*z68*z92;
	z26 = z50*z71*z92;
	z5 = z5*z8;
	z8 = z45*z59;
	z8 = z102*z45;
	z8 = z11*z3*z35*z41*z8;
	z43 = z103*z45;
	z43 = z11*z3*z34*z42*z43;
	z48 = z31*z45;
	z50 = z36*z45;
	z59 = z107*z11*z35*z41*z45;
	z11 = z108*z11*z34*z42*z45;
	z45 = -16.*z17*z46*z94;
	z68 = -16.*z105*z17*z46;
	z71 = z108*z17*z35*z41*z46;
	z17 = z107*z17*z34*z42*z46;
	z78 = z37*z46*z64;
	z37 = z37*z46*z61;
	z31 = z31*z47;
	z36 = z36*z47;
	z46 = z12*z3*z35*z41*z47*z53*z93;
	z35 = z12*z3*z35*z41*z47*z61;
	z41 = z104*z12*z3*z34*z42*z47*z91;
	z3 = z12*z3*z34*z42*z47*z64;
	z12 = z16*z21;
	z16 = z12*z15*z44;
	z12 = z12*z5*z62;
	z21 = z54*z99;
	z34 = 2.*z55*z99;
	z42 = z56*z99;
	z28 = z28*z57;
	z6 = z15*z28*z6;
	z5 = z28*z5*z63;
	z15 = 2.*flambda;
	z20 = z20 + z31 + z48 + z70 + z75 + z76 + z78 + z83 + z90;
	z27 = z27 + z36 + z37 + z39 + z50 + z69 + z74 + z87 + z88;
	z8 = z17 + z21 + z22 + z24 + z32 + z34 + z41 + z42 + z43 + z46 + z71 + z73 + z77 + z79 + z8;
	z17 = z100*z15*z49;
	z15 = z101*z15*z51;
	z16 = (z10 < frcut) ? z16 + z52 : 0;
	z10 = (z10 < frcut) ? z12 + z66 : 0;
	z6 = (z14 < frcut) ? z6 + z67 : 0;
	z5 = (z14 < frcut) ? z5 + z58 : 0;
	z12 = 2.*fA;
	z14 = 2.*z95;
	z21 = 2.*z96;
	z22 = z16*z96;
	z24 = z6*z95;
	z4 = 170.6666666666666*z4;
	z20 = z20*z4;
	z27 = z27*z4;
	z8 = z4*z8;
	z10 = z10*z12*z38;
	z1 = z1*z12*z5;
	z5 = z14*z16*z49;
	z12 = z21*z51*z6;
	z28 = z22*z54;
	z31 = z24*z54;
	z6 = z14*z55*z6;
	z14 = z16*z21*z55;
	z16 = z22*z56;
	z21 = z24*z56;
	z3 = z10 + z11 + z14 + z16 + z17 + z18 + z19 + z23 + z28 + z29 +
	     z3 + z33 + z5 + z60 + z68 + z7 + z72 + z82 + z84 + z9;
	z1 = z1 + z12 + z13 + z15 + z2 + z21 + z25 + z26 + z30 + z31 +
	     z35 + z40 + z45 + z59 + z6 + z65 + z81 + z85 + z86 + z89;
	z2 = z3*z4;
	z1 = z1*z4;

	//z1 = List(z2,z1,z80,z27,z20,z8);

	moduli(0,0) = z2;
	moduli(1,1) = z1;
	moduli(2,2) = z80;
	moduli(1,2) = z27;
	moduli(0,2) = z20;
	moduli(0,1) = z8;

	/* symmetric */
	moduli.CopySymmetric();

}

/* symmetric 2nd Piola-Kirchhoff reduced index vector */
void SWDiamond110::ComputePK2(const dSymMatrixT& E, dSymMatrixT& PK2)
{
	double z1, z2, z3, z4, z5, z6, z7, z8, z9, z10, z11, z12;
	double z13, z14, z15, z16, z17, z18, z19, z20, z21, z22, z23, z24;
	double z25, z26, z27, z28, z29, z30, z31, z32, z33, z34, z35, z36;
	double z37, z38, z39, z40, z41;

	double	E11 = E[0];
	double	E22 = E[1];
	double	E12 = E[2];

	double fa = 1.0 + ThermalElongation();

	z1 = pow(2.,-0.6666666666666666);
	z2 = pow(2.,0.3333333333333333);
	z3 = pow(2.,0.6666666666666666);
	z4 = pow(3.,-0.5);
	z5 = pow(fa,-1.);
	z6 = pow(fa,4.);
	z7 = -4.*E11;
	z8 = 4.*E11;
	z9 = -4.*E12;
	z10 = 4.*E12;
	z11 = -4.*E22;
	z12 = 4.*E22;
	z13 = -1.*frcut;
	z14 = -1. + z10;
	z10 = 1. + z10;
	z11 = -1. + z11;
	z12 = 3. + z12;
	z7 = -1. + z7;
	z8 = 3. + z8;
	z9 = -1. + z9;
	z15 = pow(z12,-3.);
	z16 = pow(z12,-2.);
	z17 = pow(z12,-1.);
	z18 = pow(z8,-3.);
	z19 = pow(z8,-2.);
	z20 = pow(z8,-1.);
	z12 = z12*z2;
	z8 = z2*z8;
	z21 = -8.*z11*z16;
	z22 = -8.*z17;
	z11 = z11*z17;
	z17 = -8.*z19*z7;
	z23 = -8.*z20;
	z7 = z20*z7;
	z20 = pow(z12,-1.5);
	z24 = pow(z12,-0.5);
	z12 = pow(z12,0.5);
	z25 = pow(z8,-1.5);
	z26 = pow(z8,-0.5);
	z8 = pow(z8,0.5);
	z1 = 9.*fB*z1*z6;
	z16 = z1*z16;
	z1 = z1*z19;
	z19 = z21 + z22;
	z11 = 0.3333333333333333 + z11;
	z17 = z17 + z23;
	z7 = 0.3333333333333333 + z7;
	z21 = z2*z24*z26;
	z12 = z12*z4;
	z8 = z4*z8;
	z16 = -1. + z16;
	z1 = -1. + z1;
	z22 = z14*z21;
	z23 = -1.*z10*z21;
	z27 = z21*z9;
	z28 = z4*z5;
	z29 = z5*z8;
	z30 = pow(z11,2.);
	z31 = pow(z7,2.);
	z5 = z12*z5;
	z22 = 0.3333333333333333 + z22;
	z23 = 0.3333333333333333 + z23;
	z27 = 0.3333333333333333 + z27;
	z29 = z13 + z29;
	z5 = z13 + z5;
	z13 = pow(z22,2.);
	z32 = pow(z23,2.);
	z33 = pow(z27,2.);
	z34 = pow(z29,-2.);
	z29 = pow(z29,-1.);
	z35 = pow(z5,-2.);
	z5 = pow(z5,-1.);
	z36 = fdelta*z29;
	z29 = fgamma*z29;
	z36 = exp(z36);
	z29 = exp(z29);
	z37 = fdelta*z5;
	z5 = fgamma*z5;
	z38 = (z8 < frcut) ? z36 : 0;
	z39 = (z8 < frcut) ? z29 : 0;
	z37 = exp(z37);
	z5 = exp(z5);
	z18 = -72.00000000000001*fA*fB*z18*z2*z38*z6;
	z38 = pow(z39,2.);
	z40 = (z12 < frcut) ? z37 : 0;
	z41 = (z12 < frcut) ? z5 : 0;
	z36 = (z8 < frcut) ? -2.*fdelta*z2*z26*z28*z34*z36 : 0;
	z8 = (z8 < frcut) ? -2.*fgamma*z2*z26*z28*z29*z34 : 0;
	z29 = flambda*z39;
	z7 = flambda*z17*z38*z7;
	z6 = -72.00000000000001*fA*fB*z15*z2*z40*z6;
	z15 = flambda*z41;
	z17 = pow(z41,2.);
	z1 = 2.*fA*z1*z36;
	z34 = flambda*z8;
	z34 = -4.*z29*z41;
	z36 = z29*z41;
	z31 = 2.*z29*z31*z8;
	z8 = z15*z8;
	z38 = 2.*z13*z8;
	z39 = z32*z8;
	z8 = z33*z8;
	z40 = z21*z34;
	z41 = z23*z40;
	z40 = z27*z40;
	z27 = z27*z34;
	z9 = z27*z3*z9;
	z27 = z24*z25*z9;
	z9 = z20*z26*z9;
	z3 = z3*z36;
	z25 = z24*z25*z3;
	z3 = z20*z26*z3;
	z20 = -8.*z14*z22*z25;
	z14 = -8.*z14*z22*z3;
	z25 = 4.*z10*z23*z25;
	z3 = 4.*z10*z23*z3;
	z10 = z21*z36;
	z10 = 8.*z10*z22;
	z21 = (z12 < frcut) ? -2.*fdelta*z2*z24*z28*z35*z37 : 0;
	z2 = (z12 < frcut) ? -2.*fgamma*z2*z24*z28*z35*z5 : 0;
	z5 = z10 + z40 + z41;
	z10 = 2.*z2;
	z12 = 2.*fA*z16*z21;
	z2 = z2*z29;
	z16 = z2*z32;
	z2 = z2*z33;
	z13 = z10*z13*z29;
	z10 = z10*z15*z30;
	z11 = flambda*z11*z17*z19;
	z1 = z1 + z18 + z20 + z25 + z27 + z31 + z38 + z39 + z7 + z8;
	z4 = 170.6666666666666*z4;
	z5 = z4*z5;
	z1 = z1*z4;
	z2 = z10 + z11 + z12 + z13 + z14 + z16 + z2 + z3 + z6 + z9;

	//{z1, z2 z4, z5}

	PK2[0] = z1;
	PK2[1] = z2*z4;
	PK2[2] = z5;
}

/* strain energy density */
double SWDiamond110::ComputeEnergyDensity(const dSymMatrixT& E)
{
	double z1, z2, z3, z4, z5, z6, z7, z8, z9, z10, z11, z12, z13;
	double z14, z15, z16, z17;

	double E11 = E[0];
	double E22 = E[1];
	double E12 = E[2];

	double fa = 1.0 + ThermalElongation();

	z1 = pow(2.,-0.6666666666666666);
	z2 = pow(2.,0.3333333333333333);
	z3 = pow(3.,-0.5);
	z4 = pow(fa,-1.);
	z5 = pow(fa,4.);
	z6 = -4.*E11;
	z7 = 4.*E11;
	z8 = -4.*E12;
	z9 = 4.*E12;
	z10 = -4.*E22;
	z11 = 4.*E22;
	z12 = -1.*frcut;
	z10 = -1. + z10;
	z11 = 3. + z11;
	z6 = -1. + z6;
	z7 = 3. + z7;
	z8 = -1. + z8;
	z13 = -1. + z9;
	z9 = 1. + z9;
	z14 = pow(z11,-2.);
	z15 = pow(z11,-1.);
	z16 = pow(z7,-2.);
	z17 = pow(z7,-1.);
	z11 = z11*z2;
	z7 = z2*z7;
	z10 = z10*z15;
	z6 = z17*z6;
	z15 = pow(z11,-0.5);
	z11 = pow(z11,0.5);
	z17 = pow(z7,-0.5);
	z7 = pow(z7,0.5);
	z14 = 9.*fB*z1*z14*z5;
	z1 = 9.*fB*z1*z16*z5;
	z5 = 0.3333333333333333 + z10;
	z6 = 0.3333333333333333 + z6;
	z10 = z11*z3;
	z8 = z15*z17*z2*z8;
	z11 = z13*z15*z17*z2;
	z2 = -1.*z15*z17*z2*z9;
	z7 = z3*z7;
	z9 = -1. + z14;
	z1 = -1. + z1;
	z13 = z4*z7;
	z5 = pow(z5,2.);
	z6 = pow(z6,2.);
	z8 = 0.3333333333333333 + z8;
	z11 = 0.3333333333333333 + z11;
	z2 = 0.3333333333333333 + z2;
	z4 = z10*z4;
	z13 = z12 + z13;
	z8 = pow(z8,2.);
	z11 = pow(z11,2.);
	z2 = pow(z2,2.);
	z4 = z12 + z4;
	z12 = pow(z13,-1.);
	z4 = pow(z4,-1.);
	z13 = fdelta*z12;
	z12 = fgamma*z12;
	z14 = fdelta*z4;
	z4 = fgamma*z4;
	z13 = (z7 < frcut) ? exp(z13) : 0;
	z7 = (z7 < frcut) ? exp(z12) : 0;
	z1 = 2.*fA*z1*z13;
	z12 = pow(z7,2.);
	z13 = (z10 < frcut) ? exp(z14) : 0;
	z4 = (z10 < frcut) ? exp(z4) : 0;
	z6 = flambda*z12*z6;
	z9 = 2.*fA*z13*z9;
	z8 = flambda*z4*z7*z8;
	z10 = 2.*flambda*z11*z4*z7;
	z2 = flambda*z2*z4*z7;
	z4 = pow(z4,2.);
	z4 = flambda*z4*z5;
	z1 = z1 + z10 + z2 + z4 + z6 + z8 + z9;
	z1 = 512*z1*z3/3;
	
	// z1

	return(z1);
}
