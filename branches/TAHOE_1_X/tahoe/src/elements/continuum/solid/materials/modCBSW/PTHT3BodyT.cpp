/* $Id: PTHT3BodyT.cpp,v 1.3 2004-06-17 07:41:03 paklein Exp $ */
/* created: paklein (10/11/1997)                                          */

#include "PTHT3BodyT.h"

#include <math.h>

#include "ifstreamT.h"
#include "dMatrixT.h"
#include "iArray2DT.h"
#include "ThermalDilatationT.h"


/* parameters */

using namespace Tahoe;

const int kNumVars = 3; //number of arguments in Phi

/* constructor */
PTHT3BodyT::PTHT3BodyT(const dArrayT& lengths,
	const dArrayT& angles, const iArray2DT& bondpairs,
	const ThermalDilatationT* thermal, ifstreamT& in):
	ThreeBodyT(lengths, angles, bondpairs, thermal)
{
	/* read potential patameters */
	in >> fB;
	in >> fZ;
}

/* triggers recomputation */
void PTHT3BodyT::Set(void)
{
	/* shallow workspace */
	dMatrixT sh_mat;
	dArrayT  sh_vec;

	/* expansion factor */
	double a = 1.0;
	if (fThermal) a *= (1.0 + fThermal->PercentElongation());
		
	for (int i = 0; i < fPairs.MajorDim(); i++)
	{
		double r1  = fLengths[fPairs(i,0)];
		double r2  = fLengths[fPairs(i,1)];
		double c12 = fAngles[i];
	
		/* phi */
		fPhi[i] = ComputePhi(r1, r2, c12, a);
		
		/* gradient */		
		sh_vec.Set(kNumVars, fdPhi(i));
		ComputeGradient(sh_vec, r1, r2, c12, a);

		/* hessian */		
		sh_mat.Set(kNumVars, kNumVars, fddPhi(i));
		ComputeHessian(sh_mat, r1, r2, c12, a);
	}
}

/**********************************************************************
* Private
**********************************************************************/

/* return the 3-body potential */	
double PTHT3BodyT::ComputePhi(double r1, double r2, double c12, double a) const
{
	r1 /= a;
	r2 /= a;

	return( fB*fZ*(1 + 3*c12*(c12*r1 - r2)*(-r1 + c12*r2)/(pow(r1,2) -
	        2*c12*r1*r2 + pow(r2,2)))/
(pow(r1,3)*pow(r2,3)*pow(pow(r1,2) -
2*c12*r1*r2 + pow(r2,2),3/2)) );
}

/* compute Gradient of the 3-body potential */	
void PTHT3BodyT::ComputeGradient(dArrayT& grad, double r1, double r2,
	double c12, double a)
{
	double z1, z2, z3, z4, z5, z6, z7, z8, z9, z10, z11, z12;
	double z13, z14, z15, z16, z17, z18, z19, z20, z21, z22, z23, z24;
	
	r1 /= a;
	r2 /= a;
	
	z1 = pow(c12,2.);
	z2 = pow(r1,-4.);
	z3 = pow(r1,-3.);
	z4 = pow(r1,-2.);
	z5 = -1.*r1;
	z6 = 2.*r1;
	z7 = -2.*c12*r1;
	z8 = c12*r1;
	z9 = pow(r1,2.);
	z10 = pow(r2,-4.);
	z11 = pow(r2,-3.);
	z12 = pow(r2,-2.);
	z13 = -1.*r2;
	z14 = 2.*r2;
	z15 = -2.*c12*r2;
	z16 = c12*r2;
	z17 = pow(r2,2.);
	z18 = r2*z7;
	z5 = z16 + z5;
	z6 = z15 + z6;
	z7 = z14 + z7;
	z13 = z13 + z8;
	z9 = z17 + z18 + z9;
	z14 = pow(z9,-2.5);
	z15 = pow(z9,-2.);
	z17 = pow(z9,-1.5);
	z9 = pow(z9,-1.);
	z18 = -3.*c12*z5;
	z19 = 3.*z1*z5*z9;
	z20 = -3.*c12*z13*z9;
	z1 = 3.*z1*z13*z9;
	z16 = 3.*z13*z16*z9;
	z21 = 3.*z13*z5*z9;
	z22 = z13*z15*z18;
	z18 = z18*z9;
	z23 = c12*z21;
	z24 = z22*z6;
	z22 = z22*z7;
	z13 = 6.*r2*z13*z15*z5*z8;
	z5 = 3.*z5*z8*z9;
	z8 = 1. + z23;
	z9 = z19 + z20 + z24;
	z1 = z1 + z18 + z22;
	z5 = z13 + z16 + z21 + z5;
	z13 = fB*fZ;
	z15 = z11*z13*z17;
	z11 = z11*z13*z3;
	z16 = z13*z17*z3;
	z13 = z13*z8;
	z2 = -3.*z15*z2*z8;
	z9 = z15*z3*z9;
	z1 = z1*z15*z3;
	z3 = z15*z3*z5;
	z5 = -1.5*z11*z14*z8;
	z8 = -3.*z10*z16*z8;
	z10 = z13*z14;
	z6 = z5*z6;
	z5 = z5*z7;
	z4 = 3.*z10*z12*z4;
	z2 = z2 + z6 + z9;
	z1 = z1 + z5 + z8;
	z3 = z3 + z4;

	//returned	
	//List(z2,z1,z3);
	
	grad[0] = z2;
	grad[1] = z1;
	grad[2] = z3;
}

/* compute Hessian of the 3-body potential */
void PTHT3BodyT::ComputeHessian(dMatrixT& hessian, double r1, double r2,
	double c12, double a)
{

	double z1, z2, z3, z4, z5, z6, z7, z8, z9, z10, z11, z12;
	double z13, z14, z15, z16, z17, z18, z19, z20, z21, z22, z23, z24;
	double z25, z26, z27, z28, z29, z30, z31, z32, z33, z34, z35, z36;
	double z37, z38, z39, z40, z41, z42, z43, z44, z45, z46, z47, z48;
	double z49, z50, z51, z52, z53, z54, z55, z56, z57, z58, z59, z60;
	double z61, z62, z63, z64, z65, z66, z67, z68, z69, z70, z71, z72, z73;
	
	r1 /= a;
	r2 /= a;
	
	z1 = -2.*c12;
	z2 = pow(c12,2.);
	z3 = pow(c12,3.);
	z4 = pow(r1,-5.);
	z5 = pow(r1,-4.);
	z6 = pow(r1,-3.);
	z7 = pow(r1,-2.);
	z8 = pow(r1,-1.);
	z9 = -1.*r1;
	z10 = 2.*r1;
	z11 = c12*r1;
	z12 = pow(r1,2.);
	z13 = pow(r2,-5.);
	z14 = pow(r2,-4.);
	z15 = pow(r2,-3.);
	z16 = pow(r2,-2.);
	z17 = pow(r2,-1.);
	z18 = -1.*r2;
	z19 = 2.*r2;
	z20 = c12*r2;
	z21 = pow(r2,2.);
	z22 = r1*z1;
	z1 = r2*z1;
	z23 = r2*z22;
	z18 = z11 + z18;
	z19 = z19 + z22;
	z1 = z1 + z10;
	z10 = z12 + z21 + z23;
	z9 = z20 + z9;
	z22 = pow(z19,2.);
	z23 = pow(z1,2.);
	z24 = pow(z10,-3.5);
	z25 = pow(z10,-3.);
	z26 = pow(z10,-2.5);
	z27 = pow(z10,-2.);
	z28 = pow(z10,-1.5);
	z10 = pow(z10,-1.);
	z25 = z18*z25*z9;
	z29 = -6.*r2*z11*z18*z27;
	z30 = 6.*r1*r2*z18*z2*z27;
	z31 = 12.*z11*z18*z21*z27;
	z32 = 3.*c12*z18*z19*z27;
	z33 = -6.*z18*z19*z2*z27;
	z34 = -3.*z18*z19*z20*z27;
	z35 = 6.*c12*z1*z18*z27;
	z36 = -3.*z1*z18*z2*z27;
	z37 = -3.*z1*z18*z20*z27;
	z38 = -6.*r2*z11*z27*z9;
	z39 = 6.*r1*r2*z2*z27*z9;
	z40 = 12.*z12*z20*z27*z9;
	z41 = -6.*c12*z18*z27*z9;
	z42 = 12.*r1*r2*z18*z27*z9;
	z43 = 6.*z11*z18*z27*z9;
	z44 = 6.*z18*z2*z27*z9;
	z45 = 6.*z18*z20*z27*z9;
	z46 = 6.*c12*z19*z27*z9;
	z47 = -3.*z11*z19*z27*z9;
	z48 = -3.*z19*z2*z27*z9;
	z49 = -3.*z18*z19*z27*z9;
	z50 = 3.*c12*z1*z27*z9;
	z51 = -3.*z1*z11*z27*z9;
	z52 = -6.*z1*z2*z27*z9;
	z27 = -3.*z1*z18*z27*z9;
	z53 = 3.*c12*z10;
	z54 = -3.*z10*z11;
	z55 = 6.*r2*z10*z11;
	z56 = -6.*z10*z2;
	z57 = 3.*r1*z10*z2;
	z58 = 3.*r2*z10*z2;
	z59 = -3.*z10*z20;
	z60 = -3.*z10*z18;
	z61 = 6.*c12*z10*z18;
	z62 = 6.*r2*z10*z18;
	z63 = 3.*z10*z18*z2;
	z20 = 3.*z10*z18*z20;
	z64 = -3.*z10*z9;
	z65 = 6.*c12*z10*z9;
	z66 = 6.*r1*z10*z9;
	z67 = 3.*z10*z11*z9;
	z2 = 3.*z10*z2*z9;
	z68 = 3.*z10*z18*z9;
	z3 = 3.*z10*z3;
	z10 = c12*z25;
	z69 = -12.*r2*z11*z19*z25;
	z11 = -12.*r2*z1*z11*z25;
	z25 = r2*z43;
	z70 = c12*z49;
	z71 = c12*z27;
	z9 = z18*z53*z9;
	z18 = c12*z60;
	z72 = c12*z64;
	z12 = 24.*z10*z12*z21;
	z21 = 6.*z1*z10*z19;
	z73 = 6.*z10*z22;
	z10 = 6.*z10*z23;
	z30 = z30 + z34 + z38 + z43 + z47 + z49 + z57 + z59 + z61 + z64;
	z27 = z27 + z29 + z37 + z39 + z45 + z51 + z54 + z58 + z60 + z65;
	z20 = z20 + z25 + z67 + z68;
	z9 = 1. + z9;
	z2 = z18 + z2 + z71;
	z18 = z63 + z70 + z72;
	z12 = z12 + z31 + z40 + z42 + z55 + z62 + z66;
	z3 = z21 + z3 + z32 + z36 + z44 + z48 + z50 + z53;
	z21 = z33 + z41 + z46 + z56 + z73;
	z10 = z10 + z35 + z41 + z52 + z56;
	z25 = z30 + z69;
	z11 = z11 + z27;
	z27 = fB*fZ;
	z2 = z2*z27;
	z29 = z15*z27;
	z15 = z15*z2;
	z30 = -6.*z27*z28;
	z31 = z2*z28;
	z32 = z28*z29;
	z33 = z32*z5;
	z34 = -3.*z18*z33;
	z35 = -6.*z15*z28*z5;
	z12 = z12*z32*z6;
	z3 = z3*z32*z6;
	z21 = z21*z32*z6;
	z10 = z10*z32*z6;
	z25 = z25*z32*z6;
	z11 = z11*z32*z6;
	z30 = z14*z18*z30*z6;
	z31 = -3.*z14*z31*z6;
	z36 = -3.*z18*z19*z26*z29*z6;
	z37 = -1.5*z15*z19*z26*z6;
	z38 = -1.5*z1*z18*z26*z29*z6;
	z15 = -3.*z1*z15*z26*z6;
	z18 = 3.*z16*z18*z26*z27*z7;
	z2 = 3.*z16*z2*z26*z7;
	z33 = -3.*z20*z33;
	z39 = -1.5*z19*z20*z26*z29*z6;
	z40 = -1.5*z1*z20*z26*z29*z6;
	z41 = -3.*z14*z20*z27*z28*z6;
	z20 = 6.*z16*z20*z26*z27*z7;
	z4 = 12.*z32*z4*z9;
	z32 = 4.5*z19*z26*z29*z5*z9;
	z42 = 9.*z1*z26*z29*z5*z9;
	z5 = 9.*z14*z27*z28*z5*z9;
	z43 = 3.75*z1*z19*z24*z29*z6*z9;
	z22 = 3.75*z22*z24*z29*z6*z9;
	z23 = 3.75*z23*z24*z29*z6*z9;
	z44 = -3.*z26*z29*z6*z9;
	z45 = 3.*c12*z26*z29*z6*z9;
	z46 = -6.*z16*z26*z27*z6*z9;
	z47 = 9.*z14*z19*z26*z27*z6*z9;
	z14 = 4.5*z1*z14*z26*z27*z6*z9;
	z6 = 12.*z13*z27*z28*z6*z9;
	z13 = -7.5*z16*z19*z24*z27*z7*z9;
	z1 = -7.5*z1*z16*z24*z27*z7*z9;
	z7 = -6.*z26*z29*z7*z9;
	z8 = 15.*z17*z24*z27*z8*z9;
	z4 = z10 + z15 + z23 + z35 + z4 + z42 + z44;
	z3 = z14 + z3 + z31 + z32 + z34 + z37 + z38 + z43 + z45 + z5;
	z5 = z21 + z22 + z30 + z36 + z44 + z47 + z6;
	z1 = z1 + z11 + z2 + z33 + z40 + z46;
	z2 = z13 + z18 + z25 + z39 + z41 + z7;
	z6 = z12 + z20 + z8;

	//return
	//{{z4, z3, z1},
	// {z3, z5, z2},
	// {z1, z2, z6}}

	hessian(0,0) = z4;
	hessian(1,1) = z5;
	hessian(2,2) = z6;
	hessian(2,1) = hessian(1,2) = z2;
	hessian(2,0) = hessian(0,2) = z1;
	hessian(1,0) = hessian(0,1) = z3;
}
