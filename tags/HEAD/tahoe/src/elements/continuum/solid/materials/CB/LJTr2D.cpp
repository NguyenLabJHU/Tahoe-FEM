/* $Id: LJTr2D.cpp,v 1.1.1.1 2001-01-29 08:20:23 paklein Exp $ */
/* created: paklein (07/01/1996)                                          */

#include "LJTr2D.h"

#include <math.h>
#include <iostream.h>

#include "fstreamT.h"

const double sqrt3 = sqrt(3.0);

/* constructor */
LJTr2D::LJTr2D(ifstreamT& in, const ElasticT& element):
	NL_E_RotMat2DT(in, element, kPlaneStress),
	fBondVectors(3,2)
{
	in >> fScale;	if (fScale < 0.0) throw eBadInputValue;

	/* undeformed bond vectors (natural axes) */
	double bonds[3][2] = {{1.0      , 0.0},
	                      {sqrt3/2.0, 0.5},
	                      {sqrt3/2.0,-0.5}};

	dArrayT shbond_nat;
	dArrayT shbond;
	for (int i = 0; i < fBondVectors.MajorDim(); i++)
	{
		shbond_nat.Set(2, bonds[i]);
		fBondVectors.RowAlias(i,shbond);
	
		shbond = TransformOut(shbond_nat);
	}
}

/* I/O functions */
void LJTr2D::Print(ostream& out) const
{
	/* inherited */
	NL_E_RotMat2DT::Print(out);

	out << " Lennard-Jones scaling constant. . . . . . . . . = " << fScale << '\n';
}

/*************************************************************************
* Protected
*************************************************************************/

void LJTr2D::PrintName(ostream& out) const
{
	NL_E_RotMat2DT::PrintName(out);

	out << "    LJ triangular 2D\n";
}


void LJTr2D::ComputeModuli(const dSymMatrixT& E, dMatrixT& moduli)
{

	double	z1, z2, z3, z4, z5, z6, z7, z8, z9, z10, z11, z12, z13, z14, z15;
	double  z16, z17, z18, z19, z20, z21, z22, z23, z24, z25;
	
	double E11 = E[0];
	double E22 = E[1];
	double E12 = E[2];
	
	double a = 1.0 + ThermalElongation();
	
	z1 = 1.0/sqrt3;
	z2 = sqrt3;
	z3 = pow(a,6);
	z4 = pow(a,12);
	z5 = 2*E11;
	z6 = 2*E22;
	z7 = -(E12*z2);
	z8 = E12*z2;
	z5 = 1 + z5;
	z6 = 1 + z6;
	z9 = -(z2*z6)/2;
	z10 = z2*z6/2;
	z6 = z2*z6;
	z11 = pow(z5,-8);
	z12 = pow(z5,-5);
	z5 = z5/2;
	z9 = E12 + z9;
	z10 = E12 + z10;
	z12 = -48*z12*z3;
	z11 = 84*z11*z4;
	z7 = z5 + z7;
	z5 = z5 + z8;
	z8 = -(z2*z9)/2;
	z9 = z10*z2/2;
	z10 = z11 + z12;
	z7 = z7/2;
	z5 = z5/2;
	z10 = 2*fScale*z10;
	z7 = z7 + z8;
	z5 = z5 + z9;
	z8 = pow(z7,-8);
	z7 = pow(z7,-5);
	z9 = pow(z5,-8);
	z5 = pow(z5,-5);
	z11 = -27*z3*z7;
	z12 = -9*z3*z7;
	z13 = -3*z3*z7;
	z7 = z3*z7;
	z14 = -27*z3*z5;
	z15 = -9*z3*z5;
	z16 = -3*z3*z5;
	z3 = z3*z5;
	z5 = z2*z7;
	z7 = z4*z8;
	z4 = z4*z9;
	z8 = z15*z2;
	z9 = z16*z2;
	z17 = 3*z5;
	z5 = 9*z5;
	z18 = 21*z7/4;
	z19 = 63*z7/4;
	z20 = 189*z7/4;
	z21 = -63*z2*z7/4;
	z22 = -21*z2*z7/4;
	z7 = z2*z7;
	z23 = 21*z4/4;
	z24 = 63*z4/4;
	z4 = 189*z4/4;
	z25 = z2*z23;
	z2 = z2*z24;
	z13 = z13 + z18;
	z12 = z12 + z19;
	z11 = z11 + z20;
	z5 = z21 + z5;
	z17 = z17 + z22;
	z16 = z16 + z23;
	z15 = z15 + z24;
	z4 = z14 + z4;
	z9 = z25 + z9;
	z2 = z2 + z8;
	z8 = 2*fScale;
	z13 = z13*z8;
	z12 = z12*z8;
	z11 = z11*z8;
	z5 = z5*z8;
	z14 = z17*z8;
	z16 = z16*z8;
	z15 = z15*z8;
	z4 = z4*z8;
	z9 = z8*z9;
	z2 = z2*z8;
	z8 = z10 + z13 + z16;
	z10 = z12 + z15;
	z4 = z11 + z4;
	z9 = z14 + z9;
	z2 = z2 + z5;
	z5 = z1*z8;
	z8 = z1*z10;
	z4 = z1*z4;
	z9 = z1*z9;
	z1 = z1*z2;
	
	/* returns reduced index 3 x 3 */
	/* Also C_33 = C_12 for LJTr2Dt */

	moduli(0,0) = z5;
	moduli(1,1) = z4;
	moduli(2,2) = z8;
	moduli(1,2) = z1;
	moduli(0,2) = z9;
	moduli(0,1) = z8;
	
	/* symmetric */
	moduli.CopySymmetric();
}

/* 2nd Piola-Kirchhoff stress vector */
void LJTr2D::ComputePK2(const dSymMatrixT& E, dSymMatrixT& PK2)
{
	double	z1, z2, z3, z4, z5, z6, z7, z8, z9, z10, z11, z12, z13, z14;//, z15;
	double	z16, z17, z18, z19, z20;
	
	double E11 = E[0];
	double E22 = E[1];
	double E12 = E[2];
	
	double a = 1.0 + ThermalElongation();

	z1 = 1/sqrt3;
	z2 = sqrt3;
	z3 = pow(a,6);
	z4 = pow(a,12);
	z5 = 2*E11;
	z6 = 2*E22;
	z7 = -(E12*z2);
	z8 = E12*z2;
	z5 = 1 + z5;
	z6 = 1 + z6;
	z9 = -(z2*z6)/2;
	z10 = z2*z6/2;
	z6 = z2*z6;
	z11 = pow(z5,-7);
	z12 = pow(z5,-4);
	z5 = z5/2;
	z9 = E12 + z9;
	z10 = E12 + z10;
	z12 = 6*z12*z3;
	z11 = -6*z11*z4;
	z7 = z5 + z7;
	z5 = z5 + z8;
	z8 = -(z2*z9)/2;
	z9 = z10*z2/2;
	z10 = z11 + z12;
	z7 = z7/2;
	z5 = z5/2;
	z10 = 2*fScale*z10;
	z7 = z7 + z8;
	z5 = z5 + z9;
	z8 = pow(z7,-7);
	z7 = pow(z7,-4);
	z9 = pow(z5,-7);
	z5 = pow(z5,-4);
	z11 = -3*z2/2;
	z12 = 3*z2/2;
	z13 = 3*z3/2;
	z14 = 9*z3/2;
	// z15 = z2*z3; not used again
	z16 = z3*z7;
	z3 = z3*z5;
	z17 = z13*z7;
	z13 = z13*z5;
	z7 = z14*z7;
	z5 = z14*z5;
	z14 = -9*z4/2;
	z18 = -3*z4/2;
	z2 = z2*z4;
	z19 = z4*z8;
	z4 = z4*z9;
	z16 = z11*z16;
	z3 = z12*z3;
	z20 = z14*z8;
	z14 = z14*z9;
	z8 = z18*z8;
	z9 = z18*z9;
	z12 = z12*z19;
	z4 = z11*z4;
	z7 = z20 + z7;
	z5 = z14 + z5;
	z8 = z17 + z8;
	z9 = z13 + z9;
	z11 = z12 + z16;
	z3 = z3 + z4;
	z4 = 2*fScale;
	z7 = z4*z7;
	z5 = z4*z5;
	z8 = z4*z8;
	z9 = z4*z9;
	z11 = z11*z4;
	z3 = z3*z4;
	z4 = z5 + z7;
	z5 = z10 + z8 + z9;
	z3 = z11 + z3;
	z4 = z1*z4;
	z5 = z1*z5;
	z1 = z1*z3;

	PK2[0] = z5;
	PK2[1] = z4;
	PK2[2] = z1;
}

/* strain energy density */
double LJTr2D::ComputeEnergyDensity(const dSymMatrixT& E)
{
	double z1, z2, z3, z4, z5, z6, z7, z8, z9, z10, z11, z12;
	
	double E11 = E[0];
	double E22 = E[1];
	double E12 = E[2];
	
	double a = 1.0 + ThermalElongation();

	z1 = 1.0/sqrt3;
	z2 = sqrt3;
	z3 = pow(a,6);
	z4 = pow(a,12);
	z5 = 2*E11;
	z6 = 2*E22;
	z7 = -(E12*z2);
	z8 = E12*z2;
	z5 = 1 + z5;
	z6 = 1 + z6;
	z9 = -(z2*z6)/2;
	z10 = z2*z6/2;
	z6 = z2*z6;
	z11 = pow(z5,-6);
	z12 = pow(z5,-3);
	z5 = z5/2;
	z9 = E12 + z9;
	z10 = E12 + z10;
	z12 = -(z12*z3);
	z11 = z11*z4/2;
	z7 = z5 + z7;
	z5 = z5 + z8;
	z8 = -(z2*z9)/2;
	z2 = z10*z2/2;
	z9 = z11 + z12;
	z7 = z7/2;
	z5 = z5/2;
	z9 = 2*fScale*z9;
	z7 = z7 + z8;
	z2 = z2 + z5;
	z5 = pow(z7,-6);
	z7 = pow(z7,-3);
	z8 = pow(z2,-6);
	z2 = pow(z2,-3);
	z3 = -z3;
	z7 = z3*z7;
	z2 = z2*z3;
	z3 = z4/2;
	z4 = z3*z5;
	z3 = z3*z8;
	z4 = z4 + z7;
	z2 = z2 + z3;
	z3 = 2*fScale;
	z4 = z3*z4;
	z2 = z2*z3;
	z2 = z2 + z4 + z9;
	z1 = z1*z2;

	return z1;
}

/*************************************************************************
* Private
*************************************************************************/

/* second derivative of the Lennard-Jones 6/12 potential */
double LJTr2D::ddU(double l) const
{
	double a = 1.0 + ThermalElongation();
	return fScale*((78.0*pow(a,12))/pow(l,14) -
	               (42.0*pow(a, 6))/pow(l, 8));
}
