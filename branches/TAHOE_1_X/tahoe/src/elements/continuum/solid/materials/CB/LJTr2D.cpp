/* $Id: LJTr2D.cpp,v 1.10 2004-06-17 07:40:15 paklein Exp $ */
/* created: paklein (07/01/1996) */
#include "LJTr2D.h"

#include <math.h>
#include <iostream.h>

#include "ifstreamT.h"

using namespace Tahoe;

const double sqrt3 = sqrt(3.0);

/* constructor */
LJTr2D::LJTr2D(ifstreamT& in, const FSMatSupportT& support):
	NL_E_Mat2DT(in, support, kPlaneStress),
	CBLatticeT(2,2,3)
{
	in >> feps;	if (feps < 0.0) throw ExceptionT::kBadInputValue;
}

void LJTr2D::Initialize(void) 
{
  NL_E_Mat2DT::Initialize();

  CBLatticeT::Initialize();
}

/* I/O functions */
void LJTr2D::Print(ostream& out) const
{
	/* inherited */
	NL_E_Mat2DT::Print(out);

	out << " Lennard-Jones energy scaling constant . . . . . = " << feps << '\n';

}

/*************************************************************************
* Protected
*************************************************************************/

void LJTr2D::PrintName(ostream& out) const
{
	NL_E_Mat2DT::PrintName(out);

	out << "    LJ triangular 2D\n";
}


void LJTr2D::ComputeModuli(const dSymMatrixT& E, dMatrixT& moduli)
{	
	ComputeDeformedLengths(E);
	
	dMatrixT fBondTensor4(3);
	fBondTensor4 = 0;
	moduli = 0.; 
	for (int i = 0; i < fNumBonds; i++)
	{
		double ri = fDefLength[i];
		double coeff = 2./sqrt3*(ddU(ri)-dUlj(ri)/ri)/ri/ri;
		BondComponentTensor4(i,fBondTensor4);
		moduli.AddScaled(coeff,fBondTensor4);
	}
	
	/* symmetric */
	moduli.CopySymmetric();
}

/* 2nd Piola-Kirchhoff stress vector */
void LJTr2D::ComputePK2(const dSymMatrixT& E, dSymMatrixT& PK2)
{
	
	ComputeDeformedLengths(E);
	
	dArrayT fBondTensor2(3);
	PK2 = 0.; 
	for (int i = 0; i < fNumBonds; i++)
	{
		double ri = fDefLength[i];
		double coeff = dUlj(ri)*2./sqrt3/ri;
		BondComponentTensor2(i,fBondTensor2);
		PK2.AddScaled(coeff,fBondTensor2);
	}
}

/* strain energy density */
double LJTr2D::ComputeEnergyDensity(const dSymMatrixT& E)
{
	ComputeDeformedLengths(E);
	
	double tmpSum  = 0.;
	
	for (int i = 0; i < fNumBonds; i++)
	{
		double r = fDefLength[i];
		tmpSum += Ulj(r);
	}
	tmpSum *= 2./sqrt3;
	
	return tmpSum;
}

void LJTr2D::LoadBondTable(void)
{
  	fBondCounts = 1;

  	fDefLength = 0.0; 
  
  	double bonddata[3][2] = {
  		{1.0, 0.0},
  		{0.5,-sqrt3/2.0},
  		{0.5,sqrt3/2.0}};
 
  	if (fBonds.MajorDim() != fNumBonds ||
     	fBonds.MinorDim() != 2) throw ExceptionT::kGeneralFail;

	for (int i = 0; i < fNumBonds; i++)
    	for (int j = 0; j < 2; j++)
      		fBonds(i,j) = bonddata[i][j];

  	/*  Multiply by lattice constant */
  	//  fBonds *= fa;
}

/*************************************************************************
* Private
*************************************************************************/

/* second derivative of the Lennard-Jones 6/12 potential */
double LJTr2D::ddU(double l) const
{
	double a = pow(1.0 + ThermalElongation(),6);
	return feps*a*(78.0*a/pow(l,14) - 42.0/pow(l, 8));
}

double LJTr2D::dUlj(double r) const
{
	double a = pow(1.0 + ThermalElongation(),6);
	return feps*a*(-6.*a/pow(r,13)+6./pow(r,7));
}

double LJTr2D::Ulj(double r) const
{
	double a = pow(1.0 + ThermalElongation(),6);
	return feps*a*(a*pow(r,-12)/2.-pow(r,-6));
}


