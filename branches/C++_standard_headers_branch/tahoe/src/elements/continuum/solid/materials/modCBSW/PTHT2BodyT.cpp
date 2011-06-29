/* $Id: PTHT2BodyT.cpp,v 1.5 2004-07-15 08:28:36 paklein Exp $ */
/* created: paklein (10/11/1997) */
#include "PTHT2BodyT.h"

#include <math.h>
#include "ThermalDilatationT.h"

using namespace Tahoe;

/* constructor */
PTHT2BodyT::PTHT2BodyT(const dArrayT& lengths,
	const ThermalDilatationT* thermal, double A, double A1, double A2):
	TwoBodyT(lengths, thermal),
	fA(A), fA1(A1), fA2(A2)
{

}

/* set free dof - triggers recomputation */
void PTHT2BodyT::Set(void)
{	
	double* pPhi   = fPhi.Pointer();
	double* pdPhi  = fdPhi.Pointer();
	double* pddPhi = fddPhi.Pointer();
	
	const double* pl = fLengths.Pointer();
	
	/* expansion factor */
	double a = 1.0;
	if (fThermal) a *= (1.0 + fThermal->PercentElongation());
	
	for (int i = 0; i < fLengths.Length(); i++)
	{
		*pPhi++   = U2body(*pl,a);
		*pdPhi++  = DU2body(*pl,a);
		*pddPhi++ = DDU2body(*pl,a);
	
		pl++;
	}
}

/**********************************************************************
* Private
**********************************************************************/

	/* 2 body potential and derivatives */
double PTHT2BodyT::U2body(double r, double a) const
{
	return( fA*(fA1/pow(r/a,12) - fA2/pow(r/a,6)) );
}

double PTHT2BodyT::DU2body(double r, double a) const
{
	return( 6*fA*(-2*fA1 + fA2*pow(r/a,6))/pow(r/a,13) );
}

double PTHT2BodyT::DDU2body(double r, double a) const
{
	return( 6*fA*(26*fA1 - 7*fA2*pow(r/a,6))/pow(r/a,14) );
}
