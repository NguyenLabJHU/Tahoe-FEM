/* $Id: HarmonicPairT.cpp,v 1.1 2002-11-25 07:19:46 paklein Exp $ */
#include "HarmonicPairT.h"

using namespace Tahoe;

/* initialize static parameters */
double HarmonicPairT::sR0 = 0.0;
double HarmonicPairT::sK = 0.0;

/* constructor */
HarmonicPairT::HarmonicPairT(double R0, double K):
	fR0(R0),
	fK(K)
{
	/* assume nearest neighbor - 10% of equilibrium spacing */
	SetRange(1.1*fR0);
}

/* return a pointer to the energy function */
PairPropertyT::EnergyFunction HarmonicPairT::getEnergyFunction(void)
{
	/* copy my data to static */
	sR0 = fR0;
	sK = fK;

	/* return function pointer */
	return HarmonicPairT::Energy;
}

PairPropertyT::ForceFunction HarmonicPairT::getForceFunction(void)
{
	/* copy my data to static */
	sR0 = fR0;
	sK = fK;

	/* return function pointer */
	return HarmonicPairT::Force;
}

PairPropertyT::StiffnessFunction HarmonicPairT::getStiffnessFunction(void)
{
	/* copy my data to static */
	sR0 = fR0;
	sK = fK;

	/* return function pointer */
	return HarmonicPairT::Stiffness;
}

/***********************************************************************
 * Private
 ***********************************************************************/

double HarmonicPairT::Energy(double r_ab, double* data_a, double* data_b)
{
#pragma unused(data_a)
#pragma unused(data_b)

	double dr = r_ab - sR0;
	return 0.5*sK*dr*dr;
}

double HarmonicPairT::Force(double r_ab, double* data_a, double* data_b)
{
#pragma unused(data_a)
#pragma unused(data_b)

	double dr = r_ab - sR0;
	return sK*dr;
}

double HarmonicPairT::Stiffness(double r_ab, double* data_a, double* data_b)
{
#pragma unused(data_a)
#pragma unused(data_b)
#pragma unused(r_ab)

	return sK;
}
