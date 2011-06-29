/* $Id: HarmonicPairT.cpp,v 1.5 2004-07-15 08:29:49 paklein Exp $ */
#include "HarmonicPairT.h"
#include <iostream.h>

using namespace Tahoe;

/* initialize static parameters */
double HarmonicPairT::sR0 = 0.0;
double HarmonicPairT::sK = 0.0;

/* constructor */
HarmonicPairT::HarmonicPairT(double mass, double R0, double K):
	fR0(R0),
	fK(K)
{
	SetName("harmonic");

	/* assume nearest neighbor - 10x of equilibrium spacing */
	SetRange(10.0*fR0);
	SetMass(mass);
	SetNearestNeighbor(fR0);
}

HarmonicPairT::HarmonicPairT(void):
	fR0(0.0),
	fK(0.0)
{
	SetName("harmonic");
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

/* describe the parameters needed by the interface */
void HarmonicPairT::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	PairPropertyT::DefineParameters(list);

	ParameterT rest_length(fR0, "rest_length");
	rest_length.AddLimit(0.0, LimitT::Lower);
	list.AddParameter(rest_length);

	ParameterT stiffness(fK, "stiffness");
	stiffness.AddLimit(0.0, LimitT::LowerInclusive);
	list.AddParameter(stiffness);
}

/* accept parameter list */
void HarmonicPairT::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	PairPropertyT::TakeParameterList(list);

	fR0 = list.GetParameter("rest_length");
	fK = list.GetParameter("stiffness");
	SetRange(10.0*fR0);
	SetNearestNeighbor(fR0);
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
