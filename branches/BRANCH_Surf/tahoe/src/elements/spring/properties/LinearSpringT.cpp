/* $Id: LinearSpringT.cpp,v 1.6 2004-07-15 08:30:22 paklein Exp $ */
/* created: paklein (05/28/1996) */
#include "LinearSpringT.h"

#include <iostream.h>

#include "ExceptionT.h"
#include "ifstreamT.h"
#include "ThermalDilatationT.h"

using namespace Tahoe;

/* constructor */
LinearSpringT::LinearSpringT(ifstreamT& in): RodMaterialT(in)
{
	in >> fSpringConstant;	if (fSpringConstant < 0.0) throw ExceptionT::kBadInputValue;
}

/* potential function and derivatives */
double LinearSpringT::Potential(double rmag, double Rmag) const
{
	/* imposed strain */
	if (fThermal->IsActive()) Rmag *= (1.0 + fThermal->PercentElongation());
	double delta = rmag - Rmag;	
	return 0.5*fSpringConstant*delta*delta;
}

double LinearSpringT::DPotential(double rmag, double Rmag) const
{
	/* imposed strain */
	if (fThermal->IsActive()) Rmag *= (1.0 + fThermal->PercentElongation());
	return fSpringConstant*(rmag - Rmag);
}

double LinearSpringT::DDPotential(double rmag, double Rmag) const
{
#pragma unused(rmag)
#pragma unused(Rmag)

	return fSpringConstant;
}
