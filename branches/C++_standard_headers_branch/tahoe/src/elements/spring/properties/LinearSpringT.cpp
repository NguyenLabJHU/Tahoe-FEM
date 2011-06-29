/* $Id: LinearSpringT.cpp,v 1.7 2005-11-06 00:37:58 paklein Exp $ */
/* created: paklein (05/28/1996) */
#include "LinearSpringT.h"

#include <iostream.h>
#include "ExceptionT.h"
//#include "ThermalDilatationT.h"

using namespace Tahoe;

/* constructor */
LinearSpringT::LinearSpringT(double mass, double k):
	RodMaterialT(mass),
	fSpringConstant(k)
{
	if (fSpringConstant < 0.0) ExceptionT::BadInputValue("LinearSpringT::LinearSpringT");
}

/* potential function and derivatives */
double LinearSpringT::Potential(double rmag, double Rmag) const
{
	/* imposed strain */
//	if (fThermal->IsActive()) Rmag *= (1.0 + fThermal->PercentElongation());
	double delta = rmag - Rmag;	
	return 0.5*fSpringConstant*delta*delta;
}

double LinearSpringT::DPotential(double rmag, double Rmag) const
{
	/* imposed strain */
//	if (fThermal->IsActive()) Rmag *= (1.0 + fThermal->PercentElongation());
	return fSpringConstant*(rmag - Rmag);
}

double LinearSpringT::DDPotential(double rmag, double Rmag) const
{
#pragma unused(rmag)
#pragma unused(Rmag)

	return fSpringConstant;
}
