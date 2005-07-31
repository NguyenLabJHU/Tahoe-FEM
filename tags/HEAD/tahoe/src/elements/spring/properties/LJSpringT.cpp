/* $Id: LJSpringT.cpp,v 1.1.1.1 2001-01-29 08:20:25 paklein Exp $ */
/* created: paklein (5/28/1996)                                           */

#include "LJSpringT.h"

#include <iostream.h>
#include <math.h>

#include "Environment.h"
#include "ExceptionCodes.h"

#include "fstreamT.h"
#include "ThermalDilatationT.h"


/*
* constructor
*/
LJSpringT::LJSpringT(ifstreamT& in): RodMaterialT(in)
{
	in >> fLJConstant;	if (fLJConstant < 0.0) throw eBadInputValue;
}

/*
* Returns true if the material has internal forces in the unloaded
* configuration, ie thermal strains.
*/
int LJSpringT::HasInternalStrain(void) const
{
	return 1;
	
	/*
	 * Always potentially has internal strain since the unstressed
	 * configurations is determined by the LJ potential and not the
	 * initial geometry.
	 */
}

/*
* I/O functions.
*/
void LJSpringT::Print(ostream& out) const
{
	/* inherited */
	RodMaterialT::Print(out);

	out << " Lennard-Jones scaling constant . . . . . . . . .= " << fLJConstant << '\n';	
}
	
/* potential function and derivatives */
double LJSpringT::Potential(double rmag, double Rmag) const
{
#pragma unused(Rmag)

	double a = 1.0 + fThermal->PercentElongation();
	
	double r = a/rmag;
	
	return fLJConstant*(-pow(r,6) + 0.5*pow(r,12));
}

double LJSpringT::DPotential(double rmag, double Rmag) const
{
#pragma unused(Rmag)

	double a = 1.0 + fThermal->PercentElongation();

	double r = a/rmag;

	return (6.0*fLJConstant/a)*(pow(r,7) - pow(r,13));
}

double LJSpringT::DDPotential(double rmag, double Rmag) const
{
#pragma unused(Rmag)

	double a = 1.0 + fThermal->PercentElongation();

	double r = a/rmag;

	return (fLJConstant/(a*a))*(-42.0*pow(r,8) + 78.0*pow(r,14));
}

/*************************************************************************
* Protected
*************************************************************************/

void LJSpringT::PrintName(ostream& out) const
{
	out << "    Lennard-Jones 6/12 spring\n";
}
