/* $Id: DPPrimitiveT.cpp,v 1.1.1.1 2001-01-29 08:20:30 paklein Exp $ */
/* created: myip (06/01/1999)                                             */
/* Base class for a Druker Prager pressure dependent plastic material     */
/* with linear isotropic hardening.                                       */

#include "DPPrimitiveT.h"

#include <iostream.h>
#include <math.h>

#include "fstreamT.h"
#include "dSymMatrixT.h"

const double sqrt23 = sqrt(2.0/3.0);
const double sqrt32 = sqrt(3.0/2.0);

/* constructor */
DPPrimitiveT::DPPrimitiveT(ifstreamT& in):
falpha_bar(-1.0),  //**mien**//
	ffriction(-1.0),
	fdilation(-1.0),
	fH_prime(0.0),
	fK_prime(0.0),
	fH_delta(1.0),
	fK_delta(1.0)
{
	/* read parameters */
in >> falpha_bar;  if (falpha_bar < 0.0 ) throw eBadInputValue;  //**mien**//
	in >> ffriction;   if (ffriction < 0.0 ) throw eBadInputValue;
	in >> fdilation;   if (fdilation < 0.0 ) throw eBadInputValue;
	in >> fH_prime;
	in >> fK_prime;
	in >> fH_delta;    if (fH_delta > 0.0 ) throw eBadInputValue;
	in >> fK_delta;    if (fK_delta > 0.0 ) throw eBadInputValue;
}

/* destructor */
DPPrimitiveT::~DPPrimitiveT(void) { }

/* write parameters */
void DPPrimitiveT::Print(ostream& out) const
{
out << " Cohesion-like strength parameter. . . . . . . . = " << falpha_bar << '\n';  //**mien**//
out << " Friction-like parameter . . . . . . . . . . . . = " << ffriction  << '\n';
out << " Dilation parameter. . . . . . . . . . . . . . . = " << fdilation  << '\n';
out << " Deviatoric hardening parameter. . . . . . . . . = " << fH_prime   << '\n';
out << " Volumetric hardening parameter. . . . . . . . . = " << fK_prime   << '\n';
out << " Localized deviatoric hardening parameter. . . . = " << fH_delta   << '\n';
out << " Localized volumetric hardening parameter. . . . = " << fK_delta   << '\n';
}

/***********************************************************************
* Protected
***********************************************************************/

void DPPrimitiveT::PrintName(ostream& out) const
{
out << "    Drucker-Prager Pressure-Dependent Plasticity Model\n";  //**mien**//
out << "    J2 Isotropic/Kinematic\n";
out << "    Hardening with Radial Return\n";
}

/*
* Returns the value value of the yield function given the
* relative stress vector and state variables, where  alpha
* represents isotropic hardening.  NOTE: the relative stress
* should already contain the correction for any kinematic
* hardening.
*/
double DPPrimitiveT::YieldCondition(const dSymMatrixT& devstress, const double meanstress,
			double alpha_dev, double alpha_vol) const  //**mien**//
{
double kTemp;

kTemp  = sqrt32*sqrt(devstress.ScalarProduct());
kTemp += sqrt(3.0)*(-falpha_bar + ffriction*meanstress);
kTemp += alpha_dev;
kTemp += sqrt(3.0)*fdilation*alpha_vol;
return   kTemp;
}
