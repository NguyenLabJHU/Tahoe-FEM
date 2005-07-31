/* $Id: SSCubic2DT.cpp,v 1.1.1.1 2001-01-29 08:20:30 paklein Exp $ */
/* created: paklein (06/11/1997)                                          */

#include "SSCubic2DT.h"
#include "ThermalDilatationT.h"

/* constructor */
SSCubic2DT::SSCubic2DT(ifstreamT& in, const ElasticT& element):
	SSHookeanMatT(in, element),
	Cubic2DT(in, fModulus, fDensity)
{
	/* transform moduli into global coords */
	TransformOut(fModulus);
}

/* print parameters */
void SSCubic2DT::Print(ostream& out) const
{
	/* inherited */
	SSHookeanMatT::Print(out);
	Cubic2DT::Print(out);
}

/* print name */
void SSCubic2DT::PrintName(ostream& out) const
{
	/* inherited */
	SSHookeanMatT::PrintName(out);
	Cubic2DT::PrintName(out);
}

/*************************************************************************
* Private
*************************************************************************/

/* set the internal thermal strain */
bool SSCubic2DT::SetThermalStrain(dSymMatrixT& thermal_strain)
{
	thermal_strain = 0.0;
	if (fThermal->IsActive())
	{
		double factor = DilatationFactor();
		thermal_strain.PlusIdentity(factor*fThermal->PercentElongation());
		return true;
	}
	else
		return false;
}
