/* $Id: FDCubic2DT.cpp,v 1.1.1.1 2001-01-29 08:20:30 paklein Exp $ */
/* created: paklein (06/11/1997)                                          */

#include "FDCubic2DT.h"
#include "ThermalDilatationT.h"

/* constructor */
FDCubic2DT::FDCubic2DT(ifstreamT& in, const ElasticT& element):
	FDHookeanMatT(in, element),
	Cubic2DT(in, fModulus, fDensity)
{
	/* transform modulus into global coords */
	TransformOut(fModulus);
}

/* print parameters */
void FDCubic2DT::Print(ostream& out) const
{
	/* inherited */
	FDHookeanMatT::Print(out);
	Cubic2DT::Print(out);
}

/* print name */
void FDCubic2DT::PrintName(ostream& out) const
{
	/* inherited */
	FDHookeanMatT::PrintName(out);
	Cubic2DT::PrintName(out);
}

/*************************************************************************
* Private
*************************************************************************/

/* set inverse of thermal transformation - return true if active */
bool FDCubic2DT::SetInverseThermalTransformation(dMatrixT& F_trans_inv)
{
	if (fThermal->IsActive())
	{
		double factor = DilatationFactor();

		/* assuming isotropic expansion */
		double Fii_inv = 1.0/(1.0 + factor*fThermal->PercentElongation());
		F_trans_inv.Identity(Fii_inv);
		return true;
	}
	else
	{
		F_trans_inv.Identity(1.0);
		return false;
	}
}
