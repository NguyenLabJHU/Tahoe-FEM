/* $Id: FDKStV2D.cpp,v 1.1.1.1 2001-01-29 08:20:30 paklein Exp $ */
/* created: paklein (06/10/1997)                                          */

#include "FDKStV2D.h"
#include "ThermalDilatationT.h"

/* constructor */
FDKStV2D::FDKStV2D(ifstreamT& in, const ElasticT& element):
	FDHookeanMatT(in, element),
	KStV2D(in, fModulus, fDensity)
{

}

/* print parameters */
void FDKStV2D::Print(ostream& out) const
{
	/* inherited */
	FDHookeanMatT::Print(out);
	KStV2D::Print(out);
}

/* print name */
void FDKStV2D::PrintName(ostream& out) const
{
	/* inherited */
	FDHookeanMatT::PrintName(out);
	KStV2D::PrintName(out);
}

/*************************************************************************
* Private
*************************************************************************/

/* set inverse of thermal transformation - return true if active */
bool FDKStV2D::SetInverseThermalTransformation(dMatrixT& F_trans_inv)
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
