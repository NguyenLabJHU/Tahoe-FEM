/* $Id: FDKStV2D.cpp,v 1.1.1.1.2.1 2001-06-06 16:20:45 paklein Exp $ */
/* created: paklein (06/10/1997)                                          */

#include "FDKStV2D.h"
#include "ThermalDilatationT.h"

/* constructor */
FDKStV2D::FDKStV2D(ifstreamT& in, const ElasticT& element):
	FDKStV(in, element),
	Material2DT(in)
{
	/* account for thickness */
	fDensity *= fThickness;
}

/* print parameters */
void FDKStV2D::Print(ostream& out) const
{
	/* inherited */
	FDHookeanMatT::Print(out);
	Material2DT::Print(out);
}

/*************************************************************************
* Protected
*************************************************************************/

/* set (material) tangent modulus */
void FDKStV2D::SetModulus(dMatrixT& modulus)
{
	IsotropicT::ComputeModuli2D(modulus, fConstraintOption);
	modulus *= fThickness;
}

/*************************************************************************
* Private
*************************************************************************/

/* set inverse of thermal transformation - return true if active */
bool FDKStV2D::SetInverseThermalTransformation(dMatrixT& F_trans_inv)
{
	if (fThermal->IsActive())
	{
		double factor = IsotropicT::DilatationFactor2D(fConstraintOption);

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
