/* $Id: FDKStV2D.cpp,v 1.7.30.3 2004-03-03 16:14:57 paklein Exp $ */
/* created: paklein (06/10/1997) */
#include "FDKStV2D.h"
#include "ThermalDilatationT.h"

using namespace Tahoe;

/* constructor */
FDKStV2D::FDKStV2D(ifstreamT& in, const FSMatSupportT& support):
	ParameterInterfaceT("large_strain_StVenant_2D"),
	FDKStV(in, support)
{

}

FDKStV2D::FDKStV2D(void):
	ParameterInterfaceT("large_strain_StVenant_2D")
{

}

/*************************************************************************
 * Protected
 *************************************************************************/

/* set (material) tangent modulus */
void FDKStV2D::SetModulus(dMatrixT& modulus)
{
	IsotropicT::ComputeModuli2D(modulus, Constraint());
}

/*************************************************************************
 * Private
 *************************************************************************/

/* set inverse of thermal transformation - return true if active */
bool FDKStV2D::SetInverseThermalTransformation(dMatrixT& F_trans_inv)
{
	if (fThermal->IsActive())
	{
		/* note - this is approximate at finite strains */
		double factor = IsotropicT::DilatationFactor2D(Constraint());

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
