/* $Id: FDCubic2DT.cpp,v 1.9 2004-07-15 08:27:09 paklein Exp $ */
/* created: paklein (06/11/1997) */
#include "FDCubic2DT.h"
#include "ThermalDilatationT.h"

using namespace Tahoe;

/* constructor */
FDCubic2DT::FDCubic2DT(void):
	ParameterInterfaceT("large_strain_cubic_2D")
{

}

double FDCubic2DT::Pressure(void) const
{
	if (Constraint() != kPlaneStress)
		ExceptionT::GeneralFail("FDCubic2DT::Pressure", "not implemented for plane strain");
	return FDCubicT::Pressure();
}

/*************************************************************************
 * Protected
 *************************************************************************/

/* set (material) tangent modulus */
void FDCubic2DT::SetModulus(dMatrixT& modulus)
{
	/* compute modulus in crystal coordinates */
	CubicT::ComputeModuli2D(modulus, Constraint());
	
	/* transform modulus into global coords */
	TransformOut(modulus);
}

/*************************************************************************
 * Private
 *************************************************************************/

/* set inverse of thermal transformation - return true if active */
bool FDCubic2DT::SetInverseThermalTransformation(dMatrixT& F_trans_inv)
{
	if (fThermal->IsActive())
	{
		/* note - this is approximate at finite strains */
		double factor = CubicT::DilatationFactor2D(Constraint());

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
