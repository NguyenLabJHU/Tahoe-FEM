/* $Id: FDCubic2DT.cpp,v 1.8.46.1 2004-04-08 07:32:48 paklein Exp $ */
/* created: paklein (06/11/1997) */
#include "FDCubic2DT.h"
#include "ThermalDilatationT.h"

using namespace Tahoe;

/* constructor */
FDCubic2DT::FDCubic2DT(ifstreamT& in, const FSMatSupportT& support):
	ParameterInterfaceT("large_strain_cubic_2D"),
	FDCubicT(in, support),
	Anisotropic2DT(in)
{

}

FDCubic2DT::FDCubic2DT(void):
	ParameterInterfaceT("large_strain_cubic_2D")
{

}

/* print parameters */
void FDCubic2DT::Print(ostream& out) const
{
	/* inherited */
	FDCubicT::Print(out);
	Anisotropic2DT::Print(out);
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
