/* $Id: SSCubic2DT.cpp,v 1.6.32.2 2004-03-02 17:46:14 paklein Exp $ */
/* created: paklein (06/11/1997) */
#include "SSCubic2DT.h"
#include "ThermalDilatationT.h"

using namespace Tahoe;

/* constructor */
SSCubic2DT::SSCubic2DT(ifstreamT& in, const SSMatSupportT& support):
	ParameterInterfaceT("small_strain_cubic_2D"),
	SSCubicT(in, support),
	Anisotropic2DT(in)
{

}

/* constructor */
SSCubic2DT::SSCubic2DT(void):
	ParameterInterfaceT("small_strain_cubic_2D")
{

}

/* print parameters */
void SSCubic2DT::Print(ostream& out) const
{
	/* inherited */
	SSCubicT::Print(out);
	Anisotropic2DT::Print(out);
}

double SSCubic2DT::Pressure(void) const
{
	if (Constraint() != kPlaneStress)
		ExceptionT::GeneralFail("SSCubic2DT::Pressure", "not implemented for plane strain");
	return SSCubicT::Pressure();
}

/*************************************************************************
 * Protected
 *************************************************************************/

/* set (material) tangent modulus */
void SSCubic2DT::SetModulus(dMatrixT& modulus)
{
	/* compute modulus in crystal coordinates */
	CubicT::ComputeModuli2D(modulus, Constraint());
	
	/* transform modulus into global coords */
	TransformOut(modulus);
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
		double factor = CubicT::DilatationFactor2D(Constraint());
		thermal_strain.PlusIdentity(factor*fThermal->PercentElongation());
		return true;
	}
	else
		return false;
}
