/* $Id: SSCubic2DT.cpp,v 1.7 2004-07-15 08:27:09 paklein Exp $ */
/* created: paklein (06/11/1997) */
#include "SSCubic2DT.h"
#include "ThermalDilatationT.h"

using namespace Tahoe;

/* constructor */
SSCubic2DT::SSCubic2DT(void):
	ParameterInterfaceT("small_strain_cubic_2D")
{

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
