/* $Id: SSCubic2DT.cpp,v 1.8 2005-01-13 00:11:24 paklein Exp $ */
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

/* information about subordinate parameter lists */
void SSCubic2DT::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	SSCubicT::DefineParameters(list);

	/* in-plane rotation angle */
	ParameterT rotation(ParameterT::Double, "rotation_angle_deg");
	rotation.SetDefault(0.0);
	list.AddParameter(rotation, ParameterListT::ZeroOrOnce);
}

/* accept parameter list */
void SSCubic2DT::TakeParameterList(const ParameterListT& list)
{
	/* need to set the in-plane angle before calling the inherited method
	 * since it needs to be set by the time SSCubic2DT::SetModulus is called */
	double angle = 0.0;
	const ParameterT* rotation = list.Parameter("rotation_angle_deg");
	if (rotation) angle = *rotation;
	SetRotation(angle);

	/* inherited */
	SSCubicT::TakeParameterList(list);
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
