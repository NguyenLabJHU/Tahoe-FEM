/* $Id: FDCubic2DT.cpp,v 1.10 2005-01-13 00:11:24 paklein Exp $ */
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

/* information about subordinate parameter lists */
void FDCubic2DT::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	FDCubicT::DefineParameters(list);

	/* in-plane rotation angle */
	ParameterT rotation(ParameterT::Double, "rotation_angle_deg");
	rotation.SetDefault(0.0);
	list.AddParameter(rotation, ParameterListT::ZeroOrOnce);
}

/* accept parameter list */
void FDCubic2DT::TakeParameterList(const ParameterListT& list)
{
	/* need to set the in-plane angle before calling the inherited method
	 * since it needs to be set by the time SSCubic2DT::SetModulus is called */
	double angle = 0.0;
	const ParameterT* rotation = list.Parameter("rotation_angle_deg");
	if (rotation) angle = *rotation;
	SetRotation(angle);

	/* inherited */
	FDCubicT::TakeParameterList(list);
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
