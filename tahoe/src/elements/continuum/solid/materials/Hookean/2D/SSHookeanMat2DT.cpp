/* $Id: SSHookeanMat2DT.cpp,v 1.1 2004-07-22 21:09:37 paklein Exp $ */
#include "SSHookeanMat2DT.h"
#include "ThermalDilatationT.h"

using namespace Tahoe;

/* constructor */
SSHookeanMat2DT::SSHookeanMat2DT(void):
	ParameterInterfaceT("small_strain_Hookean_2D")
{

}

/* describe the parameters needed by the interface */
void SSHookeanMat2DT::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	SSHookeanMatT::DefineParameters(list);
	
	/* 2D option must be plain stress */
	ParameterT& constraint = list.GetParameter("constraint_2D");
	constraint.SetDefault(kPlaneStrain);
}

/* accept parameter list */
void SSHookeanMat2DT::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	SSHookeanMatT::TakeParameterList(list);

	//TEMP
	int constraint = list.GetParameter("constraint_2D");
	if (constraint == kPlaneStress)
		ExceptionT::GeneralFail("SSHookeanMat2DT::TakeParameterList", 
			"plain stress not supported");
}

/*************************************************************************
 * Private
 *************************************************************************/

/* set the internal thermal strain */
bool SSHookeanMat2DT::SetThermalStrain(dSymMatrixT& thermal_strain)
{
	thermal_strain = 0.0;
	if (fThermal->IsActive())
	{
		//TEMP
		ExceptionT::GeneralFail("SSHookeanMat2DT::SetInverseThermalTransformation", 
			"not implemented");
		return true;
	}
	else
		return false;
}
