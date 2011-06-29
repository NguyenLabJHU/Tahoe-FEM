/* $Id: FDHookeanMat2DT.cpp,v 1.2 2004-09-10 22:38:57 paklein Exp $ */
#include "FDHookeanMat2DT.h"
#include "ThermalDilatationT.h"

using namespace Tahoe;

/* constructor */
FDHookeanMat2DT::FDHookeanMat2DT(void):
	ParameterInterfaceT("large_strain_Hookean_2D")
{
	/* reset default value */
	fConstraint = kPlaneStrain;
}

/* accept parameter list */
void FDHookeanMat2DT::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	FDHookeanMatT::TakeParameterList(list);

	//TEMP
	int constraint = list.GetParameter("constraint_2D");
	if (constraint == kPlaneStress)
		ExceptionT::GeneralFail("FDHookeanMat2DT::TakeParameterList", 
			"plain stress not supported");
}

/*************************************************************************
 * Private
 *************************************************************************/

/* set inverse of thermal transformation - return true if active */
bool FDHookeanMat2DT::SetInverseThermalTransformation(dMatrixT& F_trans_inv)
{
	if (fThermal->IsActive())
	{
		//TEMP
		ExceptionT::GeneralFail("FDHookeanMat2DT::SetInverseThermalTransformation", 
			"not implemented");
		return true;
	}
	else
	{
		F_trans_inv.Identity(1.0);
		return false;
	}
}
