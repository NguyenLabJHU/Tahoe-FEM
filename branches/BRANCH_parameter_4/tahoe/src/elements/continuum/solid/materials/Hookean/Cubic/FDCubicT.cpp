/* $Id: FDCubicT.cpp,v 1.5.54.1 2004-07-06 06:53:31 paklein Exp $ */
/* created: paklein (06/11/1997) */
#include "FDCubicT.h"

using namespace Tahoe;

/* constructor */
FDCubicT::FDCubicT(void):
	ParameterInterfaceT("large_strain_cubic")
{

}

/* information about subordinate parameter lists */
void FDCubicT:: DefineParameters(ParameterListT& list) const
{
	/* inherited */
	FDHookeanMatT::DefineParameters(list);
	CubicT::DefineParameters(list);
}

/* information about subordinate parameter lists */
void FDCubicT:: TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	CubicT::TakeParameterList(list); /* cubic parameters must be extracted first */
	FDHookeanMatT::TakeParameterList(list);
}

/*************************************************************************
* Protected
*************************************************************************/

/* set modulus */
void FDCubicT::SetModulus(dMatrixT& modulus)
{
	CubicT::ComputeModuli(modulus);
}
