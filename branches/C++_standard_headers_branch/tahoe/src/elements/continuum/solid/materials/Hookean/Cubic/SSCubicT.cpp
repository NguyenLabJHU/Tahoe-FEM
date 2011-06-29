/* $Id: SSCubicT.cpp,v 1.5 2004-07-15 08:27:05 paklein Exp $ */
/* created: paklein (06/11/1997) */
#include "SSCubicT.h"

using namespace Tahoe;

/* constructor */
SSCubicT::SSCubicT(void):
	ParameterInterfaceT("small_strain_cubic")
{

}

/* describe the parameters needed by the interface */
void SSCubicT::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	SSHookeanMatT::DefineParameters(list);
	CubicT::DefineParameters(list);
}

/* information about subordinate parameter lists */
void SSCubicT:: TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	CubicT::TakeParameterList(list); /* cubic parameters must be extracted first */
	SSHookeanMatT::TakeParameterList(list);
}

/*************************************************************************
 * Protected
 *************************************************************************/

/* set modulus */
void SSCubicT::SetModulus(dMatrixT& modulus)
{
	CubicT::ComputeModuli(modulus);
}
