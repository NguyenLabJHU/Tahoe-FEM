/* $Id: SSCubicT.cpp,v 1.4.48.3 2004-06-09 23:17:35 paklein Exp $ */
/* created: paklein (06/11/1997) */
#include "SSCubicT.h"

using namespace Tahoe;

/* constructor */
SSCubicT::SSCubicT(ifstreamT& in, const SSMatSupportT& support):
	ParameterInterfaceT("small_strain_cubic"),
	SSHookeanMatT(in, support),
	CubicT(in)
{

}

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
