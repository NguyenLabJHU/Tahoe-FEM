/* $Id: SSCubicT.cpp,v 1.4.48.1 2004-04-08 07:32:46 paklein Exp $ */
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

/* print parameters */
void SSCubicT::Print(ostream& out) const
{
	/* inherited */
	SSHookeanMatT::Print(out);
	CubicT::Print(out);
}

/* print name */
void SSCubicT::PrintName(ostream& out) const
{
	/* inherited */
	SSHookeanMatT::PrintName(out);
	CubicT::PrintName(out);
}

/* describe the parameters needed by the interface */
void SSCubicT::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	SSHookeanMatT::DefineParameters(list);
	CubicT::DefineParameters(list);
}

/*************************************************************************
 * Protected
 *************************************************************************/

/* set modulus */
void SSCubicT::SetModulus(dMatrixT& modulus)
{
	CubicT::ComputeModuli(modulus);
}
