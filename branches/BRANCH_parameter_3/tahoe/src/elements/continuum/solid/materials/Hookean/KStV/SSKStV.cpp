/* $Id: SSKStV.cpp,v 1.4.48.2 2004-06-09 23:17:37 paklein Exp $ */
/* created: paklein (06/10/1997) */
#include "SSKStV.h"

using namespace Tahoe;

/* constructor */
SSKStV::SSKStV(ifstreamT& in, const SSMatSupportT& support):
	ParameterInterfaceT("small_strain_StVenant"),
	SSHookeanMatT(in, support),
	IsotropicT(in)
{

}

SSKStV::SSKStV(void):
	ParameterInterfaceT("small_strain_StVenant")
{

}

/* information about subordinate parameter lists */
void SSKStV::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	SSHookeanMatT::DefineSubs(sub_list);
	IsotropicT::DefineSubs(sub_list);
}

/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* SSKStV::NewSub(const StringT& list_name) const
{
	/* inherited */
	ParameterInterfaceT* params = SSHookeanMatT::NewSub(list_name);
	if (params)
		return params;
	else
		return IsotropicT::NewSub(list_name);
}

/* accept parameter list */
void SSKStV::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	IsotropicT::TakeParameterList(list); /* need moduli before SSHookeanMatT::TakeParameterList */
	SSHookeanMatT::TakeParameterList(list);
}

/*************************************************************************
 * Protected
 *************************************************************************/

/* set (material) tangent modulus */
void SSKStV::SetModulus(dMatrixT& modulus)
{
	IsotropicT::ComputeModuli(modulus);
}
