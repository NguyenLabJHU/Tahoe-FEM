/* $Id: SSKStV.cpp,v 1.5 2004-07-15 08:27:14 paklein Exp $ */
/* created: paklein (06/10/1997) */
#include "SSKStV.h"

using namespace Tahoe;

/* constructor */
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
ParameterInterfaceT* SSKStV::NewSub(const StringT& name) const
{
	/* inherited */
	ParameterInterfaceT* params = SSHookeanMatT::NewSub(name);
	if (params)
		return params;
	else
		return IsotropicT::NewSub(name);
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
