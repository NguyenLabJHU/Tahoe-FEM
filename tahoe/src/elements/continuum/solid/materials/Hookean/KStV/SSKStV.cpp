/* $Id: SSKStV.cpp,v 1.4.32.2 2004-02-18 16:33:44 paklein Exp $ */
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

/* print parameters */
void SSKStV::Print(ostream& out) const
{
	/* inherited */
	SSHookeanMatT::Print(out);
	IsotropicT::Print(out);
}

/* print name */
void SSKStV::PrintName(ostream& out) const
{
	/* inherited */
	SSHookeanMatT::PrintName(out);
	out << "    Kirchhoff-St.Venant\n";
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
	SSHookeanMatT::TakeParameterList(list);
	IsotropicT::TakeParameterList(list);
}

/*************************************************************************
 * Protected
 *************************************************************************/

/* set (material) tangent modulus */
void SSKStV::SetModulus(dMatrixT& modulus)
{
	IsotropicT::ComputeModuli(modulus);
}
