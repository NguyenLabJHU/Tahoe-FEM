/* $Id: SSHookean1D.cpp,v 1.6.18.2 2004-06-09 23:17:34 paklein Exp $ */
#include "SSHookean1D.h"

using namespace Tahoe;

/* constructor */
SSHookean1D::SSHookean1D(ifstreamT& in, const SSMatSupportT& support):
	ParameterInterfaceT("small_strain_StVenant_1D"),
	SSHookeanMatT(in, support),
	IsotropicT(in)
{

}

SSHookean1D::SSHookean1D(void):
	ParameterInterfaceT("small_strain_StVenant_1D")

{

}

/* information about subordinate parameter lists */
void SSHookean1D::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	SSHookeanMatT::DefineSubs(sub_list);
	IsotropicT::DefineSubs(sub_list);
}

/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* SSHookean1D::NewSub(const StringT& list_name) const
{
	/* inherited */
	ParameterInterfaceT* params = SSHookeanMatT::NewSub(list_name);
	if (params)
		return params;
	else
		return IsotropicT::NewSub(list_name);
}

/* accept parameter list */
void SSHookean1D::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	SSHookeanMatT::TakeParameterList(list);
	IsotropicT::TakeParameterList(list);
}

/*************************************************************************
* Protected
*************************************************************************/

/* set (material) tangent modulus */
void SSHookean1D::SetModulus(dMatrixT& modulus)
{
	if (modulus.Rows() != 1 || modulus.Cols() != 1) throw ExceptionT::kSizeMismatch;
	modulus = Young();
}
