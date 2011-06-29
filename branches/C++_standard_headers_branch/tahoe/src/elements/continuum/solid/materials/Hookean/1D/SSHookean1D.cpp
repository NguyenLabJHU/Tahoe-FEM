/* $Id: SSHookean1D.cpp,v 1.7 2004-07-15 08:27:01 paklein Exp $ */
#include "SSHookean1D.h"

using namespace Tahoe;

/* constructor */
SSHookean1D::SSHookean1D(void):
	ParameterInterfaceT("linear_material_1D")

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
ParameterInterfaceT* SSHookean1D::NewSub(const StringT& name) const
{
	/* inherited */
	ParameterInterfaceT* params = SSHookeanMatT::NewSub(name);
	if (params)
		return params;
	else
		return IsotropicT::NewSub(name);
}

/* accept parameter list */
void SSHookean1D::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	IsotropicT::TakeParameterList(list); /* need to extract moduli before initializing Hookean material */
	SSHookeanMatT::TakeParameterList(list);
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
