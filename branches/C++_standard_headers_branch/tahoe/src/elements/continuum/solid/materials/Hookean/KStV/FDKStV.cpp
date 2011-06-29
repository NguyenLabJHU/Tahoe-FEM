/* $Id: FDKStV.cpp,v 1.6 2004-07-15 08:27:14 paklein Exp $ */
/* created: paklein (06/10/1997) */
#include "FDKStV.h"

using namespace Tahoe;

/* constructor */
FDKStV::FDKStV(void):
	ParameterInterfaceT("large_strain_StVenant")
{

}

/* information about subordinate parameter lists */
void FDKStV::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	FDHookeanMatT::DefineSubs(sub_list);
	IsotropicT::DefineSubs(sub_list);
}

/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* FDKStV::NewSub(const StringT& name) const
{
	/* inherited */
	ParameterInterfaceT* params = FDHookeanMatT::NewSub(name);
	if (params)
		return params;
	else
		return IsotropicT::NewSub(name);
}

/* accept parameter list */
void FDKStV::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	IsotropicT::TakeParameterList(list); /* need moduli before FDHookeanMatT::TakeParameterList */
	FDHookeanMatT::TakeParameterList(list);
}

/*************************************************************************
 * Protected
 *************************************************************************/

/* set (material) tangent modulus */
void FDKStV::SetModulus(dMatrixT& modulus)
{
	IsotropicT::ComputeModuli(modulus);
}
