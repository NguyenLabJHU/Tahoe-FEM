/* $Id: FDKStV.cpp,v 1.5.30.4 2004-03-04 06:45:25 paklein Exp $ */
/* created: paklein (06/10/1997) */
#include "FDKStV.h"

using namespace Tahoe;

/* constructor */
FDKStV::FDKStV(ifstreamT& in, const FSMatSupportT& support):
	ParameterInterfaceT("large_strain_StVenant"),
	FDHookeanMatT(in, support),
	IsotropicT(in)
{

}

FDKStV::FDKStV(void):
	ParameterInterfaceT("large_strain_StVenant")
{

}

/* print parameters */
void FDKStV::Print(ostream& out) const
{
	/* inherited */
	FDHookeanMatT::Print(out);
	IsotropicT::Print(out);
}

/* print name */
void FDKStV::PrintName(ostream& out) const
{
	/* inherited */
	FDHookeanMatT::PrintName(out);
	out << "    Kirchhoff-St.Venant\n";
}

/* information about subordinate parameter lists */
void FDKStV::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	FDHookeanMatT::DefineSubs(sub_list);
	IsotropicT::DefineSubs(sub_list);
}

/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* FDKStV::NewSub(const StringT& list_name) const
{
	/* inherited */
	ParameterInterfaceT* params = FDHookeanMatT::NewSub(list_name);
	if (params)
		return params;
	else
		return IsotropicT::NewSub(list_name);
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
