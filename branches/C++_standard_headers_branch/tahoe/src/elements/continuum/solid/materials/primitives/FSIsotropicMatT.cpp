/* $Id: FSIsotropicMatT.cpp,v 1.2 2004-07-15 08:29:19 paklein Exp $ */
#include "FSIsotropicMatT.h"

using namespace Tahoe;

/* constructor */
FSIsotropicMatT::FSIsotropicMatT(void):
	ParameterInterfaceT("large_strain_isotropic")
{

}

/* information about subordinate parameter lists */
void FSIsotropicMatT::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	FSSolidMatT::DefineSubs(sub_list);
	IsotropicT::DefineSubs(sub_list);
}

/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* FSIsotropicMatT::NewSub(const StringT& name) const
{
	/* inherited */
	ParameterInterfaceT* params = FSSolidMatT::NewSub(name);
	if (params)
		return params;
	else
		return IsotropicT::NewSub(name);
}

/* accept parameter list */
void FSIsotropicMatT::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	FSSolidMatT::TakeParameterList(list);
	IsotropicT::TakeParameterList(list);
}
