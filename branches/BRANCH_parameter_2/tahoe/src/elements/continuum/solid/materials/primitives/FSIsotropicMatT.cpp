/* $Id: FSIsotropicMatT.cpp,v 1.1.2.1 2004-02-18 16:33:51 paklein Exp $ */
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
ParameterInterfaceT* FSIsotropicMatT::NewSub(const StringT& list_name) const
{
	/* inherited */
	ParameterInterfaceT* params = FSSolidMatT::NewSub(list_name);
	if (params)
		return params;
	else
		return IsotropicT::NewSub(list_name);
}

/* accept parameter list */
void FSIsotropicMatT::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	FSSolidMatT::TakeParameterList(list);
	IsotropicT::TakeParameterList(list);
}
