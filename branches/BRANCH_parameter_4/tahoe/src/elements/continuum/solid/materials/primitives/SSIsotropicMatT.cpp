/* $Id: SSIsotropicMatT.cpp,v 1.1.6.2 2004-07-12 16:06:26 paklein Exp $ */
#include "SSIsotropicMatT.h"

using namespace Tahoe;

/* constructor */
SSIsotropicMatT::SSIsotropicMatT(void):
	ParameterInterfaceT("small_strain_isotropic")
{

}

/* information about subordinate parameter lists */
void SSIsotropicMatT::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	SSSolidMatT::DefineSubs(sub_list);
	IsotropicT::DefineSubs(sub_list);
}

/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* SSIsotropicMatT::NewSub(const StringT& name) const
{
	/* inherited */
	ParameterInterfaceT* params = SSSolidMatT::NewSub(name);
	if (params)
		return params;
	else
		return IsotropicT::NewSub(name);
}

/* accept parameter list */
void SSIsotropicMatT::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	SSSolidMatT::TakeParameterList(list);
	IsotropicT::TakeParameterList(list);
}
