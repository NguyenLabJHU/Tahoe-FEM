/* $Id: MFGP_SSIsotropicMat.cpp,v 1.1 2005-01-06 22:51:27 kyonten Exp $ */
#include "MFGP_SSIsotropicMatT.h"

using namespace Tahoe;

/* constructor */
MFGP_SSIsotropicMatT::MFGP_SSIsotropicMatT(void):
	ParameterInterfaceT("mfgp_small_strain_isotropic")
{

}

/* information about subordinate parameter lists */
void MFGP_SSIsotropicMatT::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	MFGP_SSSolidMatT::DefineSubs(sub_list);
	IsotropicT::DefineSubs(sub_list);
}

/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* MFGP_SSIsotropicMatT::NewSub(const StringT& name) const
{
	/* inherited */
	ParameterInterfaceT* params = MFGP_SSSolidMatT::NewSub(name);
	if (params)
		return params;
	else
		return IsotropicT::NewSub(name);
}

/* accept parameter list */
void MFGP_SSIsotropicMatT::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	MFGP_SSSolidMatT::TakeParameterList(list);
	IsotropicT::TakeParameterList(list);
}
