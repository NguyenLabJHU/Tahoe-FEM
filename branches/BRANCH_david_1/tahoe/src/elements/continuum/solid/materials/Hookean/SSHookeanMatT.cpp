/* $Id: SSHookeanMatT.cpp,v 1.9 2004-07-22 21:09:32 paklein Exp $ */
/* created: paklein (06/10/1997) */
#include "SSHookeanMatT.h"

using namespace Tahoe;

/* constructor */
SSHookeanMatT::SSHookeanMatT(void):
	ParameterInterfaceT("small_strain_Hookean")
{

}

/* set the material support or pass NULL to clear */
void SSHookeanMatT::SetSSMatSupport(const SSMatSupportT* support)
{
	/* inherited */
	SSSolidMatT::SetSSMatSupport(support);
	
	HookeanMatT::Dimension(NumSD());
	fStress.Dimension(dSymMatrixT::int2DimensionT(NumSD()));
}

/* spatial description */
const dMatrixT& SSHookeanMatT::c_ijkl(void) { return Modulus(); }
const dSymMatrixT& SSHookeanMatT::s_ij(void)
{
	HookeanStress(e(), fStress);
	return fStress;
}

const dMatrixT& SSHookeanMatT::C_IJKL(void) { return Modulus(); }
const dSymMatrixT& SSHookeanMatT::S_IJ(void)
{
	HookeanStress(e(), fStress);
	return fStress;
}

/* returns the strain energy density for the specified strain */
double SSHookeanMatT::StrainEnergyDensity(void)
{
	return HookeanEnergy(e());
}

/* information about subordinate parameter lists */
void SSHookeanMatT::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	SSSolidMatT::DefineSubs(sub_list);
	HookeanMatT::DefineSubs(sub_list);
}

/* return the description of the given inline subordinate parameter list */
void SSHookeanMatT::DefineInlineSub(const StringT& name, ParameterListT::ListOrderT& order, 
	SubListT& sub_lists) const
{
	/* inherited */
	if (sub_lists.Length() == 0)
		SSSolidMatT::DefineInlineSub(name, order, sub_lists);
	if (sub_lists.Length() == 0)
		HookeanMatT::DefineInlineSub(name, order, sub_lists);
}

/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* SSHookeanMatT::NewSub(const StringT& name) const
{
	ParameterInterfaceT* sub = SSSolidMatT::NewSub(name);
	if (sub)
		return sub;
	else
		return HookeanMatT::NewSub(name);
}

/* accept parameter list */
void SSHookeanMatT::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	SSSolidMatT::TakeParameterList(list);
	HookeanMatT::TakeParameterList(list);
}
