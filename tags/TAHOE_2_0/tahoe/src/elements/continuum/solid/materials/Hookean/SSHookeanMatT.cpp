/* $Id: SSHookeanMatT.cpp,v 1.8 2004-07-15 08:26:56 paklein Exp $ */
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

/* accept parameter list */
void SSHookeanMatT::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	SSSolidMatT::TakeParameterList(list);
	
	/* set the modulus */
	HookeanMatT::Initialize();
}
