/* $Id: SSHookeanMatT.cpp,v 1.7 2004-01-10 04:41:12 paklein Exp $ */
/* created: paklein (06/10/1997) */
#include "SSHookeanMatT.h"

using namespace Tahoe;

/* constructor */
SSHookeanMatT::SSHookeanMatT(ifstreamT& in, const SSMatSupportT& support):
	SSSolidMatT(in, support),
	HookeanMatT(NumSD()),
	fStress(NumSD())
{
	SetName("small_strain_Hookean");
}

SSHookeanMatT::SSHookeanMatT(void)
{
	SetName("small_strain_Hookean");
}

/* initialization */
void SSHookeanMatT::Initialize(void)
{
	/* inherited */
	SSSolidMatT::Initialize();
	HookeanMatT::Initialize();
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
