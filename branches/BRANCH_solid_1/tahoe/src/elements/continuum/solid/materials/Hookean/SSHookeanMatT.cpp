/* $Id: SSHookeanMatT.cpp,v 1.1.1.1.2.1 2001-06-06 16:22:00 paklein Exp $ */
/* created: paklein (06/10/1997)                                          */

#include "SSHookeanMatT.h"

/* constructor */
SSHookeanMatT::SSHookeanMatT(ifstreamT& in, const ElasticT& element):
	SSStructMatT(in, element),
	HookeanMatT(NumSD()),
	fStress(NumSD())
{

}

/* initialization */
void SSHookeanMatT::Initialize(void)
{
	/* inherited */
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
