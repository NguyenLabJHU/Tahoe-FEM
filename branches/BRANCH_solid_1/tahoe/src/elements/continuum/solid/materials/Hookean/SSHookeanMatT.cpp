/* $Id: SSHookeanMatT.cpp,v 1.1.1.1 2001-01-29 08:20:30 paklein Exp $ */
/* created: paklein (06/10/1997)                                          */

#include "SSHookeanMatT.h"

/* constructor */
SSHookeanMatT::SSHookeanMatT(ifstreamT& in, const ElasticT& element):
	SSStructMatT(in, element),
	fStress(NumSD()),
	fModulus(dSymMatrixT::NumValues(NumSD()))	
{

}

/* spatial description */
const dMatrixT& SSHookeanMatT::c_ijkl(void) { return fModulus; }
const dSymMatrixT& SSHookeanMatT::s_ij(void)
{
	HookeanStress(fModulus, e(), fStress);
	return fStress;
}

const dMatrixT& SSHookeanMatT::C_IJKL(void) { return fModulus; }
const dSymMatrixT& SSHookeanMatT::S_IJ(void)
{
	HookeanStress(fModulus, e(), fStress);
	return fStress;
}

/* returns the strain energy density for the specified strain */
double SSHookeanMatT::StrainEnergyDensity(void)
{
	return HookeanEnergy(fModulus, e());
}
