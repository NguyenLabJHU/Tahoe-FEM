/* $Id: FDHookeanMatT.cpp,v 1.1.1.1 2001-01-29 08:20:30 paklein Exp $ */
/* created: paklein (06/10/1997)                                          */

#include "FDHookeanMatT.h"

/* constructor */
FDHookeanMatT::FDHookeanMatT(ifstreamT& in, const ElasticT& element):
	FDStructMatT(in, element),
	fStress(NumSD()),
	fModulus(dSymMatrixT::NumValues(NumSD()))
{

}

/* spatial description */
const dMatrixT& FDHookeanMatT::c_ijkl(void)
{
	/* set continuum */
	F();
	return C_to_c(fModulus);
}

const dSymMatrixT& FDHookeanMatT::s_ij(void)
{
	HookeanStress(fModulus, E(), fStress);
	return S_to_s(fStress);
}

/* material description */
const dMatrixT& FDHookeanMatT::C_IJKL(void)
{
	return fModulus;
}

const dSymMatrixT& FDHookeanMatT::S_IJ(void)
{
	HookeanStress(fModulus, E(), fStress);
	return fStress;
}

/* returns the strain energy density for the specified strain */
double FDHookeanMatT::StrainEnergyDensity(void)
{
	return HookeanEnergy(fModulus, E());
}
