/* $Id: FDHookeanMatT.cpp,v 1.1.1.1.2.1 2001-06-06 16:22:00 paklein Exp $ */
/* created: paklein (06/10/1997)                                          */

#include "FDHookeanMatT.h"

/* constructor */
FDHookeanMatT::FDHookeanMatT(ifstreamT& in, const ElasticT& element):
	FDStructMatT(in, element),
	HookeanMatT(NumSD()),
	fStress(NumSD()),
	fModulus(dSymMatrixT::NumValues(NumSD()))
{

}

/* initialization */
void FDHookeanMatT::Initialize(void)
{
	/* inherited */
	HookeanMatT::Initialize();
}

/* spatial description */
const dMatrixT& FDHookeanMatT::c_ijkl(void)
{
	/* push forward */
	const dMatrixT& F_mat = F();
	fModulus = PushForward(F_mat, Modulus());
	fModulus /= F_mat.Det();
	return fModulus;
}

const dSymMatrixT& FDHookeanMatT::s_ij(void)
{
	HookeanStress(E(), fStress);
	const dMatrixT& F_mat = F();
	fStress = PushForward(F_mat, fStress);
	fStress /= F_mat.Det();
	return fStress;
}

/* material description */
const dMatrixT& FDHookeanMatT::C_IJKL(void)
{ 
	return Modulus();
}

const dSymMatrixT& FDHookeanMatT::S_IJ(void)
{
	HookeanStress(E(), fStress);
	return fStress;
}

/* returns the strain energy density for the specified strain */
double FDHookeanMatT::StrainEnergyDensity(void)
{
	return HookeanEnergy(E());
}
