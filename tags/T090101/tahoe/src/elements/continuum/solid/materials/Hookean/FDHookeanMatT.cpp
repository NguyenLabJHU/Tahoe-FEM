/* $Id: FDHookeanMatT.cpp,v 1.3 2001-08-21 19:13:56 paklein Exp $ */
/* created: paklein (06/10/1997)                                          */

#include "FDHookeanMatT.h"

/* constructor */
FDHookeanMatT::FDHookeanMatT(ifstreamT& in, const FiniteStrainT& element):
	FDStructMatT(in, element),
	HookeanMatT(NumSD()),
	fE(NumSD()),
	fStress(NumSD()),
	fModulus(dSymMatrixT::NumValues(NumSD()))
{

}

/* initialization */
void FDHookeanMatT::Initialize(void)
{
	/* inherited */
	FDStructMatT::Initialize();
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
	/* strain */
	Compute_E(fE);

	/* compute stress */
	HookeanStress(fE, fStress);
	const dMatrixT& F_mat = F();

	/* push forward */
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
	/* strain */
	Compute_E(fE);

	/* compute stress */
	HookeanStress(fE, fStress);
	return fStress;
}

/* returns the strain energy density for the specified strain */
double FDHookeanMatT::StrainEnergyDensity(void)
{
	Compute_E(fE);
	return HookeanEnergy(fE);
}
