/* $Id: FDHookeanMatT.cpp,v 1.5 2002-07-02 19:55:38 cjkimme Exp $ */
/* created: paklein (06/10/1997)                                          */

#include "FDHookeanMatT.h"

/* constructor */

using namespace Tahoe;

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
	const dMatrixT& F_mech = F_mechanical();
	fModulus = PushForward(F_mech, Modulus());
	fModulus /= F_mech.Det();
	return fModulus;
}

const dSymMatrixT& FDHookeanMatT::s_ij(void)
{
	/* get mechanical part of the deformation gradient */
	const dMatrixT& F_mech = F_mechanical();

	/* strain */
	Compute_E(F_mech, fE);

	/* compute stress */
	HookeanStress(fE, fStress);

	/* push forward */
	fStress = PushForward(F_mech, fStress);
	fStress /= F_mech.Det();
	return fStress;
}

/* material description */
const dMatrixT& FDHookeanMatT::C_IJKL(void)
{
	/* has thermal strain */
	if (HasThermalStrain())
	{
		/* inverse thermal strain */
		const dMatrixT& F_t_inv = F_thermal_inverse();
	
		/* pull back */
		fModulus = PushForward(F_t_inv, Modulus());
		fModulus /= F_t_inv.Det();
		return fModulus;
	}
	else /* no thermal strain */
		return Modulus();
}

const dSymMatrixT& FDHookeanMatT::S_IJ(void)
{
	/* get mechanical part of the deformation gradient */
	const dMatrixT& F_mech = F_mechanical();

	/* strain */
	Compute_E(F_mech, fE);

	/* compute stress */
	HookeanStress(fE, fStress);

	/* has thermal strain */
	if (HasThermalStrain())
	{
		/* inverse thermal strain */
		const dMatrixT& F_t_inv = F_thermal_inverse();
	
		/* pull back */
		fStress = PushForward(F_t_inv, fStress);
		fStress /= F_t_inv.Det();
	}
	
	return fStress;
}

/* returns the strain energy density for the specified strain */
double FDHookeanMatT::StrainEnergyDensity(void)
{
	/* get mechanical part of the deformation gradient */
	const dMatrixT& F_mech = F_mechanical();

	/* strain */
	Compute_E(F_mech, fE);
	
	/* compute strain energy density */
	return HookeanEnergy(fE);
}
