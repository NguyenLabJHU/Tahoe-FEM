/* $Id: FDHookeanMatT.cpp,v 1.10 2004-07-15 08:26:56 paklein Exp $ */
/* created: paklein (06/10/1997) */
#include "FDHookeanMatT.h"

using namespace Tahoe;

/* constructor */
FDHookeanMatT::FDHookeanMatT(void):
	ParameterInterfaceT("large_strain_Hookean_material"),
	fLastCall(kNone)	
{

}

/* set the material support or pass NULL to clear */
void FDHookeanMatT::SetFSMatSupport(const FSMatSupportT* support)
{
	/* inherited */
	FSSolidMatT::SetFSMatSupport(support);
	
	/* dimension */
	int nsd = NumSD();
	HookeanMatT::Dimension(nsd);
	fE.Dimension(nsd);
	fStress.Dimension(nsd);
	fModulus.Dimension(dSymMatrixT::NumValues(nsd));
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
	fLastCall = kSpatial;
	return fStress;
}

/* return the pressure associated with the last call to s_ij */
double FDHookeanMatT::Pressure(void) const
{
	if (fLastCall != kSpatial) {
		cout << "\n FDHookeanMatT::Pressure: last call to stress must be in\n"
		     <<   "     the spatial representaion" << endl;
		throw ExceptionT::kGeneralFail;
	}
	return fStress.Trace()/3.0;
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
	
	fLastCall = kMaterial;
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

/* accept parameter list */
void FDHookeanMatT::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	FSSolidMatT::TakeParameterList(list);
	
	/* set the modulus */
	HookeanMatT::Initialize();
}
