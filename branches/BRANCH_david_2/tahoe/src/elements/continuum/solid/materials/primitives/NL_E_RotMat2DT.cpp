/* $Id: NL_E_RotMat2DT.cpp,v 1.7 2004-07-15 08:29:20 paklein Exp $ */
/* created: paklein (06/13/1997) */
#include "NL_E_RotMat2DT.h"

using namespace Tahoe;

/* constructor */
NL_E_RotMat2DT::NL_E_RotMat2DT(ifstreamT& in, const FSMatSupportT& support, ConstraintT constraint):
	ParameterInterfaceT("large_strain_E_material_2D"),
//	NL_E_MatT(in, support),
	Anisotropic2DT(in)
{
ExceptionT::Stop("NL_E_RotMat2DT::NL_E_RotMat2DT", "out of date");
}

/* modulus */
const dMatrixT& NL_E_RotMat2DT::c_ijkl(void)
{
	/* compute strain */
	Compute_E(fE);

	/* compute strain in natural coords */
	const dSymMatrixT& E_nat = TransformIn(fE);

	/* derived class function */
	ComputeModuli(E_nat, fModuli);
	
	/* natural -> spatial -> material */
	//return C_to_c(fModuli, Q());
	throw ExceptionT::kGeneralFail;
	return fModuli;
}
	
/* stresses */
const dSymMatrixT& NL_E_RotMat2DT::s_ij(void)
{
	/* compute strain */
	Compute_E(fE);

	/* compute strain in natural coords */
	const dSymMatrixT& E_nat = TransformIn(fE);

	/* derived class function */
	ComputePK2(E_nat, fPK2);

	/* natural -> spatial -> material */
	//return S_to_s(fPK2, Q());
	throw ExceptionT::kGeneralFail;
	return fPK2;
}

/* strain energy density */
double NL_E_RotMat2DT::StrainEnergyDensity(void)
{
	/* compute strain */
	Compute_E(fE);

	/* compute strain in natural coords */
	const dSymMatrixT& E_nat = TransformIn(fE);

	/* derived class function */
	return ComputeEnergyDensity(E_nat);
}
