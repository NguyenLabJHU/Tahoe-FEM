/* $Id: NL_E_RotMat2DT.cpp,v 1.6 2003-01-29 07:35:08 paklein Exp $ */
/* created: paklein (06/13/1997) */
#include "NL_E_RotMat2DT.h"

using namespace Tahoe;

/* constructor */
NL_E_RotMat2DT::NL_E_RotMat2DT(ifstreamT& in, const FSMatSupportT& support,
	ConstraintOptionT constraint):
	NL_E_Mat2DT(in, support, constraint),
	Anisotropic2DT(in)
{

}

/* print parameters */
void NL_E_RotMat2DT::Print(ostream& out) const
{
	/* inherited */
	NL_E_Mat2DT::Print(out);
	Anisotropic2DT::Print(out);
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
	
	/* account for thickness */
	fModuli *= fThickness;
	
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

	/* account for thickness */
	fPK2 *= fThickness;

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
	return fThickness*ComputeEnergyDensity(E_nat);
}
