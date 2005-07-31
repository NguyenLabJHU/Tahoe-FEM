/* $Id: NL_E_RotMat2DT.cpp,v 1.1.1.1 2001-01-29 08:20:25 paklein Exp $ */
/* created: paklein (06/13/1997)                                          */
/* Base class for materials with 2D nonlinear elastic behavior            */
/* with in-plane orientation with respect to global coordinate            */
/* axes, ie. the moduli, stress, and strain energy density functions      */
/* are formulated in the material's natural coordinates.                  */
/* (See notes in NL_E_Mat2DT.h)                                           */

#include "NL_E_RotMat2DT.h"

/* constructor */
NL_E_RotMat2DT::NL_E_RotMat2DT(ifstreamT& in, const ElasticT& element,
	ConstraintOptionT constraint):
	NL_E_Mat2DT(in, element, constraint),
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
	/* compute strain in natural coords */
	const dSymMatrixT& E_nat = TransformIn(E());

	/* derived class function */
	ComputeModuli(E_nat, fModuli);
	
	/* account for thickness */
	fModuli *= fThickness;
	
	/* natural -> spatial -> material */
	return C_to_c(fModuli, Q());
}
	
/* stresses */
const dSymMatrixT& NL_E_RotMat2DT::s_ij(void)
{
	/* compute strain in natural coords */
	const dSymMatrixT& E_nat = TransformIn(E());

	/* derived class function */
	ComputePK2(E_nat, fPK2);

	/* account for thickness */
	fPK2 *= fThickness;

	/* natural -> spatial -> material */
	 return S_to_s(fPK2, Q());
}

/* strain energy density */
double NL_E_RotMat2DT::StrainEnergyDensity(void)
{
	/* compute strain in natural coords */
	const dSymMatrixT& E_nat = TransformIn(E());

	/* derived class function */
	return fThickness*ComputeEnergyDensity(E_nat);
}
