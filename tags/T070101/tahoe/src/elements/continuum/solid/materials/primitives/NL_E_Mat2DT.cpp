/* $Id: NL_E_Mat2DT.cpp,v 1.1.1.1 2001-01-29 08:20:25 paklein Exp $ */
/* created: paklein (06/13/1997)                                          */
/* Base class for materials with 2D nonlinear elastic behavior.           */
/* (See notes in NL_E_MatT.h)                                             */

#include "NL_E_Mat2DT.h"

/* constructors */
NL_E_Mat2DT::NL_E_Mat2DT(ifstreamT& in, const ElasticT& element,
	ConstraintOptionT constraint):
	NL_E_MatT(in, element),
	Material2DT(in, constraint)
{
	fDensity *= fThickness;
}

/* print parameters */
void NL_E_Mat2DT::Print(ostream& out) const
{
	/* inherited */
	NL_E_MatT::Print(out);
	Material2DT::Print(out);
}

/* modulus */
const dMatrixT& NL_E_Mat2DT::c_ijkl(void)
{
	/* derived class function */
	ComputeModuli(E(), fModuli);
	
	/* account for thickness */
	fModuli *= fThickness;
	
	/* spatial -> material */
	return C_to_c(fModuli);
}
	
/* stress */
const dSymMatrixT& NL_E_Mat2DT::s_ij(void)
{
	/* derived class function */
	ComputePK2(E(), fPK2);

	/* account for thickness */
	fPK2 *= fThickness;

	/* spatial -> material */
	 return S_to_s(fPK2);
}

/* strain energy density */
double NL_E_Mat2DT::StrainEnergyDensity(void)
{
	/* derived class function */
	return fThickness*ComputeEnergyDensity(E());
}
