/* $Id: NL_E_MatT.cpp,v 1.1.1.1.2.2 2001-06-07 03:01:26 paklein Exp $ */
/* created: paklein (06/13/1997)                                          */
/* Base class for materials with nonlinear elastic behavior               */
/* which is computed from Langrangian coordinates (by the pure            */
/* virtual functions below).                                              */
/* Note: The material tangent moduli and 2nd PK are transformed           */
/* to the spatial quantities, spatial tangent moduli and                  */
/* Cauchy stress.                                                         */
/* Note: The moduli, stress, and stored energy are assumed to be          */
/* _at_least_ functions of the Langrangian strain, to ensure that         */
/* the finite deformation continuum is current for the                    */
/* required material->spatial transformations. Any other                  */
/* dependencies, ie. visco-elasticity must be handled by                  */
/* the derived material classes.                                          */
/* Note: The particular material orientation with respect to the          */
/* global axes is also assumed to be taken care of by the                 */
/* derived material classes.                                              */

#include "NL_E_MatT.h"

/* constructors */
NL_E_MatT::NL_E_MatT(ifstreamT& in, const ElasticT& element):
	FDStructMatT(in, element),
	fE(NumSD()),
	fPK2(NumSD()),
	fModuli(dSymMatrixT::NumValues(NumSD()))
{

}

/* spatial description */
const dMatrixT& NL_E_MatT::c_ijkl(void)
{
	/* strain */
	Compute_E(fE);

	/* derived class function */
	ComputeModuli(fE, fModuli);
	
	/* material -> spatial */
	const dMatrixT& Fmat = F();
	fModuli.SetToScaled(1.0/Fmat.Det(), PushForward(Fmat, fModuli));	
	return fModuli;
}
	
const dSymMatrixT& NL_E_MatT::s_ij(void)
{
	/* strain */
	Compute_E(fE);

	/* derived class function */
	ComputePK2(fE, fPK2);

	/* material -> spatial */
	const dMatrixT& Fmat = F();
	fPK2.SetToScaled(1.0/Fmat.Det(), PushForward(Fmat, fPK2));	
	return fPK2;
}

/* material description */
const dMatrixT& NL_E_MatT::C_IJKL(void)
{
	/* strain */
	Compute_E(fE);

	/* derived class function */
	ComputeModuli(fE, fModuli);
	return fModuli;
}
	
const dSymMatrixT& NL_E_MatT::S_IJ(void)
{
	/* strain */
	Compute_E(fE);

	/* derived class function */
	ComputePK2(fE, fPK2);
	 return fPK2;
}

/* returns the strain energy density for the specified strain */
double NL_E_MatT::StrainEnergyDensity(void)
{
	/* strain */
	Compute_E(fE);

	/* derived class function */
	return ComputeEnergyDensity(fE);
}
