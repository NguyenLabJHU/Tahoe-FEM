/* $Id: NL_E_Mat2DT.cpp,v 1.5.30.1 2004-01-21 19:10:26 paklein Exp $ */
/* created: paklein (06/13/1997) */
#include "NL_E_Mat2DT.h"

using namespace Tahoe;

/* constructors */
NL_E_Mat2DT::NL_E_Mat2DT(ifstreamT& in, const FSMatSupportT& support, ConstraintOptionT constraint):
	ParameterInterfaceT("large_strain_E_material_2D"),
	NL_E_MatT(in, support),
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
	/* strain */
	Compute_E(fE);

	/* derived class function */
	ComputeModuli(fE, fModuli);
	
	/* material -> spatial */
	const dMatrixT& Fmat = F();
	fModuli.SetToScaled(1.0/Fmat.Det(), PushForward(Fmat, fModuli));	
	return fModuli;
}
	
/* stress */
const dSymMatrixT& NL_E_Mat2DT::s_ij(void)
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

/* strain energy density */
double NL_E_Mat2DT::StrainEnergyDensity(void)
{
	/* strain */
	Compute_E(fE);

	/* derived class function */
	return ComputeEnergyDensity(fE);
}
