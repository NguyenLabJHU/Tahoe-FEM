/* $Id: FDKStV.cpp,v 1.2 2001-07-03 01:35:09 paklein Exp $ */
/* created: paklein (06/10/1997)                                          */

#include "FDKStV.h"

/* constructor */
FDKStV::FDKStV(ifstreamT& in, const FiniteStrainT& element):
	FDHookeanMatT(in, element),
	IsotropicT(in)
{

}

/* print parameters */
void FDKStV::Print(ostream& out) const
{
	/* inherited */
	FDHookeanMatT::Print(out);
	IsotropicT::Print(out);
}

/* print name */
void FDKStV::PrintName(ostream& out) const
{
	/* inherited */
	FDHookeanMatT::PrintName(out);
	out << "    Kirchhoff-St.Venant\n";
}

/*************************************************************************
* Protected
*************************************************************************/

/* set (material) tangent modulus */
void FDKStV::SetModulus(dMatrixT& modulus)
{
	IsotropicT::ComputeModuli(modulus);
}
