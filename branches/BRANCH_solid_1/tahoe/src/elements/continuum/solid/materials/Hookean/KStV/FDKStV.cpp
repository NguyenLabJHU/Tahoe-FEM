/* $Id: FDKStV.cpp,v 1.1.1.1.2.1 2001-06-06 16:20:43 paklein Exp $ */
/* created: paklein (06/10/1997)                                          */

#include "FDKStV.h"

/* constructor */
FDKStV::FDKStV(ifstreamT& in, const ElasticT& element):
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
