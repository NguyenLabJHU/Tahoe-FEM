/* $Id: FDKStV.cpp,v 1.4 2002-11-14 17:06:06 paklein Exp $ */
/* created: paklein (06/10/1997) */
#include "FDKStV.h"

using namespace Tahoe;

/* constructor */
FDKStV::FDKStV(ifstreamT& in, const FDMatSupportT& support):
	FDHookeanMatT(in, support),
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
