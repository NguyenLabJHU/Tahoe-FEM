/* $Id: FDCubicT.cpp,v 1.2 2001-07-03 01:35:06 paklein Exp $ */
/* created: paklein (06/11/1997)                                          */

#include "FDCubicT.h"

/* constructor */
FDCubicT::FDCubicT(ifstreamT& in, const FiniteStrainT& element):
	FDHookeanMatT(in, element),
	CubicT(in)
{

}

/* print parameters */
void FDCubicT::Print(ostream& out) const
{
	/* inherited */
	FDHookeanMatT::Print(out);
	CubicT::Print(out);
}

/* print name */
void FDCubicT::PrintName(ostream& out) const
{
	/* inherited */
	FDHookeanMatT::PrintName(out);
	CubicT::PrintName(out);
}

/*************************************************************************
* Protected
*************************************************************************/

/* set modulus */
void FDCubicT::SetModulus(dMatrixT& modulus)
{
	CubicT::ComputeModuli(modulus);
}
