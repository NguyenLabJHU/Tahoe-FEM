/* $Id: FDCubicT.cpp,v 1.3 2002-07-02 19:55:39 cjkimme Exp $ */
/* created: paklein (06/11/1997)                                          */

#include "FDCubicT.h"

/* constructor */

using namespace Tahoe;

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
