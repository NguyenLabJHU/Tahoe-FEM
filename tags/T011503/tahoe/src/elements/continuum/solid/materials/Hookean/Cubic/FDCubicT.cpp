/* $Id: FDCubicT.cpp,v 1.4 2002-11-14 17:06:04 paklein Exp $ */
/* created: paklein (06/11/1997) */
#include "FDCubicT.h"

using namespace Tahoe;

/* constructor */
FDCubicT::FDCubicT(ifstreamT& in, const FDMatSupportT& support):
	FDHookeanMatT(in, support),
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
