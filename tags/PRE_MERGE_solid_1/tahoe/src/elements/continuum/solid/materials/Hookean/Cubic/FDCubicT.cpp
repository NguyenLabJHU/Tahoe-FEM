/* $Id: FDCubicT.cpp,v 1.1.1.1 2001-01-29 08:20:30 paklein Exp $ */
/* created: paklein (06/11/1997)                                          */

#include "FDCubicT.h"

/* constructor */
FDCubicT::FDCubicT(ifstreamT& in, const ElasticT& element):
	FDHookeanMatT(in, element),
	CubicT(in, fModulus)
{

}

/* print parameters */
void FDCubicT::Print(ostream& out) const
{
	/* inherited */
	FDHookeanMatT::Print(out);
	CubicT::Print(out);
}

/*************************************************************************
* Protected
*************************************************************************/

/* print name */
void FDCubicT::PrintName(ostream& out) const
{
	/* inherited */
	FDHookeanMatT::PrintName(out);
	CubicT::PrintName(out);
}
