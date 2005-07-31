/* $Id: FDKStV.cpp,v 1.1.1.1 2001-01-29 08:20:30 paklein Exp $ */
/* created: paklein (06/10/1997)                                          */

#include "FDKStV.h"

/* constructor */
FDKStV::FDKStV(ifstreamT& in, const ElasticT& element):
	FDHookeanMatT(in, element),
	KStV(in, fModulus)
{

}

/* print parameters */
void FDKStV::Print(ostream& out) const
{
	/* inherited */
	FDHookeanMatT::Print(out);
	KStV::Print(out);
}

/*************************************************************************
* Protected
*************************************************************************/

/* print name */
void FDKStV::PrintName(ostream& out) const
{
	/* inherited */
	FDHookeanMatT::PrintName(out);
	KStV::PrintName(out);
}
