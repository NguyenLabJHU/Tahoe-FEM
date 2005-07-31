/* $Id: SSKStV.cpp,v 1.1.1.1 2001-01-29 08:20:30 paklein Exp $ */
/* created: paklein (06/10/1997)                                          */

#include "SSKStV.h"

/* constructor */
SSKStV::SSKStV(ifstreamT& in, const ElasticT& element):
	SSHookeanMatT(in, element),
	KStV(in, fModulus)
{

}

/* print parameters */
void SSKStV::Print(ostream& out) const
{
	/* inherited */
	SSHookeanMatT::Print(out);
	KStV::Print(out);
}

/*************************************************************************
* Protected
*************************************************************************/

/* print name */
void SSKStV::PrintName(ostream& out) const
{
	/* inherited */
	SSHookeanMatT::PrintName(out);
	KStV::PrintName(out);
}
