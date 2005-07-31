/* $Id: SSCubicT.cpp,v 1.1.1.1 2001-01-29 08:20:30 paklein Exp $ */
/* created: paklein (06/11/1997)                                          */

#include "SSCubicT.h"

/* constructor */
SSCubicT::SSCubicT(ifstreamT& in, const ElasticT& element):
	SSHookeanMatT(in, element),
	CubicT(in, fModulus)
{

}

/* print parameters */
void SSCubicT::Print(ostream& out) const
{
	/* inherited */
	SSHookeanMatT::Print(out);
	CubicT::Print(out);
}

/* print name */
void SSCubicT::PrintName(ostream& out) const
{
	/* inherited */
	SSHookeanMatT::PrintName(out);
	CubicT::PrintName(out);
}
