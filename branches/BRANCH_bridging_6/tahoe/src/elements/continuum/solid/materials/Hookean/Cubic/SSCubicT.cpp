/* $Id: SSCubicT.cpp,v 1.4 2002-11-14 17:06:04 paklein Exp $ */
/* created: paklein (06/11/1997) */
#include "SSCubicT.h"

using namespace Tahoe;

/* constructor */
SSCubicT::SSCubicT(ifstreamT& in, const SSMatSupportT& support):
	SSHookeanMatT(in, support),
	CubicT(in)
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

/*************************************************************************
* Protected
*************************************************************************/

/* set modulus */
void SSCubicT::SetModulus(dMatrixT& modulus)
{
	CubicT::ComputeModuli(modulus);
}
