/* $Id: SSCubicT.cpp,v 1.3 2002-07-02 19:55:39 cjkimme Exp $ */
/* created: paklein (06/11/1997)                                          */

#include "SSCubicT.h"

/* constructor */

using namespace Tahoe;

SSCubicT::SSCubicT(ifstreamT& in, const SmallStrainT& element):
	SSHookeanMatT(in, element),
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
