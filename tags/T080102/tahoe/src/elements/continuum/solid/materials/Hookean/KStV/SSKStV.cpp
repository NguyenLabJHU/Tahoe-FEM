/* $Id: SSKStV.cpp,v 1.3 2002-07-02 19:55:41 cjkimme Exp $ */
/* created: paklein (06/10/1997)                                          */

#include "SSKStV.h"

/* constructor */

using namespace Tahoe;

SSKStV::SSKStV(ifstreamT& in, const SmallStrainT& element):
	SSHookeanMatT(in, element),
	IsotropicT(in)
{

}

/* print parameters */
void SSKStV::Print(ostream& out) const
{
	/* inherited */
	SSHookeanMatT::Print(out);
	IsotropicT::Print(out);
}

/* print name */
void SSKStV::PrintName(ostream& out) const
{
	/* inherited */
	SSHookeanMatT::PrintName(out);
	out << "    Kirchhoff-St.Venant\n";
}

/*************************************************************************
* Protected
*************************************************************************/

/* set (material) tangent modulus */
void SSKStV::SetModulus(dMatrixT& modulus)
{
	IsotropicT::ComputeModuli(modulus);
}
