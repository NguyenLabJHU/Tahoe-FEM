/* $Id: SSKStV.cpp,v 1.2.6.1 2002-06-27 18:03:12 cjkimme Exp $ */
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
