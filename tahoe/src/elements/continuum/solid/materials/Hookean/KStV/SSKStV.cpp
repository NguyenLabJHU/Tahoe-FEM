/* $Id: SSKStV.cpp,v 1.1.1.1.2.2 2001-06-22 14:18:04 paklein Exp $ */
/* created: paklein (06/10/1997)                                          */

#include "SSKStV.h"

/* constructor */
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
