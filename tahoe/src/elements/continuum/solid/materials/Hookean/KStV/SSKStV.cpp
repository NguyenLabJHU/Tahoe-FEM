/* $Id: SSKStV.cpp,v 1.1.1.1.2.1 2001-06-06 16:20:43 paklein Exp $ */
/* created: paklein (06/10/1997)                                          */

#include "SSKStV.h"

/* constructor */
SSKStV::SSKStV(ifstreamT& in, const ElasticT& element):
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
