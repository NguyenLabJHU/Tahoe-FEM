/* $Id: SSKStV.cpp,v 1.3.8.1 2002-10-28 06:48:54 paklein Exp $ */
/* created: paklein (06/10/1997) */
#include "SSKStV.h"

using namespace Tahoe;

/* constructor */
SSKStV::SSKStV(ifstreamT& in, const SSMatSupportT& support):
	SSHookeanMatT(in, support),
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
