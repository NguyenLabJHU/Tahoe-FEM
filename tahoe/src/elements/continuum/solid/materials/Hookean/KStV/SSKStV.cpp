/* $Id: SSKStV.cpp,v 1.4.32.1 2004-01-21 19:10:08 paklein Exp $ */
/* created: paklein (06/10/1997) */
#include "SSKStV.h"

using namespace Tahoe;

/* constructor */
SSKStV::SSKStV(ifstreamT& in, const SSMatSupportT& support):
	ParameterInterfaceT("small_strain_StVenant"),
	SSHookeanMatT(in, support),
	IsotropicT(in)
{

}

SSKStV::SSKStV(void):
	ParameterInterfaceT("small_strain_StVenant")
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
