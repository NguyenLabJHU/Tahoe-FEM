/* $Id: SSHookean1D.cpp,v 1.6.2.1 2004-01-21 19:10:05 paklein Exp $ */
#include "SSHookean1D.h"

using namespace Tahoe;

/* constructor */
SSHookean1D::SSHookean1D(ifstreamT& in, const SSMatSupportT& support):
	ParameterInterfaceT("small_strain_StVenant_1D"),
	SSHookeanMatT(in, support),
	IsotropicT(in)
{

}

SSHookean1D::SSHookean1D(void):
	ParameterInterfaceT("small_strain_StVenant_1D")

{

}

/* print parameters */
void SSHookean1D::Print(ostream& out) const
{
	/* inherited */
	SSHookeanMatT::Print(out);
	IsotropicT::Print(out);
}

/* print name */
void SSHookean1D::PrintName(ostream& out) const
{
        /* inherited */
        SSHookeanMatT::PrintName(out);
        out << "    1D SS Hookean\n";
}

/*************************************************************************
* Protected
*************************************************************************/

/* set (material) tangent modulus */
void SSHookean1D::SetModulus(dMatrixT& modulus)
{
	if (modulus.Rows() != 1 || modulus.Cols() != 1) throw ExceptionT::kSizeMismatch;
	modulus = Young();
}
