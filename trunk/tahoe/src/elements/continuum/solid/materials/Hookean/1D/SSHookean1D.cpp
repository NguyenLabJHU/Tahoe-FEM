/* $Id: SSHookean1D.cpp,v 1.6 2004-01-10 04:41:14 paklein Exp $ */
#include "SSHookean1D.h"

using namespace Tahoe;

/* constructor */
SSHookean1D::SSHookean1D(ifstreamT& in, const SSMatSupportT& support):
	SSHookeanMatT(in, support),
	IsotropicT(in)
{
	SetName("small_strain_Hookean_1D");
}

SSHookean1D::SSHookean1D(void)
{
	SetName("small_strain_Hookean_1D");
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
