#include "SSHookean1D.h"

using namespace Tahoe;

/* constructor */
SSHookean1D::SSHookean1D(ifstreamT& in, const SmallStrainT& element):
	SSHookeanMatT(in, element),
	IsotropicT(in)
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
