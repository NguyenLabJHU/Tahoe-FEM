#include "SSHookean1D.h"

/* constructor */
SSHookean1D::SSHookean1D(ifstreamT& in, const SmallStrainT& element):
	IsotropicT(in),
	SSHookeanMatT(in, element)
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
        IsotropicT::ComputeModuli(modulus);
}
