/* $Id: IsotropicT.cpp,v 1.2 2001-02-20 00:23:20 paklein Exp $ */
/* created: paklein (06/10/1997)                                          */

#include "IsotropicT.h"

#include <iostream.h>

#include "dMatrixT.h"
#include "fstreamT.h"

/* constructor */
IsotropicT::IsotropicT(ifstreamT& in)
{
	double E, nu;
	in >> E >> nu;
	try { Set_E_nu(E, nu); }
	catch (int exception) { throw eBadInputValue; }
}

IsotropicT::IsotropicT(void)
{
	try { Set_E_nu(0.0, 0.0); }
	catch (int exception) { throw eBadInputValue; }
}

/* set moduli */
void IsotropicT::Set_E_nu(double E, double nu)
{
	fYoung = E;
	fPoisson = nu;

	/* checks */
	if (fYoung < 0.0) throw eGeneralFail;
	if (fPoisson > 0.5 || fPoisson < -1.0) throw eGeneralFail;
	
	/* compute remaining moduli */
	fMu     = 0.5*fYoung/(1.0 + fPoisson);
	fLambda = 2.0*fMu*fPoisson/(1.0 - 2.0*fPoisson);
	fKappa  = fLambda + 2.0/3.0*fMu;
}

void IsotropicT::Set_mu_kappa(double mu, double kappa)
{
	fMu = mu;
	fKappa = kappa;

	/* checks */
	if (fMu < 0.0 || fKappa < 0.0) throw eGeneralFail;

	/* set moduli */
	fYoung = (9.0*fKappa*fMu)/(3.0*fKappa + fMu);
	fPoisson = (3.0*fKappa - 2.0*fMu)/(6.0*fKappa + 2.0*fMu);
	fLambda = 2.0*fMu*fPoisson/(1.0 - 2.0*fPoisson);
}

/* I/O operators */
void IsotropicT::Print(ostream& out) const
{
	out << " Young's modulus . . . . . . . . . . . . . . . . = " << fYoung   << '\n';
	out << " Poisson's ratio . . . . . . . . . . . . . . . . = " << fPoisson << '\n';
	out << " Shear modulus . . . . . . . . . . . . . . . . . = " << fMu      << '\n';
	out << " Bulk modulus. . . . . . . . . . . . . . . . . . = " << fKappa   << '\n';
}

/*************************************************************************
* Protected
*************************************************************************/

/* compute the symetric Cij reduced index matrix */
void IsotropicT::ComputeModuli(dMatrixT& moduli, double mu, double lambda) const
{
	moduli = 0.0;

	if (moduli.Rows() == 3) //StressDimension in dSymMatrixT
	{
		moduli(1,1) = moduli(0,0) = lambda + 2.0*mu;
		moduli(0,1) = lambda;
		moduli(2,2) = mu;
	}
	else if (moduli.Rows() == 6)
	{
		moduli(2,2) = moduli(1,1) = moduli(0,0) = lambda + 2.0*mu;
		moduli(1,2) = moduli(0,1) = moduli(0,2) = lambda;
		moduli(5,5) = moduli(4,4) = moduli(3,3) = mu;
	}
	else throw eGeneralFail;

	/* symmetric */
	moduli.CopySymmetric();
}
