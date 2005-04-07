/* $Id: IsotropicT.cpp,v 1.12 2005-01-19 08:57:20 paklein Exp $ */
/* created: paklein (06/10/1997) */
#include "IsotropicT.h"

#include <iostream.h>

#include "dMatrixT.h"
#include "ifstreamT.h"
#include "ofstreamT.h"
#include "ParameterContainerT.h"

using namespace Tahoe;

/* constructor */
IsotropicT::IsotropicT(ifstreamT& in):
	ParameterInterfaceT("isotropic")
{
	double E, nu;
	in >> E >> nu;
	try { Set_E_nu(E, nu); }
	catch (ExceptionT::CodeT exception) { throw ExceptionT::kBadInputValue; }
}

IsotropicT::IsotropicT(void):
	ParameterInterfaceT("isotropic")
{
	try { Set_E_nu(0.0, 0.0); }
	catch (ExceptionT::CodeT exception) { throw ExceptionT::kBadInputValue; }
}

/* set moduli */
void IsotropicT::Set_E_nu(double E, double nu)
{
	fYoung = E;
	fPoisson = nu;

	/* checks */
	if (fYoung < 0.0) throw ExceptionT::kGeneralFail;
	if (fPoisson > 0.5 || fPoisson < -1.0) throw ExceptionT::kGeneralFail;
	
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
	if (fMu < 0.0 || fKappa < 0.0) throw ExceptionT::kGeneralFail;

	/* set moduli */
	fYoung = (9.0*fKappa*fMu)/(3.0*fKappa + fMu);
	fPoisson = (3.0*fKappa - 2.0*fMu)/(6.0*fKappa + 2.0*fMu);
	fLambda = 2.0*fMu*fPoisson/(1.0 - 2.0*fPoisson);
}
void IsotropicT::Set_PurePlaneStress_mu_lambda(double mu, double lambda)
{
	fMu = mu;
	fLambda = lambda;
	fKappa = mu+lambda;

	/* checks */
	if (fMu < 0.0 || fKappa < 0.0) throw ExceptionT::kGeneralFail;

	/* set moduli */
	fYoung = 4.0*mu*(lambda + mu)/(lambda + 2.0*mu);
	fPoisson = lambda/(lambda + 2.0*mu);
}

/* I/O operators */
void IsotropicT::Print(ostream& out) const
{
	out << " Young's modulus . . . . . . . . . . . . . . . . = " << fYoung   << '\n';
	out << " Poisson's ratio . . . . . . . . . . . . . . . . = " << fPoisson << '\n';
	out << " Shear modulus . . . . . . . . . . . . . . . . . = " << fMu      << '\n';
	out << " Bulk modulus. . . . . . . . . . . . . . . . . . = " << fKappa   << '\n';
	out << " Lame modulus  . . . . . . . . . . . . . . . . . = " << fLambda  << '\n';	
}

/* information about subordinate parameter lists */
void IsotropicT::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	ParameterInterfaceT::DefineSubs(sub_list);

	/* choice of how to define moduli */
	sub_list.AddSub("modulus_definition_choice", ParameterListT::Once, true);
}

/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* IsotropicT::NewSub(const StringT& name) const
{
	if (name == "modulus_definition_choice")
	{
		ParameterContainerT* choice = new ParameterContainerT(name);
		choice->SetSubSource(this);

		/* set the choices */		
		choice->SetListOrder(ParameterListT::Choice);

		ParameterContainerT E_and_nu("E_and_nu");
		ParameterT E(ParameterT::Double, "Young_modulus");
		E.AddLimit(0.0, LimitT::Lower);
		E_and_nu.AddParameter(E);
		ParameterT Poisson(ParameterT::Double, "Poisson_ratio");
		Poisson.AddLimit(-1.0, LimitT::Lower);
		Poisson.AddLimit( 0.5, LimitT::Upper);
		E_and_nu.AddParameter(Poisson);
		choice->AddSub(E_and_nu);
		
		ParameterContainerT bulk_and_shear("bulk_and_shear");		
		ParameterT kappa(ParameterT::Double, "bulk_modulus");
		kappa.AddLimit(0.0, LimitT::Lower);
		bulk_and_shear.AddParameter(kappa);
		ParameterT mu(ParameterT::Double, "shear_modulus");
		mu.AddLimit(0.0, LimitT::Lower);
		bulk_and_shear.AddParameter(mu);
		choice->AddSub(bulk_and_shear);
	
		return choice;
	}
	else /* inherited */
		return ParameterInterfaceT::NewSub(name);
}

/* accept parameter list */
void IsotropicT::TakeParameterList(const ParameterListT& list)
{
	const ParameterListT* E_and_nu = list.List("E_and_nu");
	if (E_and_nu) {
		double  E = E_and_nu->GetParameter("Young_modulus");
		double nu = E_and_nu->GetParameter("Poisson_ratio");
		Set_E_nu(E, nu);
	}
	else {
		const ParameterListT* bulk_and_shear = list.List("bulk_and_shear");
		if (!bulk_and_shear) ExceptionT::GeneralFail("IsotropicT::TakeParameterList");
		double    mu = bulk_and_shear->GetParameter("shear_modulus");
		double kappa = bulk_and_shear->GetParameter("bulk_modulus");
		Set_mu_kappa(mu, kappa);
	}
}

/*************************************************************************
 * Protected
 *************************************************************************/

/* compute the symetric Cij reduced index matrix */
void IsotropicT::ComputeModuli(dMatrixT& moduli) const
{
	if (moduli.Rows() == 6)
	{
		double mu = Mu();
		double lambda = Lambda();
		moduli = 0.0;
		moduli(2,2) = moduli(1,1) = moduli(0,0) = lambda + 2.0*mu;
		moduli(1,2) = moduli(0,1) = moduli(0,2) = lambda;
		moduli(5,5) = moduli(4,4) = moduli(3,3) = mu;

		/* symmetric */
		moduli.CopySymmetric();
	}
	else
		ExceptionT::SizeMismatch("IsotropicT::ComputeModuli", "3D only");
}

void IsotropicT::ComputeModuli2D(dMatrixT& moduli, 
	SolidMaterialT::ConstraintT constraint) const
{
	if (moduli.Rows() == 3)
	{
		double mu = Mu();
		double lambda = Lambda();

		/* plane stress correction */
		if (constraint == SolidMaterialT::kPlaneStress) {
		
			double lam_2_mu = lambda + 2.0*mu;
			if (fabs(lam_2_mu) < kSmall) {
				if (fabs(mu) > kSmall)
					ExceptionT::GeneralFail("IsotropicT::ComputeModuli2D", "bad plane stress modulus");
			}
			else
				lambda *= 2.0*mu/lam_2_mu;
		}
		
		moduli = 0.0;
		moduli(1,1) = moduli(0,0) = lambda + 2.0*mu;
		moduli(0,1) = moduli(1,0) = lambda;
		moduli(2,2) = mu;
	}
	else 
		throw ExceptionT::kSizeMismatch;
}

void IsotropicT::ComputeModuli1D(dMatrixT& moduli) const
{
	if (moduli.Rows() == 1)
	{
		double young = Young();
		moduli = 0.0;
		moduli(0,0) = young;
	}
	else
		throw ExceptionT::kSizeMismatch;
}

/* scale factor for constrained dilatation */
double IsotropicT::DilatationFactor2D(SolidMaterialT::ConstraintT constraint) const
{
	if (constraint == SolidMaterialT::kPlaneStrain)
		return 1.0 + Poisson();
	else
		return 1.0;
}
