/* $Id: RG_NeoHookean2D.cpp,v 1.5 2003-01-29 07:34:46 paklein Exp $ */
/* created:   TDN (5/31/2001) */
/* Phi(I1,J) = mu/2*(I1-3)+kappa/4*(J^2-1-2*ln(J)) */
/* I1 = trace(C); J=sqrt(det(C)) */

#include "fstreamT.h"
#include "ExceptionT.h"
#include "RG_NeoHookean2D.h"
#include <math.h>
#include <iostream.h>

using namespace Tahoe;

RG_NeoHookean2D::RG_NeoHookean2D(ifstreamT& in, const FSMatSupportT& support):
          RG_VDSplit2D(in, support),
	  fMu(2),
	  fKappa(2)
{
         /*read in material parameters*/
  
        if (fConstraintOption == Material2DT::kPlaneStress)
	{
	        cout <<"Plane stress formulation is not implemented.\n";
		throw ExceptionT::kBadInputValue;
	}

         double& mu_EQ = fMu[kEquilibrium];
	 double& kappa_EQ = fKappa[kEquilibrium];
	 double& mu_NEQ = fMu[kNonEquilibrium];
	 double& kappa_NEQ = fKappa[kNonEquilibrium];

         in >> mu_EQ;
	 in >> kappa_EQ;
	 in >> mu_NEQ;
	 in >> kappa_NEQ;

	 /*         fMu[kEquilibrium] = mu_EQ;
	 fKappa[kEquilibrium] = kappa_EQ;
	 fMu[kNonEquilibrium] = mu_NEQ; 
	 fKappa[kNonEquilibrium] = mu_NEQ; */
}

void RG_NeoHookean2D::Print(ostream& out) const
{
         RG_VDSplit2D::Print(out);
	 out << "Equilibrium Potential: \n";
	 out<<"      mu = "<<fMu[0]<<'\n';
	 out<<"      kappa = "<<fKappa[0]<<'\n';
	 out << "Non-Equilibrium Potential: \n";
	 out<<"      mu = "<<fMu[1]<<'\n';
	 out<<"      kappa = "<<fKappa[1]<<'\n';
}

void RG_NeoHookean2D::PrintName(ostream& out) const
{
         RG_VDSplit2D::PrintName(out);

         out << "2D Compressible Neo-Hookean Potential\n";
	 out<<"     Phi = mu/2*(I1_bar - 3) + kappa/4*(J^2-1-2lnJ)\n";
}

void RG_NeoHookean2D::Initialize(void)
{
         double mu = fMu[0]+fMu[1];
	 double kappa = fKappa[0]+fKappa[1];
	 IsotropicT::Set_mu_kappa(mu, kappa);
}

double RG_NeoHookean2D::Phi(const dArrayT& eigenstretch_bar, const double& J,
			    const int SpringType)
{
        double I1 = eigenstretch_bar.Sum();
	I1 += fl2_bar;
	double phi = 0.5*fMu[SpringType]*(I1-3)+
	             0.25*fKappa[SpringType]*(J*J-1-2*log(J));
	return(phi);
}
void RG_NeoHookean2D::devTau(const dArrayT& eigenstretch_bar,
			   dArrayT& eigenstress, const int SpringType)
{
	double& mu = fMu[SpringType];

	double& l0 = eigenstretch_bar[0];
	double& l1 = eigenstretch_bar[1];
	double& l2 = fl2_bar;

	eigenstress[0] = mu*fthird*(2.0*l0-l1-l2);
	eigenstress[1] = mu*fthird*(2.0*l1-l2-l0);
}

double RG_NeoHookean2D::meanTau(const double& J, const int SpringType)
{
       return(0.5*fKappa[SpringType]*(J*J-1));
}

void RG_NeoHookean2D::DdevTauDepsilon(const dArrayT& eigenstretch_bar, 
				      dSymMatrixT& eigenmodulus, 
				      const int SpringType)
{
	double& mu = fMu[SpringType];
	double ninth = fthird*fthird;

	double& l0 = eigenstretch_bar[0];
	double& l1 = eigenstretch_bar[1];
	double& l2 = fl2_bar;

	eigenmodulus[0] = 2.0*mu*ninth*(4.0*l0+l1+l2);
	eigenmodulus[1] = 2.0*mu*ninth*(4.0*l1+l2+l0);
	
	eigenmodulus[2] =2.0*mu*ninth*(-2.0*l0-2.0*l1+l2);
}

double RG_NeoHookean2D::DmeanTauDepsilon(const double& J,const int SpringType) 
{
        return(fKappa[SpringType]*J*J);
}

void RG_NeoHookean2D::OutOfPlaneStretch(const dArrayT& eigenstretch_bar, const double& J, const int SpringType)
{
#pragma unused(eigenstretch_bar)
#pragma unused(SpringType)

	double& l2_bar = fl2_bar;
	double l2 = 1.0;
	double m = pow(J,2.0*fthird);
	l2_bar = l2/m;	
}

