/* $Id: RG_NeoHookean3D.cpp,v 1.3 2002-11-14 17:06:10 paklein Exp $ */
/* created:   TDN (5/31/2001) */
/* Phi(I1,J) = mu/2*(I1-3)+kappa/4*(J^2-1-2*ln(J)) */
/* I1 = trace(C); J=sqrt(det(C)) */

#include "fstreamT.h"
#include "ExceptionT.h"
#include "RG_NeoHookean3D.h"
#include <math.h>
#include <iostream.h>

using namespace Tahoe;

RG_NeoHookean3D::RG_NeoHookean3D(ifstreamT& in, const FDMatSupportT& support):
	RG_VDSplit3D(in, support),
	fMu(2),
	fKappa(2)
{
         /*read in material parameters*/
  
         double& mu_EQ = fMu[kEquilibrium];
	 double& kappa_EQ = fKappa[kEquilibrium];
	 double& mu_NEQ = fMu[kNonEquilibrium];
	 double& kappa_NEQ = fKappa[kNonEquilibrium];

         in >> mu_EQ;
	 in >> kappa_EQ;
	 in >> mu_NEQ;
	 in >> kappa_NEQ;
}

void RG_NeoHookean3D::Print(ostream& out) const
{
         RG_VDSplit3D::Print(out);
	 out << "Equilibrium Potential: \n";
	 out<<"      mu = "<<fMu[0]<<'\n';
	 out<<"      kappa = "<<fKappa[0]<<'\n';
	 out << "Non-Equilibrium Potential: \n";
	 out<<"      mu = "<<fMu[1]<<'\n';
	 out<<"      kappa = "<<fKappa[1]<<'\n';
}

void RG_NeoHookean3D::PrintName(ostream& out) const
{
         RG_VDSplit3D::PrintName(out);
         out << "Compressible Neo-Hookean Potential\n";
	 out<<"Phi = mu/2*(I1_bar - 3) + kappa/4*(J^2-1-2lnJ)\n";
}
void RG_NeoHookean3D::Initialize(void)
{
         double mu = fMu[0]+fMu[1];
	 double kappa = fKappa[0]+fKappa[1];
	 IsotropicT::Set_mu_kappa(mu, kappa);
}

double RG_NeoHookean3D::Phi(const dArrayT& eigenstretch_bar, const double& J,
			    const int SpringType)
{
        double I1 = eigenstretch_bar.Sum();
	double phi = 0.5*fMu[SpringType]*(I1-3)+
	             0.25*fKappa[SpringType]*(J*J-1-2*log(J));
	return(phi);
}
void RG_NeoHookean3D::devTau(const dArrayT& eigenstretch_bar,
			   dArrayT& eigenstress, const int SpringType)
{
	double& mu = fMu[SpringType];

	double& l0 = eigenstretch_bar[0];
	double& l1 = eigenstretch_bar[1];
	double& l2 = eigenstretch_bar[2];

	eigenstress[0] = mu*fthird*(2.0*l0-l1-l2);
	eigenstress[1] = mu*fthird*(2.0*l1-l2-l0);
	eigenstress[2] = mu*fthird*(2.0*l2-l0-l1);
}

double RG_NeoHookean3D::meanTau(const double& J, const int SpringType)
{
       return(0.5*fKappa[SpringType]*(J*J-1));
}

void RG_NeoHookean3D::DdevTauDepsilon(const dArrayT& eigenstretch_bar, 
				      dSymMatrixT& eigenmodulus, 
				      const int SpringType)
{
	double& mu = fMu[SpringType];
	double ninth = fthird*fthird;

	double& l0 = eigenstretch_bar[0];
	double& l1 = eigenstretch_bar[1];
	double& l2 = eigenstretch_bar[2];

	eigenmodulus(0,0) = 2.0*mu*ninth*(4.0*l0+l1+l2);
	eigenmodulus(1,1) = 2.0*mu*ninth*(4.0*l1+l2+l0);
	eigenmodulus(2,2) = 2.0*mu*ninth*(4.0*l2+l0+l1);
	
	eigenmodulus(1,2) =2.0*mu*fthird*fthird*(-2.0*l1-2.0*l2+l0);
	eigenmodulus(0,2) =2.0*mu*fthird*fthird*(-2.0*l0-2.0*l2+l1);
	eigenmodulus(0,1) =2.0*mu*fthird*fthird*(-2.0*l0-2.0*l1+l2);
}

double RG_NeoHookean3D::DmeanTauDepsilon(const double& J,const int SpringType) 
{
        return(fKappa[SpringType]*J*J);
}

