/* $Id: SV_NeoHookean2D.cpp,v 1.4 2002-11-14 17:06:15 paklein Exp $ */
/* created:   TDN (5/31/2001) */
/* Phi(I1,J) = mu/2*(I1-3)+kappa/4*(J^2-1-2*ln(J)) */
/* I1 = trace(C); J=sqrt(det(C)) */

#include "ExceptionT.h"
#include "SV_NeoHookean2D.h"
#include "ifstreamT.h"

#include <math.h>
#include <iostream.h>

using namespace Tahoe;

SV_NeoHookean2D::SV_NeoHookean2D(ifstreamT& in, const FDMatSupportT& support):
         FDSimoVisco2D(in, support),
	 fMu(2),
	 fKappa(2),
	 fCbar(2)
{
         /*read in material parameters*/
         if (fConstraintOption == Material2DT::kPlaneStress)
	 {
	         cout << "Plane Stress formulation is not implemented yet\n";
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
}

void SV_NeoHookean2D::Print(ostream& out) const
{
         FDSimoVisco2D::Print(out);
	 out << "Equilibrium Potential: \n";
	 out<<"      mu = "<<fMu[0]<<'\n';
	 out<<"      kappa = "<<fKappa[0]<<'\n';
	 out << "Non-Equilibrium Potential: \n";
	 out<<"      mu = "<<fMu[1]<<'\n';
	 out<<"      kappa = "<<fKappa[1]<<'\n';
}

void SV_NeoHookean2D::PrintName(ostream& out) const
{
         FDSimoVisco2D::PrintName(out);
         out << "2D Compressible Neo-Hookean Potential\n";
	 out << "    Phi = mu/2*(I1_bar - 3) + kappa/4*(J^2-1-2lnJ)\n";
}

double SV_NeoHookean2D::Phi(const dMatrixT& Fbar, const double& J, 
			    const int SpringType)
{
         dSymMatrixT C(2);
	 C.MultATA(Fbar);

         double I1 = C.Trace()+fCbar33;
	 return(0.5*fMu[SpringType]*(I1-3)+
		0.25*fKappa[SpringType]*(J*J-1-2*log(J)));
}

void SV_NeoHookean2D::Sig_dev(const dMatrixT& Fbar, const double& J, 
			 dSymMatrixT& devsig, const int SpringType)
{
	 fCbar.MultATA(Fbar);
         dSymMatrixT s(2);
	 s.Identity();
	 s *= -fthird*(fCbar.Trace()+fCbar33);

	 devsig.MultAAT(Fbar);
	 devsig += s;
	 
	 devsig *= fMu[SpringType]/J;
}

void SV_NeoHookean2D::Devmod(const dMatrixT& Fbar, const double& J, 
			 dMatrixT& devmod, const int SpringType)
{
	 fCbar.MultATA(Fbar);
         dSymMatrixT sig(2);
	 dSymMatrixT s(2);
	 dMatrixT mod(3);

	 devmod.ReducedIndexI();
	 mod.ReducedIndexII();
	 mod *= fthird;
	 devmod -= mod;
	 devmod *= 2.0*fthird/J*fMu[SpringType]*(fCbar.Trace()+fCbar33);

	 Sig_dev(Fbar, J, sig, SpringType);
	 mod.DyadAB(sig, s.Identity());
	 mod.Symmetrize();
	 mod *= -2.0*fthird;
	 devmod += mod;
}

double SV_NeoHookean2D::dUdJ(const double& J, const int SpringType)
{
	 return(0.5*fKappa[SpringType]*(J-1/J));
}

double SV_NeoHookean2D::ddUddJ(const double& J, const int SpringType)
{
         double iJ = 1.0/J;
         return (0.5*fKappa[SpringType]*(1.0+iJ*iJ));
}

void SV_NeoHookean2D::OutOfPlaneStretch(const dMatrixT& Fbar,const double& J,
				       const int SpringType)
{
#pragma unused(Fbar)
#pragma unused(SpringType)

	double& l2_bar = fCbar33;
	double l2 = 1.0;
	double m = pow(J,2.0*fthird);
	l2_bar = l2/m; 
}	
