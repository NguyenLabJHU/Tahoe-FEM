/* $Id: SV_NeoHookean3D.cpp,v 1.1 2003-03-19 19:03:20 thao Exp $ */
/* created:   TDN (5/31/2001) */
/* Phi(I1,J) = mu/2*(I1-3)+kappa/4*(J^2-1-2*ln(J)) */
/* I1 = trace(C); J=sqrt(det(C)) */

#include "ExceptionT.h"
#include "SV_NeoHookean3D.h"
#include "ifstreamT.h"

#include <math.h>
#include <iostream.h>

using namespace Tahoe;

SV_NeoHookean3D::SV_NeoHookean3D(ifstreamT& in, const FSMatSupportT& support):
         FDSimoVisco3D(in, support),
	 fCbar(3),
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

void SV_NeoHookean3D::Print(ostream& out) const
{
         FDSimoVisco3D::Print(out);
	 out << "Equilibrium Potential: \n";
	 out<<"      mu = "<<fMu[0]<<'\n';
	 out<<"      kappa = "<<fKappa[0]<<'\n';
	 out << "Non-Equilibrium Potential: \n";
	 out<<"      mu = "<<fMu[1]<<'\n';
	 out<<"      kappa = "<<fKappa[1]<<'\n';
}

void SV_NeoHookean3D::PrintName(ostream& out) const
{
         FDSimoVisco3D::PrintName(out);
         out << "3D Compressible Neo-Hookean Potential\n";
	 out << "    Phi = mu/2*(I1_bar - 3) + kappa/4*(J^2-1-2lnJ)\n";
}

double SV_NeoHookean3D::Phi(const dMatrixT& Fbar, const double& J, 
			    const int SpringType)
{
         dSymMatrixT C(3);
	 C.MultATA(Fbar);

         double I1 = C.Trace();
	 return(0.5*fMu[SpringType]*(I1-3)+
		0.25*fKappa[SpringType]*(J*J-1-2*log(J)));
}

void SV_NeoHookean3D::Sig_dev(const dMatrixT& Fbar, const double& J, 
			 dSymMatrixT& devsig, const int SpringType)
{
	 fCbar.MultATA(Fbar);
         dSymMatrixT s(3);
	 s.Identity();
	 s *= -fthird*(fCbar.Trace());

	 devsig.MultAAT(Fbar);
	 devsig += s;
	 
	 devsig *= fMu[SpringType]/J;
}

void SV_NeoHookean3D::Devmod(const dMatrixT& Fbar, const double& J, 
			 dMatrixT& devmod, const int SpringType)
{
	 fCbar.MultATA(Fbar);
         dSymMatrixT sig(3);
	 dSymMatrixT s(3);
	 dMatrixT mod(6);

	 devmod.ReducedIndexI();
	 mod.ReducedIndexII();
	 mod *= fthird;
	 devmod -= mod;
	 devmod *= 2.0*fthird/J*fMu[SpringType]*fCbar.Trace();

	 Sig_dev(Fbar, J, sig, SpringType);
	 mod.DyadAB(sig, s.Identity());
	 mod.Symmetrize();
	 mod *= -2.0*fthird;
	 devmod += mod;
}

double SV_NeoHookean3D::dUdJ(const double& J, const int SpringType)
{
	 return(0.5*fKappa[SpringType]*(J-1/J));
}

double SV_NeoHookean3D::ddUddJ(const double& J, const int SpringType)
{
         double iJ = 1.0/J;
         return (0.5*fKappa[SpringType]*(1.0+iJ*iJ));
}
