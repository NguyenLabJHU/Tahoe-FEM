/* $Id: NeoHookean.cpp,v 1.5 2004-07-15 08:27:26 paklein Exp $ */
/* created:   TDN (5/31/2001) */
/* Phi(I1,J) = mu/2*(I1-3)+kappa/4*(J^2-1-2*ln(J)) */
/* I1 = trace(C); J=sqrt(det(C)) */
#include "NeoHookean.h"
#include "ExceptionT.h"

#include <math.h>
#include <iostream.h>


using namespace Tahoe;
const double third = 1.0/3.0;

NeoHookean::NeoHookean(void):
	fMu(0.0),
	fKappa(0.0)
{

}

/* set parameters */
void NeoHookean::SetKappaMu(double kappa, double mu)
{
	fMu = mu;
	fKappa = kappa;
}

void NeoHookean::Print(ostream& out) const
{
  out<<"      Shear Modulus = "<<fMu<<'\n';
  out<<"      Bulk Modulus = "<<fKappa<<'\n';
}

void NeoHookean::PrintName(ostream& out) const
{
  out << "Compressible Neo-Hookean Potential\n";
  out<<"        Phi = mu/2*(I1_bar - 3) + kappa/4*(J^2-1-2lnJ)\n";
}
void NeoHookean::Initialize(void)
{}

double NeoHookean::Energy(const dArrayT& lambda_bar, const double& J)
{
  double I1 = lambda_bar[0]+lambda_bar[1]+lambda_bar[2];
  double phi = 0.5*fMu*(I1-3)+0.25*fKappa*(J*J-1-2*log(J));
  return(phi);
}
void NeoHookean::DevStress(const dArrayT& lambda_bar,dArrayT& tau)
{
  int nsd = tau.Length();
  
  const double& l0 = lambda_bar[0];
  const double& l1 = lambda_bar[1];
  const double& l2 = lambda_bar[2];
  
  tau[0] = fMu*third*(2.0*l0-l1-l2);
  tau[1] = fMu*third*(2.0*l1-l0-l2);
  
  if (nsd == 3)
    tau[2] = fMu*third*(2.0*l2-l0-l1);
}

double NeoHookean::MeanStress(const double& J) {return(0.5*fKappa*(J*J-1));}

void NeoHookean::DevMod(const dArrayT& lambda_bar, dSymMatrixT& eigenmodulus)
{
  int nsd = eigenmodulus.Rows();
  double ninth = third*third;
  
  const double& l0 = lambda_bar[0];
  const double& l1 = lambda_bar[1];
  const double& l2 = lambda_bar[2];
  
  eigenmodulus[0] = 2.0*fMu*ninth*(4.0*l0+l1+l2);
  eigenmodulus[1] = 2.0*fMu*ninth*(4.0*l1+l2+l0);
  if (nsd == 2)
  {
    eigenmodulus[2] = 2.0*fMu*ninth*(-2.0*l0-2.0*l1+l2); 
  }
  else 
  {
    eigenmodulus[2] = 2.0*fMu*ninth*(4.0*l2+l0+l1);	
    eigenmodulus[3] = 2.0*fMu*ninth*(-2.0*l1-2.0*l2+l0);
    eigenmodulus[4] = 2.0*fMu*ninth*(-2.0*l0-2.0*l2+l1);
    eigenmodulus[5] = 2.0*fMu*ninth*(-2.0*l0-2.0*l1+l2);
  }
}

double NeoHookean::MeanMod(const double& J) {return(fKappa*J*J);}

