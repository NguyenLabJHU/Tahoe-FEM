/* created:   TDN (5/31/2001) */
/* Phi(I1,J) = mu/2*(I1-3)+kappa/4*(J^2-1-2*ln(J)) */
/* I1 = trace(C); J=sqrt(det(C)) */

#include "NeoHookean.h"

#include "ExceptionT.h"
#include "ifstreamT.h"
#include <math.h>
#include <iostream.h>
#include "ifstreamT.h"

using namespace Tahoe;

NeoHookean::NeoHookean(ifstreamT& in): fthird(1.0/3.0)
{
  /*read in material parameters*/
  in >> fMu;
  in >> fKappa;
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
  
  tau[0] = fMu*fthird*(2.0*l0-l1-l2);
  tau[1] = fMu*fthird*(2.0*l1-l0-l2);
  
  if (nsd == 3)
    tau[2] = fMu*fthird*(2.0*l2-l0-l1);
}

double NeoHookean::MeanStress(const double& J) {return(0.5*fKappa*(J*J-1));}

void NeoHookean::DevMod(const dArrayT& lambda_bar, dSymMatrixT& eigenmodulus)
{
  int nsd = eigenmodulus.Rows();
  double ninth = fthird*fthird;
  
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

