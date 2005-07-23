/* $Id: VWPotentialT.cpp,v 1.1 2005-07-14 19:44:38 regueiro Exp $ */
#include "VWPotentialT.h"
#include "ExceptionT.h"

#include <math.h>
#include <iostream.h>


using namespace Tahoe;
const double third = 1.0/3.0;

VWPotentialT::VWPotentialT(void):
	fMu(0.0),
	fKappa(0.0)
{

}

/* set parameters */
void VWPotentialT::SetKappaMu(double kappa, double mu)
{
	fMu = mu;
	fKappa = kappa;
}

void VWPotentialT::Print(ostream& out) const
{
  out<<"      Shear Modulus = "<<fMu<<'\n';
  out<<"      Bulk Modulus = "<<fKappa<<'\n';
}

void VWPotentialT::PrintName(ostream& out) const
{
  out << "Compressible Neo-Hookean Potential\n";
  out<<"        Phi = mu/2*(I1_bar - 3) + kappa/4*(J^2-1-2lnJ)\n";
}
void VWPotentialT::Initialize(void)
{}

double VWPotentialT::Energy(const dArrayT& lambda_bar, const double& J)
{
  double I1 = lambda_bar[0]+lambda_bar[1]+lambda_bar[2];
  double phi = 0.5*fMu*(I1-3)+0.25*fKappa*(J*J-1-2*log(J));
  return(phi);
}
void VWPotentialT::DevStress(const dArrayT& lambda_bar,dArrayT& tau)
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

double VWPotentialT::MeanStress(const double& J) {return(0.5*fKappa*(J*J-1));}

void VWPotentialT::DevMod(const dArrayT& lambda_bar, dSymMatrixT& eigenmodulus)
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

double VWPotentialT::MeanMod(const double& J) {return(fKappa*J*J);}

