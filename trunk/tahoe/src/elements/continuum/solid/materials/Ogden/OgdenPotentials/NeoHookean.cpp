/* $Id: NeoHookean.cpp,v 1.6 2007-04-09 23:33:26 tdnguye Exp $ */
/* created:   TDN (5/31/2001) */
/* Phi(I1,J) = mu/2*(I1-3)+kappa/4*(J^2-1-2*ln(J)) */
/* I1 = trace(C); J=sqrt(det(C)) */
#include "NeoHookean.h"
#include "ExceptionT.h"

#include <math.h>
#include <iostream.h>


using namespace Tahoe;
const double third = 1.0/3.0;

NeoHookean::NeoHookean(void)
{
	SetName("neo-hookean");
}

/* set parameters */
void NeoHookean::SetKappaMu(double kappa, double mu)
{
	fMu = mu;
	SetKappa(kappa);
}

void NeoHookean::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	PotentialT::DefineParameters(list);

	list.AddParameter(fMu, "mu");
	
	/* set the description */
	list.SetDescription("Psi(Cbar) = 0.3*mu(I1bar-3)");	
}

void NeoHookean::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	PotentialT::TakeParameterList(list);

	fMu = list.GetParameter("mu");

	/* check */
	if (fMu < kSmall) ExceptionT::BadInputValue("NeoHookean::TakeParameterList",
		"expecting a non-negative value mu: %d", fMu);
}

double NeoHookean::Energy(const dArrayT& lambda_bar, const double& J)
{
  double I1 = lambda_bar[0]+lambda_bar[1]+lambda_bar[2];
  double phi = 0.5*fMu*(I1-3);
  phi += MeanEnergy(J);
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


