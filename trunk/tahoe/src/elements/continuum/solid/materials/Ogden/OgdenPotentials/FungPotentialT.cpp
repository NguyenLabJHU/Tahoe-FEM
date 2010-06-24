/* $Id: FungPotentialT.cpp,v 1.2 2010-06-24 13:26:03 tdnguye Exp $ */
#include "FungPotentialT.h"
#include "ExceptionT.h"

#include <math.h>
#include <iostream.h>


using namespace Tahoe;
const double third = 1.0/3.0;
const double kBig = 1e+12;

FungPotentialT::FungPotentialT(void):
	falpha(0.0),
	fbeta(0.0),
	fmu(0.0)
{
	SetName("fung-potential");
}

/* set parameters */
void FungPotentialT::SetKappaMu(double kappa, double mu)
{
	fMu = fmu;
	SetKappa(kappa);
}

void FungPotentialT::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	PotentialT::DefineParameters(list);
		
	list.AddParameter(falpha, "alpha");
	list.AddParameter(fbeta, "beta");
	list.AddParameter(fmu, "mu");
	
	/* set the description */
	list.SetDescription("Psi(Cbar) = 0.5*mu*(I1bar - 3) + 0.5*alpha/beta*exp(0.5*beta*(I1bar-3)^2)");	
}

void FungPotentialT::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	PotentialT::TakeParameterList(list);

	falpha = list.GetParameter("alpha");
	fbeta = list.GetParameter("beta");
	fmu = list.GetParameter("mu");

	if (fmu < kSmall) ExceptionT::BadInputValue("FungPotentialT::TakeParameterList",
		"expecting a nonnegative value mu: %d", fMu);
	if (fbeta < -kSmall) ExceptionT::BadInputValue("FungPotentialT::TakeParameterList",
		"expecting a positive value beta: %d", fbeta);
	if (falpha < -kSmall) ExceptionT::BadInputValue("FungPotentialT::TakeParameterList",
		"expecting a positive value alpha: %d", falpha);
}

double FungPotentialT::Energy(const dArrayT& lambda_bar, const double& J,  double temperature)
{
  double I1 = lambda_bar[0]+lambda_bar[1]+lambda_bar[2];

  double phi = 0.5*falpha/fbeta*(exp(0.5*fbeta*(I1-3.0)*(I1-3.0))-1.0) + 0.5*fmu*(I1-3.0);
  phi += MeanEnergy(J);
  return(phi);
}
void FungPotentialT::DevStress(const dArrayT& lambda_bar,dArrayT& tau,  double temperature)
{
  int nsd = tau.Length();
  
  const double& l0 = lambda_bar[0];
  const double& l1 = lambda_bar[1];
  const double& l2 = lambda_bar[2];
	
  double I1 = lambda_bar[0]+lambda_bar[1]+lambda_bar[2];
  double coeff = falpha*(I1-3.0)*exp(0.5*fbeta*(I1-3.0)*(I1-3.0)) + fmu;
  if (coeff > kBig) ExceptionT::GeneralFail("FungPotentialT::DevStress",
		"Infinite stress.");
  
  tau[0] = coeff*third*(2.0*l0-l1-l2);
  tau[1] = coeff*third*(2.0*l1-l0-l2);
  
  if (nsd == 3)
  {
    tau[2] = coeff*third*(2.0*l2-l0-l1);
  }
}

/*dtau/dep*/
void FungPotentialT::DevMod(const dArrayT& lambda_bar, dSymMatrixT& eigenmodulus,  double temperature)
{
  int nsd = eigenmodulus.Rows();
  double ninth = third*third;
  
  const double& l0 = lambda_bar[0];
  const double& l1 = lambda_bar[1];
  const double& l2 = lambda_bar[2];
  
  double I1 = lambda_bar[0]+lambda_bar[1]+lambda_bar[2];
  double coeff = falpha*exp(0.5*fbeta*(I1-3.0)*(I1-3.0));
  if (coeff > kBig) ExceptionT::GeneralFail("VWPotentialT::DevMod",
		"Infinite modulus.");

//  cout << "\ncoef: "<<coeff;
	eigenmodulus[0] = 2.0*ninth*(coeff*(I1-3.0)+fmu)*(4.0*l0 + l1 + l2) 
		+ 2.0*ninth*coeff*(1.0+fbeta*(I1-3.0)*(I1-3.0))*(2.0*l0 - l1 - l2)*(2.0*l0 - l1 - l2);
		
	eigenmodulus[1] =2.0*ninth*(coeff*(I1-3.0)+fmu)*(4.0*l1 + l2 + l0) 
		+ 2.0*ninth*coeff*(1.0+fbeta*(I1-3.0)*(I1-3.0))*(2.0*l1 - l2 - l0)*(2.0*l1 - l2 - l0);

	if (nsd == 2)
	{
		eigenmodulus[2] = 2.0*ninth*(coeff*(I1-3.0)+fmu)*(-2.0*l0 - 2.0*l1 + l2) 
			+  2.0*ninth*coeff*(1.0+fbeta*(I1-3.0)*(I1-3.0))*(-l0 + 2.0*l1 - l2)*(2.0*l0 - l1 - l2);
	}
	else 
	{
		eigenmodulus[2] = 2.0*ninth*(coeff*(I1-3.0)+fmu)*(4.0*l2 + l0 + l1)
			 +2.0*ninth*coeff*(1.0+fbeta*(I1-3.0)*(I1-3.0))*(2.0*l2 - l0 - l1)*(2.0*l2 - l0 - l1);
			
		eigenmodulus[3] = 2.0*ninth*(coeff*(I1-3.0)+fmu)*(-2.0*l1 - 2.0*l2 + l0)
			+ 2.0*ninth*coeff*(1.0+fbeta*(I1-3.0)*(I1-3.0))*(-l2 + 2.0*l1 - l0)*(2.0*l2 - l1 - l0);
			
		eigenmodulus[4] = 2.0*ninth*(coeff*(I1-3.0)+fmu)*(-2.0*l0 - 2.0*l2 + l1)
			+ 2.0*ninth*coeff*(1.0+fbeta*(I1-3.0)*(I1-3.0))*(-l2 + 2.0*l0 - l1)*(2.0*l2 - l0 - l1);
			
		eigenmodulus[5] = 2.0*ninth*(coeff*(I1-3.0)+fmu)*(-2.0*l0 - 2.0*l1 + l2)
			+ 2.0*ninth*coeff*(1.0+fbeta*(I1-3.0)*(I1-3.0))*(-l0 + 2.0*l1 - l2)*(2.0*l0 - l1 - l2);

	}
}

