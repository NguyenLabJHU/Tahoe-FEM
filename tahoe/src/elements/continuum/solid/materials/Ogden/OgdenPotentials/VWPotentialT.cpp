/* $Id: VWPotentialT.cpp,v 1.2 2009-04-23 03:22:46 tdnguye Exp $ */
#include "VWPotentialT.h"
#include "ExceptionT.h"

#include <math.h>
#include <iostream.h>


using namespace Tahoe;
const double third = 1.0/3.0;

VWPotentialT::VWPotentialT(void):
	falpha(0.0),
	fbeta(0.0),
	fgamma(0.0)
{
	SetName("veronda-westmann");
}

/* set parameters */
void VWPotentialT::SetKappaMu(double kappa, double mu)
{
	fMu = (falpha*fbeta*fbeta-fgamma);
	SetKappa(kappa);
}

void VWPotentialT::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	PotentialT::DefineParameters(list);
		
	list.AddParameter(falpha, "alpha");
	list.AddParameter(fbeta, "beta");
	list.AddParameter(fgamma, "gamma");
	
	/* set the description */
	list.SetDescription("Psi(Cbar) = 0.5*alpha*exp(beta*(I1bar-3))-0.5*gamma(I2bar-3)");	
}

void VWPotentialT::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	PotentialT::TakeParameterList(list);

	falpha = list.GetParameter("alpha");
	fbeta = list.GetParameter("beta");
	fgamma = list.GetParameter("gamma");

	if (falpha < kSmall) ExceptionT::BadInputValue("VWPotentialT::TakeParameterList",
		"expecting a nonnegative value alpha: %d", fMu);
	if (fbeta < -kSmall) ExceptionT::BadInputValue("VWPotentialT::TakeParameterList",
		"expecting a positive value beta: %d", fMu);
	if (fgamma < -kSmall) ExceptionT::BadInputValue("VWPotentialT::TakeParameterList",
		"expecting a positive value gamma: %d", fMu);
}

double VWPotentialT::Energy(const dArrayT& lambda_bar, const double& J,  double temperature)
{
  double I1 = lambda_bar[0]+lambda_bar[1]+lambda_bar[2];
  double I2 = lambda_bar[0]*lambda_bar[1]+lambda_bar[1]*lambda_bar[2]+lambda_bar[0]*lambda_bar[2];

  double phi = 0.5*falpha*exp(fbeta*(I1-3.0) - 1.0);
  phi -= 0.5*fgamma*(I2-3.0);
  phi += MeanEnergy(J);
  return(phi);
}
void VWPotentialT::DevStress(const dArrayT& lambda_bar,dArrayT& tau,  double temperature)
{
  int nsd = tau.Length();
  
  const double& l0 = lambda_bar[0];
  const double& l1 = lambda_bar[1];
  const double& l2 = lambda_bar[2];
	
  double I1 = lambda_bar[0]+lambda_bar[1]+lambda_bar[2];
  double coeff = falpha*fbeta*exp(fbeta*(I1-3));
  
  tau[0] = coeff*third*(2.0*l0-l1-l2);
  tau[1] = coeff*third*(2.0*l1-l0-l2);
  
	/*fgamma(I2bar - 3.0) eigen kirchhoff stress*/
	tau[0] += fgamma*third*(2.0/l0 - 1.0/l1 - 1.0/l2);
	tau[1] += fgamma*third*(2.0/l1 - 1.0/l2 - 1.0/l0);

  if (nsd == 3)
  {
    tau[2] = coeff*third*(2.0*l2-l0-l1);
	tau[2] += fgamma*third*(2.0/l2 - 1.0/l0 - 1.0/l1);
  }
}

void VWPotentialT::DevMod(const dArrayT& lambda_bar, dSymMatrixT& eigenmodulus,  double temperature)
{
  int nsd = eigenmodulus.Rows();
  double ninth = third*third;
  
  const double& l0 = lambda_bar[0];
  const double& l1 = lambda_bar[1];
  const double& l2 = lambda_bar[2];
  
  double I1 = lambda_bar[0]+lambda_bar[1]+lambda_bar[2];
  double coeff = falpha*fbeta*exp(fbeta*(I1-3));

	eigenmodulus[0] = 2.0*coeff*ninth*(4.0*l0 + l1 + l2) 
		+ 2.0*fbeta*coeff*ninth*(2.0*l0 - l1 - l2)*(2.0*l0 - l1 - l2);
		
	eigenmodulus[1] = 2.0*coeff*ninth*(4.0*l1 + l2 + l0) 
		+ 2.0*fbeta*coeff*ninth*(2.0*l1 - l2 - l0)*(2.0*l1 - l2 - l0);
	
	eigenmodulus[0] += -2.0*fgamma*ninth*(4.0/l0 + 1.0/l1 + 1.0/l2);
	eigenmodulus[1] += -2.0*fgamma*ninth*(4.0/l1 + 1.0/l2 + 1.0/l0);

	if (nsd == 2)
	{
		eigenmodulus[2] = 2.0*coeff*ninth*(-2.0*l0 - 2.0*l1 + l2) 
			+ 2.0*fbeta*coeff*ninth*(-l0 + 2.0*l1 - l2)*(2.0*l0 - l1 - l2);
			
		eigenmodulus[2] += -2.0*fgamma*ninth*(-2.0/l0 - 2.0/l1 + 1.0/l2);
	}
	else 
	{
		eigenmodulus[2] = 2.0*coeff*ninth*(4.0*l2 + l0 + l1);
			+ 2.0*fbeta*coeff*ninth*(2.0*l2 - l0 - l1)*(2.0*l2 - l0 - l1);
			
		eigenmodulus[3] = 2.0*coeff*ninth*(-2.0*l1 - 2.0*l2 + l0);
			+ 2.0*fbeta*coeff*ninth*(-l1 + 2.0*l2 - l0)*(2.0*l2 - l1 - l0);
			
		eigenmodulus[4] = 2.0*coeff*ninth*(-2.0*l0 - 2.0*l2 + l1);
			+ 2.0*fbeta*coeff*ninth*(-l0 + 2.0*l2 - l1)*(2.0*l2 - l0 - l1);
			
		eigenmodulus[5] = 2.0*coeff*ninth*(-2.0*l0 - 2.0*l1 + l2);
			+ 2.0*fbeta*coeff*ninth*(-l0 + 2.0*l1 - l2)*(2.0*l0 - l1 - l2);

		eigenmodulus[2] += -2.0*fgamma*ninth*(4.0/l2 + 1.0/l0 + 1.0/l1);
		eigenmodulus[3] += -2.0*fgamma*ninth*(-2.0/l1 - 2.0/l2 + 1.0/l0);
		eigenmodulus[4] += -2.0*fgamma*ninth*(-2.0/l0 - 2.0/l2 + 1.0/l1);
		eigenmodulus[5] += -2.0*fgamma*ninth*(-2.0/l0 - 2.0/l1 + 1.0/l2);
	}
}

