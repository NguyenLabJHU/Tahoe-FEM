/* $Id: ArrudaBoyce.cpp,v 1.1 2007-07-17 20:11:49 tdnguye Exp $ */
/* created:   TDN (7/2007) */
/* Arruda, E.M and Boyce, M.C., JMPS, v41, pp. 389-412*/
#include "ArrudaBoyce.h"
#include "ExceptionT.h"

#include <math.h>
#include <iostream.h>


using namespace Tahoe;
const double third = 1.0/3.0;

ArrudaBoyce::ArrudaBoyce(void)
{
	SetName("arruda-boyce");
}

/* set parameters */
void ArrudaBoyce::SetKappaMu(double kappa, double mu)
{
	fMu = mu;
	SetKappa(kappa);
}

void ArrudaBoyce::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	PotentialT::DefineParameters(list);

	LimitT zero(0.0, LimitT::Lower);
	LimitT one(1.0, LimitT::Lower);

	ParameterT muN(ParameterT::Double, "network_stiffness");
	ParameterT lambdaL(ParameterT::Double, "locking_stretch");

	muN.AddLimit(zero);
	lambdaL.AddLimit(one);

	list.AddParameter(muN);
	list.AddParameter(lambdaL);
}

void ArrudaBoyce::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	PotentialT::TakeParameterList(list);

	fmuN = list.GetParameter("network_stiffness");
	flambdaL = list.GetParameter("locking_stretch");
	fMu = fmuN*flambdaL*fLangevin.Function(1.0/flambdaL);
}

double ArrudaBoyce::Energy(const dArrayT& lambda_bar, const double& J)
{
  double lam_eff = sqrt(1.0/3.0*(lambda_bar[0]+lambda_bar[1]+lambda_bar[2]));
  double r = lam_eff/flambdaL;
  double r0 = 1.0/flambdaL;
  double x = fLangevin.Function(r);
  double y = fLangevin.Function(r0);
  
  double phi = fmuN*flambdaL*flambdaL*(r*x + log(x/sinh(x)) - 1.0/flambdaL*y - log(y/sinh(y)));
  phi += MeanEnergy(J);

  return(phi);
}
void ArrudaBoyce::DevStress(const dArrayT& lambda_bar,dArrayT& tau)
{
  int nsd = tau.Length();
  
  double lam_eff = sqrt(1.0/3.0*(lambda_bar[0]+lambda_bar[1]+lambda_bar[2]));
  double r = lam_eff/flambdaL;
  
  fMu = fmuN/r*fLangevin.Function(r);
  
  const double& l0 = lambda_bar[0];
  const double& l1 = lambda_bar[1];
  const double& l2 = lambda_bar[2];
  
  tau[0] = fMu*third*(2.0*l0-l1-l2);
  tau[1] = fMu*third*(2.0*l1-l0-l2);
  
  if (nsd == 3)
    tau[2] = fMu*third*(2.0*l2-l0-l1);
}

void ArrudaBoyce::DevMod(const dArrayT& lambda_bar, dSymMatrixT& eigenmodulus)
{

	/*dtau_A/depsilon_B*/
  int nsd = eigenmodulus.Rows();
  double ninth = third*third;
  
  
  double lam_eff = sqrt(1.0/3.0*(lambda_bar[0]+lambda_bar[1]+lambda_bar[2]));
  double r = lam_eff/flambdaL;
  
  fMu = fmuN/r*fLangevin.Function(r);
  double dMu_dleff = fmuN/(lam_eff) * (1.0/lam_eff*fLangevin.DFunction(r) - fLangevin.Function(r)/r);

  const double& l0 = lambda_bar[0];
  const double& l1 = lambda_bar[1];
  const double& l2 = lambda_bar[2];
  
  double coeff = dMu_dleff/lam_eff * third;
  eigenmodulus[0] = 2.0*fMu*ninth*(4.0*l0+l1+l2) + coeff*ninth*(2.0*l0-l1-l2)*(2.0*l0-l1-l2);
  eigenmodulus[1] = 2.0*fMu*ninth*(4.0*l1+l2+l0) + coeff*ninth*(2.0*l1-l0-l2)*(2.0*l1-l0-l2);
  if (nsd == 2)
  {
    eigenmodulus[2] = 2.0*fMu*ninth*(-2.0*l0-2.0*l1+l2) + coeff*ninth*(2.0*l0-l1-l2)*(2.0*l1-l0-l2); 
  }
  else 
  {
    eigenmodulus[2] = 2.0*fMu*ninth*(4.0*l2+l0+l1) + coeff*ninth*(2.0*l2-l0-l1)*(2.0*l2-l0-l1);	
    eigenmodulus[3] = 2.0*fMu*ninth*(-2.0*l1-2.0*l2+l0) + coeff*ninth*(2.0*l1-l0-l2)*(2.0*l2-l0-l1);
    eigenmodulus[4] = 2.0*fMu*ninth*(-2.0*l0-2.0*l2+l1) + coeff*ninth*(2.0*l0-l1-l2)*(2.0*l2-l0-l1);
    eigenmodulus[5] = 2.0*fMu*ninth*(-2.0*l0-2.0*l1+l2) + coeff*ninth*(2.0*l0-l1-l2)*(2.0*l1-l0-l2);
  }
}

/* set parameters */
double ArrudaBoyce::GetMu(void)
{
	return(fMu);
}


