/* $Id: PotentialT.cpp,v 1.1 2007-04-09 23:33:26 tdnguye Exp $ */
/* created:   TDN (5/31/2001) */
/* Phi(I1,J) = mu/2*(I1-3)+kappa/4*(J^2-1-2*ln(J)) */
/* I1 = trace(C); J=sqrt(det(C)) */
#include "PotentialT.h"
#include "ExceptionT.h"

#include <math.h>
#include <iostream.h>


using namespace Tahoe;
const double third = 1.0/3.0;

PotentialT::PotentialT(void):
	ParameterInterfaceT("strain-energy-potential"),
	fKappa(0.0),
	fMu(0.0)
{}

PotentialT::~PotentialT(void){};

/* set parameters */
void PotentialT::SetKappa(double kappa) {fKappa = kappa; }

double PotentialT::MeanEnergy(const double& J) {0.25*fKappa*(J*J-1-2*log(J));}

double PotentialT::MeanStress(const double& J) {return(0.5*fKappa*(J*J-1));}

double PotentialT::MeanMod(const double& J) {return(fKappa*J*J);}

void PotentialT::DefineParameters(ParameterListT& list) const
{
	list.AddParameter(fKappa, "kappa");
	
	/* set the description */
	list.SetDescription("U(J) = 0.25*kappa(J^2-1-2*log(J))");	
}

void PotentialT::TakeParameterList(const ParameterListT& list)
{
	fKappa = list.GetParameter("kappa");

	/* check */
	if (fKappa < -kSmall) ExceptionT::BadInputValue("PotentialT::TakeParameterList",
		"expecting a nonnegative value kappa: %d", fKappa);
}
