/* $Id: VerondaWestmannT.cpp,v 1.2 2007-04-12 16:32:43 thao Exp $ */
/* created: paklein (02/17/2001) */
#include "VerondaWestmannT.h"

#include <math.h>

using namespace Tahoe;

/* constructor */
VerondaWestmannT::VerondaWestmannT(void):
	ParameterInterfaceT("veronda_westmann_potential"),
	fIk(3),
	fdWdIk(3),
	fddWddIk(3),
	fMat(3)
{

}

/* strain energy density */
double VerondaWestmannT::StrainEnergyDensity(void)
{
	/* stretch */
	Compute_Ik(fIk);
	
	const double& I1 = fIk[0];
	const double& I2 = fIk[1];
	const double& I3 = fIk[2];
	
	double energy= falpha*(exp(fbeta*(I1-3.0))-1.0) - 0.5*falpha*fbeta*(I2-3.0);
	energy += 0.25*fgamma*(I3-log(I3)-1.0);

	return(energy);
}

/***********************************************************************
 * Protected
 ***********************************************************************/

/* principal values given principal values of the stretch tensors,
 * i.e., the principal stretches squared */
void VerondaWestmannT::dWdE(const dArrayT& eigenstretch2, dArrayT& eigenstress)
{
	/* compute invariants and derivatives of strain energy density wrt invariants */
	Compute_Ik(eigenstretch2, fIk);
	Compute_dWdIk(fIk, fdWdIk);
	
	const double l1 = eigenstretch2[0];
	const double l2 = eigenstretch2[1];
	const double l3 = eigenstretch2[2];
	
	const double& I1 = fIk[0];
	const double& I2 = fIk[1];
	const double& I3 = fIk[2];

	const double& dWdI1 = fdWdIk[0];
	const double& dWdI2 = fdWdIk[1];
	const double& dWdI3 = fdWdIk[2];

	/* stress */
	eigenstress[0] = 2.0*(dWdI1 + I1*dWdI2 - l1*dWdI2 + 1.0/l1*I3*dWdI3);
	eigenstress[1] = 2.0*(dWdI1 + I1*dWdI2 - l2*dWdI2 + 1.0/l2*I3*dWdI3);
	eigenstress[2] = 2.0*(dWdI1 + I1*dWdI2 - l3*dWdI2 + 1.0/l3*I3*dWdI3);
}

void VerondaWestmannT::ddWddE(const dArrayT& eigenstretch2, dArrayT& eigenstress,
	dSymMatrixT& eigenmod)
{
	Compute_Ik(eigenstretch2, fIk);
	Compute_ddWddIk(fIk, fdWdIk, fddWddIk);
	
	const double l1 = eigenstretch2[0];
	const double l2 = eigenstretch2[1];
	const double l3 = eigenstretch2[2];
	
	const double& I1 = fIk[0];
	const double& I2 = fIk[1];
	const double& I3 = fIk[2];

	const double& dWdI1 = fdWdIk[0];
	const double& dWdI2 = fdWdIk[1];
	const double& dWdI3 = fdWdIk[2];

	const double& ddWddI1 = fddWddIk[0];
	const double& ddWddI2 = fddWddIk[1];
	const double& ddWddI3 = fddWddIk[2];
	const double& ddWddI23 = fddWddIk[3];
	const double& ddWddI13 = fddWddIk[4];
	const double& ddWddI12 = fddWddIk[5];

	/* stress */
	eigenstress[0] = 2.0*(dWdI1 + I1*dWdI2 - l1*dWdI2 + 1.0/l1*I3*dWdI3);
	eigenstress[1] = 2.0*(dWdI1 + I1*dWdI2 - l2*dWdI2 + 1.0/l2*I3*dWdI3);
	eigenstress[2] = 2.0*(dWdI1 + I1*dWdI2 - l3*dWdI2 + 1.0/l3*I3*dWdI3);
	
	/* moduli */
	eigenmod(0,0) = 4.0*(ddWddI1 + ddWddI12*(I1-2.0*l1) + ddWddI13*2.0/l1*I3 +
						 ddWddI2*(I1-l1)*(I1-l1) + ddWddI23*I3*2.0/l1*(I1-l1) +
						 ddWddI3*I3*I3/l1/l1 );
	eigenmod(1,1) = 4.0*(ddWddI1 + ddWddI12*(I1-2.0*l2) + ddWddI13*2.0/l2*I3 +
						ddWddI2*(I1-l2)*(I1-l2) + ddWddI23*I3*2.0/l2*(I1-l2) +
						 ddWddI3*I3*I3/l2/l2 );
	eigenmod(2,2) = 4.0*(ddWddI1 + ddWddI12*(I1-2.0*l3) + ddWddI13*2.0/l3*I3 +
						 ddWddI2*(I1-l3)*(I1-l3) + ddWddI23*I3*2.0/l3*(I1-l3) +
						 ddWddI3*I3*I3/l3/l3 );

	eigenmod(0,1) = 4.0*(ddWddI1 + ddWddI12*(I1-l1-l2) + ddWddI13*(1.0/l1 + 1.0/l2)*I3 +
						 ddWddI2*(I1-l1)*(I1-l2) + ddWddI23*I3*(1.0/l1*(I1-l2) + 1.0/l2*(I1-l1)) +
						 ddWddI3*I3*I3/l1/l2  + dWdI2 + dWdI3*I3/l1/l2 );
	eigenmod(0,2) = 4.0*(ddWddI1 + ddWddI12*(I1-l1-l3) + ddWddI13*(1.0/l1 + 1.0/l3)*I3 +
						 ddWddI2*(I1-l1)*(I1-l3) + ddWddI23*I3*(1.0/l1*(I1-l3) + 1.0/l3*(I1-l1)) +
						 ddWddI3*I3*I3/l1/l3  + dWdI2 + dWdI3*I3/l1/l3 );
	eigenmod(1,2) = 4.0*(ddWddI1 + ddWddI12*(I1-l2-l3) + ddWddI13*(1.0/l2 + 1.0/l3)*I3 +
						 ddWddI2*(I1-l2)*(I1-l3) + ddWddI23*I3*(1.0/l2*(I1-l3) + 1.0/l3*(I1-l2)) +
						 ddWddI3*I3*I3/l2/l3  + dWdI2 + dWdI3*I3/l2/l3 );
}

/*protected*/
/* describe the parameters needed by the interface */
void VerondaWestmannT::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	OgdenIsotropicT::DefineParameters(list);
	list.AddParameter(falpha, "alpha");
	list.AddParameter(fbeta, "beta");
	list.AddParameter(fgamma, "gamma");
}

/* accept parameter list */
void VerondaWestmannT::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	OgdenIsotropicT::TakeParameterList(list);
	falpha = list.GetParameter("alpha");
	fbeta = list.GetParameter("beta");
	fgamma = list.GetParameter("gamma");
}

/*private*/
void VerondaWestmannT::Compute_dWdIk(const dArrayT& Ik, dArrayT& dWdIk)
{
	const double& I1 = Ik[0];
	const double& I2 = Ik[1];
	const double& I3 = Ik[2];
	
	dWdIk[0] = falpha*fbeta*exp(fbeta*(I1-3.0));
	dWdIk[1] = -0.5*falpha*fbeta;
	dWdIk[2] = 0.25*fgamma*(1.0-1.0/I3);
}

void VerondaWestmannT::Compute_ddWddIk(const dArrayT& Ik, dArrayT& dWdIk, dSymMatrixT& ddWddIk)
{
	const double& I1 = Ik[0];
	const double& I2 = Ik[1];
	const double& I3 = Ik[2];

	dWdIk[0] = falpha*fbeta*exp(fbeta*(I1-3.0));
	dWdIk[1] = -0.5*falpha*fbeta;
	dWdIk[2] = 0.25*fgamma*(1.0-1.0/I3);

	ddWddIk[0] = falpha*fbeta*fbeta*exp(fbeta*(I1-3.0));
	ddWddIk[1] = 0.0;
	ddWddIk[2] = 0.25*fgamma/(I3*I3);
	ddWddIk[3] = 0.0;
	ddWddIk[4] = 0.0;
	ddWddIk[5] = 0.0;
}

void VerondaWestmannT::Compute_Ik(dArrayT& Ik)
{
	Compute_C(fC);

	/*C^2*/
	fMat.MultAB(fC, fC);

	/*compute invariants*/
	Ik[0] = fC.Trace();
	Ik[1] = 0.5*(Ik[0]*Ik[0]-fMat.Trace());
	Ik[2] = fC.Det();
}

void VerondaWestmannT::Compute_Ik(const dArrayT& eigenstretch2, dArrayT& Ik)
{
	const double& l1 = eigenstretch2[0];
	const double& l2 = eigenstretch2[1];
	const double& l3 = eigenstretch2[2];

	/*compute invariants*/
	Ik[0] = l1+l2+l3;
	Ik[1] = l1*l2 + l1*l3 + l2*l3;
	Ik[2] = l1*l2*l3;
}
