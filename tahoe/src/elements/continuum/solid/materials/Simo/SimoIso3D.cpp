/* $Id: SimoIso3D.cpp,v 1.1.1.1 2001-01-29 08:20:25 paklein Exp $ */
/* created: paklein (03/02/1997)                                          */

#include "SimoIso3D.h"
#include <iostream.h>
#include <math.h>
#include "ElasticT.h"

/* constructor */
SimoIso3D::SimoIso3D(ifstreamT& in, const ElasticT& element):
	FDStructMatT(in, element),
	IsotropicT(in),
	fStress(3),
	fModulus(dSymMatrixT::NumValues(3)),
	
	/* work space */
	fb_bar(3),
	fnorm(3),
	frank4(dSymMatrixT::NumValues(3)),
	
	/* fixed forms */
	fIdentity(3),
	fIcrossI(dSymMatrixT::NumValues(3)),
	fIdentity4(dSymMatrixT::NumValues(3)),
	fDevOp4(dSymMatrixT::NumValues(3))
{	
	MuAndKappa(fmu, fkappa);
	
	/* initialize work matricies */
	fIdentity.Identity();
	fIcrossI.Outer(fIdentity, fIdentity);
	fIdentity4.ReducedIndexI();	
	fDevOp4.ReducedIndexDeviatoric();
}

/* print parameters */
void SimoIso3D::Print(ostream& out) const
{
	/* inherited */
	FDStructMatT::Print(out);
	IsotropicT::Print(out);
}

/* modulus */
const dMatrixT& SimoIso3D::c_ijkl(void)
{
	/* compute b_bar */
	const dSymMatrixT& b_3D = b();
	double J = sqrt(b_3D.Det());
	fb_bar.SetToScaled(pow(J,-2.0/3.0),b_3D);

	ComputeModuli(J, fb_bar, fModulus);
	
	return fModulus;
}
	
/* stress */
const dSymMatrixT& SimoIso3D::s_ij(void)
{
	/* compute b_bar */
	const dSymMatrixT& b_3D = b();
	double J = sqrt(b_3D.Det());
	fb_bar.SetToScaled(pow(J,-2.0/3.0),b_3D);

	ComputeCauchy(J, fb_bar, fStress);
	
	return fStress;
}

/* material description */
const dMatrixT& SimoIso3D::C_IJKL(void)
{
	cout << "\n SimoIso3D::C_IJKL: use updated Lagrangian formulation" << endl;
	throw eGeneralFail;

	return fModulus; // dummy
}

const dSymMatrixT& SimoIso3D::S_IJ(void)
{
	cout << "\n SimoIso3D::S_IJ: use updated Lagrangian formulation" << endl;
	throw eGeneralFail;

	return fStress; // dummy
}

/* returns the strain energy density for the specified strain */
double SimoIso3D::StrainEnergyDensity(void)
{
	/* compute b_bar */
	const dSymMatrixT& b_3D = b();
	double J = sqrt(b_3D.Det());
	fb_bar.SetToScaled(pow(J,-2.0/3.0),b_3D);

	return ComputeEnergy(J, fb_bar);
}

/*************************************************************************
* Protected
*************************************************************************/

void SimoIso3D::PrintName(ostream& out) const
{
	/* inherited */
	FDStructMatT::PrintName(out);

	out << "    Simo Isotropic\n";
}

/* computation routines */
void SimoIso3D::ComputeModuli(double J, const dSymMatrixT& b_bar,
	dMatrixT& moduli)
{
	/* initialize */
	moduli = 0.0;

	/* volumetric */
	double du  = dU(J);
	double ddu = ddU(J);
	
	moduli.AddScaled(du + J*ddu, fIcrossI);
	moduli.AddScaled(-2.0*du,fIdentity4);
	
	/* deviatoric */
	double mu_bar = fmu*b_bar.Trace()/(J*3.0);
	moduli.AddScaled(2.0*mu_bar, fDevOp4);

	fStress.SetToScaled(fmu,b_bar);
	fStress.Deviatoric();
	
	frank4.Outer(fStress,fIdentity);
	frank4.Symmetrize();
	moduli.AddScaled(-4.0/(J*3.0), frank4);
}

void SimoIso3D::ComputeCauchy(double J, const dSymMatrixT& b_bar,
	dSymMatrixT& cauchy)
{
	/* deviatoric */
	cauchy.SetToScaled(fmu/J,b_bar);
	cauchy.Deviatoric();
	
	/* volumetric */
	cauchy.PlusIdentity(dU(J));
}

double SimoIso3D::ComputeEnergy(double J, const dSymMatrixT& b_bar)
{
	return U(J) +                         /* volumetric */
	       0.5*fmu*(b_bar.Trace() - 3.0); /* deviatoric */
}

/* Volumetric energy function and derivatives */
double SimoIso3D::U(double J) const
{
	return 0.5*fkappa*(0.5*(J*J - 1.0) - log(J));
}

double SimoIso3D::dU(double J) const
{
	return 0.5*fkappa*(J - 1.0/J);
}

double SimoIso3D::ddU(double J) const
{
	return 0.5*fkappa*(1.0 + 1.0/(J*J));
}
