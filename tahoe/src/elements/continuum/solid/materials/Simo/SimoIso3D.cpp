/* $Id: SimoIso3D.cpp,v 1.9.46.2 2004-06-09 23:17:41 paklein Exp $ */
/* created: paklein (03/02/1997) */
#include "SimoIso3D.h"
#include <iostream.h>
#include <math.h>

using namespace Tahoe;

/* constructor */
SimoIso3D::SimoIso3D(ifstreamT& in, const FSMatSupportT& support):
	ParameterInterfaceT("Simo_isotropic"),
	fStress(3),
	fModulus(dSymMatrixT::NumValues(3)),
	
	/* work space */
	fb(3),
	fb_bar(3),
	frank4(dSymMatrixT::NumValues(3)),
	
	/* fixed forms */
	fIdentity(3),
	fIcrossI(dSymMatrixT::NumValues(3)),
	fIdentity4(dSymMatrixT::NumValues(3)),
	fDevOp4(dSymMatrixT::NumValues(3))
{	
	/* initialize work matricies */
	fIdentity.Identity();
	fIcrossI.Outer(fIdentity, fIdentity);
	fIdentity4.ReducedIndexI();	
	fDevOp4.ReducedIndexDeviatoric();
}

SimoIso3D::SimoIso3D(void):
	ParameterInterfaceT("Simo_isotropic")
{	

}

/* modulus */
const dMatrixT& SimoIso3D::c_ijkl(void)
{
	/* get mechanical part of the deformation gradient */
	const dMatrixT& F_mech = F_mechanical();

	/* b */
	Compute_b(F_mech, fb);

	/* compute b_bar */
	double J = fb.Det();
	if (J <= 0.0) throw ExceptionT::kBadJacobianDet;
	J = sqrt(J);
	fb_bar.SetToScaled(pow(J,-2.0/3.0), fb);

	ComputeModuli(J, fb_bar, fModulus);
	
	return fModulus;
}
	
/* stress */
const dSymMatrixT& SimoIso3D::s_ij(void)
{
	/* get mechanical part of the deformation gradient */
	const dMatrixT& F_mech = F_mechanical();

	/* b */
	Compute_b(F_mech, fb);

	/* compute b_bar */
	double J = fb.Det();
	if (J <= 0.0) throw ExceptionT::kBadJacobianDet;
	J = sqrt(J);
	fb_bar.SetToScaled(pow(J,-2.0/3.0), fb);

	ComputeCauchy(J, fb_bar, fStress);
	
	return fStress;
}

/* material description */
const dMatrixT& SimoIso3D::C_IJKL(void)
{
	cout << "\n SimoIso3D::C_IJKL: use updated Lagrangian formulation" << endl;
	throw ExceptionT::kGeneralFail;

	return fModulus; // dummy
}

const dSymMatrixT& SimoIso3D::S_IJ(void)
{
	cout << "\n SimoIso3D::S_IJ: use updated Lagrangian formulation" << endl;
	throw ExceptionT::kGeneralFail;

	return fStress; // dummy
}

/* returns the strain energy density for the specified strain */
double SimoIso3D::StrainEnergyDensity(void)
{
	/* get mechanical part of the deformation gradient */
	const dMatrixT& F_mech = F_mechanical();

	/* b */
	Compute_b(F_mech, fb);

	/* compute b_bar */
	double J = fb.Det();
	if (J <= 0.0) throw ExceptionT::kBadJacobianDet;
	J = sqrt(J);
	fb_bar.SetToScaled(pow(J,-2.0/3.0), fb);

	return ComputeEnergy(J, fb_bar);
}

/* accept parameter list */
void SimoIso3D::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	FSIsotropicMatT::TakeParameterList(list);
	
	/* dimension work space */
	fStress.Dimension(3);
	fModulus.Dimension(dSymMatrixT::NumValues(3));
	fb.Dimension(3);
	fb_bar.Dimension(3);
	frank4.Dimension(dSymMatrixT::NumValues(3));
	fIdentity.Dimension(3);
	fIcrossI.Dimension(dSymMatrixT::NumValues(3));
	fIdentity4.Dimension(dSymMatrixT::NumValues(3));
	fDevOp4.Dimension(dSymMatrixT::NumValues(3));

	/* initialize work matricies */
	fIdentity.Identity();
	fIcrossI.Outer(fIdentity, fIdentity);
	fIdentity4.ReducedIndexI();	
	fDevOp4.ReducedIndexDeviatoric();
}

/*************************************************************************
 * Protected
 *************************************************************************/

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
	double mu_bar = Mu()*b_bar.Trace()/(J*3.0);
	moduli.AddScaled(2.0*mu_bar, fDevOp4);

	fStress.SetToScaled(Mu(), b_bar);
	fStress.Deviatoric();
	
	frank4.Outer(fStress,fIdentity);
	frank4.Symmetrize();
	moduli.AddScaled(-4.0/(J*3.0), frank4);
}

void SimoIso3D::ComputeCauchy(double J, const dSymMatrixT& b_bar,
	dSymMatrixT& cauchy)
{
	/* deviatoric */
	cauchy.SetToScaled(Mu()/J,b_bar);
	cauchy.Deviatoric();
	
	/* volumetric */
	cauchy.PlusIdentity(dU(J));
}

double SimoIso3D::ComputeEnergy(double J, const dSymMatrixT& b_bar)
{
	return U(J) +                          /* volumetric */
	       0.5*Mu()*(b_bar.Trace() - 3.0); /* deviatoric */
}
