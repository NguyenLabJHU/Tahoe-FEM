/* $Id: QuadLog3D.cpp,v 1.8 2003-01-29 07:34:48 paklein Exp $ */
/* created: paklein (06/27/1997) */
#include "QuadLog3D.h"

#include <iostream.h>
#include <math.h>

using namespace Tahoe;

/* constructor */
QuadLog3D::QuadLog3D(ifstreamT& in, const FSMatSupportT& support):
	FSSolidMatT(in, support), //in principal stress space
	IsotropicT(in),
	fSpectral(3),
	fb(3),
	fStress(3),
	fModulus(dSymMatrixT::NumValues(3)),
	fDevOp3(3),
		
	/* spectral decomposition */
	fEigs(3),
	floge(3),
	fBeta(3),
	fEigMod(3)
{			
	/* elastic modulis in principal stress space */
	fDevOp3 = -1.0/3.0;
	fDevOp3.PlusIdentity();
	fEigMod = Kappa();
	fEigMod.AddScaled(2.0*Mu(), fDevOp3);
}

/* print parameters */
void QuadLog3D::Print(ostream& out) const
{
	/* inherited */
	FSSolidMatT::Print(out);
	IsotropicT::Print(out);
}

void QuadLog3D::PrintName(ostream& out) const
{
	/* inherited */
	FSSolidMatT::PrintName(out);

	out << "    Quadratic logarithmic isotropic model\n";
}

/* modulus */
const dMatrixT& QuadLog3D::c_ijkl(void)
{
	Compute_b(fb);
	ComputeModuli(fb, fModulus);	
	return fModulus;
}
	
/* stresses */
const dSymMatrixT& QuadLog3D::s_ij(void)
{
	Compute_b(fb);
	ComputeCauchy(fb, fStress);	
	return fStress;
}

/* material description */
const dMatrixT& QuadLog3D::C_IJKL(void)
{
	cout << "\n QuadLog3D::C_IJKL: use updated Lagrangian formulation" << endl;
	throw ExceptionT::kGeneralFail;

	return fModulus; // dummy
}

const dSymMatrixT& QuadLog3D::S_IJ(void)
{
	cout << "\n QuadLog3D::S_IJk: use updated Lagrangian formulation" << endl;
	throw ExceptionT::kGeneralFail;

	return fStress; // dummy
}

/* strain energy density for the specified strain */
double QuadLog3D::StrainEnergyDensity(void)
{
	Compute_b(fb);

	/* principal values */
	fb.PrincipalValues(fEigs);

	/* logarithmic stretches */
	LogStretches(fEigs);

	return ComputeEnergy(floge);
}

/*************************************************************************
* Protected
*************************************************************************/

/* computation routines */
void QuadLog3D::ComputeModuli(const dSymMatrixT& b, dMatrixT& moduli)
{
	/* principal values */
	b.PrincipalValues(fEigs);
	
	/* treat undeformed state separately */
	if ( fabs(fEigs[0] - 1.0) < kSmall &&
	     fabs(fEigs[1] - 1.0) < kSmall &&
	     fabs(fEigs[2] - 1.0) < kSmall )
	{
		IsotropicT::ComputeModuli(moduli);
	}
	/* compute moduli */
	else
	{
		/* full spectral decomposition (and perturb) */
		fSpectral.DecompAndModPrep(b, true);

		/* logarithmic stretches */
		fEigs = fSpectral.Eigenvalues();
		LogStretches(fEigs);
	
		/* construct modulus */
		moduli = fSpectral.EigsToRank4(fEigMod);
		
		/* principal stresses */
		fEigMod.Multx(floge, fBeta);

		/* stress part */
		moduli.AddScaled(2.0*fBeta[0], fSpectral.SpatialTensor(b, 0));
		moduli.AddScaled(2.0*fBeta[1], fSpectral.SpatialTensor(b, 1));
		moduli.AddScaled(2.0*fBeta[2], fSpectral.SpatialTensor(b, 2));
		
		/* factor of J */
		moduli /= sqrt(fEigs[0]*fEigs[1]*fEigs[2]);
	}
}

void QuadLog3D::ComputeCauchy(const dSymMatrixT& b, dSymMatrixT& cauchy)
{
	/* spectral decomposition (don't perturb roots) */
	fSpectral.SpectralDecomp_new(b, false);

	/* logarithmic stretches */
	fEigs = fSpectral.Eigenvalues();
	LogStretches(fEigs);
	
	/* principal stresses */
	fEigMod.Multx(floge, fBeta);

	/* Kirchhoff -> Cauchy */
	fBeta /= sqrt(fEigs[0]*fEigs[1]*fEigs[2]);

	/* QuadLog3D stress */
	cauchy = fSpectral.EigsToRank2(fBeta);
}

double QuadLog3D::ComputeEnergy(const dArrayT& loge)
{
	return 0.5*Lambda()*pow(loge.Sum(), 2.0) +
	           Mu()*dArrayT::Dot(loge, loge);
}

/* logarithmic stretches from the given eigenvalues */
void QuadLog3D::LogStretches(const dArrayT& eigs)
{
	/* logarithmic stretches */
	floge[0] = 0.5*log(eigs[0]);
	floge[1] = 0.5*log(eigs[1]);
	floge[2] = 0.5*log(eigs[2]);
}
