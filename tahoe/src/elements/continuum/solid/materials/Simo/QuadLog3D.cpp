/* $Id: QuadLog3D.cpp,v 1.1.1.1 2001-01-29 08:20:25 paklein Exp $ */
/* created: paklein (06/27/1997)                                          */
/* Hyperelastic material governed by quadratic logarithmic potential.     */

#include "QuadLog3D.h"

#include <iostream.h>
#include <math.h>

#include "ElasticT.h"

const double kpert  = 1.0e-08; //repeated root perturbation

/* constructor */
QuadLog3D::QuadLog3D(ifstreamT& in, const ElasticT& element):
	FDStructMatT(in, element), //in principal stress space
	IsotropicT(in),
	fStress(3),
	fModulus(dSymMatrixT::NumValues(3)),
		
	/* fixed forms - principal stress space */
	fIdentity3(3),
	f1x1(3),
	fDevOp3(3),
	
	/* spectral decomposition */
	fISym(3),
	fEigs(3),
	floge(3),
	fBeta(3),
	fEigMod(3),
	fm(3),
	fn0xn0(fm[0]),
	fn1xn1(fm[1]),
	fn2xn2(fm[2]),
	
	/* decomp work space */
	fm1(3),
	fm2(3),
	fSpatTensor(dSymMatrixT::NumValues(3)),
	fc_b(dSymMatrixT::NumValues(3)),
	fRank4(dSymMatrixT::NumValues(3)),
	fRank2(3),
	fIdentity4(dSymMatrixT::NumValues(3))
{	
	/* initialize fixed forms */
	fISym.Identity();
	fIdentity3.Identity();	
	f1x1 = 1.0;
	fDevOp3.SetToCombination(1.0, fIdentity3, -1.0/3.0, f1x1);
	
	/* bulk and shear moduli */
	Lame(fmu, flambda);
	double kappa = flambda + 2.0/3.0*fmu;
	
	/* elastic modulis in principal stress space */
	fEigMod.SetToCombination(kappa, f1x1, 2.0*fmu, fDevOp3);	
	
	/* dimension rank 1 matrices */
	fn0xn0.Allocate(3);
	fn1xn1.Allocate(3);
	fn2xn2.Allocate(3);
	
	/* spatial tensor work space */
	fIdentity4.ReducedIndexI();
}

/* print parameters */
void QuadLog3D::Print(ostream& out) const
{
	/* inherited */
	FDStructMatT::Print(out);
	IsotropicT::Print(out);
}

void QuadLog3D::PrintName(ostream& out) const
{
	/* inherited */
	FDStructMatT::PrintName(out);

	out << "    Quadratic logarithmic isotropic model\n";
}

/* modulus */
const dMatrixT& QuadLog3D::c_ijkl(void)
{
	ComputeModuli(b(), fModulus);	
	return fModulus;
}
	
/* stresses */
const dSymMatrixT& QuadLog3D::s_ij(void)
{
	ComputeCauchy(b(), fStress);	
	return fStress;
}

/* material description */
const dMatrixT& QuadLog3D::C_IJKL(void)
{
	cout << "\n QuadLog3D::C_IJKL: use updated Lagrangian formulation" << endl;
	throw eGeneralFail;

	return fModulus; // dummy
}

const dSymMatrixT& QuadLog3D::S_IJ(void)
{
	cout << "\n QuadLog3D::S_IJk: use updated Lagrangian formulation" << endl;
	throw eGeneralFail;

	return fStress; // dummy
}

/* strain energy density for the specified strain */
double QuadLog3D::StrainEnergyDensity(void)
{
	/* principal values */
	b().PrincipalValues(fEigs);

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
	//Assume spectral decomp occurs only during s_ij???
	
	/* principal values */
	b.PrincipalValues(fEigs);
	
	/* treat undeformed state separately */
	if ( fabs(fEigs[0] - 1.0) < kSmall &&
	     fabs(fEigs[1] - 1.0) < kSmall &&
	     fabs(fEigs[2] - 1.0) < kSmall )
	{
		double mu, lambda;
		Lame(mu, lambda);
		
		IsotropicT::ComputeModuli(moduli, mu, lambda);
	}
	/* compute moduli */
	else
	{
		/* spectral decomposition (and perturb) */
		SpectralDecomp(b, fEigs, true);

		/* logarithmic stretches */
		LogStretches(fEigs);
	
		/* initialize */
		moduli = 0.0;

		/* using symmetry in A and B */
		for (int B = 0; B < 3; B++)
			for (int A = B; A < 3; A++)
			{
				fRank4.Outer(fm[A],fm[B]);
				
				double gamma = 1.0;
				if (A != B)
				{
					gamma = 2.0;
					fRank4.Symmetrize();
				}
				
				moduli.AddScaled(gamma*fEigMod(A,B), fRank4);
			}
		
		/* principal stresses */
		fEigMod.Multx(floge,fBeta);

		/* spatial tensor component */
		Set_b_Tensor(b);
		
		/* principal spatial tensor */
		for (int A = 0; A < 3; A++)
			moduli.AddScaled(2.0*fBeta[A], SpatialTensor(b,A) );
		
		/* factor of J */
		moduli /= sqrt( fEigs[0]*fEigs[1]*fEigs[2] );
	}
}

void QuadLog3D::ComputeCauchy(const dSymMatrixT& b, dSymMatrixT& cauchy)
{
	/* principal values */
	b.PrincipalValues(fEigs);

	/* spectral decomposition (don't perturb roots) */
	SpectralDecomp(b, fEigs, false);

	/* logarithmic stretches */
	LogStretches(fEigs);
	
	/* principal stresses */
	fEigMod.Multx(floge,fBeta);

	/* Kirchhoff -> Cauchy */
	fBeta /= sqrt(fEigs[0]*fEigs[1]*fEigs[2]);

	/* QuadLog3D stress */
	cauchy.SetToCombination(fBeta[0], fn0xn0,
	                        fBeta[1], fn1xn1,
	                        fBeta[2], fn2xn2);
}

double QuadLog3D::ComputeEnergy(const dArrayT& loge)
{
	return 0.5*flambda*pow(loge.Sum(),2.0) +
	           fmu*dArrayT::Dot(loge,floge);
}

/* compute spectral decomposition of b
*
* NOTE: Repeated eigenvalues are perturbed */
void QuadLog3D::SpectralDecomp(const dSymMatrixT& b, dArrayT& eigs, bool perturb_repeated)
{
	/* check for repeated roots */
	double d01 = eigs[0] - eigs[1];
	double d12 = eigs[1] - eigs[2];
	double d02 = eigs[0] - eigs[2];

	/* 3 distinct roots */
	if (Min(fabs(d01), fabs(d12), fabs(d02)) > kSmall)
	{
		double k01 = 1.0/d01;
		double k12 = 1.0/d12;
		double k02 = 1.0/d02;

		fm1.SetToCombination(-k01, b, k01*eigs[1], fISym);
		fm2.SetToCombination(-k02, b, k02*eigs[2], fISym);
		fn0xn0.MultAB(fm1,fm2); /* m1 and m2 commute */

		fm1.SetToCombination(-k12, b, k12*eigs[2], fISym);
		fm2.SetToCombination( k01, b,-k01*eigs[0], fISym);
		fn1xn1.MultAB(fm1,fm2); /* m1 and m2 commute */

		fm1.SetToCombination( k02, b,-k02*eigs[0], fISym);
		fm2.SetToCombination( k12, b,-k12*eigs[1], fISym);
		fn2xn2.MultAB(fm1,fm2); /* m1 and m2 commute */
	}
	// Perturb repeated roots: no perturbations are applied to
	// to the out-of-plane root, to preserve possible plane
	// strain constraints.

	/* 3 identical roots */
	else if (fabs(d01) < kSmall &&
	         fabs(d12) < kSmall &&
	         fabs(d02) < kSmall)
	{
		fn0xn0 = 0.0;
		fn1xn1 = 0.0;
		fn2xn2 = 0.0;

		fn0xn0(0,0) = 1.0;
		fn1xn1(1,1) = 1.0;
		fn2xn2(2,2) = 1.0;
		
		/* perturb roots */
		if (perturb_repeated)
		{
			eigs[0] -= kpert;
			eigs[1] += kpert;
		}
	}
	else /* 2 repeated roots */
	{
		if (fabs(d01) < kSmall)
		{
			double k2x =-1.0/d02;
		
			fn2xn2.SetToCombination(k2x, b, -eigs[0]*k2x, fISym);
			fn0xn0.SetToCombination(1.0, fISym,-1.0, fn2xn2);
			//fn1xn1 = fn0xn0;
			fn1xn1 = 0.0;
			
			/* perturb roots */
			if (perturb_repeated)
			{
				if (eigs[2] > eigs[0])
				eigs[0] -= kpert;
					else
				eigs[0] += kpert;
			}
		}
		else if (fabs(d12) < kSmall)
		{
			double k0x = 1.0/d01;
		
			fn0xn0.SetToCombination(k0x, b, -eigs[1]*k0x, fISym);
			fn1xn1.SetToCombination(1.0, fISym,-1.0, fn0xn0);
			//fn2xn2 = fn1xn1;
			fn2xn2 = 0.0;
		
			/* perturb roots */
			if (perturb_repeated)
			{
				if (eigs[0] > eigs[1])
					eigs[1] -= kpert;
				else
					eigs[1] += kpert;
			}
		}
		else if (fabs(d02) < kSmall)
		{
			double k1x = 1.0/d12;
		
			fn1xn1.SetToCombination(k1x, b, -eigs[2]*k1x, fISym);
			fn2xn2.SetToCombination(1.0, fISym,-1.0, fn1xn1);
			//fn0xn0 = fn2xn2;
			fn0xn0 = 0.0;
		
			/* perturb roots */
			if (perturb_repeated)
			{
				if (eigs[1] > eigs[2])
					eigs[0] -= kpert;
				else
					eigs[0] += kpert;
			}
		}
	}
}

/* logarithmic stretches from the given eigenvalues */
void QuadLog3D::LogStretches(dArrayT& eigs)
{
	/* logarithmic stretches */
	floge[0] = 0.5*log(eigs[0]);
	floge[1] = 0.5*log(eigs[1]);
	floge[2] = 0.5*log(eigs[2]);
}

/* compute contribution to spatial tensor that depends only on b */
void QuadLog3D::Set_b_Tensor(const dSymMatrixT& b)
{
	double z1, z2, z3, z4, z5, z6, z7, z8, z9, z10, z11, z12;
	double z13, z14, z15, z16, z17, z18, z19, z20, z21;
	
	z1 = b[0];
	z2 = b[1];
	z3 = b[2];
	z4 = b[3];
	z5 = b[4];
	z6 = b[5];
	z7 = z1*z1;
	z8 = z1*z2;
	z9 = z2*z2;
	z10 = z1*z3;
	z11 = z2*z3;
	z12 = z3*z3;
	z13 = z1*z4;
	z14 = z2*z4;
	z15 = z3*z4;
	z16 = z4*z4;
	z17 = z1*z5;
	z18 = z2*z5;
	z19 = z3*z5;
	z20 = z4*z5;
	z21 = z5*z5;
	z1 = z1*z6;
	z2 = z2*z6;
	z3 = z3*z6;
	z4 = z4*z6;
	z5 = z5*z6;
	z6 = z6*z6;
	z11 = z11 + z16;
	z10 = z10 + z21;
	z3 = z20 + z3;
	z18 = z18 + z4;
	z13 = z13 + z5;
	z8 = z6 + z8;
	z11 = 0.5*z11;
	z10 = 0.5*z10;
	z3 = 0.5*z3;
	z18 = 0.5*z18;
	z13 = 0.5*z13;
	z8 = 0.5*z8;

	//{{z7, z6, z21, z5, z17, z1},
	// {z6, z9, z16, z14, z4, z2},
	// {z21, z16, z12, z15, z19, z20},
	// {z5, z14, z15, z11, z3, z18},
	// {z17, z4, z19, z3, z10, z13},
// {z1, z2, z20, z18, z13, z8}}

fc_b(0,0) = z7;
fc_b(0,1) = fc_b(1,0) = z6;
fc_b(0,2) = fc_b(2,0) = z21;
fc_b(0,3) = fc_b(3,0) = z5;
fc_b(0,4) = fc_b(4,0) = z17;
fc_b(0,5) = fc_b(5,0) = z1;
fc_b(1,1) = z9;
fc_b(1,2) = fc_b(2,1) = z16;
fc_b(1,3) = fc_b(3,1) = z14;
fc_b(1,4) = fc_b(4,1) = z4;
fc_b(1,5) = fc_b(5,1) = z2;
fc_b(2,2) = z12;
fc_b(2,3) = fc_b(3,2) = z15;
fc_b(2,4) = fc_b(4,2) = z19;
fc_b(2,5) = fc_b(5,2) = z20;
fc_b(3,3) = z11;
fc_b(3,4) = fc_b(4,3) = z3;
fc_b(3,5) = fc_b(5,3) = z18;
fc_b(4,4) = z10;
fc_b(4,5) = fc_b(5,4) = z13;
fc_b(5,5) = z8;

	/* b (x) b */
	fRank4.Outer(b,b);
	
	/* assemble */
	fc_b -= fRank4;
}

/* compute the principal spatial tensor associated with the Ath
* eigenvalue */
const dMatrixT& QuadLog3D::SpatialTensor(const dSymMatrixT& b, int A)
{
	/* initialize */
	fSpatTensor = 0.0;

	/* cyclic permutation */
	int B = int(fmod(A + 1, 3));
	int C = int(fmod(B + 1, 3));
	double dA = (fEigs[A] - fEigs[B])*(fEigs[A] - fEigs[C]);
			
	/* I_b - b (x) b */
	fSpatTensor.AddScaled(1.0/dA, fc_b);

	/* I - (1-m_A) (x) (1-m_A) */	
	fRank2.SetToCombination(1.0, fISym, -1.0, fm[A]);
	fRank4.Outer(fRank2,fRank2);
	double k = b.Det()/(dA*fEigs[A]);
	fSpatTensor.AddCombination(-k, fIdentity4, k, fRank4);
	
	/* [ b (x) m_A ]^s */
	fRank4.Outer(b,fm[A]);
	fRank4.Symmetrize();
	fSpatTensor.AddScaled(2.0*fEigs[A]/dA, fRank4);
		
	/* m_A (x) m_A */
	fRank4.Outer(fm[A],fm[A]);
	fSpatTensor.AddScaled((fEigs[A]/dA)*(b.Trace() - 4.0*fEigs[A]), fRank4);

	return fSpatTensor;
}
