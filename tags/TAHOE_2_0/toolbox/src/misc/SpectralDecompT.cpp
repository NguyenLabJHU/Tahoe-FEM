/* $Id: SpectralDecompT.cpp,v 1.11 2003-11-21 22:41:57 paklein Exp $ */
/* created: paklein (11/09/1997)                                          */
/* Spectral decomposition solver                                          */

#include "SpectralDecompT.h"
#include <iostream.h>
#include <math.h>
#include <float.h>


using namespace Tahoe;

const double kpert     = 1.0e-08; // repeated root perturbation
const double sqrt1by3  = 1.0/sqrt(3.0);
const double sqrt1by6  = 1.0/sqrt(6.0);
const double kSameRoot = 1.0e-8; // kSmall is too small

/* constructor */
SpectralDecompT::SpectralDecompT(int nsd):
	/* fixed forms - principal stress space */
	f_I_Rank2(nsd),
	f_I_Rank4(dSymMatrixT::NumValues(nsd)),

	/* spectral decomp */
	fEigs(nsd),
	fm(nsd),
	
	/* spectral decomp work space */
	fm1(nsd),
	fm2(nsd),
	fEvecMatrix(nsd),
	fEvecs(nsd),

	/* polar decomp work space */
	fInvEigs(nsd),
	fUInv(nsd),
	
	/* spatial tensor work space */
	fSpatTensor(dSymMatrixT::NumValues(nsd)),
	fc_b(dSymMatrixT::NumValues(nsd)),
	fRank4(dSymMatrixT::NumValues(nsd)),
	fRank2(nsd)	
{	
	/* initialize fixed forms */
	f_I_Rank2.Identity();
	f_I_Rank4.ReducedIndexI();
		
	/* dimension rank 1 matrices */
	for (int i = 0; i < nsd; i++)
		fm[i].Dimension(nsd);

	/* eigenvectors in columns */
	for (int j = 0; j < nsd; j++)
		fEvecs[j].Set(nsd, fEvecMatrix(j));
}

/* compute spectral decomposition of rank2
*
* NOTE: Repeated eigenvalues are NOT perturbed */
void SpectralDecompT::SpectralDecomp(const dSymMatrixT& rank2, const dArrayT& eigs,
	bool perturb_repeated)
{
	/* copy in */
	fEigs = eigs;

	/* dispatch */
	if (rank2.Rows() == 2)
	{
		/* check for repeated roots */
		double deig = fEigs[0] - fEigs[1];
	
		/* 2 distinct roots */
		if (fabs(deig) > kSmall)
		{
			fm[0].SetToCombination(1.0/deig, rank2, -fEigs[1]/deig, f_I_Rank2);
			fm[1].SetToCombination(1.0, f_I_Rank2,-1.0, fm[0]);
		}
		else
		{
			fm[0] = 0.0;
			fm[0](0,0) = 1.0;
	
			fm[1] = 0.0;
			fm[1](1,1) = 1.0;
		}
	}
	else
	{
		/* shift */
		double shift = (fEigs[0] + fEigs[1] + fEigs[2])/3.0;
		double* tmp = (double*) rank2.Pointer();
		tmp[0] -= shift;
		tmp[1] -= shift;
		tmp[2] -= shift;
		fEigs[0] -= shift;
		fEigs[1] -= shift;
		fEigs[2] -= shift;
	
		/* do decomp */
		SpectralDecomp3D(rank2, fEigs, perturb_repeated);	
		
		/* shift back */
		tmp[0] += shift;
		tmp[1] += shift;
		tmp[2] += shift;
		fEigs[0] += shift;
		fEigs[1] += shift;
		fEigs[2] += shift;
	}
}

void SpectralDecompT::SpectralDecomp(const dSymMatrixT& rank2, bool perturb_repeated)
{
	/* compute principal values */
	rank2.PrincipalValues(fEigs);
	
	SpectralDecomp(rank2, fEigs, perturb_repeated);
}

/* compute spectral decomposition using Jacobi iterations */
void SpectralDecompT::SpectralDecomp_Jacobi(const dSymMatrixT& rank2, bool perturb_repeated)
{
	bool sort_descending = false;
	rank2.Eigensystem(fEigs, fEvecMatrix, sort_descending);

	/* set rank 1 tensors */
	int nsd = rank2.Rows();
	if (nsd == 2)
	{
		double* pvec = fEvecMatrix.Pointer();
		dSymMatrixT& n0n0 = fm[0];
		n0n0[0] = pvec[0]*pvec[0];
		n0n0[1] = pvec[1]*pvec[1];
		n0n0[2] = pvec[0]*pvec[1];

		pvec += fEvecMatrix.Rows();
		dSymMatrixT& n1n1 = fm[1];
		n1n1[0] = pvec[0]*pvec[0];
		n1n1[1] = pvec[1]*pvec[1];
		n1n1[2] = pvec[0]*pvec[1];		
	}
	else if (nsd == 3)
	{
		double* pvec = fEvecMatrix.Pointer();
		dSymMatrixT& n0n0 = fm[0];
		n0n0[0] = pvec[0]*pvec[0];
		n0n0[1] = pvec[1]*pvec[1];
		n0n0[2] = pvec[2]*pvec[2];
		n0n0[3] = pvec[1]*pvec[2];
		n0n0[4] = pvec[0]*pvec[2];
		n0n0[5] = pvec[0]*pvec[1];

		pvec += fEvecMatrix.Rows();
		dSymMatrixT& n1n1 = fm[1];
		n1n1[0] = pvec[0]*pvec[0];
		n1n1[1] = pvec[1]*pvec[1];
		n1n1[2] = pvec[2]*pvec[2];
		n1n1[3] = pvec[1]*pvec[2];
		n1n1[4] = pvec[0]*pvec[2];
		n1n1[5] = pvec[0]*pvec[1];

		pvec += fEvecMatrix.Rows();
		dSymMatrixT& n2n2 = fm[2];
		n2n2[0] = pvec[0]*pvec[0];
		n2n2[1] = pvec[1]*pvec[1];
		n2n2[2] = pvec[2]*pvec[2];
		n2n2[3] = pvec[1]*pvec[2];
		n2n2[4] = pvec[0]*pvec[2];
		n2n2[5] = pvec[0]*pvec[1];
	}
	else
		throw ExceptionT::kGeneralFail;
		
	/* perturb repeated roots */
	if (perturb_repeated) PerturbRepeated(fEigs);
}

void SpectralDecompT::ModulusPrep(const dSymMatrixT& rank2)
{
	/* modulus tensors */
	if (rank2.Rows() == 2)
		fc_b = 0.0;
		//TEMP
		//Set_I4_Tensor2D(rank2, fc_b);
		//fc_b.ReducedIndexI();
//DEV - what's going on here?
	else
	{
		Set_I4_Tensor3D(rank2, fc_b);
	
		/* b (x) b */
		fRank4.Outer(rank2,rank2);
	
		/* assemble */
		fc_b -= fRank4;
	}
}

/* return the tensor for the given eigenvalues */
const dSymMatrixT& SpectralDecompT::EigsToRank2(const dArrayT& eigs)
{
	if (eigs.Length() == 2)
		fRank2.SetToCombination(eigs[0], fm[0],
	                            eigs[1], fm[1]);
	else
		fRank2.SetToCombination(eigs[0], fm[0],
	                            eigs[1], fm[1],
	                            eigs[2], fm[2]);
	return fRank2;
}

/* compute the polar decomposition of F = RU, where R^T R = 1 and U = U^T */
void SpectralDecompT::PolarDecomp(const dMatrixT& F, dMatrixT& R, dSymMatrixT& U,
	bool perturb_repeated)
{
#if __option(extended_errorcheck)
	if (fEigs.Length() != F.Rows() ||
	          F.Rows() != R.Rows() ||
	          R.Rows() != U.Rows()) throw ExceptionT::kSizeMismatch;
#endif

	/* construct stretch C */
	U.MultATA(F);
	
	/* spectral decomposition */
	//SpectralDecomp(U, perturb_repeated);
	SpectralDecomp_Jacobi(U, perturb_repeated);

	/* eigenvalues to stretches */
	if (fEigs[0] <= 0.0) throw ExceptionT::kBadJacobianDet; fEigs[0] = sqrt(fEigs[0]);
	if (fEigs[1] <= 0.0) throw ExceptionT::kBadJacobianDet; fEigs[1] = sqrt(fEigs[1]);
	if (fEigs.Length() == 3)
	{
		if (fEigs[2] <= 0.0) throw ExceptionT::kBadJacobianDet;
		fEigs[2] = sqrt(fEigs[2]);
	}
	
	/* inverse */
	fInvEigs[0] = 1.0/fEigs[0];
	fInvEigs[1] = 1.0/fEigs[1];
	if (fEigs.Length() == 3) fInvEigs[2] = 1.0/fEigs[2];
	
	/* U^-1 */
	EigsToRank2(fInvEigs).ToMatrix(fUInv);
	
	/* R */
	R.MultAB(F, fUInv);
	
	/* U */
	U = EigsToRank2(fEigs);
}

/* function to perturb repeated roots */
bool SpectralDecompT::PerturbRepeated(dArrayT& values) const
{
	/* check */
	if (values.Length() != 2 && values.Length() != 3)
	{
		cout << "\n SpectralDecompT::PerturbRepeated: expecting array length 2 or 3" << endl;
		throw ExceptionT::kGeneralFail;
	}

	/* perturb repeated values */
	bool perturbed = false;
	if (values.Length() == 2)
	{
		if (fabs(values[0] - values[1]) < kSameRoot)
		{
			values[0] *= (1.0 + kpert);
			values[1] /= (1.0 + kpert);
			perturbed = true;
		}
	}
	else
	{
		if (fabs(values[0] - values[1]) < kSameRoot)
		{
			values[0] *= (1.0 + kpert);
			values[1] /= (1.0 + kpert);
			if (fabs(values[0] - values[2]) < kSameRoot)
				values[0] /= (1.0 + kpert);
			else if (fabs(values[1] - values[2]) < kSameRoot)
				values[1] *= (1.0 + kpert);
			perturbed = true;
		}
		else if (fabs(values[0] - values[2]) < kSameRoot)
		{
			values[0] *= (1.0 + kpert);
			values[2] /= (1.0 + kpert);
			if (fabs(values[0] - values[1]) < kSameRoot)
				values[0] /= (1.0 + kpert);
			else if (fabs(values[2] - values[1]) < kSameRoot)
				values[2] *= (1.0 + kpert);
			perturbed = true;
		}
		else if (fabs(values[1] - values[2]) < kSameRoot)
		{
			values[1] *= (1.0 + kpert);
			values[2] /= (1.0 + kpert);
			if (fabs(values[1] - values[0]) < kSameRoot)
				values[1] /= (1.0 + kpert);
			else if (fabs(values[2] - values[0]) < kSameRoot)
				values[2] *= (1.0 + kpert);
			perturbed = true;
		}
	}
	return perturbed;
}

const dMatrixT& SpectralDecompT::EigsToRank4(const dSymMatrixT& eigs)
{
	/* initialize */
	fSpatTensor = 0.0;

	/* stress derivative term */
	if (eigs.Rows() == 2)
	{		
		fRank4.Outer(fm[0],fm[0]);
		fSpatTensor.AddScaled(eigs(0,0),fRank4);
		
		fRank4.Outer(fm[1],fm[1]);
		fSpatTensor.AddScaled(eigs(1,1),fRank4);
		
		fRank4.Outer(fm[0],fm[1]);
		fRank4.Symmetrize();
		fSpatTensor.AddScaled(2.0*eigs(0,1),fRank4);
	}
	else
		for (int B = 0; B < 3; B++)
			for (int A = B; A < 3; A++) /* using symmetry in A and B */
			{
				fRank4.Outer(fm[A],fm[B]);
				double gamma = 1.0;
				if (A != B)
				{
					gamma = 2.0;
					fRank4.Symmetrize();
				}
				fSpatTensor.AddScaled(gamma*eigs(A,B), fRank4);
			}

	return fSpatTensor;
}
/*Created by TDN: 03/05/2001
 Forms rank 4 tensor from nonsymmetric matrix of derivative of principal stress*/
const dMatrixT& SpectralDecompT::NonSymEigsToRank4(const dMatrixT& eigs)
{
	/*initialize*/
	fSpatTensor = 0.0;
	
	/*stress derivative term*/
	int nsd = eigs.Rows();
	for (int B = 0; B < nsd; B++)
		for (int A = 0; A < nsd; A++)
		{
			fRank4.Outer(fm[A],fm[B]);
			fSpatTensor.AddScaled(eigs(A,B),fRank4);
		}
	return fSpatTensor;
}
	
/* compute the principal spatial tensor associated with the Ath
* eigenvalue */
const dMatrixT& SpectralDecompT::SpatialTensor(const dSymMatrixT& b, int A)
{
	/* dispatch */
	if (b.Rows() == 2)
		return(SpatialTensor2D(b,A));
	else
		return(SpatialTensor3D(b,A));
}

/***********************************************************************
* Private
***********************************************************************/

//TEMP
/* returns max */
double SpectralDecompT::Min(double d1, double d2, double d3)
{
	return ( (d1 < d2) ?
			((d1 < d3) ? d1 : d3) :
			((d2 < d3) ? d2 : d3) );
}

/* compute contribution to spatial tensor that depends only on mat */
void SpectralDecompT::Set_I4_Tensor3D(const dSymMatrixT& mat,
	dMatrixT& rank4)
{
	double z1, z2, z3, z4, z5, z6, z7, z8, z9, z10, z11, z12;
	double z13, z14, z15, z16, z17, z18, z19, z20, z21;
	
	z1 = mat[0];
	z2 = mat[1];
	z3 = mat[2];
	z4 = mat[3];
	z5 = mat[4];
	z6 = mat[5];
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

rank4(0,0) = z7;
rank4(0,1) = rank4(1,0) = z6;
rank4(0,2) = rank4(2,0) = z21;
rank4(0,3) = rank4(3,0) = z5;
rank4(0,4) = rank4(4,0) = z17;
rank4(0,5) = rank4(5,0) = z1;
rank4(1,1) = z9;
rank4(1,2) = rank4(2,1) = z16;
rank4(1,3) = rank4(3,1) = z14;
rank4(1,4) = rank4(4,1) = z4;
rank4(1,5) = rank4(5,1) = z2;
rank4(2,2) = z12;
rank4(2,3) = rank4(3,2) = z15;
rank4(2,4) = rank4(4,2) = z19;
rank4(2,5) = rank4(5,2) = z20;
rank4(3,3) = z11;
rank4(3,4) = rank4(4,3) = z3;
rank4(3,5) = rank4(5,3) = z18;
rank4(4,4) = z10;
rank4(4,5) = rank4(5,4) = z13;
rank4(5,5) = z8;
}

void SpectralDecompT::Set_I4_Tensor2D(const dSymMatrixT& mat,
	dMatrixT& rank4)
{
	/* pull values */
	double m11 = mat(0,0);
	double m22 = mat(1,1);
	double m12 = mat(0,1);

rank4(0,0) = m11*m11;
rank4(0,1) = rank4(1,0) = m12*m12;
rank4(0,2) = rank4(2,0) = m11*m12;

rank4(1,1) = m22*m22;
rank4(1,2) = rank4(2,1) = m22*m12;

rank4(2,2) = 0.5*(m11*m22 + m12*m12);
}

/* work routines */
void SpectralDecompT::SpectralDecomp3D(const dSymMatrixT& rank2, dArrayT& eigs,
	bool perturb_repeated)
{
	/* aliases */
	dSymMatrixT& n0xn0 = fm[0];
	dSymMatrixT& n1xn1 = fm[1];
	dSymMatrixT& n2xn2 = fm[2];

	/* check for repeated roots */
	double d01 = eigs[0] - eigs[1];
	double d12 = eigs[1] - eigs[2];
	double d02 = eigs[0] - eigs[2];
	
	/* 3 identical roots or rank2 is axial */
	if ((fabs(d01) < kSameRoot && fabs(d12) < kSameRoot && fabs(d02) < kSameRoot) ||
	(fabs(rank2[3]) < kSameRoot && fabs(rank2[4]) < kSameRoot && fabs(rank2[5]) < kSameRoot))
	{
		n0xn0 = 0.0;
		n1xn1 = 0.0;
		n2xn2 = 0.0;

		/* spherical */
		n0xn0(0,0) = 1.0;
		n1xn1(1,1) = 1.0;
		n2xn2(2,2) = 1.0;
		
		/* remove repeated */
		//if (perturb_repeated) PerturbRepeated(eigs);
	}
	else
	{
		/* (symmetrically) perturb repeated roots */
		if (perturb_repeated)
		{
			if (fabs(d01) < kSmall)
			{
				double pert = 0.5*(eigs[0] + eigs[1])*kpert;
				eigs[0] += pert;				
				eigs[1] -= pert;				
				d01 = eigs[0] - eigs[1];
			}
			else if (fabs(d12) < kSmall)
			{
				double pert = 0.5*(eigs[1] + eigs[2])*kpert;
				eigs[1] += pert;				
				eigs[2] -= pert;				
				d12 = eigs[1] - eigs[2];
			}
			else if (fabs(d02) < kSmall)
			{
				double pert = 0.5*(eigs[0] + eigs[2])*kpert;
				eigs[0] += pert;				
				eigs[2] -= pert;				
				d02 = eigs[0] - eigs[2];
			}
		}

		/* 3 distinct roots */
	//	if (Min(fabs(d01), fabs(d12), fabs(d02)) > kSmall) // too loose
		if (d01 != 0.0 && d12 != 0.0 && d02 != 0.0)
		{
			double k01 = 1.0/d01;
			double k12 = 1.0/d12;
			double k02 = 1.0/d02;
	
			fm1.SetToCombination(-k01, rank2, k01*eigs[1], f_I_Rank2);
			fm2.SetToCombination(-k02, rank2, k02*eigs[2], f_I_Rank2);
			n0xn0.MultAB(fm1,fm2);

			fm1.SetToCombination(-k12, rank2, k12*eigs[2], f_I_Rank2);
			fm2.SetToCombination( k01, rank2,-k01*eigs[0], f_I_Rank2);
			n1xn1.MultAB(fm1,fm2);

			fm1.SetToCombination( k02, rank2,-k02*eigs[0], f_I_Rank2);
			fm2.SetToCombination( k12, rank2,-k12*eigs[1], f_I_Rank2);
			n2xn2.MultAB(fm1,fm2);
		}
		else /* 2 repeated roots - construct orthogonal basis */
		{
			if (fabs(d01) < kSmall)
				SchmidtDecompose(rank2, eigs[2], n2xn2, eigs[0], n0xn0, n1xn1);
			else if (fabs(d12) < kSmall)
				SchmidtDecompose(rank2, eigs[0], n0xn0, eigs[1], n1xn1, n2xn2);
			else if (fabs(d02) < kSmall)
				SchmidtDecompose(rank2, eigs[1], n1xn1, eigs[2], n2xn2, n0xn0);
		}
	}
}

/* compute the principal spatial tensor associated with the Ath
* eigenvalue - from Simo, CMAME (1992) */
const dMatrixT& SpectralDecompT::SpatialTensor3D(const dSymMatrixT& b, int A)
{
	/* initialize */
	fSpatTensor = 0.0;

	/* cyclic permutation */
	double dA;

	int B = int(fmod(A + 1.0, 3));
	int C = int(fmod(B + 1.0, 3));
	dA = (fEigs[A] - fEigs[B])*(fEigs[A] - fEigs[C]);

	/* checks */
	if (fabs(dA) < kSmall)
	{
//		cout << "\n SpectralDecompT::SpatialTensor3D: detected repeated roots:\n"
//		     << fEigs << endl;
		fSpatTensor.Identity();
		return fSpatTensor;
	}
	else if (fEigs[A] <= 0.0) throw ExceptionT::kBadJacobianDet;

	/* I_b - b (x) b */
	double k1 = 1.0/dA;
	fSpatTensor.AddScaled(k1, fc_b);

	/* I - (1-m_A) (x) (1-m_A) */	
	fRank2.SetToCombination(1.0, f_I_Rank2, -1.0, fm[A]);
	fRank4.Outer(fRank2,fRank2);
	double k2 = b.Det()/(dA*fEigs[A]);
	fSpatTensor.AddCombination(-k2,f_I_Rank4,k2,fRank4);

	/* [ b (x) m_A ]^s */
	fRank4.Outer(b,fm[A]);
	fRank4.Symmetrize();
	double k3 = 2.0*fEigs[A]/dA;
	fSpatTensor.AddScaled(k3,fRank4);

	/* m_A (x) m_A */
	fRank4.Outer(fm[A],fm[A]);
	double k4 = (fEigs[A]/dA)*(b.Trace() - 4.0*fEigs[A]);
	fSpatTensor.AddScaled(k4,fRank4);

	return fSpatTensor;
}

/* compute the principal spatial tensor associated with the Ath
* eigenvalue */
const dMatrixT& SpectralDecompT::SpatialTensor2D(const dSymMatrixT& b, int A)
{
#pragma unused(b)

	/* initialize */
	fSpatTensor = 0.0;

	/* cyclic permutation */
	int B     = (A == 1) ? 0 : 1;
	double dA = fEigs[A] - fEigs[B];
	double k  = fEigs[B]/dA;

	/* I - (1-m_A) (x) (1-m_A) */	
	fRank2.SetToCombination(1.0, f_I_Rank2, -1.0, fm[A]);
	fRank4.Outer(fRank2,fRank2);
	fSpatTensor.AddCombination(k,f_I_Rank4,-k,fRank4);

	/* m_A (x) m_A */
	fRank4.Outer(fm[A],fm[A]);
	fSpatTensor.AddScaled(-fEigs[A]/dA, fRank4);
							
	return(fSpatTensor);
}

/* construct orthogonal basis for 2 repeated roots */
void SpectralDecompT::SchmidtDecompose(const dSymMatrixT& rank2,
	double l2, dSymMatrixT& n2xn2,
	double l , dSymMatrixT& n0xn0, dSymMatrixT& n1xn1)
{
	double d2 = 1.0/(l2 - l);		
	n2xn2.SetToCombination(d2, rank2, -l*d2, f_I_Rank2);

	/* check values */
	if ((n2xn2[0] < 0.0 && fabs(n2xn2[0]) > kSmall) ||
	    (n2xn2[1] < 0.0 && fabs(n2xn2[1]) > kSmall) ||
	    (n2xn2[2] < 0.0 && fabs(n2xn2[2]) > kSmall))
	{
		int prec_old = cout.precision();
		cout.precision(DBL_DIG);
		cout << "\n SpectralDecompT::SchmidtDecompose: invalid decomposition\n";
		cout << " tensor =\n" << rank2 << '\n';
		cout << " 2x eig = " << l  << '\n';
		cout << " 1x eig = " << l2 << endl;
		cout << " 2x nxn = (diagonals must be > 0)\n" << n2xn2 << '\n';
		cout.precision(prec_old);
		throw ExceptionT::kGeneralFail;
	}
			
	/* components of unique eigenvector - chop noise */
	double a0 = (fabs(n2xn2[0]) < kSmall) ? 0.0 : sqrt(n2xn2[0]);
	double a1 = (fabs(n2xn2[1]) < kSmall) ? 0.0 : sqrt(n2xn2[1]);
	double a2 = (fabs(n2xn2[2]) < kSmall) ? 0.0 : sqrt(n2xn2[2]);
			
	/* construct orthogonal vector (Gramm-Schmidt) */
	double b0, b1, b2;
	if (fabs(a0 - a1) < kSmall && fabs(a1 - a2) < kSmall)
	{
		b0 = -2.0*sqrt1by6;
		b1 = b2 = sqrt1by6;		
	}
	else
	{
		double dot = sqrt1by3*(a0 + a1 + a2);
		b0 = sqrt1by3 - dot*a0;
		b1 = sqrt1by3 - dot*a1;
		b2 = sqrt1by3 - dot*a2;

		double lb = 1.0/sqrt(b0*b0 + b1*b1 + b2*b2);
		b0 *= lb;
		b1 *= lb;
		b2 *= lb;
	}
			
	/* first rank 1 tensor */
	n0xn0[0] = b0*b0;
	n0xn0[1] = b1*b1;
	n0xn0[2] = b2*b2;
	n0xn0[3] = b2*b1;
	n0xn0[4] = b2*b0;
	n0xn0[5] = b1*b0;
			
	/* second rank 1 tensor */
	n1xn1.SetToCombination(1.0/l, rank2,
                           -l2/l, n2xn2,
                            -1.0, n0xn0);
}	
