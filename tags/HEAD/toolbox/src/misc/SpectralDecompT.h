/* $Id: SpectralDecompT.h,v 1.1.1.1 2001-01-25 20:56:25 paklein Exp $ */
/* created: paklein (11/09/1997)                                          */
/* Spectral decomposition solver                                          */

#ifndef _SPECTRAL_DECOMP_T_H_
#define _SPECTRAL_DECOMP_T_H_

/* direct members */
#include "dMatrixT.h"
#include "dSymMatrixT.h"

/* forward declarations */
class dArrayT;

class SpectralDecompT
{
public:

	/* constructor */
	SpectralDecompT(int nsd);
		
	/* compute spectral decomposition of rank2
* NOTE: Repeated eigenvalues are NOT perturbed */
	void SpectralDecomp(const dSymMatrixT& rank2, bool perturb_repeated);
	void SpectralDecomp(const dSymMatrixT& rank2, const dArrayT& eigs, bool perturb_repeated);
	void DecompAndModPrep(const dSymMatrixT& rank2, bool perturb_repeated); //before forming rank4's

	/* compute spectral decomposition (using Jacobi iterations in 3D) */
	void SpectralDecomp_new(const dSymMatrixT& rank2);

	/* compute the polar decomposition of F = RU, where R^T R = 1 and U = U^T */
	void PolarDecomp(const dMatrixT& F, dMatrixT& R, dSymMatrixT& U, bool perturb_repeated);

	/* from last decomposition */
	const dArrayT& Eigenvalues(void) const;
	const ArrayT<dArrayT>& Eigenvectors(void) const;

	/* rank-1 tensors of the principal stretches from last decomp */
	const dSymMatrixT& Rank1_Principal(int A) const;

	/* return the tensor for the given eigenvalues */
	const dSymMatrixT& EigsToRank2(const dArrayT& eigs);
	const dMatrixT& EigsToRank4(const dSymMatrixT& eigs);

	/* principal spatial tensor for eigenvalue A */
	const dMatrixT& SpatialTensor(const dSymMatrixT& b, int A);

	/* access to fixed forms */
	const dSymMatrixT& I_rank2() const;
	const dMatrixT&    I_rank4() const;
	const dMatrixT&   Ib_rank4() const;
//NOTE: Ib_rank4 is a misleading name. Is really (I_b - b(x)b)

	/* 4th rank mixed index tensor */
	void Set_I4_Tensor3D(const dSymMatrixT& mat, dMatrixT& rank4);
	void Set_I4_Tensor2D(const dSymMatrixT& mat, dMatrixT& rank4);
//TEMP - new IsoVIB needs this, but this is the only place the
//       class displays nsd, so should probably write something
//       different.

	/* find an eigenvalue based using the Rayleigh quotient iteration */
	void RayleighSolve(const dSymMatrixT& rank2, double& eig, dArrayT& vec);	

private:

	/* returns min */
	static double Min(double d1, double d2, double d3);

	/* work routines */
	void SpectralDecomp3D(const dSymMatrixT& rank2, dArrayT& eigs, bool perturb_repeated);

	/* find eigenvectors/values using Jacobi iterations - returns eigenvectors
	 * in the columns of evecs */
	int EigenSystem3D(const dSymMatrixT& matrix, dArrayT& evals, dMatrixT& evecs);

	/* principal spatial tensor for eigenvalue A */
	const dMatrixT& SpatialTensor2D(const dSymMatrixT& b, int A);
	const dMatrixT& SpatialTensor3D(const dSymMatrixT& b, int A);

	/* construct orthogonal basis for 2 repeated roots */
	void SchmidtDecompose(const dSymMatrixT& rank2,
	                 double l3, dSymMatrixT& n2xn2,
	                 double l , dSymMatrixT& n0xn0, dSymMatrixT& n1xn1);
protected:
	
	/* fixed forms */
	dSymMatrixT f_I_Rank2;
	dMatrixT    f_I_Rank4;

	/* spectral decomposition */
	dArrayT fEigs;
	ArrayT<dSymMatrixT> fm; //array of rank 1 matrices	
	
	/* spectral decomp work space */
	dSymMatrixT fm1;
	dSymMatrixT fm2;
	dMatrixT    fEvecMatrix;
	ArrayT<dArrayT> fEvecs;

	/* polar decomp work space */
	dArrayT  fInvEigs;
	dMatrixT fUInv;
	
	/* spatial tensor work space */
	dMatrixT   fSpatTensor;
	dMatrixT   fc_b; //part of spatial tensor dependent only on b
	dMatrixT   fRank4;
	dSymMatrixT fRank2;
};

/* inlines */

/* access to fixed forms */
inline const dSymMatrixT& SpectralDecompT::I_rank2() const { return f_I_Rank2; }
inline const dMatrixT& SpectralDecompT::I_rank4() const { return f_I_Rank4; }
inline const dMatrixT& SpectralDecompT::Ib_rank4() const { return fc_b; }

/* rank-1 tensors of the principal stretches from last decomp */
inline const dSymMatrixT& SpectralDecompT::Rank1_Principal(int A) const { return fm[A]; }

/* eigenvalues from last SpectralDecomp */
inline const dArrayT& SpectralDecompT::Eigenvalues(void) const { return fEigs; }
inline const ArrayT<dArrayT>& SpectralDecompT::Eigenvectors(void) const { return fEvecs; }

#endif /* _SPECTRAL_DECOMP_T_H_ */
