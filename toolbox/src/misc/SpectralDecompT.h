/* $Id: SpectralDecompT.h,v 1.8 2002-07-05 22:26:31 paklein Exp $ */
/* created: paklein (11/09/1997)                                          */
/* Spectral decomposition solver                                          */

#ifndef _SPECTRAL_DECOMP_T_H_
#define _SPECTRAL_DECOMP_T_H_

/* direct members */
#include "dMatrixT.h"
#include "dSymMatrixT.h"

namespace Tahoe {

/* forward declarations */
class dArrayT;

class SpectralDecompT
{
public:

//DEV - soon to be gone
void SpectralDecomp_new(const dSymMatrixT& rank2, bool perturb_repeated) { SpectralDecomp_Jacobi(rank2, perturb_repeated); };
void DecompAndModPrep(const dSymMatrixT& rank2, bool perturb_repeated)
{
	SpectralDecomp_Jacobi(rank2, perturb_repeated);
	ModulusPrep(rank2);
};
//DEV

	/** constructor. \param nsd number of spatial dimensions */
	SpectralDecompT(int nsd);
		
	/** spectral decomposition.
	 * decomposition computed in closed form using Cardano's equations
     * \param rank2 tensor to decompose 
     * \param perturb_repeated flag to control perturbation of repeated roots */
	void SpectralDecomp(const dSymMatrixT& rank2, bool perturb_repeated);

	/** spectral decomposition.
	 * decomposition computed in closed form using Cardano's equations
     * \param rank2 tensor to decompose
     * \param eigs given eigenvalues of rank2
     * \param perturb_repeated flag to control perturbation of repeated roots */
	void SpectralDecomp(const dSymMatrixT& rank2, const dArrayT& eigs, bool perturb_repeated);

	/** spectral decomposition.
	 * decomposition computed in using Jacobi iterations 
     * \param rank2 tensor to decompose
     * \param perturb_repeated flag to control perturbation of repeated roots
     * after computing the decomposition "exactly" */
	void SpectralDecomp_Jacobi(const dSymMatrixT& rank2, bool perturb_repeated);

	/** compute the polar decomposition of F = RU.
	 * \param F given deformation gradient
	 * \param R rotation tensor satisfying R^T R = 1
	 * \param U stretching part of the decomposition
     * \param perturb_repeated flag to control perturbation of repeated roots */
	void PolarDecomp(const dMatrixT& F, dMatrixT& R, dSymMatrixT& U, bool perturb_repeated);

	/** eigenvalues from the most recent decomposition */
	const dArrayT& Eigenvalues(void) const;

	/** eigenvectors from the most recent decomposition */
	const ArrayT<dArrayT>& Eigenvectors(void) const;

	/** eigenmatrix from the most recent decomposition.
	 * matrix has eigen vectors in columns */
	const dMatrixT& Eigenmatrix(void) const;

	/** perturb repeated values.
	 * \param roots source array must be length 2 or 3
	 * \return true if any values were perturbed, false otherwise */
	bool PerturbRepeated(dArrayT& values) const;
	
	/** perturb internal roots.
	 * \return true if roots were perturbed, false otherwise */
	bool PerturbRoots(void) { return PerturbRepeated(fEigs); };

	/* rank-1 tensors of the principal stretches from last decomp */
	const dSymMatrixT& Rank1_Principal(int A) const;

	/** tensor reconstruction.
	 * using the eigenvectors from the most recent decomposition
	 * \param eigs eigenvalues for the reconstruction */
	const dSymMatrixT& EigsToRank2(const dArrayT& eigs);

	/** tensor reconstruction.
	 * using the eigenvectors from the most recent decomposition
	 * \param eigs derivative of eigenvalues for the reconstruction */
	const dMatrixT& EigsToRank4(const dSymMatrixT& eigs);
	const dMatrixT& NonSymEigsToRank4(const dMatrixT& eigs);

	/** set rank 4 spatial tensor. needed for computing moduli using 
	 * Simo's spectral formulation. must be called before any call to 
	 * SpatialTensor and after a spectral decomposition of rank2.
	 * \param rank2 source tensor */
	void ModulusPrep(const dSymMatrixT& rank2);

	/** principal spatial tensor.
	 * \param b source stretch tensor
	 * \param A eigenvalue number
	 * \note Derivation of this tensor assumes the roots are distinct.
	 * Repeated roots must be perturbed to avoid division by zero. */
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
inline const dMatrixT& SpectralDecompT::Eigenmatrix(void) const { return fEvecMatrix; }

} // namespace Tahoe 
#endif /* _SPECTRAL_DECOMP_T_H_ */
