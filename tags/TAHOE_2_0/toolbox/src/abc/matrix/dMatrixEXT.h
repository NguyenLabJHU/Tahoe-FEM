/* $Id: dMatrixEXT.h,v 1.11 2003-11-21 22:41:36 paklein Exp $ */
/* created: paklein (03/06/1998) */

#ifndef _DMATRIXEX_T_H_
#define _DMATRIXEX_T_H_

/* base class */
#include "dMatrixT.h"

/* direct members */
#include "dArrayT.h"
#include "iArrayT.h"
#include "dTensor4DT.h"

namespace Tahoe {

/** interface for dMatrixT plus special matrix functions */
class dMatrixEXT: public dMatrixT
{
public:

	/* constructor */
	dMatrixEXT(void);
	explicit dMatrixEXT(int squaredim);
	dMatrixEXT(int squaredim, const double* p);

	/** dimension to a square matrix */
	void Dimension(int squaredim);

	/** \deprecated replaced by dMatrixEXT::Dimension on 02/13/2002 */
	void Allocate(int squaredim) { Dimension(squaredim); };

	/* diagonalize (using symmetric QR algorithm).
	 * assumes the matrix is symmetric. Returns the
	 * number of QR iterations needed to diagonalize */
	int Diagonalize(dArrayT& eigs);

	/* return the {eigenvalue,eigenvector} pair corresponding
	 * to the approximate eigenvalue that is passed in */
	void Eigenvector(double& eig_guess, dArrayT& eigenvector) const;

	/* returns eigenvalues of a general matrix by first putting it in 
	 * Hessian form */
	void eigvalfinder (dMatrixEXT& matrix, dArrayT& realev, dArrayT& imev);
    void eigenvalue3x3(dMatrixEXT& J, dArrayT& reroot, dArrayT& imroot);
	void eigenvector3x3(dMatrixEXT& J, double value, int numvector,  dArrayT& vector, dArrayT& vector2, dArrayT& vector3);
	/*forms acoustic tensor from rank 4 tangent modulus, normal */
	//void formacoustictensor(dMatrixEXT& A, double C [3] [3] [3] [3], dArrayT& normal);
	void formacoustictensor(dMatrixEXT& A, dTensor4DT& C, dArrayT& normal);


	/** generate singular value decomposition of *this = U*W*V^T. 
	 * \param U return matrix: [n_rows] x [n_cols]
	 * \param W diagonal matrix of singular values: [n_cols]
	 * \param V square return matrix: [n_cols] x [n_cols] 
	 * \param threshold for singular values allows to be nonzero, i.e., upper bound relative
	 *        to the maximum singular value
	 * \param max_its maximum number of iterations, 30 by default */
	void Compute_SVD(dMatrixT& U, dArrayT& W, dMatrixT& V, double threshold, int max_its = 30) const;

	/** back substitute given a decomposition computed with dMatrixEXT::Compute_SVD.
	 * \param U return matrix: [n_rows] x [n_cols]
	 * \param W diagonal matrix of singular values: [n_cols]
	 * \param V square return matrix: [n_cols] x [n_cols] 
	 * \param RHS goes in as the RHS vector and returns as the solution
 	 * \note Aside from dimension checks, this function doesn't really involve
 	 *       the matrix stored in *this. */
	void BackSubstitute_SVD(const dMatrixT& U, const dArrayT& W, const dMatrixT& V, dArrayT& RHS) const;

	/* assignment operator. Operator will re-dimension matrix as needed.
	 * \param RHS source
	 * \return reference to *this */
	dMatrixT& operator=(const dMatrixT& RHS) { return dMatrixT::operator=(RHS); };

	/* assignment operator. Set all entries in the matrix to value
	 * \return reference to *this */
	dMatrixT& operator=(const double value) { return dMatrixT::operator=(value); };

private:

	/* compute tridiagonal decomposition */
	void TriDiagonalForm(void);
	
	/* peform single implicit, symmetric QR step (with
	 * Wilkinson shift) upto the maxkk diagonal element.
	 * assumes (this) is symmetric and tridiagonal */
	void SymmetricQRStep(int maxkk);
	
	/* perform Givens rotation. assumes (this) is tridiagonal
	 * and symmetric, ie. apply TriDiagonalForm(), first.
	 *
	 *     k  : start point (0...dimension-2)
	 *     c,s: Givens scalars
	 *     z  : previous -> next non-zero off-diagonal value
	 *          (previous ignored for k = 0)
	 */
	void ApplyGivens(int k, double c, double s, double& z, int restart = 0);
	
	/* returns the Householder vector (5.1.2-3) in x2v, where the
	 * Householder reflection is given by:
	 *
	 *           P = I - 2/(v.v) v (x) v
	 *
	 * where beta is then defined as 2/(v.v)
	 *
	 * such that P.x = sqrt(x.x) e_1
	 */
	void HouseholderVector(const double* x, double* v, double& beta, int length) const;
	
	/* returns the Givens scalars { cos(t),sin(t) } to transform
	 * {a,b}^T into {r,0} */
	void GivensScalars(double a, double b, double& c, double& s) const;

	/* private */
	double Dot(const double* v1, const double* v2, int length) const;

	/* pull up tridiagonal matrix overwriting the diagonal at i*/
	void PullUpTri(int i, int length);

	/* Numerical Recipies */
	double pythag(double a, double b) const;
	int tqli(double d[], double e[], int n);

	/* Numerical recipe, puts general matrix in Hessenberg form */
	void elmhes(dMatrixEXT& a,int n);

	/* numerical recipe, finds eigenvalues of Hessenberg matrix
	 *real part of each value is in wr, imaginary in wi */
	void hqr(dMatrixEXT& a, int n, dArrayT& wr, dArrayT& wi);

	/** compute SVD */
	int svdcmp(double* a, int m, int n, double* w, double* v, double* rv1, int max_its) const;

	/** back substitute SVD */
	void svbksb(const double* u, const double* w, const double* v, int m, int n, 
		const double* b, double* x, double* tmp) const;

private:

	/* work vectors length fRows (== fCols) */
	double* v1;
	double* v2;

	/* workspace */
	dArrayT	fworkspace;  	
};

/* inlines */

/* private */
inline double dMatrixEXT::Dot(const double* v1, const double* v2, int length) const
{
	double prod = 0.0;
	for (int i = 0; i < length; i++)
		prod += (*v1++)*(*v2++);
	return prod;
}

} // namespace Tahoe 
#endif /* _DMATRIXEX_T_H_ */
