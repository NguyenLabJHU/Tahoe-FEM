/* $Id: dMatrixEXT.h,v 1.1.1.1 2001-01-25 20:56:23 paklein Exp $ */
/* created: paklein (03/06/1998)                                          */
/* dMatrixT plus special matrix functions                                 */

#ifndef _DMATRIXEX_T_H_
#define _DMATRIXEX_T_H_

/* base class */
#include "dMatrixT.h"

/* direct members */
#include "dArrayT.h"
#include "iArrayT.h"

class dMatrixEXT: public dMatrixT
{
public:

	/* constructor */
	dMatrixEXT(void);
	dMatrixEXT(int squaredim);
	dMatrixEXT(int squaredim, double* p);

	/* post constructor (re-)dimensioning */
	void Allocate(int squaredim);

	/* diagonalize (using symmetric QR algorithm).
	 * assumes the matrix is symmetric. Returns the
	 * number of QR iterations needed to diagonalize */
	int Diagonalize(dArrayT& eigs);

	/* return the {eigenvalue,eigenvector} pair corresponding
	 * to the approximate eigenvalue that is passed in */
	void Eigenvector(double& eig_guess, dArrayT& eigenvector) const;

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
	double pythag(double a, double b);
	int tqli(double d[], double e[], int n);

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

#endif /* _DMATRIXEX_T_H_ */
