/* $Id: dSymMatrixT.h,v 1.18 2004-03-01 02:41:01 paklein Exp $ */
/* created: paklein (05/24/1996) */
#ifndef _DSYM_MATRIX_T_H_
#define _DSYM_MATRIX_T_H_

/* base class */
#include "dArrayT.h"

namespace Tahoe {

/* forward declarations */
class dMatrixT;

/** interface for a 1D/2D/3D reduced index symmetric matrix stored as 
 * a reduced index vector */
class dSymMatrixT: public dArrayT
{
public:

	/** enumerated type to indicate matrix dimension */
	enum DimensionT {
		kNone = 0,
		k1D = 1,
		k2D = 2,
		k3D = 3,
		k3D_plane = 4 /**< 2D plane components, plus axial out-of-plane component */
	};

	/** integer to dSymMatrixT::DimensionT type conversion operator. Needed
	 * for backward compatibility with code written before the enum was defined */
	static DimensionT int2DimensionT(int i) {
#if __option(extended_errorcheck)
		if (i < kNone || i > k3D_plane) ExceptionT::OutOfRange("dSymMatrixT::int2DimensionT");
#endif
		DimensionT dims[5] = {kNone, k1D, k2D, k3D, k3D_plane};
		return dims[i];
	};

	/** \name constructors */
	/*@{*/
	dSymMatrixT(void);
	explicit dSymMatrixT(DimensionT nsd);
	dSymMatrixT(const dSymMatrixT& source);

	/** construct alias */
	dSymMatrixT(DimensionT nsd, const double* array);
	/*@}*/

	/** dimension the matrix for the number of spatial dimensions. No change 
	 * occurs if the array is already the specified length. The previous contents 
	 * of the array is not preserved. */
	void Dimension(DimensionT nsd);

	/** create shallow matrix. Explicitly set the memory used the matrix.
	 * \param nsd matrix dimension
	 * \param array pointer to memory at least length dSymMatrixT::NumValues(nsd) */
	void Alias(DimensionT nsd, const double* array);

	/** \name assignment operators */
	/*@{*/
	dSymMatrixT& operator=(const dSymMatrixT& RHS);
	dSymMatrixT& operator=(const double value);
	/*@}*/

	/** \name element accessor */
	/*@{*/
	double& operator()(int row, int col);
	const double& operator()(int row, int col) const;
	/*@}*/

	/* returns the number of dimensions in the reduced
	 * index symmetric matrix vectors given the number of spatial
	 * dimensions.  nsd can only be 2 or 3. */
	static int NumValues(DimensionT nsd) {
#if __option(extended_errorcheck)
		if (nsd < kNone || nsd > k3D_plane) ExceptionT::OutOfRange("dSymMatrixT::NumValues");
#endif
		int map[5] = {0, 1, 3, 6, 4};
		return map[nsd];	
	};
	static void ExpandIndex(DimensionT nsd, int dex, int& dex_1, int& dex_2);

	/** \name dimensions */
	/*@{*/
	int Rows(void) const;
	int Cols(void) const;
	/*@}*/
	
	/* I/O operators */
	friend ostream& operator<<(ostream& out, const dSymMatrixT& array);
	friend istream& operator>>(istream& in, dSymMatrixT& array);

	/* return eigenvalues and eigenvectors (in columns) */
	void PrincipalValues(dArrayT& val) const; // will get phased out
	void Eigenvalues(dArrayT& val, bool sort_descending) const;
	void Eigensystem(dArrayT& val, dMatrixT& vec, bool sort_descending) const;

	/* returns the scalar inner product of a tensor with itself
	 * define as:
	 *
	 *             T:T = T_ij T_ij
	 */
	double ScalarProduct(void) const;
	/* 
	 *             A:B = A_ij B_ij
	 */
	double ScalarProduct(const dSymMatrixT& matrix) const;
	
	/* the second invariant */
	double Invariant2(void) const;

	/* matrix determinant */
	double Det(void) const;

	/* returns the trace of *this.  For 2D, the out-of-plane
	 * component is assumed to be zero, ie if *this represents stress,
	 * then plane stress is assumed and if *this is strain then plane
	 * strain is assumed. */
	double Trace(void) const;

	/* identity operations */
	void PlusIdentity(double value = 1.0);
	dSymMatrixT& Identity(double value = 1.0);
	
	/* returns the deviatoric part of *this.
	 *
	 * Note: This function should only be called when working in
	 * 3D, ie. a vector length 6, since there's no way to enforce
	 * Trace(Deviatoric{vector)) = 0, in 2D. */
	dSymMatrixT& Deviatoric(const dSymMatrixT& tensor);
	dSymMatrixT& Deviatoric(void);

	/* matrix inverse */
	dSymMatrixT& Inverse(const dSymMatrixT& matrix);
	dSymMatrixT& Inverse(void);

	/* take the symmetric part of the matrix */
	dSymMatrixT& Symmetrize(const dMatrixT& matrix);

	/* reduced index <-> matrix transformations */	
	void ToMatrix(dMatrixT& matrix) const;
	dSymMatrixT& FromMatrix(const dMatrixT& matrix);

	/* 2D <-> 3D translations.  Return references to *this */
	dSymMatrixT& ExpandFrom2D(const dSymMatrixT& vec2D); /* assumed plane strain */
	dSymMatrixT& ReduceFrom3D(const dSymMatrixT& vec3D);

	/** ? <-> ? translations */
	dSymMatrixT& Translate(const dSymMatrixT& matrix);

	/** outer product */
	void Outer(const dArrayT& v, double scale = 1.0);

	/* matrix-matrix multiplication */
	void MultAB(const dSymMatrixT& A, const dSymMatrixT& B); /* A and B must commute */
	void MultATA(const dMatrixT& A);
	void MultAAT(const dMatrixT& A);

	/* matrix-matrix-matrix operations, ie. tensor basis transformations */
	void MultQBQT(const dMatrixT& Q, const dSymMatrixT& B);
	void MultQTBQ(const dMatrixT& Q, const dSymMatrixT& B);

	/* matrix-vector multiplication */
	void Multx(const dArrayT& x, dArrayT& b) const;

	/* vector-matrix-vector product */
	double MultmBn(const dArrayT& m, const dArrayT& n) const;

	/* symmetric rank-4 - rank-2 contraction */
	void A_ijkl_B_kl(const dMatrixT& A, const dSymMatrixT& B);
	void A_ijkl_B_ij(const dMatrixT& A, const dSymMatrixT& B);

	/* symmetric rank-(2-4-2) contraction with this */
	double B_ij_A_ijkl_B_kl(const dMatrixT& A) const;

	/** scale off-diagonal value by the given factor */
	void ScaleOffDiagonal(double factor);

	/** \name deprecated methods */
	/*@{*/
	explicit dSymMatrixT(int nsd) ;
	dSymMatrixT(int nsd, double* array);
	void Dimension(int nsd) { Dimension(int2DimensionT(nsd)); };
	void Allocate(int nsd) { Dimension(nsd); };
	void Set(int nsd, double* array);
	void Set(DimensionT nsd, double* array);
	static int NumValues(int nsd) { return NumValues(int2DimensionT(nsd)); };
	static void ExpandIndex(int nsd, int dex, int& dex_1, int& dex_2) { 
	  ExpandIndex(int2DimensionT(nsd), dex, dex_1, dex_2); };
	/*@}*/

private:

	/* iterative eigensystem routines */
	int Eigenvalues3D(dArrayT& evals, bool sort_descending, int max_iterations = 15) const;
	int Eigensystem3D(dArrayT& evals, dMatrixT& evecs, bool sort_descending,
		int max_iterations = 15) const; // append algorithm name
		
	/* closed form eigenvalues - unstable for repeated roots */
	void Eigenvalues3D_Cardano(dArrayT& evals) const;
	
private:

	DimensionT fNumSD;	
};

/* inlines */

/* assigment operators */
inline dSymMatrixT& dSymMatrixT::operator=(const dSymMatrixT& RHS)
{
	/* must be same dimension or nsd must be zero */
	if (fNumSD != 0 && fNumSD != RHS.fNumSD)
		ExceptionT::SizeMismatch("dSymMatrixT::operator=");
	else
		fNumSD = RHS.fNumSD;

	/* inherited */
	dArrayT::operator=(RHS);
	return *this;
}

inline dSymMatrixT& dSymMatrixT::operator=(const double value)
{
	/* inherited */
	dArrayT::operator=(value);
	return *this;
}

/* dimensions */
inline int dSymMatrixT::Rows(void) const { 
  return (fNumSD<3) ? fNumSD:3; 
}
inline int dSymMatrixT::Cols(void) const { 
return (fNumSD<3) ? fNumSD:3; 
}
inline dSymMatrixT& dSymMatrixT::Deviatoric(void) { return Deviatoric(*this); }
inline dSymMatrixT& dSymMatrixT::Inverse(void) { return Inverse(*this); }

/* scale off-diagonal value by the given factor */
inline void dSymMatrixT::ScaleOffDiagonal(double factor)
{
	if (fNumSD == 2 || fNumSD == dSymMatrixT::k3D_plane)
		fArray[2] *= factor;
	else if (fNumSD == 3) {
		fArray[3] *= factor;
		fArray[4] *= factor;	
		fArray[5] *= factor;	
	}
}

/* set fields */
inline void dSymMatrixT::Alias(DimensionT nsd, const double* array)
{
	fNumSD = int2DimensionT(nsd);
#if __option(extended_errorcheck)
	if (fNumSD < 1 || fNumSD > 4) 
		ExceptionT::GeneralFail("dSymMatrixT::Alias", "invalid dimension %d", nsd);
#endif
	/* inherited */
	dArrayT::Alias(NumValues(fNumSD), array);
}

inline void dSymMatrixT::Set(DimensionT nsd, double* array)
{
	fNumSD = int2DimensionT(nsd);
#if __option(extended_errorcheck)
	if (fNumSD < 1 || fNumSD > 4) 
		ExceptionT::GeneralFail("dSymMatrixT::Set", "invalid dimension %d", nsd);
#endif
	/* inherited */
	dArrayT::Set(NumValues(fNumSD), array);
}

inline void dSymMatrixT::Set(int nsd, double* array)
{
	Set(int2DimensionT(nsd), array); 
};

/* const element accessor */
inline double& dSymMatrixT::operator()(int row, int col) {
	const double& a = ((const dSymMatrixT*) this)->operator()(row,col);
	return const_cast<double&>(a);
}

} /* namespace Tahoe */

#endif /* _DSYM_MATRIX_T_H_ */
