/* $Id: TriDiagdMatrixT.h,v 1.4 2002-10-20 22:38:54 paklein Exp $ */
/* created: paklein (01/15/1998)                                          */
/* Triadiagonal matrix with Gauss elimination.                            */

#ifndef _TRIDIAG_DMATRIX_T_H_
#define _TRIDIAG_DMATRIX_T_H_

/* Environmental */
#include "Environment.h"

/* base class */
#include "nArrayT.h"

namespace Tahoe {

/* forward declarations */
class dArrayT;

class TriDiagdMatrixT: public nArrayT<double>
{
public:

	/* constructor */
	TriDiagdMatrixT(int rows);

	/* assigning values to a row - first and last row have overrun */
	void SetRow(int row, double L, double D, double R);
	void AddToRow(int row, double L, double D, double R);

	/* dimensions */
	int Rows(void) const;

	/* Gaussian elimination with the given RHS vector or RHS matrix */
	void LinearSolve(dArrayT& RHS);
		
private:
	
	int	fRows;

	/* pointers to each of the bands */
	double *pL;
	double *pD;
	double *pR;
};

/* in-lines */

/* row accessor */
inline void TriDiagdMatrixT::SetRow(int row, double L, double D, double R)
{
/* range checking */
#if __option (extended_errorcheck)
	if (row < 0 || row >= fRows) throw ExceptionT::kOutOfRange;
#endif

	pL[row] = L;	
	pD[row] = D;	
	pR[row] = R;	
};

inline void TriDiagdMatrixT::AddToRow(int row, double L, double D, double R)
{
/* range checking */
#if __option (extended_errorcheck)
	if (row < 0 || row >= fRows) throw ExceptionT::kOutOfRange;
#endif

	pL[row] += L;	
	pD[row] += D;	
	pR[row] += R;	
};

/* dimensions */
inline int TriDiagdMatrixT::Rows(void) const { return (fRows); }

} // namespace Tahoe 
#endif /* _TRIDIAG_DMATRIX_T_H_ */
