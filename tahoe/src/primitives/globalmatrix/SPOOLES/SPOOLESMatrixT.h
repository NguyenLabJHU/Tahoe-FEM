/* $Id: SPOOLESMatrixT.h,v 1.18 2005-01-07 21:23:02 paklein Exp $ */
/* created: paklein (09/13/2000) */
#ifndef _SPOOLES_MATRIX_T_H_
#define _SPOOLES_MATRIX_T_H_

/* base class */
#include "MSRMatrixT.h"

/* library support options */
#ifdef __SPOOLES__

/* direct members */
#include "AutoArrayT.h"
#include "iArrayT.h"
#include "dArrayT.h"
#include "dMatrixT.h"
#include "nVariMatrixT.h"

namespace Tahoe {

/* forward declarations */
class MSRBuilderT;

/** interface to SPOOLES sparse, direct linear solver. Solution driver supports
 * multiple solves with the same factorized matrix. Sparsity and symbolic factorization
 * is not preserved between matricies. All previous information about the previous
 * matrix is delete with the call to SPOOLESMatrixT::Clear. */
class SPOOLESMatrixT: public MSRMatrixT
{
public:

	/* constuctor */
	SPOOLESMatrixT(ostream& out, int check_code, bool symmetric,
		bool pivoting);

	/* copy constructor */
	SPOOLESMatrixT(const SPOOLESMatrixT& source);

	/** destructor */
	virtual ~SPOOLESMatrixT(void);

	/** SPOOLESMatrixT::Solve does preserve the data in the matrix */
	virtual bool SolvePreservesData(void) const { return true; };	  

	/** clear values for next assembly */
	virtual void Clear(void);

	/** return the form of the matrix */
	virtual GlobalT::SystemTypeT MatrixType(void) const { return GlobalT::kNonSymmetric; };

	/** assignment operatpr */
	SPOOLESMatrixT& operator=(const SPOOLESMatrixT& rhs);

	/** return a clone of self. Caller is responsible for disposing of the matrix */
	virtual GlobalMatrixT* Clone(void) const;

protected:

	/** precondition matrix */
	virtual void Factorize(void);

	/** determine new search direction and put the results in result */
	virtual void BackSubstitute(dArrayT& result);

protected:

	/** parameters */
	bool fPivoting;	

	/** data for repeated solves with the same matrix (not for MT) */
	void* pLU_dat;

	/** runtime flag */
	bool fIsFactorized;
};

} // namespace Tahoe 
#endif /*__SPOOLES__ */
#endif /* _SPOOLES_MATRIX_T_H_ */




