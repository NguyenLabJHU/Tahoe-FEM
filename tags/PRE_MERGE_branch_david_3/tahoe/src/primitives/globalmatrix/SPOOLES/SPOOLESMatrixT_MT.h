/* $Id: SPOOLESMatrixT_MT.h,v 1.1 2005-04-05 16:08:59 paklein Exp $ */
#ifndef _SPOOLES_MATRIX_T_MT_H_
#define _SPOOLES_MATRIX_T_MT_H_

/* base class */
#include "MSRMatrixT.h"

/* library support options */
#ifdef __SPOOLES_MT__

/* direct members */
#include "AutoArrayT.h"
#include "iArrayT.h"
#include "dArrayT.h"
#include "dMatrixT.h"
#include "nVariMatrixT.h"

namespace Tahoe {

/* forward declarations */
class MSRBuilderT;

/** interface to the multi-threaded SPOOLES sparse, direct linear solver. Solution driver 
 * supports multiple solves with the same factorized matrix. Sparsity and symbolic factorization
 * is not preserved between matricies. All previous information about the previous
 * matrix is delete with the call to SPOOLESMatrixT_MT::Clear. */
class SPOOLESMatrixT_MT: public MSRMatrixT
{
public:

	/* constuctor */
	SPOOLESMatrixT_MT(ostream& out, int check_code, bool symmetric,
		bool pivoting, int num_threads);

	/* copy constructor */
	SPOOLESMatrixT_MT(const SPOOLESMatrixT_MT& source);

	/** destructor */
	virtual ~SPOOLESMatrixT_MT(void);

	/** SPOOLESMatrixT_MT::Solve does preserve the data in the matrix */
	virtual bool SolvePreservesData(void) const { return true; };	  

	/** clear values for next assembly */
	virtual void Clear(void);

	/** return the form of the matrix */
	virtual GlobalT::SystemTypeT MatrixType(void) const { return GlobalT::kNonSymmetric; };

	/** assignment operatpr */
	SPOOLESMatrixT_MT& operator=(const SPOOLESMatrixT_MT& rhs);

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
	
	/** number of threads to use in calls to SPOOLES_MT library */
	int fNumThreads;
};

} // namespace Tahoe 
#endif /*__SPOOLES_MT__ */
#endif /* _SPOOLES_MATRIX_T_MT_H_ */




