/* $Id: SPOOLESMatrixT.h,v 1.15 2004-03-14 02:51:35 paklein Exp $ */
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

/** interface to SPOOLES sparse, direct linear solver */
class SPOOLESMatrixT: public MSRMatrixT
{
public:

	/* constuctor */
	SPOOLESMatrixT(ostream& out, int check_code, bool symmetric,
		bool pivoting);

	/* copy constructor */
	SPOOLESMatrixT(const SPOOLESMatrixT& source);

	/** SPOOLESMatrixT::Solve does preserve the data in the matrix */
	virtual bool SolvePreservesData(void) const { return true; };	  

	/** return the form of the matrix */
	virtual GlobalT::SystemTypeT MatrixType(void) const { return GlobalT::kNonSymmetric; };
	
	/** assignment operator */
	virtual GlobalMatrixT& operator=(const SPOOLESMatrixT& rhs);

	/** assignment operator */
	virtual GlobalMatrixT& operator=(const GlobalMatrixT& rhs);
	
	/** return a clone of self. Caller is responsible for disposing of the matrix */
	virtual GlobalMatrixT* Clone(void) const;

protected:
	
	/* determine new search direction and put the results in result */
	virtual void BackSubstitute(dArrayT& result);

protected:

	/** parameters */
	bool fPivoting;	
};

} // namespace Tahoe 
#endif /*__SPOOLES__ */
#endif /* _SPOOLES_MATRIX_T_H_ */




