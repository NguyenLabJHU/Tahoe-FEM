/* $Id: FSFiberOptimize_MatT.h,v 1.1 2009-04-23 03:03:50 thao Exp $ */
/* created: paklein (06/09/1997) */
#ifndef _FSFiberOptimize_MatT_H_
#define _FSFiberOptimize_MatT_H_

/* direct members */
#include "dSymMatrixT.h"
#include "StringT.h"
#include "dArrayT.h"
#include "dArray2DT.h"

namespace Tahoe {

/* forward declarations */
class FSFiberMatSupportT;

/** defines the interface for small strain continuum materials */
class FSFiberOptimize_MatT
{
public:

	/** constructor */
	FSFiberOptimize_MatT(void);

	/*destructors*/
	~FSFiberOptimize_MatT(void);
	
	/** set the material support or pass NULL to clear */
	void SetSupport(FSFiberMatSupportT* support);
		
	/*computes the derivative of the stress tensor with respect to the parameters*/
	/*Implementation in the base class uses numerical finite difference */
	virtual const double constraint(void) = 0;
	
	virtual const dArrayT& constraint_grad(void) = 0;
	
	virtual const dArray2DT& ds_ij_dlambda_q(void) = 0;

	virtual const ArrayT<StringT>& ParamLabels(void) const {return flabels;};
			
protected:

	/** material support */
	FSFiberMatSupportT* fSupport;
	dArray2DT fParamGrads;
	dArrayT fconstraint_grad;
	
	dArrayT fparams;
	ArrayT<StringT> flabels;
	
	
//	const SSMatSupportT* fSSMatSupport;
	
};


} // namespace Tahoe 
#endif /* _SS_STRUCT_MAT_T_H_ */
