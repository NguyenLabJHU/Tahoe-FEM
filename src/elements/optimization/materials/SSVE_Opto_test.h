/* $Id: SSVE_Opto_test.h,v 1.1 2009-04-23 03:03:52 thao Exp $ */
/* created: TDN (5/31/2001) */
#ifndef _SSVE_OPTO_TEST_H_
#define _SSVE_OPTO_TEST_H_
 
#include "SSOptimize_MatT.h"
#include "SSVE_test.h"

#include "dSymMatrixT.h"

namespace Tahoe {

/* forward declarations */

/** small strain linear viscoelastic constitutive law */
class SSVE_Opto_test: public SSOptimize_MatT, public SSVE_test
{
	public:

	/** constructor */
	SSVE_Opto_test(void);
	
	virtual bool Need_Strain_last(void) const {return true;};

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** information about subordinate parameter lists */
	virtual void DefineParameters(ParameterListT& list) const;

	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/
		
	virtual const double constraint(void);

	virtual const dArrayT& constraint_grad(void);
		
	virtual const dArray2DT& ds_ij_dlambda_q(void);
	
	private:
	int fnum_params;
};

} // namespace Tahoe 
#endif /*_SSVE_OPTO_TEST_H_*/
