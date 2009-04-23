/* $Id: SSKStV_Optimize.h,v 1.1 2009-04-23 03:03:51 thao Exp $ */
/* created: paklein (06/10/1997) */
#ifndef _SS_KSTV_OPTIMIZE_H_
#define _SS_KSTV_OPTIMIZE_H_

/* base classes */
#include "SSOptimize_MatT.h"
#include "IsotropicT.h"
#include "SSKStV.h"

namespace Tahoe {

class SSKStV_Optimize: public SSOptimize_MatT, public SSKStV
{
public:

	/* constructor */
	SSKStV_Optimize(void);

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** information about subordinate parameter lists */
	virtual void DefineSubs(SubListT& sub_list) const;

	/** a pointer to the ParameterInterfaceT of the given subordinate */
	virtual ParameterInterfaceT* NewSub(const StringT& name) const;

	/** information about subordinate parameter lists */
	virtual void DefineParameters(ParameterListT& list) const;

	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/
	
	virtual const double constraint(void){return 0.0;};
	
	virtual const dArrayT& constraint_grad(void){fconstraint_grad = 0.0; return(fconstraint_grad);}

	virtual const dArray2DT& ds_ij_dlambda_q(void);
			
private:
	double fE, fnu;	

};

} // namespace Tahoe 
#endif /* _SS_KSTV_H_ */
