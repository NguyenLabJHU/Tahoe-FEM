/* $Id: SSVisco_Optimize.h,v 1.1 2009-04-23 03:03:52 thao Exp $ */
/* created: TDN (5/31/2001) */
#ifndef _SS_VISCO_OPTO_H_
#define _SS_VISCO_OPTO_H_
 
#include "SSOptimize_MatT.h"
#include "SSViscoelasticityT.h"

#include "dSymMatrixT.h"

namespace Tahoe {

/* forward declarations */

/** small strain linear viscoelastic constitutive law */
class SSVisco_Optimize: public SSOptimize_MatT, public SSViscoelasticityT
{
	public:

	/** constructor */
	SSVisco_Optimize(void);
	
	virtual bool Need_Strain_last(void) const {return true;};

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
	
	virtual const dArrayT& constraint_grad(void){return(fconstraint_grad);}

	virtual const dArray2DT& ds_ij_dlambda_q(void);
	
	private:
	int fnum_params;
};

} // namespace Tahoe 
#endif /*_SS_VISCO_H_*/
