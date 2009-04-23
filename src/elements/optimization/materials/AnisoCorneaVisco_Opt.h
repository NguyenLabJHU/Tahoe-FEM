/* $Id: AnisoCorneaVisco_Opt.h,v 1.1 2009-04-23 03:03:50 thao Exp $ */
/* created: TDN (5/31/2001) */
#ifndef _SS_VISCO_OPTO_H_
#define _SS_VISCO_OPTO_H_
 
#include "FSFiberOptimize_MatT.h"
#include "AnisoCorneaVisco.h"

#include "dSymMatrixT.h"

namespace Tahoe {

/* forward declarations */

/** small strain linear viscoelastic constitutive law */
class AnisoCorneaVisco_Opt: public FSFiberOptimize_MatT, public AnisoCorneaVisco
{
	public:

	/** constructor */
	AnisoCorneaVisco_Opt(void);
	
	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** information about subordinate parameter lists */
	virtual void DefineParameters(ParameterListT& list) const;

	/** information about subordinate parameter lists */
	virtual void DefineSubs(SubListT& sub_list) const;

	/** a pointer to the ParameterInterfaceT of the given subordinate */
	virtual ParameterInterfaceT* NewSub(const StringT& name) const;

	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/
		 
	virtual const double constraint(void){return 0.0;};
	
	virtual const dArrayT& constraint_grad(void){return(fconstraint_grad);}

	virtual const dArray2DT& ds_ij_dlambda_q(void);
	
	private:
	void Construct(void);
		
	private:
	int fnum_params;
	dArrayT fI4e;
	
	dSymMatrixT fb;
	
	ArrayT<dArrayT> fjac_dk; // for an inhomogeneous material
	ArrayT<dArrayT> fjac_dphi; // for an inhomogeneous material
	
	/*workspace*/
	dSymMatrixT fdS;
	
	dSymMatrixT fdsf_dalpha;
	dSymMatrixT fdsf_dbeta;
	dSymMatrixT fdsf_deta0;
	dSymMatrixT fdsf_dtau0;
	dSymMatrixT fdsf_dk;
	dSymMatrixT fdsf_dphi;

	
};

} // namespace Tahoe 
#endif /*_SS_VISCO_H_*/
