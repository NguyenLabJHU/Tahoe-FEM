/* $Header: /home/regueiro/tahoe_cloudforge_repo_snapshots/development/src/elements/fluid_element/materials/FluidMatListT.h,v 1.1 2006-07-13 17:57:28 a-kopacz Exp $ */
/* created: tdnguye (07/12/2006) */
#ifndef _FLUID_MAT_LIST_T_H_
#define _FLUID_MAT_LIST_T_H_

/* base class */
#include "MaterialListT.h"

namespace Tahoe {

/* forward declarations */
class FluidMatSupportT;
class FluidMaterialT;

/** list of materials for diffusion analysis */
class FluidMatListT: public MaterialListT
{
public:

	/** enum defining material types */
	enum TypeT {
        kLinear = 1,
     kNonLinear = 2};

	/** constructors */
	FluidMatListT(int length, const FluidMatSupportT& support);
	FluidMatListT(void);

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** information about subordinate parameter lists */
	virtual void DefineSubs(SubListT& sub_list) const;

	/** return the description of the given inline subordinate parameter list */
	virtual void DefineInlineSub(const StringT& name, ParameterListT::ListOrderT& order, 
		SubListT& sub_lists) const;

	/** a pointer to the ParameterInterfaceT of the given subordinate */
	virtual ParameterInterfaceT* NewSub(const StringT& name) const;

	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/

	/** construct the specified material or NULL if the request cannot be completed */
	FluidMaterialT* NewFluidMaterial(const StringT& name) const;
	
private:

	/** support for diffusion materials */
	const FluidMatSupportT* fFluidMatSupport;
};

} // namespace Tahoe 
#endif /* _FLUID_MAT_LIST_T_H_ */
