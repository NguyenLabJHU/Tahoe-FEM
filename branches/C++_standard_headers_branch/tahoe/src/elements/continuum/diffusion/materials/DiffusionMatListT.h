/* $Id: DiffusionMatListT.h,v 1.9 2004-07-15 08:26:22 paklein Exp $ */
/* created: paklein (10/02/1999) */
#ifndef _DIFFUSE_MAT_LIST_T_H_
#define _DIFFUSE_MAT_LIST_T_H_

/* base class */
#include "MaterialListT.h"

namespace Tahoe {

/* forward declarations */
class DiffusionMatSupportT;
class DiffusionMaterialT;

/** list of materials for diffusion analysis */
class DiffusionMatListT: public MaterialListT
{
public:

	/** enum defining material types */
	enum TypeT {
        kLinear = 1,
     kNonLinear = 2};

	/** constructors */
	DiffusionMatListT(int length, const DiffusionMatSupportT& support);
	DiffusionMatListT(void);

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
	DiffusionMaterialT* NewDiffusionMaterial(const StringT& name) const;
	
private:

	/** support for diffusion materials */
	const DiffusionMatSupportT* fDiffusionMatSupport;
};

} // namespace Tahoe 
#endif /* _DIFFUSE_MAT_LIST_T_H_ */
