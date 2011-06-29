/* $Id: CSEIsoT.h,v 1.11 2004-07-15 08:25:57 paklein Exp $ */
/* created: paklein (11/19/1997) */

#ifndef _CSE_ISO_T_H_
#define _CSE_ISO_T_H_

/* base class */
#include "CSEBaseT.h"

/* direct members */
#include "pArrayT.h"

namespace Tahoe {

/* forward declarations */
class C1FunctionT;

/** Cohesive surface elements with scalar traction potentials,
 * i.e., the traction potential is a function of the gap magnitude,
 * or effective gap magnitude only. */
class CSEIsoT: public CSEBaseT
{
public:

#ifndef _FRACTURE_INTERFACE_LIBRARY_
	/** constructor */
	CSEIsoT(const ElementSupportT& support);
#else
	CSEIsoT(ElementSupportT& support);
#endif

	/** form of tangent matrix */
	virtual GlobalT::SystemTypeT TangentType(void) const;

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** information about subordinate parameter lists */
	virtual void DefineSubs(SubListT& sub_list) const;

	/** return the description of the given inline subordinate parameter list */
	virtual void DefineInlineSub(const StringT& name, ParameterListT::ListOrderT& order, 
		SubListT& sub_lists) const;

	/** a pointer to the ParameterInterfaceT */
	virtual ParameterInterfaceT* NewSub(const StringT& name) const;

	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/

protected:

	/** tangent matrix */
	virtual void LHSDriver(GlobalT::SystemTypeT sys_type);

	/** force vector */
	virtual void RHSDriver(void);

	/** compute output values */
	virtual void ComputeOutput(const iArrayT& n_codes, dArray2DT& n_values,
		const iArrayT& e_codes, dArray2DT& e_values);
	
protected:

	/** cohesive surface potentials */
	pArrayT<C1FunctionT*> fSurfPots;
};

} /* namespace Tahoe */

#endif /* _CSE_ISO_T_H_ */
