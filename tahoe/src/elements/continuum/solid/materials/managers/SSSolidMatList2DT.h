/* $Id: SSSolidMatList2DT.h,v 1.1.2.1 2004-01-21 19:10:18 paklein Exp $ */
/* created: paklein (02/14/1997) */
#ifndef _SS_MATLIST_2D_T_H_
#define _SS_MATLIST_2D_T_H_

/* base classes */
#include "SolidMatListT.h"
#include "SolidT.h"

namespace Tahoe {

/** small strain materials list for 2D structural analysis */
class SSSolidMatList2DT: public SolidMatListT, public SolidT
{
public:

	/** constructor */
	SSSolidMatList2DT(int length, const SSMatSupportT& support);
	SSSolidMatList2DT(void);

	/** read material data from the input stream */
	virtual void ReadMaterialData(ifstreamT& in);

	/** return true if the list contains plane stress models */
	virtual bool HasPlaneStress(void) const;

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** information about subordinate parameter lists */
	virtual void DefineSubs(SubListT& sub_list) const;

	/** return the description of the given inline subordinate parameter list */
	virtual void DefineInlineSub(const StringT& sub, ParameterListT::ListOrderT& order, 
		SubListT& sub_sub_list) const;

	/** a pointer to the ParameterInterfaceT of the given subordinate */
	virtual ParameterInterfaceT* NewSub(const StringT& list_name) const;
	/*@}*/

private:

	/** support for small strain materials */
	const SSMatSupportT* fSSMatSupport;

	/** support for gradient enhanced small strain materials */
	const GradSSMatSupportT* fGradSSMatSupport;
};

} /* namespace Tahoe  */

#endif /* _MATLIST_2D_T_H_ */
