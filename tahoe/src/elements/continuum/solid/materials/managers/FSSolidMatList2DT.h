/* $Id: FSSolidMatList2DT.h,v 1.1.6.1 2004-07-06 06:53:52 paklein Exp $ */
/* created: paklein (02/14/1997) */
#ifndef _MATLIST_2D_T_H_
#define _MATLIST_2D_T_H_

/* base classes */
#include "SolidMatListT.h"
#include "SolidT.h"

namespace Tahoe {

/* forward declaration */
class FSSolidMatT;

/** materials list for 2D structural analysis */
class FSSolidMatList2DT: public SolidMatListT, public SolidT
{
public:

	/** constructor */
	FSSolidMatList2DT(int length, const FSMatSupportT& support);
	FSSolidMatList2DT(void);

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

	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/

	/** construct the specified material or NULL if the request cannot be completed */
	FSSolidMatT* NewFSSolidMat(const StringT& name) const;

private:

	/** support for finite strain materials */
	const FSMatSupportT* fFSMatSupport;
};

} // namespace Tahoe 
#endif /* _MATLIST_2D_T_H_ */
