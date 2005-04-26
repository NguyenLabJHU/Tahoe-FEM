/* $Id: MFGPSSSolidMatList3DT.h,v 1.1 2005-04-26 22:27:01 kyonten Exp $ */
#ifndef _MFGP_SS_MATLIST_3D_T_H_
#define _MFGP_SS_MATLIST_3D_T_H_

/* base class */
#include "MFGPSolidMatListT.h"
#include "SolidT.h"

#include "MFGPMatSupportT.h"

namespace Tahoe {

/* forward declarations */
class MFGPSSSolidMatT;

/** materials list for 3D structural analysis */
class MFGPSSSolidMatList3DT: public MFGPSolidMatListT, public SolidT
{
public:

	/** constructors */
	MFGPSSSolidMatList3DT(int length, const MFGPMatSupportT& support);
	MFGPSSSolidMatList3DT(void);

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
	MFGPSSSolidMatT* NewMFGPSSSolidMat(const StringT& name) const;

private:

	/** support for finite strain materials */
	const MFGPMatSupportT* fMFGPMatSupport;
	
};

} // namespace Tahoe 
#endif /* _MFGP_SS_MATLIST_3D_T_H_ */
