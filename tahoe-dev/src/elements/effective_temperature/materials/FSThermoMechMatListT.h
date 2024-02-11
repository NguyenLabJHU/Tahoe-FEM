/* $Header: /home/regueiro/tahoe_cloudforge_repo_snapshots/development/src/elements/effective_temperature/materials/FSThermoMechMatListT.h,v 1.1 2015-09-20 03:48:45 tahoe.vickynguyen Exp $ */
/* created: rxiao (01/01/2014) */
#ifndef _FSThermo_MAT_LIST_T_H_
#define _FSThermo_MAT_LIST_T_H_

/* base class */
#include "SolidMatListT.h"
#include "SolidT.h"

namespace Tahoe {

/* forward declarations */
class FSSolidMatT;
//class SMP_coupled;
class FSThermoMechSupportT;
class FSThermoMechMatT;

/** list of materials for fluid analysis */
class FSThermoMechMatListT: public SolidMatListT, public SolidT
{
public:

	/** enum defining material types */
/*	enum TypeT {
		kLinear = 1,
		kNonLinear = 2
	}; */

	/** constructors */
	FSThermoMechMatListT(int length, const FSThermoMechSupportT& support);
	FSThermoMechMatListT(void);

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
	FSThermoMechMatT* NewFSThermoMechMat(const StringT& name) const;

private:

	/** support for thermomech materials */
	const FSThermoMechSupportT* fFSThermoMechSupport;

	
};

} // namespace Tahoe 
#endif /* _FLUID_MAT_LIST_T_H_ */
