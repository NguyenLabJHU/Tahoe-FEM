/* $Id: ParameterContainerT.h,v 1.1 2003-11-01 21:10:37 paklein Exp $ */
#ifndef _PARAMETER_CONTAINER_T_H_
#define _PARAMETER_CONTAINER_T_H_

/* base classes */
#include "ParameterInterfaceT.h"
#include "ParameterListT.h"

namespace Tahoe {

/** standalone container of parameters. Class that both describes and contains a
 * user-definable list of parameters. Parameters are defined by the user through
 * the ParameterListT and then the class internally takes care of calls to the
 * ParameterInterfaceT interface. The only sublists which are supported are
 * those predefined for ParameterInterfaceT. */
class ParameterContainerT: public ParameterListT, public ParameterInterfaceT
{
public:

	/** constructor */
	ParameterContainerT(const StringT& name);

	/** \name identifier */
	/*@{*/
	const StringT& Name(void) const { return ParameterInterfaceT::Name(); };
	void SetName(const StringT& name);
	/*@}*/

	/** \name add a sublist */
	/*@{*/
	void AddSub(const StringT& name, 
		ParameterListT::OccurrenceT occur = ParameterListT::Once, 
		bool is_inline = false); 
	void AddSub(const SubListDescriptionT& sub);
	/*@}*/

protected:

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** describe the parameters needed by the interface. */
	virtual void DefineParameters(ParameterListT& list) const;

	/** information about subordinate parameter lists
	 * \param sub_lists description of subordinate parameter lists */
	virtual void DefineSubs(SubListT& sub_list) const;
	/*@}*/

protected:

	/** sublists registered by ParameterContainerT::AddSub */
	SubListT fSubs;
};

} /* namespace Tahoe */

#endif /* _PARAMETER_CONTAINER_T_H_ */
