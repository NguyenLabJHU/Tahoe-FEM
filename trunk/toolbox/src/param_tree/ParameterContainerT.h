/* $Id: ParameterContainerT.h,v 1.2 2004-01-21 17:05:40 paklein Exp $ */
#ifndef _PARAMETER_CONTAINER_T_H_
#define _PARAMETER_CONTAINER_T_H_

/* base classes */
#include "ParameterInterfaceT.h"
#include "ParameterListT.h"

namespace Tahoe {

/** standalone container of parameters. Class that both describes and contains a
 * user-definable list of parameters. Parameters are defined by the user through
 * the ParameterListT and then the class internally takes care of calls to the
 * ParameterInterfaceT interface. The only sublists which the container can
 * construct are those predefined for ParameterInterfaceT. All others will be
 * requested from the ParameterContainerT::fSubSource, if it has been set with
 * ParameterContainerT::SetSubSource. */
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

	/** set source for subs not defined by the container. The source must remain
	 * in scope at least as long as the container since the container may ask the
	 * source to construct any subs it does not define. Pass NULL to clear the
	 * sub source. */
	void SetSubSource(const ParameterInterfaceT* sub_source);

	/** a pointer to the ParameterInterfaceT of the given subordinate. If the container
	 * does not define the given sub, it will attempt to get it from ParameterContainerT::fSubSource,
	 * if the source has been defined. */
	virtual ParameterInterfaceT* NewSub(const StringT& list_name) const;

	/** return the description of the given inline subordinate parameter list.
     * If the container does not define the given sub, it will attempt to get 
     * it from ParameterContainerT::fSubSource, if the source has been defined. */
     virtual void DefineInlineSub(const StringT& sub, ParameterListT::ListOrderT& order, 
		SubListT& sub_sub_list) const;

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
	
	/** source for subs that are not defined by the container */
	const ParameterInterfaceT* fSubSource;
};

} /* namespace Tahoe */

#endif /* _PARAMETER_CONTAINER_T_H_ */
