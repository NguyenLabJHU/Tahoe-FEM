/* $Id: ParameterInterfaceT.h,v 1.2.2.3 2003-04-28 17:05:30 paklein Exp $ */
#ifndef _PARAMETER_INTERFACE_T_H_
#define _PARAMETER_INTERFACE_T_H_

/* direct members */
#include "ParameterListT.h"

namespace Tahoe {

/* forward declarations */
class StringT;

/** abstract interface for classes which define and use parameters */
class ParameterInterfaceT
{
public:

	/** constructor */
	ParameterInterfaceT(const StringT& name);

	/** identifier */
	const StringT& Name(void) const { return fName; };

	/** \name parameters */
	/*@{*/
	/** accept completed parameter list */
	virtual void SetParameters(const ParameterListT& list);
	
	/** build parameter list description.
	 * \param list destination for parameter description. The list will have the
	 *        name either of ParameterInterfaceT::Name or of any sub-list, returned
	 *        by ParameterInterfaceT::SubNames that is defined as inline. */
	virtual void DefineParameters(ParameterListT& list) const;
	/*@}*/

	/** \name subordinates that define parameters lists
	 * There are two types of subordinate lists, inlined and non-inlined. Non-inlined
	 * subordinates are returned with a call to ParameterInterfaceT::NewSub. Inlined
	 * subordinates are defined by the call to ParameterInterfaceT::DefineInlineSub. */
	/*@{*/
	/** information about subordinate parameter lists
	 * \param names list of subordinate list names
	 * \param occur occurrence specifier of subordinate list names 
	 * \param is_inline flag indicating if list is inline */
	virtual void SubNames(ArrayT<StringT>& names, ArrayT<ParameterListT::OccurrenceT>& occur,
		ArrayT<bool>& is_inline) const;

	/** return the description of the given inline subordinate parameter list.
	 * Method will be called for each subordinate defined as inline by ParameterInterfaceT::SubNames
	 * or defined recursively by ParameterInterfaceT::DefineInlineSub.
	 * \param sub name of the inlined subordinate list
	 * \param order defines whether list is a sequence or choice.
	 * \param names list of subordinate list names
	 * \param occur occurrence specifier of subordinate list names 
	 * \param is_inline flag indicating if list is inline */
	virtual void DefineInlineSub(const StringT& sub, ParameterListT::ListOrderT& order, ArrayT<StringT>& names, 
		ArrayT<ParameterListT::OccurrenceT>& occur, ArrayT<bool>& is_inline) const;

	/** a pointer to the ParameterInterfaceT of the given subordinate
	 * or NULL if the name is invalid. Responsibility for deleteting instantiations
	 * resides with the client who requested them. */
	virtual ParameterInterfaceT* NewSub(const StringT& list_name) const;
	/*@}*/

private:

	StringT fName;
};

} /* namespace Tahoe */

#endif /* _PARAMETER_SOURCE_T_H_ */
