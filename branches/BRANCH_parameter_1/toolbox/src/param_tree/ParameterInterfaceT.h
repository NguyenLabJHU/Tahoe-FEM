/* $Id: ParameterInterfaceT.h,v 1.2.2.4 2003-05-03 09:06:52 paklein Exp $ */
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
	/** build parameter list description.
	 * \param list destination for parameter description. The list will have the
	 *        name either of ParameterInterfaceT::Name or of any sub-list, returned
	 *        by ParameterInterfaceT::SubNames that is defined as inline. */
	virtual void DefineParameters(ParameterListT& list) const;

	/** extract validated parameters. Take a raw list of parameters and produce 
	 * a validated parameter list. If the validated list cannot be 
	 * produced for any reason, the class throws a ExceptionT::kBadInputValue 
	 * \param raw_list source list in which all parameters are stored as
	 *        strings, as read from a source file. 
	 * \param valid_list returns as a validated list witb values of the appropriate data 
	 *        type, validating against constraints and applying any unspecified default values. */
	virtual void ValidateParameters(const ParameterListT& raw_list, ParameterListT& valid_list) const;

	/** accept completed parameter list */
	virtual void SetParameters(const ParameterListT& list);
	/*@}*/

	/** \name subordinates that define parameters lists
	 * There are two types of subordinate lists, inlined and non-inlined. Non-inlined
	 * subordinates are returned with a call to ParameterInterfaceT::NewSub. Inlined
	 * subordinates are defined by the call to ParameterInterfaceT::DefineInlineSub. */
	/*@{*/
	/** the order of subordinate lists */
	virtual ParameterListT::ListOrderT ListOrder(void) const;
	
	/** information about subordinate parameter lists
	 * \param order defines whether list is a sequence or choice
	 * \param names list of subordinate list names
	 * \param occur occurrence specifier of subordinate list names 
	 * \param is_inline flag indicating if list is inline */
	virtual void SubNames(ArrayT<StringT>& names, ArrayT<ParameterListT::OccurrenceT>& occur,
		ArrayT<bool>& is_inline) const;

	/** return the description of the given inline subordinate parameter list.
	 * Method will be called for each subordinate defined as inline by ParameterInterfaceT::SubNames
	 * or defined recursively by ParameterInterfaceT::DefineInlineSub. Nested inlines are
	 * not supported.
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
