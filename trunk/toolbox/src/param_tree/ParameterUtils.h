/* $Id: ParameterUtils.h,v 1.1 2003-08-14 01:22:03 paklein Exp $ */
#ifndef _PARAMETER_UTILS_H_
#define _PARAMETER_UTILS_H_

/* direct members */
#include "ParameterInterfaceT.h"

namespace Tahoe {

/** list of integers. Defines a list of integers with default ParameterInterfaceT name 
 * "IntegerList" which contains zero or more "Integer" entries. */
class IntegerListT: public ParameterInterfaceT
{
public:
	/** \name constructors */
	/*@{*/
	IntegerListT(const StringT& name);
	IntegerListT(void);
	/*@}*/

	/** the list name */
	const StringT& ListName(void) const { return fListName; };

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** describe the parameters needed by the interface */
	virtual void DefineParameters(ParameterListT& list) const;

	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);

	/** information about subordinate parameters */
	virtual void DefineSubs(SubListT& sub_list) const;

	/** a pointer to the ParameterInterfaceT of the given subordinate */
	virtual ParameterInterfaceT* NewSub(const StringT& list_name) const;
	/*@}*/

private:

	/** list name */
	StringT fListName;
};

/** an integer with default ParameterInterfaceT name "Integer" */
class IntegerT: public ParameterInterfaceT
{
public:

	/** \name constructors */
	/*@{*/
	IntegerT(void);
	IntegerT(const StringT& name);
	/*@{*/

	/** the value name */
	const StringT& ValueName(void) const { return fValueName; };

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** describe the parameters needed by the interface */
	virtual void DefineParameters(ParameterListT& list) const;

	/** accept parameter list.
	 * \param list input parameter list, which should be validated using ParameterInterfaceT::ValidateParameterList
	 *        to ensure the list conforms to the description defined by the interface. */
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/

private:

	/** integer name */
	StringT fValueName;

	/** integer value */
	int fValue;
};

} /* namespace Tahoe */

#endif /* _PARAMETER_UTILS_H_ */
